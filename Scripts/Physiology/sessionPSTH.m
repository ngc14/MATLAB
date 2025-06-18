date = '06_03_2019';
monkey = 'Gilligan';
binSize = .01;
sigma = 10; 
MIN_BLOCKS_FOR_UNIT = 15;
conds = ["Extra Small Sphere", "Large Sphere", "Photocell"];
alignSegsConds = {["StartReach"],["StartReach"],["StartReach"]};
alignLims = {[-.5, 2.5]};
saveDir = 'S:\Lab\ngc14\Working\EMG_UNITS\Stim_Triggered';
sessionPath = ['S:\Lab\',monkey,'\All Data\',monkey,'_',date,'\'];
allmuscles = [{'Deltoid'}, {'Biceps'},{'Triceps'},{'Wrist extensor'},{'Wrist flexor'},{'Digit extensor'},{'Digit flexor'}];
EMGFile = dir([sessionPath, '\EMG\*.txt']);
sessionPhys = dir([sessionPath,'Physiology\Results']);
sessionPhys = sessionPhys(~[sessionPhys.isdir]);
[~,sortInd] = natsort({sessionPhys.name});
sessionPhys = sessionPhys(sortInd);
%%
fileID = fopen([EMGFile.folder,'\',EMGFile.name]);
muscles = {};
channelNum = [];
while(~feof(fileID))
    line = fgetl(fileID);
    if(contains(line,'Channel'))
        if(contains(line, 'Channel B','IgnoreCase',true))
            channelNumStart = regexp(line, 'Channel B','end');
        elseif(contains(line, 'Channel A', 'IgnoreCase', true))
            channelNumStart = regexp(line, 'Channel A','end');
        end
        channelNumEnd = regexp(line, ':');
        channelNum(end+1) = str2double(line(channelNumStart+1:channelNumEnd-1)) + ...
            (128*double(contains(line, 'Channel B','IgnoreCase',true)));
        muscles{end+1} = line(channelNumEnd+2:end);
    end
end
fclose(fileID);
muscles = muscles(~cellfun(@isempty,muscles) & cellfun(@(m) any(isletter(m)),muscles));
channelNum = unique(channelNum);
channelNum = channelNum(~cellfun(@isempty,muscles) & cellfun(@(m) any(isletter(m)),muscles));
[~,present]= ismember(muscles,allmuscles);
[~,present] = sort(present);
muscles = muscles(present);
if(~exist([sessionPath,'\EMG\All_Trials\'],'dir'))
    EMG_SortRawData(date,monkey);
end
%%
[allUnits,avgSegs] = deal(cell(1,length(conds)));
for f = 1:length(sessionPhys)
    m = matfile_m([sessionPhys(f).folder, '\', sessionPhys(f).name]);
    m = m.sortedSpikeData;
    if(size(m.SegTimes,1)>0)
        load([sessionPhys(f).folder,'\',sessionPhys(f).name]);
        allGoodTrials = ~(cellfun(@(a,b) length(a) <=(b(end)-b(1)) ...
            | length(a)>200*(b(end)-b(1)) | any(isnan(a)), sortedSpikeData.SpikeTimes, sortedSpikeData.SegTimes));
        blockInds = cumsum(mod(1:length(sortedSpikeData.ArduinoData),length(conds))==1);
        unitTrials = num2cell(allGoodTrials.*blockInds,2);
        goodUnitsOnChannel = find(cellfun(@(tc) length(unique(tc))-1, unitTrials)>MIN_BLOCKS_FOR_UNIT);
        for u = 1:length(goodUnitsOnChannel)
            [alignedSig,avgSegs] = condSignalAligned(conds,alignSegsConds,...
                sortedSpikeData,'SpikeTimes',goodUnitsOnChannel(u));
            trialHistsSmooth = cellfun(@(a) cellfun(@(ha,al)cell2mat(cellfun(@(h) ...
                conv(histcounts(h,al(1):binSize:al(end))./binSize,gausswin(sigma)/sum(gausswin(sigma)),'same'),...
                ha,'UniformOutput',false)'),a,alignLims,'UniformOutput',false),alignedSig, 'UniformOutput', false);
            if(all(cellfun(@isempty,allUnits)))
                allUnits = trialHistsSmooth;
            else
                allUnits = cellfun(@(a,ts) [a,ts],allUnits,trialHistsSmooth,'UniformOutput',false);
            end
        end
    end
end
%%
allEMGs = cell(1,length(conds));
for m = 1:length(muscles)
    muscleEMG = load([sessionPath, '\EMG\All_Trials\',monkey,'_',date,'_',muscles{m},'.mat'],'sortedEMGData');
    muscleEMG = muscleEMG.sortedEMGData;
    Fs = muscleEMG.SampleRate;
    smoothKernel = sigma./(1000/Fs);
    [alignedTimes,allSegTimes] = condSignalAligned(conds,alignSegsConds,muscleEMG,'SegTimes',1);
    [alignedSig,allSegs] = condSignalAligned(conds,cell(1,length(conds)),muscleEMG,'EMGData',1);
    alignedTimes = cellfun(@(c) {cellfun(@(a) cellfun(@(s) s(1):(1/Fs):s(end),a,'UniformOutput',false),...
        c,'UniformOutput',false)},alignedTimes);
    alignedTimes = cellfun(@(s,a,at) {at{1}(cellfun(@(t) find(cellfun(@(s2) isequal(t,s2), a)), s))},...
        allSegs,allSegTimes,alignedTimes,'UniformOutput',false);
    allSegs = cellfun(@(c) cellfun(@transpose,c,'UniformOutput',false),allSegs,'UniformOutput',false);
    goodTrials = cellfun(@(s,t) getBadTrials(s{1},cellfun(@(tt) uint64(1000*(tt-tt(1))./(1000/Fs)),t,'UniformOutput',false),Fs)==0,...
        alignedSig, allSegs, 'UniformOutput', false);
     alignedEMG = cellfun(@(a,at,gt) cellfun(@(ha,as,al) cellfun(@(h,l) conv(abs(h(...
         (l>=al(1) & l<=al(end)))),gausswin(smoothKernel)/sum(gausswin(smoothKernel)),...
         'same'),ha(gt),as(gt),'UniformOutput',false),a,at,alignLims,'UniformOutput',false),alignedSig,alignedTimes,goodTrials,'UniformOutput',false);
     voltData = cellfun(@(a) cellfun(@(ha,al) cell2mat(cellfun(@(h) [h,NaN(1,length(al(1):(1/Fs):al(end))-1-length(h))],...
         ha,'UniformOutput',false)'),a,alignLims,'UniformOutput',false),alignedEMG,'UniformOutput',false);
     if(all(cellfun(@isempty,allEMGs)))
         allEMGs = voltData;
     else
         allEMGs = cellfun(@(a,ts) [a,ts],allEMGs,voltData,'UniformOutput',false);
     end
end
%%
gapWind = 0.10;
close all; 
unitsPerPlot = 3;
allColors = cellfun(@(u) distinguishable_colors(length(u)+length(muscles)),allUnits,'UniformOutput',false);
unitColors = cellfun(@(c) c(1:(end-length(muscles)),:), allColors, 'UniformOutput',false);
muscleColors = cellfun(@(c) c((size(allUnits{1},2)+1):end,:), allColors, 'UniformOutput',false);
unitCols = ceil(size(allUnits{1},2)/unitsPerPlot);
set(0, 'DefaultFigureRenderer', 'painters');
figure();
f=tiledlayout(length(alignSegsConds),unitCols+2,"TileIndexing","columnmajor");
%normVals = max(cell2mat(cellfun(@(cp) cellfun(@(p) cellfun(@(n) max(mean(n,1,'omitnan')),p), cp(~isempty(cp)),'UniformOutput',false), allUnits)'))';
for a = 1:unitCols
    startInd = (unitsPerPlot*(a-1))+1;
    uInds = startInd:min(size(allUnits{1},2),startInd+(unitsPerPlot-1));
    currUnits = cellfun(@(al) {al(uInds)}, allUnits, 'UniformOutput',false);
    colors = cellfun(@(ac) ac(uInds,:),unitColors, 'UniformOutput',false);
    f=plotPSTH(currUnits,alignSegsConds,avgSegs,sortedSpikeData.ConditionSegments,alignLims,binSize,colors,gapWind);
    ch = flipud(get(f.Children,'Children'));
    if(a==1)
        arrayfun(@(n) ylabel(ch(n),conds(n)), 1:length(conds));
    end
    cols = find(arrayfun(@(m) mod(m,unitsPerPlot)-1,1:length(ch))==0);
    ch(cols(a)).Title.String = "Units "+num2str(uInds);
end
%%
armMuscles = cellfun(@(s) ismember(s,allmuscles(1:3)),muscles);
f=plotPSTH(cellfun(@(m) {m(armMuscles)},allEMGs,'UniformOutput',false),alignSegsConds,...
   allSegs,muscleEMG.ConditionSegments,alignLims,1/Fs,cellfun(@(m) m(armMuscles,:),muscleColors,'UniformOutput',false),gapWind);
lh = legend;
lh.String = lh.String(end-sum(armMuscles):end-1);
lh.String = allmuscles(armMuscles);
f=plotPSTH(cellfun(@(m) {m(~armMuscles)},allEMGs,'UniformOutput',false),alignSegsConds,...
    allSegs,muscleEMG.ConditionSegments,alignLims,1/Fs,cellfun(@(m) m(~armMuscles,:),muscleColors,'UniformOutput',false),gapWind);
lh = legend;
lh.String = lh.String(end-sum(~armMuscles):end-1);
lh.String = allmuscles(~armMuscles);
allChildren = flipud(get(f.Children,'Children'));
allChildren = allChildren(cellfun(@(r) strcmp(r,'axes'),get(allChildren,'Type')));
cols = find(arrayfun(@(m) mod(m,unitsPerPlot)-1,1:length(allChildren))==0);
allChildren(cols(end-1)).Title.String = "Arm";
allChildren(cols(end)).Title.String = "Hand";
saveFigures(gcf,saveDir+"\"+string(date)+"\","PSTH+EMGs",[])



function [alignedSig,avgSegs] = condSignalAligned(conds,alignSegsConds,sigStruct,fieldName,ind)
signals = sigStruct.(fieldName);
avgSegs = cell(1,length(conds));
alignedSig = cellfun(@(c) cell(1,length(c)),alignSegsConds,'UniformOutput',false);
for c = 1:length(conds)
    alignSegs =  alignSegsConds{c};
    currTrialInds = find(cellfun(@(a) contains(a, conds{c}),sigStruct.ArduinoData(:,1)));
    currTrialInds = currTrialInds(cellfun(@(a) ~any(isnan(a)), signals(ind,currTrialInds)));
    currSig = signals(ind,currTrialInds);
    currSegs = sigStruct.SegTimes(ind,currTrialInds);
    avgSegs{c} = currSegs;
    alignedSig{c}{1} = currSig;
    for a = 1:length(alignSegs)
        alignInd = find(strcmp(sigStruct.ConditionSegments{c}, alignSegs(a)));
        alignTimes = cellfun(@(s) getAlignedTimes(s,alignInd), currSegs, 'UniformOutput', false);
        alignedSig{c}{a} = cellfun(@minus, currSig, alignTimes, 'UniformOutput', false);
    end
end
end
function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end