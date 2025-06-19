date = '06_06_2019';
monkey = 'Gilligan';
binSize = .01;
sigma = 10; 
MIN_BLOCKS_FOR_UNIT = 15;
conds = ["Extra Small Sphere", "Large Sphere", "Photocell"];
alignSegsConds = {["StartReach"],["StartReach"],["StartReach"]};
alignLims = {[-.5, 2.5]};
saveDir = 'S:\Lab\ngc14\Working\EMG_UNITS\Stim_Triggered\';
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
[allUnits,allSegs] = deal(cell(1,length(conds)));
taskAlign = containers.Map(conds,{{["GoSignal" "StartLift"]},{["GoSignal","StartLift"]},...
    {["GoSignal","StartHold"]}});
params = PhysRecording(conds,binSize,sigma/1000,-6,5);
tUnits = [];
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
                sortedSpikeData,'SpikeTimes');
            trialHistsSmooth = cellfun(@(a) cellfun(@(ha,al)cell2mat(cellfun(@(h) ...
                conv(histcounts(h,al(1):binSize:al(end))./binSize,gausswin(sigma)/sum(gausswin(sigma)),'same'),...
                ha,'UniformOutput',false)'),a,alignLims,'UniformOutput',false),alignedSig, 'UniformOutput', false);

            [taskBaseline,taskFR] = calculatePhases(params,taskAlign,{[0.2 0]},...
                cellfun(@(c) cellfun(@(t) {t-t(2)},c,'UniformOutput',false), avgSegs,'UniformOutput',false),...
                cellfun(@(a)cellfun(@(ha,al)cellfun(@(h) ...
                {conv(histcounts(h,[params.bins,max(params.bins)+binSize])./binSize,gausswin(sigma)/sum(gausswin(sigma)),'same')},...
                ha,'uniformoutput',false),a,alignLims','UniformOutput',false),alignedSig),...
                false,true);
            [vals,tUnit] = cellfun(@(b,cn) ttestTrials({{cellfun(@cell2mat,[b{:}])}},...
                {{cellfun(@cell2mat,[cn{:}])}},1,true,0.05),...
                taskBaseline,taskFR,'UniformOutput',false);
            tUnits(end+1) = any(cell2mat(tUnit));
            if(all(cellfun(@isempty,allUnits)))
                allUnits = trialHistsSmooth;
                allSegs = avgSegs;
            else
                allUnits = cellfun(@(a,ts) [a,ts],allUnits,trialHistsSmooth,'UniformOutput',false);
                allSegs = cellfun(@(a,ts) [a,ts],allSegs,avgSegs,'UniformOutput',false);
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
    [alignedTimes,allSegTimes] = condSignalAligned(conds,alignSegsConds,muscleEMG,'SegTimes');
    [alignedSig,allSegs] = condSignalAligned(conds,cell(1,length(conds)),muscleEMG,'EMGData');
    alignedTimes = cellfun(@(c) {cellfun(@(a) cellfun(@(s) s(1):(1/Fs):s(end),a,'UniformOutput',false),...
        c,'UniformOutput',false)},alignedTimes);
    alignedTimes = cellfun(@(s,a,at) {at{1}(cellfun(@(t) find(cellfun(@(s2) isequal(t,s2), a)), s))},...
        allSegs,allSegTimes,alignedTimes,'UniformOutput',false);
    allSegs = cellfun(@(c) cellfun(@transpose,c,'UniformOutput',false),allSegs,'UniformOutput',false);
    goodTrials = cellfun(@(s,t) getBadTrials(s{1},cellfun(@(tt) uint64(1000*(tt-tt(1))./(1000/Fs)),t,'UniformOutput',false),Fs)==0,...
        alignedSig, allSegs, 'UniformOutput', false);
     alignedEMG = cellfun(@(a,at,gt) cellfun(@(ha,as,al) cellfun(@(h,l) conv(abs(h(...
         (l>=al(1) & l<=al(end)))),gausswin(smoothKernel)/sum(gausswin(smoothKernel)),...
         'same'),ha(),as(),'UniformOutput',false),a,at,alignLims,'UniformOutput',false),alignedSig,alignedTimes,goodTrials,'UniformOutput',false);
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
f=tiledlayout(length(alignSegsConds),unitCols+3,"TileIndexing","columnmajor");
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
normVals = max(cell2mat(cellfun(@(cp) cellfun(@(n) max(mean(n,1,'omitnan')), cp(~cellfun(@isempty,cp))),allUnits, 'UniformOutput',false)'),[],1);
f=plotPSTH(cellfun(@(m)cellfun(@(n,v) {mean(n,1,'omitnan')./v},m,num2cell(normVals),'UniformOutput',false),allUnits,'UniformOutput',false),alignSegsConds,...
   avgSegs,sortedSpikeData.ConditionSegments,alignLims,binSize,{[1 0 0],[224,144,38]./255, [0 0 1]},gapWind);
ch = flipud(get(f.Children,'Children'));
cols = find(arrayfun(@(m) mod(m,unitsPerPlot)-1,1:length(ch))==0);
ch(cols(end)).Title.String = "Average of Units";
%%
armMuscles = cellfun(@(s) ismember(s,allmuscles(1:3)),muscles);
normEMG = max(cell2mat(cellfun(@(cp) cellfun(@(n) max(mean(n,1,'omitnan')), cp(~cellfun(@isempty,cp))),allEMGs, 'UniformOutput',false)'),[],1);
allEMGs = cellfun(@(m) cellfun(@(a,v) a./v,m,num2cell(normEMG),'UniformOutput',false), allEMGs,'UniformOutput',false);
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
%%
figure();
tiledlayout(1,length(conds))
for c = 1:length(conds)
    ax{c} = nexttile();
    emg = allEMGs{c};
    psth = allUnits{c};
    hold on;
    [lagMap,corrMap] = deal([]);
    title(conds(c));
    for p = 1:size(psth,2)
        for m = 1:size(emg,2)
            % [rvals,lags] = xcorr(conv(interp1(1:size(psth{p},2),mean(psth{p},1,'omitnan'),linspace(1,size(psth{p},2),size(emg{m},2)),'pchip'),...
            %     gausswin(sigma)/sum(gausswin(sigma)),'same'),mean(emg{m},1,'omitnan'),400,'normalized');
            [rvals,lags] = cellfun(@(a,b) xcorr(conv(interp1(1:size(a,2),a,linspace(1,size(a,2)+1,size(b,2)),'pchip'),...
                gausswin(sigma)/sum(gausswin(sigma)),'same'),b,400,'normalized'),...
                num2cell(psth{p},2),num2cell(emg{m},2),'UniformOutput',false);
            maxVals = cellfun(@max,rvals);
            corrMap(p,m)=mean(maxVals,'omitnan');
        end
    end
    imagesc(corrMap);
    xticks(1:m);
    yticks(1:p);
    xlim([.5,m+.5]);
    ylim([.5,p+.5]);
    xticklabels(muscles);
    if(c==length(conds))
    colorbar;
    end
    colormap('jet');
    clim([0.5 1]);
    title(conds(c));
end


function [alignedSig,avgSegs] = condSignalAligned(conds,alignSegsConds,sigStruct,fieldName)
signals = sigStruct.(fieldName);
avgSegs = cell(1,length(conds));
alignedSig = cellfun(@(c) cell(1,length(c)),alignSegsConds,'UniformOutput',false);
for c = 1:length(conds)
    alignSegs =  alignSegsConds{c};
    currTrialInds = find(cellfun(@(a) contains(a, conds{c}),sigStruct.ArduinoData(:,1)));
    currTrialInds = currTrialInds(cellfun(@(a) ~any(isnan(a)), signals(:,currTrialInds)));
    currTrialInds = true(1,length(currTrialInds));
    currSig = signals(:,currTrialInds);
    currSegs = sigStruct.SegTimes(:,currTrialInds);
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