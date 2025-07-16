monkey ='Gilligan';
sessionDate ='06_03_2019';
msWind = [-30 50];
trigs = -40:1:40;
testWindow = [5 15];
preWindow = [-20 -5];
smoothKernel = 5;
minSpikes = 50;
spikeORstim = "stim";
plotISA = true;
saveDir = 'S:\Lab\ngc14\Working\EMG_UNITS\Stim_Triggered\';
allmuscles = [{'Deltoid'}, {'Biceps'},{'Triceps'},{'Wrist extensor'},{'Wrist flexor'},{'Digit extensor'},{'Digit flexor'}];
if(strcmp(monkey, 'Gilligan'))
    sessionDateF = char(datetime(sessionDate,'InputFormat','MM_dd_yyyy','Format','MM_dd_yyyy'));
else
    sessionDateF = char(datetime(sessionDate,'InputFormat','yyyy_MM_dd','Format','yyyy_MM_dd'));
end
sessionPath = ['S:\Lab\',monkey,'\All Data\', monkey,'_', sessionDateF];
EMGFile = dir([sessionPath, '\EMG\*.txt']);
physDir = [sessionPath, '\Physiology'];
%% get spike and trial times
if(strcmpi(spikeORstim,"stim"))
    nFiles = dir([sessionPath, '\EMG\']);
    nevFilesInd = find(cellfun(@(a) contains(a,'stim') & any(cellfun(@(s)strcmp(s,'nev'),...
        regexp(a, '\.','split'))), {nFiles.name}));
    nevFiles = nFiles(nevFilesInd);
    [~,mostRecentOrder] = sort(cellfun(@datenum,{nevFiles.date}),'descend');
    nevInd = 4;
    [~, hFile] = ns_OpenFile([physDir,'\',nevFiles(mostRecentOrder(nevInd)).name]);
    eventEntityID = find(cellfun(@(a) strcmpi(a,'Event'),{hFile.Entity.EntityType}));
    if(isempty(eventEntityID))
        eventEntityID = find(cellfun(@(a) a>5121 & a<10241,{hFile.Entity.ElectrodeID}));
    end
    stimElectrodes = length(eventEntityID);
    allLabels = {hFile.Entity.Label};
    allLabels = allLabels(~cellfun(@isempty,allLabels));
    channelNames = unique(allLabels(cellfun(@(a) contains(a,'Nano2_stim'), allLabels)));
    [spikeSegs,spikeChannels] = deal(cell(1,length(channelNames)));
    for e = 1:stimElectrodes
        [~, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(e));
        numCount = eventEntityInfo.ItemCount;
        ch = eventEntityInfo.EntityLabel;
        [eventTimeStamps]= deal(NaN(1,numCount));
        eventData = cell(1,numCount);
        for i = 1:numCount
            if(strcmpi(hFile.Entity(eventEntityID(e)).EntityType,"Segment"))
                [~,eventTimeStamps(i),~,~,~] = ns_GetSegmentData(hFile,eventEntityID(e),i);
            else
                [~, eventTimeStamps(i),eventData{i}, ~] = ns_GetEventData(hFile, eventEntityID(e), i);
            end
        end
        if(strcmpi(hFile.Entity(eventEntityID(e)).EntityType,"Segment"))
            trainStart = find([1,diff(eventTimeStamps)>.05]==1);
        else
            eventTimeStamps = eventTimeStamps(cell2mat(eventData)>0);
            trainStart = find([1,diff(eventTimeStamps)>.2]==1);
        end
        spikeChannels{str2double(ch(end-1:end))} = eventTimeStamps(trainStart);
        spikeSegs{str2double(ch(end-1:end))} = eventTimeStamps;
    end
    ns_CloseFile(hFile);
else
    resPath = [physDir,'\All_Trials\'];
    if(~exist(resPath,'dir'))
        Spike_SortRawData(sessionDateF,monkey);
    end
    sortedDir = dir(resPath);
    sortedDir = sortedDir(~[sortedDir.isdir]);
    [spikeChannels, spikeSegs] = deal(cell(1,length(sortedDir)));
    [~,~,alltrials,ardTrials,~,channel,~,labs] = getSessionInfo2(resPath,'single');
    if(isempty(labs))
        labelSingleUnits(monkey,datetime(sessionDateF,'InputFormat','MM_dd_yyyy'));%
        [~,~,trials,~,~,channel,~,labs] = getSessionInfo2(resPath,'single');
    end
    for c = 1:length(channel)
        currChannel = [monkey+"_"+sessionDateF+"_"+num2str(channel(c))+".mat"];
        sortedSpks = load([resPath+currChannel],'sortedSpikeData');
        moveConds = ~cellfun(@(a) strcmp(a,"Rest"), sortedSpks.sortedSpikeData.ArduinoData(:,1));
        for u = 1:sum(channel==channel(c))
            spikeChannels{u,channel(c)}= cellfun(@(ss,st) ss(ss>st(1) & ss<st(max(end-2,1))),...
                sortedSpks.sortedSpikeData.SpikeTimes(u,moveConds), sortedSpks.sortedSpikeData.SegTimes(u,moveConds),"UniformOutput", false);
            spikeSegs{u,channel(c)} = sortedSpks.sortedSpikeData.SegTimes(u,moveConds);
        end
    end
end
currDir = dir(physDir);
ns5FilesInd = find(cellfun(@(a) ~contains(a,"stim") & any(cellfun(@(s)strcmp(s,'ns5'),...
        regexp(a, '\.','split'))), {currDir.name}));
if(~isempty(ns5FilesInd))
    [~, hFileRaw] = ns_OpenFile([currDir(ns5FilesInd).folder, '\', currDir(ns5FilesInd).name]);
    rawSigInds = find([hFileRaw.Entity.FileType] == find(strcmp({hFileRaw.FileInfo.Type},'ns2')) & ...
                cellfun(@(s) strcmp(deblank(s),'uV'), {hFileRaw.Entity.Units}));
    [~,channelMap] = natsort({hFileRaw.Entity(rawSigInds).Label});
    if(strcmp(monkey,"Gilligan"))
        if(channelMap~=[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32])
            disp("Check channel map");
        end
    else
        if(channelMap~=[32:-1:1])
            disp("Check channel map");
        end
    end
    ns_CloseFile(hFileRaw);
end
spikeChannels = spikeChannels(:,channelMap);
spikeSegs = spikeSegs(:,channelMap);
[unitNum,spkChannel] = find(~cellfun(@isempty,spikeChannels));
%% clean spike data and align trials
spikeTimes = {};
trialSpikes = {};
for f = 1:length(spkChannel)
    if(strcmpi(spikeORstim,"spike"))
        tInds = ~cellfun(@(m,n) size(m,1)==0 || all(isnan(n)),...
            spikeChannels{unitNum(f),spkChannel(f)},spikeSegs{unitNum(f),spkChannel(f)});
        spikes = spikeChannels{unitNum(f),spkChannel(f)};
        spikeTimes{f} = spikes;%cellfun(@(n) all(isnan(n)),spikes));
        trialSpikes{f} = cell(1,length(tInds));
        trialSpikes{f}(tInds) = cellfun(@(a,b) a-b(2), [spikeChannels{unitNum(f),spkChannel(f)}(tInds)], ...
            [spikeSegs{unitNum(f),spkChannel(f)}(tInds)], 'UniformOutput', false);
    else
        spikeTimes{f} = num2cell(spikeChannels{unitNum(f),spkChannel(f)});
        trialSpikes{f} = num2cell(spikeSegs{unitNum(f),spkChannel(f)});
    end
end
nSpikes = cellfun(@(t) sum(cellfun(@length,t)), spikeTimes);
spikeTimes = spikeTimes(nSpikes>minSpikes);
trialSpikes = trialSpikes(nSpikes>minSpikes);
unitNum = unitNum(nSpikes>minSpikes);
spkChannel = spkChannel(nSpikes>minSpikes);
%% get EMG channels from txt file
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
channelNum = channelNum(present);
if(~exist(saveDir+"\"+sessionDateF,'dir'))
    mkdir(saveDir+"\"+sessionDateF);
    if(plotISA)
        mkdir(saveDir+"\"+sessionDateF+"\TRACE\");
        %mkdir(saveDir+"\"+sessionDateF+"\ISA\");
        %mkdir(saveDir+"\"+sessionDateF+"\SUB\");
    end
end
%% get EMG voltage trace and trial times
unitXmuscle = cell(length(spikeTimes),length(channelNum));
for s=1:length(spikeTimes)
    aSpikes = spikeTimes{s};
    aSegs = trialSpikes{s};
    if(plotISA)
        close all;
        f1 = figure("Units","normalized","Position",[0 0 1 1]);
        %f2 = figure("Units","normalized","Position",[0 0 1 1]);
        %f3 = figure("Units","normalized","Position",[0 0 1 1]);
    end
    tic
    for m = 1:length(channelNum)
        muscle = muscles{m};
        tSpikes = cell(1,length(aSpikes));
        if(strcmpi(spikeORstim,"stim"))
            nFiles = dir(['S:\Lab\',monkey,'\All Data\',monkey,'_', sessionDateF, '\EMG\']);
            nf3FilesInd = find(cellfun(@(a) any(cellfun(@(s)strcmp(s,'nf3'),...
                regexp(a, '\.','split'))), {nFiles.name}));
            nf3FilesInd = nf3FilesInd(contains(cellfun(@(e) extractBefore(extractAfter(e,"stim_"),characterListPattern(".")),{nFiles(nf3FilesInd).name},'UniformOutput',false),...
                extract(extractAfter(hFile.Name,"stim_"),digitsPattern)));
            if(~isempty(nf3FilesInd))
                nf3DotInd = regexp(nFiles(nf3FilesInd).name, '[_.]','split');
                nf3FileNumber = str2double(nf3DotInd{end-1});
                [voltData, dataTimes, Fs, ~, ~ ] = loadRippleData([nFiles(nf3FilesInd).folder,'\',nFiles(nf3FilesInd).name],channelNum(m));
                voltData = {voltData};
                cSpikes = aSpikes;
                nFrags = length(cSpikes);
                trigInd = 0;
                trigs = 0;
                startInd = 1;
                for p = 1:length(aSpikes)-1
                    [~,nextInd] = find(cell2mat(aSegs)==aSpikes{p+1});
                    tSpikes{p} = [aSegs{startInd:nextInd-1}]-aSpikes{p};
                    startInd = nextInd;
                end
                tSpikes{p+1} = [aSegs{startInd:end}]-aSpikes{end};
            end
        else
            if(~exist(['S:\Lab\',monkey,'\All Data\',monkey,'_', sessionDateF, '\EMG\All_Trials\'],'dir'))
                EMG_SortRawData(sessionDateF,monkey);
            end
            muscleEMG = load(['S:\Lab\',monkey,'\All Data\',monkey,'_', sessionDateF, '\EMG\All_Trials\',...
                monkey,'_',sessionDateF,'_',muscle,'.mat'],'sortedEMGData');
            muscleEMG = muscleEMG.sortedEMGData;
            Fs = muscleEMG.SampleRate;
            moveConds = ~cellfun(@(a) strcmp(a,"Rest"), muscleEMG.ArduinoData(:,1));
            voltData = muscleEMG.EMGData(moveConds);
            segTimes = muscleEMG.SegTimes(moveConds);
            badTrials = getBadTrials(voltData,cellfun(@(t) int64(1000*((t-t(1))./(1000/Fs))), segTimes, 'UniformOutput', false),Fs);
            dataTimes = cell2mat(cellfun(@(t) t(1):(1/Fs):t(end), segTimes(badTrials==0),'UniformOutput',false));
            voltData = (voltData(badTrials==0));
            cSpikes = aSpikes(badTrials==0);
            tSpikes = trialSpikes{s}(badTrials==0);
            nFrags = ceil(sqrt(length(cSpikes)));
            trigInd = trigs./(1000/Fs);
        end
        %% spike-triggered averaging
        cSpikes = cell2mat(reshape(cSpikes,[],1));
        voltData = cell2mat(cellfun(@(v)conv(abs(v),gausswin(smoothKernel./(1000/Fs))/...
            sum(gausswin(smoothKernel/(1000/Fs))),'same'),voltData,'UniformOutput',false));
        indexWindow = msWind./(1000/Fs);
        [~,~,testInd] = intersect(testWindow./(1000/Fs),indexWindow(1):indexWindow(end));
        [~,~,preInd] = intersect(preWindow./(1000/Fs),indexWindow(1):indexWindow(end));
        [isaAvg,muscleWinds] = deal(NaN(length(cSpikes),1+range(indexWindow)));
        %% spike segments
        parfor u = 1:length(cSpikes)
            isa = NaN(length(trigInd),1+(range(indexWindow)));
            [val,winInd] = min(abs(dataTimes -cSpikes(u)));
            if(indexWindow(1)+winInd+trigInd(1)>0 && indexWindow(end)+winInd+trigInd(end)<length(voltData) && val<(1/Fs))
                if(dataTimes(indexWindow(1)+winInd+trigInd(1))==dataTimes(winInd)+(trigs(1)+msWind(1))/1000 && ...
                        dataTimes(indexWindow(end)+winInd+trigInd(end))==dataTimes(winInd)+(trigs(end)+msWind(end))/1000)
                    muscleWinds(u,:) = voltData(winInd+indexWindow(1):winInd+indexWindow(end));
                    for t = 1:length(trigInd)
                        currTrig = winInd + trigInd(t);
                        currRange = currTrig+indexWindow(1):currTrig+indexWindow(end);
                        isa(t,:) = voltData(currRange);
                    end
                end
            end
            isaAvg(u,:)= mean(isa,'omitnan');
        end
        fragInds = round(linspace(1,length(cSpikes),nFrags));
        [isaFrags,muscleFrags] = deal(NaN(nFrags,1+range(indexWindow)));
        %% fragments for stats
        parfor f = 1:nFrags
            if(f==nFrags)
                muscleFrags(f,:) = mean(muscleWinds(fragInds(f):size(muscleWinds,1),:),1,'omitnan');
                isaFrags(f,:) = mean(isaAvg(fragInds(f):size(isaAvg,1),:),1,'omitnan');
            else
                muscleFrags(f,:) = mean(muscleWinds(fragInds(f):fragInds(f+1),:),1,'omitnan');
                isaFrags(f,:) = mean(isaAvg(fragInds(f):fragInds(f+1),:),1,'omitnan');
            end
        end
        subFrags = muscleFrags-isaFrags;
        dVals = mean(subFrags(:,testInd(1):testInd(end)),2,'omitnan') - mean(...
            subFrags(:,preInd(1):preInd(end)),2,'omitnan');
        [~,pval,~,stats] = ttest(dVals);
        unitXmuscle{s,m} = {muscleFrags, isaFrags};
        if(plotISA)
            set(0,'CurrentFigure',f1);plotISAFigs(m,spkChannel(s),tSpikes,muscle,muscleFrags,msWind);
            % set(0,'CurrentFigure',f2);plotISAFigs(m,spkChannel(s),tSpikes,muscle,isaFrags,msWind);
            % set(0,'CurrentFigure',f3);plotISAFigs(m,spkChannel(s),tSpikes,muscle+" p = "+...
            %     num2str(pval,'%0.3f'),subFrags,msWind);
            %plot(preInd,repmat(mean(subFrags(:,preInd(1):preInd(end)),'all','omitnan'),length(preInd)),'y','LineWidth',4);
            %plot(testInd,repmat(mean(subFrags(:,testInd(1):testInd(end)),'all','omitnan'),length(testInd)),'b','LineWidth',4);
        end
    end
    if(plotISA)
        saveFigures(f1,saveDir+string(sessionDateF)+"\TRACE\",sessionDateF+"_"+num2str(nf3FileNumber)+"_Cropped",[]);
        % saveFigures(f2,saveDir+string(sessionDateF)+"\ISA\",sessionDateF+"_"+num2str(spkChannel(s))+"_"+num2str(unitNum(s)),[]);
        % saveFigures(f3,saveDir+string(sessionDateF)+"\SUB\",sessionDateF+"_"+num2str(spkChannel(s))+"_"+num2str+"+(unitNum(s)),[]);
    end
    disp(toc);
end
%save([saveDir,sessionDateF,'\unitXmuscle.mat'],'unitXmuscle','spkChannel','unitNum','muscles','indexWindow');

function plotISAFigs(mind,unitN,spks,muscleTitle,frags,xLab)
if(mind==4)
    subplot(2,4,4);hold on;
    title(["Channel "+num2str(mind)+" Unit #"+num2str(unitN)+": "+ num2str(length(spks)/length(unique(cellfun(@length,spks)))) + " trials,"+...
        num2str(unique(cellfun(@length, spks)))+ " spikes"]);
    cellfun(@(x,y) arrayfun(@(xx)plot([xx,xx],[y-1,y],'k'), x(1:min(length(x),100))),...
        reshape(cellfun(@(m) round(m,3), spks, 'UniformOutput',false),1,[]),num2cell(1:length(spks)));
    xlim([-.005 .1]);
end
subplot(2,4,(5*double(mind>3))+mod(mind,4));hold on;
title(muscleTitle);
plot(frags');
plot(mean(frags',2,'omitnan'),'k','LineWidth',2);
ylim([min(0,min(mean(frags,'omitnan'))) max(30,max(mean(frags,1,'omitnan')))]);
xlim([0 size(frags,2)]);
xticks([0,abs(xLab(1)/range(xLab))*(size(frags,2)-1),size(frags,2)]);
xticklabels([xLab(1),0,xLab(2)]);
end