function done = labelSingleUnits(processDates,monkey)
subDirName = '\Results_New';
writeLabels = true;
plotWaveforms = false;
if(~exist('monkey','var'))
    monkey='Gilligan';
end
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'MM_dd_uuuu';
else
    dateFormat = 'uuuu_MM_dd';
end
if(~exist('processDates','var'))
    dateArray = datetime(2019,05,20):datetime(2019,12,12);
else
    dateArray = datetime(processDates,'InputFormat',dateFormat);
end
%%
for d = 1:length(dateArray)
    currName = ['S:\Lab\',monkey,'\All Data\', monkey,'_', char(datetime(dateArray(d),'InputFormat' ,dateFormat,'Format',dateFormat)),'\Physiology'];
    if(exist(currName, 'dir'))
        currDir = dir(currName);
        names = {currDir.name};
        sortedFile = find((contains(names, '-sorted') | contains(names, 'RASort'))...
            & ~contains(names, 'stim') &  (contains(names, '.nev') | contains(names, '.plx')));
        ns5FilesInd = find(cellfun(@(a) ~contains(a,"stim") & ...
            contains(a, '.ns5'), names));
        if(~isempty(ns5FilesInd))
            [~,mostRecent]= sort({currDir(ns5FilesInd).date});
            ns5FilesInd = ns5FilesInd(mostRecent);
            ns5FilesInd = ns5FilesInd(end);
            rawFile = [currDir(ns5FilesInd).folder, '\', currDir(ns5FilesInd).name];
            [~, hFileRaw] = ns_OpenFile(rawFile);
            rawSigInds = find([hFileRaw.Entity.FileType] == find(strcmp({hFileRaw.FileInfo.Type},'ns5')) & ...
                cellfun(@(s) strcmp(deblank(s),'uV'), {hFileRaw.Entity.Units}));
            rawSigInds = rawSigInds(sort(cellfun(@(r) str2double(extract(...
                extractAfter(r,characterListPattern('. ')+wildcardPattern),digitsPattern)),{hFileRaw.Entity(rawSigInds).Label})));
            if(any(sortedFile))
                if(any(contains(names(sortedFile),'.plx')))
                    NEVFile = 0;
                    [~,mostRecent]=sort([currDir(cellfun(@(s) contains(s,'.plx'),...
                        {currDir(sortedFile).name})).datenum]);
                    loadFile = currDir(sortedFile(mostRecent(1)));
                    [~,channels] = plx_chan_names([loadFile.folder, '\', loadFile.name]);
                else
                    NEVFile = 1;
                    [~,mostRecent] = sort([currDir(sortedFile).datenum],'descend');
                    loadFile = currDir(sortedFile(mostRecent(1)));
                    spkInds = find(cellfun(@(s) strcmp(string(s), "Segment"),...
                        {hFileRaw.Entity(:).EntityType}, 'UniformOutput', true));
                end
                if(~isempty(dir([currName, subDirName,'\*.mat'])))
                    channelMats = dir([currName,subDirName,'\*.mat']);
                    [channelMats,~] = natsort({channelMats.name});
                    unlabeledChannels =cellfun(@(c) str2double(regexp(c,'(?<=\_)\d*(?=\.)',...
                        'match')), channelMats);
                    tic
                    if(isempty(rawSigInds))
                        rawSigInds = find(cellfun(@(s) strcmp(string(s), "Analog"),...
                            {hFileRaw.Entity(:).EntityType}, 'UniformOutput', true) & ...
                            ismember([hFileRaw.Entity(:).ElectrodeID],unlabeledChannels));
                        for c = 1:length(unlabeledChannels)
                            [~,~,data] = ns_GetAnalogData(hFileRaw, rawSigInds(c),1,Inf);
                            avgVs(c) = mean(data,'omitnan');
                        end
                    else
                        avgVs = loadAvgVoltage(hFileRaw,rawSigInds(unlabeledChannels));
                    end
                    toc;
                    meanBG = 5.*abs(avgVs);
                    dirTotal = length(channelMats);
                    if(dirTotal~=32)
                        disp([char(dateArray(d)), ' less than 32 ch']);
                    end
                    [~,channelMap] = natsort({hFileRaw.Entity(rawSigInds).Label});
                    %cellfun(@(r) str2double(r(end-1:end)),{hFileRaw.Entity(cellfun(@(hs) strcmp(hs,'Analog'),{hFileRaw.Entity.EntityType}) & [hFileRaw.Entity.FileType]==3).Label});
                    if(plotWaveforms)
                        figure();
                        hold on;
                    end
                    if(NEVFile)
                        [~,nevFile] = ns_OpenFile([loadFile.folder, '\', loadFile.name]);
                    end
                    %%
                    for e = 1:dirTotal
                        sChInd = channelMap(e);
                        savedStruct = matfile([currName, subDirName,'\', channelMats{e}], 'Writable', true);
                        savedFields = fieldnames(savedStruct);
                        if(~any(contains(savedFields, 'label')) | writeLabels | plotWaveforms)
                            if(NEVFile)
                                [dataTime, ids] = loadSpikeData(nevFile, sChInd);
                                sortedIDs = 1:nevFile.Entity(sChInd).nUnits;
                                %sortedIDs = 1:sum(unique(ids) > 0 & unique(ids)<255);
                            else
                                info = plx_info([loadFile.folder, '\', loadFile.name],1);
                                sortedIDs = 1:sum(info(2:end,sChInd+1)>0);
                            end
                            label = string([]);
                            for a = 1:length(sortedIDs)
                                if(NEVFile)
                                    unitIndx = find(ids==sortedIDs(a));
                                    unitIndx = datasample(unitIndx,min(350,length(unitIndx)));
                                    currUnitTime = dataTime(ids==sortedIDs(a));
                                    data = [];
                                    for s = 1:length(unitIndx)
                                        [~,~,data(:,s),~,~] = ns_GetSegmentData(hFileRaw, ...
                                           find([hFileRaw.Entity(spkInds).ElectrodeID]== ...
                                           nevFile.Entity(sChInd).ElectrodeID,1), unitIndx(s));
                                    end
                                else
                                    [~,~,currUnitTime,data] = plx_waves_v(...
                                        [loadFile.folder, '\', loadFile.name],sChInd, sortedIDs(a));
                                    data = 1000*data';
                                end
                                meanWF = mean(data,2,'omitnan');
                                isi = 1000*(currUnitTime(2:end) - currUnitTime(1:end-1));
                                if(abs(max(meanWF)-min(meanWF)) > meanBG(sChInd) &...
                                        100*(sum(isi<=1.8)/length(isi)) <= 0.5)
                                    label(a) = 's';
                                else
                                    label(a) = 'm';
                                end
                                isi =[];
                                currUnitTime = [];
                                if(plotWaveforms)
                                    subplot(7,6,e);
                                    hold on;
                                    colors = [1 .5 .8; 1 1 0; 1 0 0;];
                                    shadedErrorBar(1:52, mean(data,2)', std(data,0,2)', 'lineprops', {'color',colors(a,:)},'patchSaturation',0.3)
                                    plot(1:52,zeros(52,1),'Color', 'k');
                                    xlim([0 52]);
                                    yL = ylim;
                                    if(abs(yL(end) - yL(1))>50)
                                        yticks(yL(1):50:yL(end));
                                        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)),yL(1):50:yL(end)));
                                    else
                                        yticks(yL(1):25:yL(end));
                                        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)), yL(1):25:yL(end)));
                                    end
                                    title(['CH:', num2str(e), ',E:',num2str(sChInd),', Unit(s):', num2str(a)]);
                                    grid off
                                    ylim(yL)
                                    xticks(30.6)
                                    xticklabels(1);
                                end
                            end
                            numUnitsTrials = size(getfield(savedStruct.sortedSpikeData, 'SpikeTimes'));
                            if(length(label)<numUnitsTrials(1))
                                oldUnits = ~ismember(1:numUnitsTrials(1), sortedIDs);
                                savedStruct = load(savedStruct.Properties.Source);
                                sortedSpikeData = savedStruct.sortedSpikeData;
                                sortedSpikeData.SegTimes(oldUnits,:) = [];
                                sortedSpikeData.SpikeTimes(oldUnits,:) = [];
                                savedStruct.sortedSpikeData = sortedSpikeData;
                            end
                            if(writeLabels)
                                if(isfield(savedStruct.sortedSpikeData,'label'))
                                    savedStruct.label = [];
                                end
                                if(isfield(savedStruct.sortedSpikeData,'labels'))
                                    savedStruct.labels = [];
                                end
                                parSave([currName, subDirName,'\', channelMats{e}],...
                                    savedStruct,label);
                            end
                        else
                            disp(['Labels already created, make overwriteLabels true for: ',...
                                char(dateArray(d))]);
                            break;
                        end
                        if(e==length(channelMats))
                            disp(['Processed and labeled: ', char(dateArray(d))]);
                            if(plotWaveforms)
                                saveFigures(gcf,[currName,'\'],"Waveforms",[]);
                                close all;
                            end
                        end
                    end
                else
                    disp(['No sorted file for: ', char(dateArray(d))]);
                end
                if(NEVFile)
                    ns_CloseFile(hFileRaw);
                    ns_CloseFile(nevFile);
                end
            end
            fclose('all');
        else
            disp(['No raw signal for: ', char(dateArray(d))]);
        end
    end
end
end
function parSave(fileName, savedStruct, lab)
label = lab;
sortedSpikeData = savedStruct.sortedSpikeData;
save(fileName,'sortedSpikeData','label','-append');
end