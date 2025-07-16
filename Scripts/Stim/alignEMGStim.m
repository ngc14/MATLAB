clc
clear all
close all

monkey = 'Gilligan';
date = '02_25_2021';
run = '0001';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chanMap = [1:2:32, 2:2:32];

sigma = 5; % smoothing window in ms

file = ['S:\Lab\Gilligan\All Data\Gilligan_', ...
    date, '\EMG\Gilligan_', date, '_EMG.txt'];
fileID = fopen(file);
fileInd = 1;
muscles = {};
while(fopen(fileID))
    line = fgetl(fileID);
    if(contains(line, 'Channel B','IgnoreCase',true))
        channelNumStart = regexp(line, 'B');
        channelNumEnd = regexp(line, ':');
        muscles{end+1} = line(channelNumEnd+2:end);
        channelNum(fileInd) = 128 + str2double(line(channelNumStart+1:channelNumEnd-1));
        fileInd = fileInd + 1;
    end
end
% fclose(fileID);
channelNum = unique(channelNum);
time_before_event = -.1;
time_after_event = .1;

for m = 1:length(muscles)
    filePath = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date,'\EMG\Gilligan_',date,'_stim_', run,'.nf3'];
    %% get basic file info
    [ns_status, hFile] = ns_OpenFile(filePath);
    
    %% Determine correct entityID for desired datastream
    
    % find index of specified electrode
    EntityIndices = find([hFile.Entity(:).ElectrodeID] == channelNum(m));
    fileTypeNums = [hFile.Entity(EntityIndices).FileType];
    fileTypes = {hFile.FileInfo(fileTypeNums).Type};
    entityID = EntityIndices(find(contains(fileTypes,filePath(end-2:end))));
    %% get analog info contains things like range and sampling rate
    [ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID);
    
    %% extract data and data time
    TimeStamps = hFile.FileInfo(hFile.Entity(entityID).FileType).TimeStamps;
    numSamples = sum(TimeStamps(:,end));
    analogInputData = zeros(1,numSamples);
    startIndex = 1;
    indexCount = TimeStamps(2,1);
    for i = 1:size(TimeStamps,2)
        [a, b, tempData] = ns_GetAnalogData(hFile, entityID, startIndex, indexCount);
        dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
        
        % data matrix
        analogInputData(dataRange) = tempData';
        
        clear tempData
        if i ~= size(TimeStamps,2)
            startIndex = startIndex + TimeStamps(2,i);
            indexCount = TimeStamps(2,i+1);
        end
    end
    
    % data time matrix
    analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
    
    
    %% filter out noise spikes (abnormally large data values)
    
    % detect noise spike and replace it with NaN
    noiseSpike = find(analogInputData > 3E38);
    analogInputData(noiseSpike)=NaN;
    
    % make sure data and data time matrix are an even number in length
    if rem(length(analogInputData),2)==1
        analogInputData = analogInputData(1:end-1);
        analogInputDataTime_s = analogInputDataTime_s(1:end-1);
    end
    EMGSig{m} = abs(analogInputData);
    EMGTimes{m} = analogInputDataTime_s;
    ns_CloseFile(hFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['\\univ.pitt.edu\sni\Gharbawie\Lab\' monkey, '\All Data\', monkey,'_', date,'\Physiology\', monkey, '_', date, '_stim_',run];
[n,names] = plx_chan_names([filename, '-sorted.plx']);
fileID = fopen([filename, '.txt']);
fileInfo = textscan(fileID, '%s%f%f%f%f%f%f%f','Delimiter', '\tn', 'headerlines', 1);
fclose(fileID);
fileInfo{1} = cellfun(@(a)  str2num(string(a)), fileInfo{1}, 'UniformOutput', false);
trialInfo = [cellfun(@(a) a(1), fileInfo{1}), fileInfo{3}, fileInfo{4}, fileInfo{5}, fileInfo{6}, fileInfo{7}, fileInfo{8}];
condCombs = unique(trialInfo, 'rows');
alignedTrial = cell(size(EMGSig));
%%
stimElectrodes = unique(trialInfo(:,1));
stim_ts = cell(size(stimElectrodes));

for e = 1:length(stimElectrodes)
    stimE = stimElectrodes(e);
    name = sprintf('Nano2_stim.e%02.f',stimE);
    stim_chan = find(contains(string(names), name));
    [n, stim_ts_raw] = plx_ts([filename, '-sorted.plx'], stim_chan, 1);
    currTrialsInd = trialInfo(:,1)==stimE;
    pulses = trialInfo(currTrialsInd,2);
    tInd = 1;
    stim_ts{e}(end+1,:) = stim_ts_raw(1);
    for t = 1:length(pulses)-1
        stim_ts{e}(end+1,:) = stim_ts_raw(tInd+pulses(t));
        tInd = tInd + pulses(t);
    end
end
for e = 1:size(condCombs,1)
    cond = condCombs(e,:);
    name = sprintf('Nano2_stim.e%02.f',cond(1));
    stim_chan = find(contains(string(names), name));
    currTrials = trialInfo(trialInfo(:,1)==cond(1),:);
    currTrialInds = currTrials(:,2)==cond(2) &  currTrials(:,3)==cond(3) & ...
        currTrials(:,4)==cond(4) & currTrials(:,5)==cond(5) & currTrials(:,6)==cond(6) & ...
        currTrials(:,7)==cond(7);
    stim_ts_cond = stim_ts{find(cond(1)==stimElectrodes)}(currTrialInds);
    for s = 1:length(stim_ts_cond)
        startInds = cellfun(@(a) getEMGInds(a, stim_ts_cond(s)+time_before_event), EMGTimes, 'UniformOutput', false);
        endInds = cellfun(@(a) getEMGInds(a, stim_ts_cond(s)+time_after_event), EMGTimes, 'UniformOutput', false);
        centeredSig = cellfun(@(a,b,c) a(b:c),EMGSig,startInds, endInds,'UniformOutput', false);
        [row, col] = find(~cellfun(@isempty, centeredSig));
        for i = 1:length(row)
            alignedTrial{row(i), col(i)}{end+1,:} = centeredSig{row(i), col(i)};
        end
    end

    for m = 1:length(alignedTrial)
        s = subplot(size(condCombs,1),length(muscles),length(muscles)*(e-1) + m);
        if(m==1)
            ylabel(cond);
        end
        hold on;
        trials = alignedTrial{m};
        yInd = 0;
        for t = 1:length(trials)
            plot(s, -.1:1/2000:.1,  yInd + trials{t});
        end
        plot(s,-.1:1/2000:.1,nanmean(cell2mat(trials),1), 'k', 'LineWidth', 2);
        title(s,muscles{m});
    end

end
function minInd = getEMGInds(sig, stim)
    [~,minInd] = min(abs(sig-stim));
end