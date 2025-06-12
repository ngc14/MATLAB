clc
clear all
close all

monkey = 'Gilligan';
date = '06_13_2019';
run = '0004';
eventChannelName = 'SMA 2';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 5; % smoothing window in ms

file = ['S:\Lab\Gilligan\All Data\Gilligan_', ...
    date, '\EMG\Gilligan_', date, '_NotesEMG.txt'];
fileID = fopen(file);
fileInd = 1;
muscles = {};
while(~feof(fileID))
    line = fgetl(fileID);
    if(contains(line, 'Channel B','IgnoreCase',true))
        channelNumStart = regexp(line, 'B');
        channelNumEnd = regexp(line, ':');
        muscles{end+1} = line(channelNumEnd+2:end);
        channelNum(fileInd) = 128 + str2double(line(channelNumStart+1:channelNumEnd-1));
        fileInd = fileInd + 1;
    end
end
fclose(fileID);
channelNum = unique(channelNum);
time_before_event = -.1;
time_after_event = .1;
%%
for m = 1:length(muscles)
    filePath = ['S:\Lab\Gilligan\All Data\Gilligan_',date,'\EMG\Gilligan_',date,'_stim_', run,'.nf3'];
    % get basic file info
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
filename = ['S:\Lab\' monkey, '\All Data\', monkey,'_', date,'\Physiology\', monkey, '_', date, '_stim_',run];
[ns_status, hFile] = ns_OpenFile([filename, '.nev']);
eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({eventChannelName}, size({hFile.Entity.Reason}))));

% get event entity information
[ns_RESULT, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));

% event data
numCount = eventEntityInfo.ItemCount;
eventData = NaN(1, numCount);
eventTimeStamps = NaN(1, numCount);
dataSize = NaN(1, numCount);

% get actual event time stamps
dispstat('','init');
for i = 1:numCount
    [~, eventTimeStamps(i), eventData(i), dataSize(i)] = ns_GetEventData(hFile, eventEntityID, i);
    dispstat(['Loading Event Data... ',num2str(round(i/numCount*100)),'% Done.']);
end
dispstat('Loading Event Data... 100% Done.');

% Only consider time stamps for rising edge (eventData = 0 on falling edge)
% This may not work with parallel input
eventTimeStamps_fallingEdge = eventTimeStamps(eventData == 0);
%% NOTE: THIS LINE IS CONTINGENT ON THE TRAINS BEING LONGER THAN .3 SECONDS
%%APART
trainStart = find(diff(eventTimeStamps_fallingEdge)>.3);
alignedTrial = cell(size(EMGSig));
%%

stim_ts_cond = eventTimeStamps_fallingEdge(trainStart);
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
    s = subplot(1,length(muscles), m);
    hold on;
    trials = alignedTrial{m};
    yInd = 0;
    for t = 1:length(trials)
        plot(s, -.1:1/2000:.1,  yInd + trials{t});
    end
    plot(s,-.1:1/2000:.1,nanmean(cell2mat(trials),1), 'k', 'LineWidth', 2);
    title(s,muscles{m});
end

function minInd = getEMGInds(sig, stim)
[~,minInd] = min(abs(sig-stim));
end