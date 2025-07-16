function [analogInputData, analogInputDataTime_s, Fs, eventTimeStamps_risingEdge, eventTimeStamps_fallingEdge] = loadRippleData(fName,channelNum,eventChannel)
[~,hFile] = ns_OpenFile(fName);
dataChannelID = find(channelNum == [hFile.Entity.ElectrodeID]);
%% get analog info contains things like range and sampling rate
[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, dataChannelID);

%% extract data and data time
TimeStamps = hFile.FileInfo(hFile.Entity(dataChannelID(1)).FileType).TimeStamps;
numSamples = sum(TimeStamps(:,end));

analogInputData = zeros(length(dataChannelID),numSamples);
startIndex = 1;
% handling data pauses during session
indexCount = TimeStamps(2,1);
for i = 1:size(TimeStamps,2)
    if(length(dataChannelID)>1)
        [~, tempData] = ns_GetAnalogDataBlock(hFile, dataChannelID, startIndex, indexCount);
    else
        [~, ~, tempData] = ns_GetAnalogData(hFile, dataChannelID, startIndex, indexCount);
    end
    dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));

    % data matrix
    analogInputData(:,dataRange) = tempData';

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
analogInputData(analogInputData > 3E38)=NaN;

% make sure data and data time matrix are an even number in length
if rem(length(analogInputData),2)==1
    analogInputData = analogInputData(1:end-1);
    analogInputDataTime_s = analogInputDataTime_s(1:end-1);
end

%% filtering
Fs = analogInfo.SampleRate;
% Notch filter at 60 Hz
% Wo = 60/(Fs/2);  BW = Wo/35;
% [b,a] = iirnotch(Wo, BW);
% analogInputData = filtfilt(b, a, analogInputData);   %60Hz filtering of signal

%% digital event timestamps
eventTimeStamps_risingEdge = [];
eventTimeStamps_fallingEdge = [];
if exist('eventChannel','var')
    if ~isempty(eventChannel)
        % Retrieve the ID number for the chosen event channel
        eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({eventChannel}, size({hFile.Entity.Reason}))));

        % get event entity information
        [ns_RESULT, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));

        % event data
        numCount = eventEntityInfo.ItemCount;
        eventData = NaN(1, numCount);
        eventTimeStamps = NaN(1, numCount);

        % get actual event time stamps
        dispstat('','init');
        for i = 1:numCount
            [~, eventTimeStamps(i), eventData(i),~] = ns_GetEventData(hFile, eventEntityID, i);
            dispstat(['Loading Event Data... ',num2str(round(i/numCount*100)),'% Done.']);
        end
        dispstat('Loading Event Data... 100% Done.');

        % Only consider time stamps for rising edge (eventData = 0 on falling edge)
        % This may not work with parallel input
        eventTimeStamps_risingEdge = eventTimeStamps(find(eventData == 32767));

        eventTimeStamps_fallingEdge = eventTimeStamps(find(eventData == 0));
    else
        eventTimeStamps_risingEdge = [];
        eventTimeStamps_fallingEdge = [];
    end
end
end