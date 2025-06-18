clear all
sessionDate = '';

events = repmat({{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}}, 2,1);
events(end+1) = {{'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}};
events(end+1) = {{'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}};
events = events';
conds = {'Extra Small Sphere','Large Sphere','Photocell', 'Rest'};
eventChannel = 'SMA 1';
nFiles = dir(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate, '\Physiology\']);
nevFilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.'):end), '.nev'), {nFiles.name});
nevFiles = nFiles(nevFilesInd);
[~, nevInd] = find(contains({nevFiles.name}, 'sort', 'IgnoreCase',true));
filePath = ['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate, '\Physiology\',nevFiles(nevInd).name];
[ns_status, hFile] = ns_OpenFile(filePath);
eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({eventChannel}, size({hFile.Entity.Reason}))));

%% digital event timestamps
GET_EVENTS = false;
eventTimeStamps = [];
if exist('eventChannel','var')
    if ~isempty(eventChannel)
        GET_EVENTS = true;
    end
end

if GET_EVENTS
    % Retrieve the ID number for the chosen event channel
    eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({eventChannel}, size({hFile.Entity.Reason}))));
    
    % get event entity information
    [ns_RESULT, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));
    
    % event data
    numCount = eventEntityInfo.ItemCount;
    eventTimes = NaN(1, numCount);
    eventTimeStamps = NaN(1, numCount);
    dataSize = NaN(1, numCount);
    
    % get actual event time stamps
    dispstat('','init');
    for i = 1:numCount
        [~, eventTimeStamps(i), eventTimes(i), dataSize(i)] = ns_GetEventData(hFile, eventEntityID, i);
        dispstat(['Loading Event Data... ',num2str(round(i/numCount*100)),'% Done.']);
    end
    dispstat('Loading Event Data... 100% Done.');
    
    % Only consider time stamps for rising edge (eventData = 0 on falling edge)
    tempValidTimeStampInds = find(eventTimes == 32767);
    eventTimes_risingEdge = eventTimeStamps(tempValidTimeStampInds);
    
    tempValidTimeStampInds = find(eventTimes == 0);
    eventTimes_fallingEdge = eventTimeStamps(tempValidTimeStampInds);
end
ns_status = ns_CloseFile(hFile);
% Find start of trial
pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;

% Find start of correct trials
startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);

for i = 1:min(length(startEventIdx),length(endEventIdx))-1
    if(startEventIdx(i+1)<endEventIdx(i))
        startEventIdx(i+1) = [];
        break;
    end
end
for i = 1:min(length(startEventIdx),length(endEventIdx))-1
    if(endEventIdx(i)<startEventIdx(i))
        endEventIdx(i) = [];
    end
end
if(length(startEventIdx)~=length(endEventIdx))
    disp('Mistmatch lengths');
end
fclose('all');

for f = 1:32
    if(any(ismember([hFile.Entity.ElectrodeID],f)))
        segmentTimes = cell(1,length(startEventIdx));
        
        for n = 1:length(startEventIdx) % step through each event time
            % find index of event n
            [~, eventTimeInd] = min(abs(LFPTIMESTAMPS-eventTimes_risingEdge(startEventIdx(n))));
            [~, eventTimeEndInd] = min(abs(LFPTIMESTAMPS-eventTimes_risingEdge(endEventIdx(n))));
            % put LFP data into trials
           eventData{n} = data(eventTimeInd:eventTimeEndInd);
           % times of segments in trial
           [~, segmentTimes{n}] = min(abs(LFPTIMESTAMPS-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
           segmentTimes{n} = segmentTimes{n}-eventTimeInd;
        end
    end
end