sessionDate = '04_19_2019';
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

events = repmat({{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}}, 2,1);
events(end+1) = {{'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}};
events(end+1) = {{'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}};
events = events';
conds = {'Extra Small Sphere','Large Sphere','Photocell', 'Rest'};
eventChannel = 'SMA 1';
relevantTrials = [];
channelMap = [1,17,2,18,3,19,4,20,5,21,6,22,7,23,8,24,9,25,10,26,11,27,12,28,...
    13,29,14,30,15,31,16,32];


%% File info
nFiles = dir(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\']);
nevFilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.'):end), '.nev'), {nFiles.name});
nevFiles = nFiles(nevFilesInd);
[~, nevInd] = find(contains({nevFiles.name}, 'sorted'));

filePath = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\',nevFiles(nevInd).name];
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
    % This may not work with parallel input
    tempValidTimeStampInds = find(eventTimes == 32767);
    eventTimes_risingEdge = eventTimeStamps(tempValidTimeStampInds);
    
    tempValidTimeStampInds = find(eventTimes == 0);
    eventTimes_fallingEdge = eventTimeStamps(tempValidTimeStampInds);
end
ns_status = ns_CloseFile(hFile);
% Find start of trial
pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;

pulseLength(2173) = 0.01
pulseLength(2602) = 0.01

% Find start of correct trials
startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
if(length(startEventIdx)~=length(endEventIdx))
    disp('Mistmatch lengths');
    return
end
%%
spikeTimestamps = double(readNPY(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\spike_times.npy']))./30000;
spikeClusters = readNPY(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\spike_clusters.npy']);
clusterIDs = unique(spikeClusters);
clusterInfo = tdfread(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\cluster_info.tsv']);
goodClusters = arrayfun(@(a) contains(a, 'good'), string(clusterInfo.group));
clusterIDs = clusterIDs(goodClusters);
%% Process session
for f = 1:length(clusterIDs)
    fullName = ['Gilligan_', sessionDate,'_', num2str(f)];
    spikeTimes = cell(1,length(startEventIdx));
    segmentTimes = cell(1,length(startEventIdx));
    
    dataUnits = spikeTimestamps(spikeClusters==clusterIDs(f));
    for n = 1:length(startEventIdx) % step through each event time
        % find index of event n
        [closestTimeStart, eventTimeInd] = min(abs(dataUnits-eventTimes_risingEdge(startEventIdx(n))));
        [closestTimeEnd, eventTimeEndInd] = min(abs(dataUnits-eventTimes_risingEdge(endEventIdx(n))));
        if(min(closestTimeStart,closestTimeEnd)>eventTimes_risingEdge(endEventIdx(n))-eventTimes_risingEdge(startEventIdx(n)))
            spikeTimes{1,n} = NaN;
            segmentTimes{1,n} = NaN;
        else
            spikeTimes{1,n} = dataUnits(eventTimeInd:eventTimeEndInd);
            segmentTimes{1,n} = eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1);
        end
    end
    %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
    arduinoFilename = ['\\pitt\sni\gharbawie\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',sessionDate,'\Gilligan_',sessionDate,'.txt'];
    
    % Open the text file.
    fileID = fopen(arduinoFilename,'r');
    
    % Read columns of data according to the format.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
    
    % Close the text file.
    fclose(fileID);
    
    % Create output variable
    recordedTrial = [dataArray{1:end-1}];
    
    % remove lines that are not task related (pause, headers)
    recordedTrial = recordedTrial(ismember(recordedTrial(:,1), conds),:);
    if ~isempty(relevantTrials)
        recordedTrial=recordedTrial(relevantTrials(1):relevantTrials(end),:);
    end
    
    %[row,~]=find(cellfun(@(x) removeShortSegs(x), recordedTrial(:,2:6)));
    
    % Remove failed trials.
    falseTrialIdx=cellfun(@str2num,recordedTrial(:,8),'UniformOutput',false);
    correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
    %correctTrialIdx = setdiff(correctTrialIdx, row);
    
    successfulTrial=recordedTrial(correctTrialIdx,[1:7,9]);
    spikeTimes = spikeTimes(:,correctTrialIdx);
    segmentTimes = segmentTimes(:,correctTrialIdx);
    
    sortedSpikeData.SpikeTimes= spikeTimes;
    sortedSpikeData.SegTimes = segmentTimes;
    
    sortedSpikeData.Date=sessionDate;
    sortedSpikeData.DataChannel=f;
    sortedSpikeData.EventChannel=eventChannel;
    
    sortedSpikeData.ArduinoData = successfulTrial;
    sortedSpikeData.Conditions = conds;
    sortedSpikeData.ConditionSegments = events;
    sortedSpikeData.Locations = unique(successfulTrial(:,end));
    
    if(~exist(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results2\'], 'dir'))
        mkdir(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results2\']);
    end
    save(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results2\',fullName],'sortedSpikeData');
    clear sortedSpikeData spikeTimes segmentTimes
    
end
function x = removeShortSegs(x)
if(isnumeric(str2double(x)))
    x = str2double(x)<= 10;
else
    x = 0;
end
end