function Spike_SortRawData(date)
sessionDate = '05_22_2019';
if(exist('date', 'var'))
    sessionDate = date;
end
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
monkey = 'Gilligan';
channelMap = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32];
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'mm_dd_yyyy';
else
    dateFormat = 'yyyy_mm_dd';
end

%% File info
nFiles = dir(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\']);
nevFilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.'):end), '.nev'), {nFiles.name});
nevFiles = nFiles(nevFilesInd);
nevInd = 1;
filePath = ['S:\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\',nevFiles(nevInd).name];
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
switch(datestr(sessionDate,dateFormat))
    case('05_14_2019')
        pulseLength(1157) = 0.01;
        pulseLength(1161) = 0.01;
    case('05_15_2019')
        pulseLength(2199) = 0.01;
    case('05_22_2019')
        pulseLength(1319)= 0.01;
        pulseLength(1107) = 0.01;
end
% Find start of correct trials
startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
for i = 1:min(length(startEventIdx),length(endEventIdx))-1
    if(startEventIdx(i+1)<endEventIdx(i))
        startEventIdx(i+1) = [];
    end
end
for i = 1:min(length(startEventIdx),length(endEventIdx))-1
    if(endEventIdx(i)<startEventIdx(i))
        endEventIdx(i) = [];
    end
end
if(length(startEventIdx)~=length(endEventIdx))
    disp('Mistmatch lengths');
    return
end
%% Process session
plxInd = contains({nFiles.name}, '-sorted.plx');
plxName = ['S:\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\',nFiles(plxInd).name];
[~,channels] = plx_chan_names(plxName);
info = plx_info(plxName,1);
channelNums = find(info(1,:)>0)-1;

for f = 1:length(channelNums)
    fullName = ['Gilligan_', sessionDate,'_', num2str(f)];
    sortedIDs = sum(info(2:end,channelNums(f)+1)>0);
    spikeTimes = cell(length(sortedIDs),length(startEventIdx));
    segmentTimes = cell(length(sortedIDs),length(startEventIdx));

    for u = 1:length(sortedIDs)
        [~,dataUnits] = plx_ts(plxName,channels(channelNums(f),:),u);
        for n = 1:length(startEventIdx) % step through each event time
            % find index of event n
            [closestTimeStart, eventTimeInd] = min(abs(dataUnits-eventTimes_risingEdge(startEventIdx(n))));
            [closestTimeEnd, eventTimeEndInd] = min(abs(dataUnits-eventTimes_risingEdge(endEventIdx(n))));
            if(min(closestTimeStart,closestTimeEnd)>eventTimes_risingEdge(endEventIdx(n))-eventTimes_risingEdge(startEventIdx(n)))
                spikeTimes{u,n} = NaN;
                segmentTimes{u,n} = NaN;
            else
                spikeTimes{u,n} = dataUnits(eventTimeInd:eventTimeEndInd);
                segmentTimes{u,n} = eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1);
            end
        end
    end
    %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
    if(~exist(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Arduino\Gilligan_',sessionDate, '.txt'],'file'))
        copyfile(['S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',sessionDate,'\Gilligan_',sessionDate,'.txt'],...
            ['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Arduino\Gilligan_',sessionDate, '.txt']);
    end
    arduinoFilename = ['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Arduino\Gilligan_',sessionDate, '.txt'];
    
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
    if(length(recordedTrial)~=max(size((segmentTimes))))
        disp('Arduino file mistmatch trials');
    end
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
    
    if(~exist(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\'], 'dir'))
        mkdir(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\']);
    end
    save(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\',fullName],'sortedSpikeData');
    
    clear sortedSpikeData spikeTimes segmentTimes
    
end
end
function x = removeShortSegs(x)
if(isnumeric(str2double(x)))
    x = str2double(x)<= 10;
else
    x = 0;
end
end