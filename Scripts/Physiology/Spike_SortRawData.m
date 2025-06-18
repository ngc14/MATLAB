function Spike_SortRawData(date, monkeyName)
sessionDate = '11_19_2019';
monkey = 'Gilligan';
if(exist('date', 'var'))
    sessionDate = date;
end
if(exist('monkeyName', 'var'))
    monkey =monkeyName;
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

if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'MM_dd_yyyy';
else
    dateFormat = 'yyyy_MM_dd';
end
sessionDate = char(datetime(sessionDate,"InputFormat",dateFormat,"Format",dateFormat));
%% File info
fileDir = ['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate, '\Physiology\'];
if(~exist(fileDir,'dir'))
    fprintf('%s not exist', datetime(date,dateFormat));
    return;
end
nFiles = dir(fileDir);
nevFilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.'):end), '.nev'), {nFiles.name});
nevFiles = nFiles(nevFilesInd);
[~,mostRecentOrder] = sort(cellfun(@datenum,{nevFiles.date}),'descend');
[~, nevInd] = find(contains({nevFiles(mostRecentOrder).name}, 'sort', 'IgnoreCase',true),1);
filePath = ['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate, '\Physiology\',nevFiles(mostRecentOrder(nevInd)).name];
if(isempty(nevInd))
    disp('Not sorted');
    disp(sessionDate);
    return;
end
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
% corrections to specific recording sessions that output 
% spurious errors of TTL pulse encoding of event times
switch(char(datetime(sessionDate,'InputFormat',dateFormat,'Format',dateFormat)))
    case('2021_09_23')
        eventTimes_fallingEdge(629) = eventTimes_risingEdge(629)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_09_24')
        eventTimes_fallingEdge(713) = eventTimes_risingEdge(713)+0.02;
        eventTimes_fallingEdge(2425) = eventTimes_risingEdge(2425)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_09_27')
        eventTimes_fallingEdge(1121) = eventTimes_risingEdge(1121)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_09_28')
        eventTimes_fallingEdge(1613) = eventTimes_risingEdge(1613)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_09_29')
        eventTimes_fallingEdge(382) = eventTimes_risingEdge(382)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_09_30')
        eventTimes_fallingEdge(1405) = eventTimes_risingEdge(1405)+0.02;
        eventTimes_fallingEdge(1731) = eventTimes_risingEdge(1731)+0.02;
        eventTimes_fallingEdge(1754) = eventTimes_risingEdge(1754)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_10_01')
        eventTimes_fallingEdge(660) = eventTimes_risingEdge(660)+0.02;
        eventTimes_fallingEdge(701) = eventTimes_risingEdge(701)+0.02;
        eventTimes_fallingEdge(1124) = eventTimes_risingEdge(1124)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_10_04')
        eventTimes_fallingEdge(928) = eventTimes_risingEdge(928)+0.02;
        eventTimes_fallingEdge(1148) = eventTimes_risingEdge(1148)+0.02;
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    case('2021_10_13')
        pulseLength(1172) = 0.02;
    case('2021_11_09')
        pulseLength(1310) = 0.02;
    case('2021_11_15')
        pulseLength(1612) = 0.02;
    case('2021_11_17')
        pulseLength(839) = 0.02;
    case('2021_11_19')
        pulseLength(1) = 0.02;
    case('2021_12_16')
        pulseLength(2824) = 0.02;
        pulseLength(2866) = 0.02;
        pulseLength(3074) = 0.02;
    case('2022_06_20')
        pulseLength();
    case('2022_07_03')
        pulseLength(1) = 0;
        pulseLength(2) = 0.02;
    case('2022_07_04')
        pulseLength(1510:1512) = 0;
    case('2022_12_17')
        pulseLength(746)  = 0;
    case('2019_04_29')
        pulseLength(365)=0.01;
         pulseLength(1425) = 0.01;
    case('06_10_2019')
        pulseLength(2329)= 0.01;
    case('02_07_2020')
        pulseLength(2184)=0.01;
end
startEventIdx=find(pulseLength>0.0175 & pulseLength<0.03);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);

passed = false;
while(~passed)
    for i = 1:min(length(startEventIdx),length(endEventIdx))-1
        if(startEventIdx(i)>endEventIdx(i+1))
            startEventIdx(i) = [];
            removed = true;
            break;
        elseif(startEventIdx(i+1)<endEventIdx(i))
            startEventIdx(i+1) = [];
            removed = true;
            break;
        elseif(endEventIdx(i)<startEventIdx(i))
            endEventIdx(i) = [];
            removed = true;
            break;

        else
            removed = false;
        end
    end
    passed = ~removed;
end
switch(datestr(sessionDate,'yyyy_mm_dd'))
    case ('2021_09_22')
        startEventIdx(end) = [];
end
if(length(startEventIdx)==length(endEventIdx)+1)
    startEventIdx=startEventIdx(1:end-1);
end
if(length(startEventIdx)~=length(endEventIdx))
    disp('Mistmatch lengths');
    disp(sessionDate)
    return;
end
spikeChannels = cellfun(@(a) strcmp(a, 'Segment'), {hFile.Entity.EntityType});
%% Process session
for f = 1:sum(spikeChannels)
    fullName = [monkey,'_', sessionDate,'_', num2str(f)];
    [dataTime, ids] = loadSpikeData(hFile,f);
    if(~isempty(dataTime))
        sortedIDs = unique(ids);
        sortedIDs = sortedIDs(sortedIDs>0 & sortedIDs<255);
        spikeTimes = cell(length(sortedIDs),length(startEventIdx));
        segmentTimes = cell(length(sortedIDs),length(startEventIdx));
        for u = 1:length(sortedIDs)
            dataUnits = dataTime(ids==sortedIDs(u));
            for n = 1:length(startEventIdx) % step through each event time
                % find index of event n
                [closestTimeStart, eventTimeInd] = min(abs(dataUnits-eventTimes_risingEdge(startEventIdx(n)+1)));
                [closestTimeEnd, eventTimeEndInd] = min(abs(dataUnits-eventTimes_risingEdge(endEventIdx(n)-1)));
                if(min(closestTimeStart,closestTimeEnd)>eventTimes_risingEdge(endEventIdx(n)-1)-eventTimes_risingEdge(startEventIdx(n)+1))
                    spikeTimes{u,n} = NaN;
                    segmentTimes{u,n} = NaN;
                else
                    spikeTimes{u,n} = dataUnits(eventTimeInd:eventTimeEndInd);
                    segmentTimes{u,n} = eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1);
                end
            end
        end
        %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
        if(~exist(['S:\Lab\', monkey, '\All Data\', monkey,'_', ...
                sessionDate,'\Arduino\', monkey,'_',sessionDate, '.txt'],'file'))
            if(~exist(['S:\Lab\', monkey, '\All Data\', monkey, '_', sessionDate, '\Arduino\'], 'dir'))
                mkdir(['S:\Lab\', monkey, '\All Data\', monkey, '_', sessionDate, '\Arduino\']);
            end
            if(strcmp(monkey, 'Gilligan'))
                copyfile(['S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',...
                    sessionDate,'\Gilligan_',sessionDate,'.txt'],...
                    ['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Arduino\Gilligan_',sessionDate, '.txt']);
            else
                copyfile(['S:\Lab\BehaveBackup\Monkey_Training\Skipper_Macaque_152_17\Skipper_',...
                    datestr(sessionDate, 'mm_dd_yyyy'),'\Skipper_',datestr(sessionDate, 'mm_dd_yyyy'),'.txt'],...
                    ['S:\Lab\Skipper\All Data\Skipper_', sessionDate, '\Arduino\Skipper_',sessionDate, '.txt']);
            end
        end
        arduinoFilename = ['S:\Lab\', monkey, '\All Data\', monkey,'_',...
            sessionDate,'\Arduino\', monkey,'_',sessionDate, '.txt'];
        % Open the text file.
        fileID = fopen(arduinoFilename,'r');
        % Read columns of data according to the format.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
        % Close the text file.
        fclose(fileID);
        fclose('all');
        
        % Create output variable
        recordedTrial = [dataArray{1:end-1}];
        % remove lines that are not task related (pause, headers)
        recordedTrial = recordedTrial(ismember(recordedTrial(:,1), conds),:);
        if ~isempty(relevantTrials)
            recordedTrial=recordedTrial(relevantTrials(1):relevantTrials(end),:);
        end
        % Remove failed trials.
        if(length(recordedTrial)~=max(size((segmentTimes))))
            disp('Arduino file mistmatch trials');
            return;
        end
        falseTrialIdx=cellfun(@(s) strcmp(s,"False Start") | strcmp(s,"Failed-to-Reach"),recordedTrial(:,8),'UniformOutput',false);
        correctTrialIdx = find(~cell2mat(falseTrialIdx));
        % falseTrialIdx=cellfun(@num2str,recordedTrial(:,8),'UniformOutput',false);
        % correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
        
        %[row,~]=find(cellfun(@(x) removeShortSegs(x), recordedTrial(:,2:6)));
        
        
        successfulTrial=recordedTrial(correctTrialIdx,[1:7,9]);
        sortedSpikeData.SpikeTimes= spikeTimes(:,correctTrialIdx);
        sortedSpikeData.SegTimes= segmentTimes(:,correctTrialIdx);
        sortedSpikeData.Date= sessionDate;
        sortedSpikeData.DataChannel= f;
        sortedSpikeData.EventChannel= eventChannel;
        sortedSpikeData.ArduinoData= successfulTrial;
        sortedSpikeData.Conditions= conds;
        sortedSpikeData.ConditionSegments= events;
        sortedSpikeData.Locations= unique(successfulTrial(:,end));
        
        if(~exist(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
                '\Physiology\Results_New\'], 'dir'))
            mkdir(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
                '\Physiology\Results_New\']);
        end
        save(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
            '\Physiology\Results_New\',fullName],'sortedSpikeData');
        clear sortedSpikeData spikeTimes segmentTimes
    end
end
fclose('all');
end