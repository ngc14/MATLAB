function Spike_SortRawData(date, monkeyName)
sessionDate = '03_17_2020';
monkey = 'Gilligan';
if(exist('date', 'var'))
    sessionDate = date;
end
if(exist('monkeyName', 'var'))
    monkey =monkeyName;
end
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

events = repmat({{'StartTrial','GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward','EndTrial'}}, 2,1);
events(end+1) = {{'StartTrial','GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward','EndTrial'}};
events(end+1) = {{'StartTrial','GoSignal','StartHold', 'StartReplaceHold','StartReplaceSuccess','StartReward','EndTrial'}};
events = events';
conds = {'Extra Small Sphere','Large Sphere','Photocell', 'Rest'};
eventChannel = 'SMA 1';
relevantTrials = [];

if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'MM_dd_uuuu';
else
    dateFormat = 'uuuu_MM_dd';
end
sessionDate = char(datetime(sessionDate,"InputFormat",dateFormat,"Format",dateFormat));
%% File info
fileDir = ['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate, '\Physiology\'];
if(~exist(fileDir,'dir'))
    fprintf('%s not exist', datetime(date,dateFormat));
    return;
end
nFiles = dir(fileDir);
nevFilesInd = cellfun(@(a) contains(a, '.nev'), {nFiles.name});
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
        pulseLength(1415)=0.01;
        pulseLength(2184)=0.01;
    case('03_17_2020')
        pulseLength(1664:end) = 0.01;

end
startEventIdx=find(pulseLength>0.0175 & pulseLength<0.03);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);

passed=false;
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
                [closestTimeStart, eventTimeInd] = min(abs(dataUnits-eventTimes_risingEdge(startEventIdx(n))));
                [closestTimeEnd, eventTimeEndInd] = min(abs(dataUnits-eventTimes_risingEdge(endEventIdx(n))));
                %if(min(closestTimeStart,closestTimeEnd)>eventTimes_risingEdge(endEventIdx(n)-1)-eventTimes_risingEdge(startEventIdx(n)+1))
                %    spikeTimes{u,n} = NaN;
                %    segmentTimes{u,n} = NaN;
                %else
                    spikeTimes{u,n} = dataUnits(eventTimeInd:eventTimeEndInd);
                    segmentTimes{u,n} = eventTimes_risingEdge(startEventIdx(n):endEventIdx(n));
                %end
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
            disp('Arduino file mistmatch trials, using digital lines');
            for n = 1:length(startEventIdx)
                segmentTimes{n} = eventTimes_risingEdge(startEventIdx(n):endEventIdx(n));
            end
            trialLength = endEventIdx-startEventIdx;
            errorKeys = [3 4 5 6 7 8 9 10];
            errorVals = {'False Start', 'Failed-to-Reach','Failed-to-Contact','Failed-to-Lift','Failed-to-Lift'...
                'Failed-to-Hold','Failed-to-Replace','Failed-to-ReplaceHold'};
            [~,mostSegments] = max(cellfun(@length, events));
            mostSegments = events{mostSegments}(1:end-1);
            errorCondSegs = cellfun(@(e) find(~contains(mostSegments,e)), events, 'UniformOutput',false);
            errorDiffs = cellfun(@(a) -1.*ismember(1:length(mostSegments),a),errorCondSegs,'UniformOutput',false);
            eKeys = cellfun(@(e) cumsum(ones(size(mostSegments)) + e), errorDiffs,'UniformOutput',false);
            condKeys = containers.Map(conds,cellfun(@(k) unique(k(ismember(k,errorKeys))), eKeys, 'UniformOutput',false));
            condVals = containers.Map(conds,cellfun(@(f) string(errorVals(~ismember(1:length(errorKeys),find(f~=0)))), errorDiffs,'UniformOutput',false));
            
            correctTrialIdx = pulseLength(endEventIdx)<0.05;
            [~, hFile] = ns_OpenFile(filePath);
            eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({'SMA 2'}, size({hFile.Entity.Reason}))));
            [~, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));
            numCount = eventEntityInfo.ItemCount;
            pulseVal = NaN(1, numCount);
            SMA2Times = NaN(1, numCount);
            for i = 1:numCount
                [~, SMA2Times(i), pulseVal(i), ~] = ns_GetEventData(hFile, eventEntityID, i);
            end
            SMA2Times = SMA2Times(pulseVal == 32767);
            eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({'SMA 3'}, size({hFile.Entity.Reason}))));
            [~, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));
            numCount = eventEntityInfo.ItemCount;
            pulseVal = NaN(1, numCount);
            SMA3Times = NaN(1, numCount);
            for i = 1:numCount
                [~, SMA3Times(i), pulseVal(i), ~] = ns_GetEventData(hFile, eventEntityID, i);
            end
            SMA3Times = SMA3Times(pulseVal == 32767);
            eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({'SMA 4'}, size({hFile.Entity.Reason}))));
            [~, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));
            numCount = eventEntityInfo.ItemCount;
            pulseVal = NaN(1, numCount);
            SMA4Times = NaN(1, numCount);
            for i = 1:numCount
                [~, SMA4Times(i), pulseVal(i), ~] = ns_GetEventData(hFile, eventEntityID, i);
            end
            SMA4Times = SMA4Times(pulseVal == 32767);
            SMA2Times = SMA2Times(correctTrialIdx(trialLength>3 & trialLength~=5 | correctTrialIdx~=1 & trialLength==5));
            [minVal,minInd]=arrayfun(@(a) min([min(abs(a-SMA3Times)),min(abs(a-SMA4Times))]),SMA2Times);
  

            photocellInds = minVal>1;
            essInds = minInd==1 & ~photocellInds;
            succesfulTrials = cumsum(correctTrialIdx);
            arduinoPseudoTable = cell(length(correctTrialIdx),length(errorKeys));
            restInds = find(correctTrialIdx & trialLength==5);
            moveConds = find(correctTrialIdx & trialLength~=5);
            condInds = {moveConds(essInds),moveConds(~essInds & ~photocellInds),moveConds(photocellInds),restInds};
            [~,maxSegConds] = max(cellfun(@length,events));
            maxSegConds = events{maxSegConds};
            for c = 1:length(conds)
                arduinoPseudoTable(condInds{c},1) = conds(c);
                segCondTimes = NaN(length(condInds{c}),length(maxSegConds));
                matchedTrials = cellfun(@(t) [t(1:2), NaN(length(events{c})-length(t),1),t(3:end)], ...
                    segmentTimes(condInds{c}),'UniformOutput',false);
                segCondTimes(:,ismember(maxSegConds,events{c})) = round(1000*...
                    [zeros(length(condInds{c}),1),diff(cell2mat(matchedTrials'),1,2)],0);
                arduinoPseudoTable(condInds{c},2:end-1) = num2cell(segCondTimes(:,3:size(arduinoPseudoTable,2)));
                arduinoPseudoTable(condInds{c},end) = num2cell(succesfulTrials(condInds{c}));
            end
            failedTrials = ~correctTrialIdx;
            startErrorRep = find([isempty(arduinoPseudoTable{1,1}) diff(failedTrials)]==1);
            endErrorRep = find(diff(correctTrialIdx)==1);
            if(length(startErrorRep)>length(endErrorRep))
                startErrorRep = startErrorRep(1:end-1);
                arduinoPseudoTable = arduinoPseudoTable(1:size(arduinoPseudoTable,1)-1,:);
            end
            for r = 1:length(startErrorRep)
                trialRange = startErrorRep(r):endErrorRep(r);
                arduinoPseudoTable(trialRange,1) = arduinoPseudoTable(endErrorRep(r)+1,1);
                for t = 1:length(trialRange)
                    arduinoPseudoTable(trialRange(t),3:length(segmentTimes{trialRange(t)})) = ...
                        num2cell(1000*diff(segmentTimes{trialRange(t)}(2:end)-eventTimes_risingEdge(startEventIdx(trialRange(t)))));
                    currErrors = cell2mat(condVals.values(arduinoPseudoTable(trialRange(t),1)));
                    errorInd = cell2mat(condKeys.values(arduinoPseudoTable(trialRange(t),1)))==length(segmentTimes{trialRange(t)})+1;
                    arduinoPseudoTable(trialRange(t),find(errorInd,1)+1:end) = cellstr(currErrors(errorInd));
                end
            end
            arduinoPseudoTable(cellfun(@(n) all(isnan(n)), arduinoPseudoTable)) = {''};
            arduinoPseudoTable = cellfun(@char, cellfun(@string, arduinoPseudoTable,'UniformOutput',false), 'UniformOutput', false);
            falseTrialIdx=false(size(arduinoPseudoTable,1),1); %%cellfun(@(s) strcmp(s,"False Start") | strcmp(s,"Failed-to-Reach") | isempty(s),arduinoPseudoTable(:,8),'UniformOutput',false);
            correctTrialIdx = find(~falseTrialIdx);
            sortedSpikeData.SpikeTimes= spikeTimes(:,correctTrialIdx);
            sortedSpikeData.SegTimes= segmentTimes(:,correctTrialIdx);
            sortedSpikeData.Date= sessionDate;
            sortedSpikeData.DataChannel= f;
            sortedSpikeData.EventChannel= eventChannel;
            sortedSpikeData.ArduinoData = arduinoPseudoTable(correctTrialIdx,:);
            sortedSpikeData.Conditions= conds;
            sortedSpikeData.ConditionSegments= events;
            sortedSpikeData.Locations= "ARDUINO_OUTPUT_ERROR";
        else
            falseTrialIdx=false(size(recordedTrial,1),1);
            correctTrialIdx = find(~falseTrialIdx);
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
        end
        if(~exist(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
                '\Physiology\Results_All\'], 'dir'))
            mkdir(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
                '\Physiology\Results_All\']);
        end
        save(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDate,...
            '\Physiology\Results_All\',fullName,'.mat'],'sortedSpikeData','-mat','-append');
        clear sortedSpikeData spikeTimes segmentTimes
    end
end
fclose('all');
end