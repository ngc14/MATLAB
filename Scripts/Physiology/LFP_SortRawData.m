clear all
sessionDate = '07_25_2019';
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

events = repmat({{'TrialStart','GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}}, 3,1);
events(end+1) = {{'TrialStart','GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}};
events = events';
conds = {'Extra Small Sphere','Large Sphere','Photocell', 'Rest'};
eventChannel = 'SMA 1';
relevantTrials = [];
%% Process session
for f = 1:32
    fullName = ['Gilligan_', sessionDate,'_', num2str(f)];
    nFiles = dir(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate, '\Physiology\']);
    ns2FilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.'):end), '.ns2'), {nFiles.name});
    if(~isempty(ns2FilesInd))
        ns2Files = nFiles(ns2FilesInd);
        [~, nevInd] = max([ns2Files.bytes]);
        ns2DotInd = regexp(ns2Files(nevInd).name, '\.');
        ns2FileNumber = ns2Files(nevInd).name(ns2DotInd-1);
        filePath = ['S:\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\Gilligan_',sessionDate,'_000',num2str(ns2FileNumber),'.ns2'];
        

        [data, dataTime, Fs, eventTimes_risingEdge, eventTimes_fallingEdge ] = ...
            loadRippleData(filePath,f,eventChannel);
        
        % Find start of trial
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;

        % Find start of correct trials
        startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
        endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
        
        for i = 1:min(length(startEventIdx),length(endEventIdx))-1
            if(startEventIdx(i+1)<endEventIdx(i))
                startEventIdx(i+1) = NaN;
            end
        end
        for i = 1:min(length(startEventIdx),length(endEventIdx))-1
            if(endEventIdx(i)<startEventIdx(i))
                endEventIdx(i) = NaN;
            end
        end
        if(length(startEventIdx)~=length(endEventIdx))
            disp('Mistmatch lengths');
            return
        end
        if(length(startEventIdx)==length(endEventIdx))
            % segment into trials
            parfor n = 1:length(startEventIdx)
                [~, eventTimeInd] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n))));
                [~, eventTimeEndInd] = min(abs(dataTime-eventTimes_risingEdge(endEventIdx(n))));
                eventData{n} = data(eventTimeInd:eventTimeEndInd);
                [~, segmentTimes{n}] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
                segmentTimes{n} = segmentTimes{n}-eventTimeInd;
            end
            %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
            if(~exist(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Arduino\'],'dir'))
                mkdir(['S:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Arduino\']);
            end
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
            falseTrialIdx=cellfun(@str2num,recordedTrial(:,8),'UniformOutput',false);
            correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
            %correctTrialIdx = setdiff(correctTrialIdx, row);
            
            successfulTrial=recordedTrial(correctTrialIdx,1:7);
            eventData = eventData(correctTrialIdx);
            segmentTimes = segmentTimes(correctTrialIdx);
            
            sortedLFPData.LFPData= eventData;
            sortedLFPData.SegTimes = segmentTimes;
            
            sortedLFPData.Date=sessionDate;
            sortedLFPData.DataChannel=f;
            sortedLFPData.EventChannel=eventChannel;
            sortedLFPData.SampleRate=Fs;
            
            sortedLFPData.ArduinoData = successfulTrial;
            sortedLFPData.Conditions = conds;
            sortedLFPData.ConditionSegments = events;
            
            if(~exist(['U:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\'], 'dir'))
                mkdir(['U:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\']);
            end
            save(['U:\Lab\Gilligan\All Data\Gilligan_', sessionDate,'\Physiology\Results\',fullName],'sortedLFPData');
            clear sortedLFPData;
        else
            disp('Mismatch lengths');
        end
    end
end
function x = removeShortSegs(x)
if(isnumeric(str2double(x)))
    x = str2double(x)<= 10;
else
    x = 0;
end
end