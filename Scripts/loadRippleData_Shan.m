for f = 1:32
    %load raw data and event times
    [data, dataTime, Fs, eventTimes_risingEdge, eventTimes_fallingEdge ] = ...
        loadRippleData(filePath,f,eventChannel);
    
    % Find pulse lengths
    pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    % Find start and end of trials
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
    if(length(startEventIdx)==length(endEventIdx))
        % segment into individual trials
        parfor n = 1:length(startEventIdx)
            [~, eventTimeInd] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n))));
            [~, eventTimeEndInd] = min(abs(dataTime-eventTimes_risingEdge(endEventIdx(n))));
            eventData{n} = data(eventTimeInd:eventTimeEndInd);
            [~, segmentTimes{n}] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
            segmentTimes{n} = segmentTimes{n}-eventTimeInd;
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
        recordedTrial = recordedTrial(ismember(recordedTrial(:,1), conds) & ~cellfun(@isempty, recordedTrial(:,8)),:);
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
    else
        disp('Mistmatch trial lengths');
    end
end