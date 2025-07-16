function EMG_SortRawData(dateIn,monkeyIn)
if(~exist('dateIn','var'))
    sessionDate = '08_05_2019';
else
    sessionDate = dateIn;
end
if(~exist('monkey','var'))
    monkey = 'Gilligan';
else
    monkey = monkeyIn;
end
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

events = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}});
conds = {'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'};
eventChannel = 'SMA 1';
relevantTrials = [];
%% Process session
file = dir(['S:\Lab\',monkey,'\All Data\', monkey,'_', ...
    sessionDate, '\EMG\*.txt']);
fileID = fopen([file.folder,'\',file.name]);
fileInd = 1;
muscles = {};
while(~feof(fileID))
    line = fgetl(fileID);
    if(contains(line, 'Channel B','IgnoreCase',true))
        channelNumStart = regexp(line, 'Channel B','end');
        channelNumEnd = regexp(line, ':');
        channelNum(fileInd) = 128 + str2double(line(channelNumStart+1:channelNumEnd-1));
        muscles{end+1} = line(channelNumEnd+2:end);
        fileInd = fileInd + 1;
    elseif(contains(line, 'Channel A', 'IgnoreCase', true))
        channelNumStart = regexp(line, 'Channel A','end');
        channelNumEnd = regexp(line, ':');
        channelNum(fileInd) = str2double(line(channelNumStart+1:channelNumEnd-1));
        muscles{end+1} = line(channelNumEnd+2:end);
        fileInd = fileInd + 1;
    end
end
fclose(fileID);
channelNum = unique(channelNum);
channelNum = channelNum(~cellfun(@isempty,muscles) & cellfun(@(m) any(isletter(m)),muscles));
for f = 1:length(channelNum)
    muscle = muscles{f};
    fullName = [monkey,'_', sessionDate,'_', muscle];
    nFiles = dir(['S:\Lab\',monkey,'\All Data\',monkey,'_', sessionDate, '\EMG\']);
    nf3FilesInd = cellfun(@(a) strcmp(a(regexp(a, '\.','once'):length(a)), '.nf3'), {nFiles.name});
    if(~isempty(nf3FilesInd))
        nf3Files = nFiles(nf3FilesInd);
        [~, nf3Ind] = max([nf3Files.bytes]);
        nf3DotInd = regexp(nf3Files(nf3Ind).name, '\.');
        nf3FileNumber = nf3Files(nf3Ind).name(nf3DotInd-1);
        filePath = dir(['S:\Lab\',monkey,'\All Data\',monkey,'_',sessionDate,'\EMG\*000',num2str(nf3FileNumber),'.nf3']);
        filePath = [filePath(nf3Ind).folder, '\' filePath(nf3Ind).name];
        if(f>1)
            eventChannel = [];
            [data, dataTime, Fs, ~, ~ ] = loadRippleData(char(filePath),channelNum(f),eventChannel);
        else
            [data, dataTime, Fs, eventTimes_risingEdge, eventTimes_fallingEdge ] = ...
                loadRippleData(char(filePath),channelNum(f),eventChannel);
        end
        % Find start of trial
        pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
        % Find start of correct trials
        switch(datestr(sessionDate,'yyyy_mm_dd'))
            case('2021_09_23')
                eventTimes_fallingEdge(629) = eventTimes_risingEdge(629)+0.02;
            case('2021_09_24')
                eventTimes_fallingEdge(713) = eventTimes_risingEdge(713)+0.02;
                eventTimes_fallingEdge(2425) = eventTimes_risingEdge(2425)+0.02;
            case('2021_09_27')
                eventTimes_fallingEdge(1121) = eventTimes_risingEdge(1121)+0.02;
            case('2021_09_28')
                eventTimes_fallingEdge(1613) = eventTimes_risingEdge(1613)+0.02;
            case('2021_09_29')
                eventTimes_fallingEdge(382) = eventTimes_risingEdge(382)+0.02;
            case('2021_09_30')
                eventTimes_fallingEdge(1405) = eventTimes_risingEdge(1405)+0.02;
                eventTimes_fallingEdge(1731) = eventTimes_risingEdge(1731)+0.02;
                eventTimes_fallingEdge(1754) = eventTimes_risingEdge(1754)+0.02;
            case('2021_10_01')
                eventTimes_fallingEdge(660) = eventTimes_risingEdge(660)+0.02;
                eventTimes_fallingEdge(701) = eventTimes_risingEdge(701)+0.02;
                eventTimes_fallingEdge(1124) = eventTimes_risingEdge(1124)+0.02;
            case('2021_10_04')
                eventTimes_fallingEdge(928) = eventTimes_risingEdge(928)+0.02;
                eventTimes_fallingEdge(1148) = eventTimes_risingEdge(1148)+0.02;
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
            case('2021_10_26')
                pulseLength(1560) = 0.02;
            case('2021_10_27')
                pulseLength(396) = 0.02;
            case('2019_04_16')
                pulseLength(1792) = 0.01;
                pulseLength(2655) = 0.01;
            case('2019_04_17')
                pulseLength(1573) = 0.01;
            case('2019_04_19')
                pulseLength(2165) = 0.01;
                pulseLength(2596) = 0.01;
            case('2019_04_22')
                pulseLength(652) = 0.01;
                pulseLength(112) = 0.01;
                pulseLength(2287) = 0.01;
                pulseLength(1823) = 0.01;
                pulseLength(2238) = 0.01;
            case('2019_04_23')
                pulseLength(2018) = 0.01;
                pulseLength(2432) = 0.01;
                pulseLength(2756) = 0.01;
            case('2019_04_26')
                pulseLength(619) = 0.01;
                pulseLength(1091) = 0.01;
                pulseLength(2018) = 0.04;
            case('2019_04_29')
                pulseLength(1425) = 0.01;
            case('2020_06_20')
                pulseLength(1425) = 0.01;
            case('2020_06_03')
                pulseLength(2062:end) = 0.01;
            case('2020_06_05')
                pulseLength(3092) = 0.01;
            case('2019_06_10')
                pulseLength(2329:2337) = [];
            case('2020_07_02')
                pulseLength(110) = 0.01;
                pulseLength(659) = 0.01;
                pulseLength(668) = 0.01;
                pulseLength(942) = 0.01;
                pulseLength(2742) = 0.02;
            case('2020_07_09')
                pulseLength(295) = 0.04;
                pulseLength(989) = 0.01;
                pulseLength(1429) = 0.04;
                pulseLength(1456) = 0.01;
            case('2023_03_28')
                pulseLength(331) = 0.06;
                pulseLength(351) = 0.06;
                pulseLength(352) = 0.02;
                pulseLength(388) = 0.02;
                pulseLength(2175) = 0.01;
            case('2023_04_06')
                pulseLength(240:243) = 0.01;
                pulseLength(1995) = 0.01;
            case('2019_09_05')
                pulseLength(824)= 0.01;
            case('2019_11_26')
                pulseLength(1835)= 0.01;
            case('2019_04_18')
                pulseLength(1272) = 0.01;
            case('2019_06_04')
                pulseLength(308) = 0.01;
            case('2019_06_07')
                pulseLength(1997) = 0.01;
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
        if(length(startEventIdx)~=length(endEventIdx))
            disp('Mistmatch lengths');
            disp(sessionDate)
            return;
        end
        if(length(startEventIdx)==length(endEventIdx))
            parfor n = 1:length(startEventIdx) % step through each event time
                % find index of event n
                [~, eventTimeInd] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1)));
                [~, eventTimeEndInd] = min(abs(dataTime-eventTimes_risingEdge(endEventIdx(n)-1)));
                eventData{n} = data(eventTimeInd:eventTimeEndInd);
                [~, segmentTimes{n}] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
                %segmentTimes{n} = segmentTimes{n}- eventTimeInd;
                segmentTimes{n} = dataTime(segmentTimes{n});
            end
            %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
            if(strcmp(monkey,'Gilligan'))
                if(~exist(strjoin(["S:\Lab\Gilligan\All Data\Gilligan_", sessionDate,"\Arduino\Gilligan_",sessionDate, ".txt"],''),'file'))
                    mkdir(strjoin(["S:\Lab\Gilligan\All Data\Gilligan", sessionDate,"\Arduino\"],''));
                    copyfile(strjoin(["S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_",sessionDate,"\Gilligan_",sessionDate,".txt"],''),...
                        strjoin(["S:\Lab\Gilligan\All Data\Gilligan_", sessionDate, "\Arduino\Gilligan_",sessionDate, ".txt"],''));
                end
            else
                if(~exist(strjoin(["S:\Lab\Skipper\All Data\Skipper_", sessionDate,"\Arduino\Skipper_",sessionDate, ".txt"],''),'file'))
                    mkdir(strjoin(["S:\Lab\Skipper\All Data\Skipper_", sessionDate,"\Arduino\"],''));
                    copyfile(strjoin(["S:\Lab\BehaveBackup\Monkey_Training\Skipper_Macaque_152_17\Skipper_",datestr(sessionDate,'yyyy_mm_dd'),"\Skipper_",datestr(sessionDate,'yyyy_mm_dd'),".txt"],''),...
                        strjoin(["S:\Lab\Skipper\All Data\Skipper_", sessionDate, "\Arduino\Skipper_",sessionDate, ".txt"],''));
                end
            end
            arduinoFilename = strjoin(["S:\Lab\",monkey,"\All Data\",monkey,"_", sessionDate, "\Arduino\",monkey,"_",sessionDate, ".txt"],'');

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
            falseTrialIdx=cellfun(@(s) strcmp(s,"False Start") | strcmp(s,"Failed-to-Reach"),recordedTrial(:,8),'UniformOutput',false);
            correctTrialIdx = find(~cell2mat(falseTrialIdx));
            correctTrialIdx = ~cellfun(@isnan,cellfun(@str2double,recordedTrial(:,8),'UniformOutput',false));

            successfulTrial=recordedTrial(correctTrialIdx,1:end);
            eventData = eventData(correctTrialIdx);
            segmentTimes = segmentTimes(correctTrialIdx);

            sortedEMGData.EMGData= eventData;
            sortedEMGData.SegTimes = segmentTimes;

            sortedEMGData.Date=sessionDate;
            sortedEMGData.Muscle=muscle;
            sortedEMGData.DataChannel=channelNum(f);
            sortedEMGData.EventChannel=eventChannel;
            sortedEMGData.SampleRate=Fs;

            sortedEMGData.ArduinoData = successfulTrial;
            sortedEMGData.Conditions = unique(successfulTrial(:,1));
            sortedEMGData.ConditionSegments = values(events,sortedEMGData.Conditions);

            if(~exist(strjoin(["S:\Lab\",monkey,"\All Data\",monkey,"_", sessionDate,"\EMG\All_Trials\"],''), 'dir'));
                mkdir(strjoin(["S:\Lab\",monkey,"\All Data\",monkey,"_", sessionDate,"\EMG\All_Trials\"],''));
            end
            save(strjoin(["S:\Lab\",monkey,"\All Data\",monkey,"_", sessionDate,"\EMG\All_Trials\",fullName],''),'sortedEMGData');
            clear sortedEMGData;
        else
            disp('Mismatch lengths');
        end
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