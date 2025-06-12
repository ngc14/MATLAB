sessionDate='04_09_2018';
nf3FileNumber=1; %generally 1 or 2
dataChannel = 130;          % EMG data channel. Must be the same length as the muscle variable.
eventChannel = 'SMA 1';     % digital event channel name
notch = 60;
%%
graspConditions = {'Small Sphere','Large Sphere','Extra Small Sphere'}; % spacing and capitalization must match Arduino text file
nonGraspConditions = {'Photocell'};
conditions={graspConditions{:},nonGraspConditions{:}};
filePath = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\EMG\Gilligan_',sessionDate,'_000',num2str(nf3FileNumber),'.nf3'];
arduinoFilename = ['\\130.49.229.252\gharbawie\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',sessionDate,'\Gilligan_',sessionDate,'.txt'];
[data, dataTime, Fs, eventTimes_risingEdge, eventTimes_fallingEdge ] = loadRippleData(filePath,dataChannel,eventChannel);

% Find start of trial
pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;

% Find start of correct trials
startEventIdx=find(pulseLength>0.017 & pulseLength<0.025);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
if(length(startEventIdx)~=length(endEventIdx))
    disp('Mismatch with pulses, check beginning of pulseLengths');
    return;
end
%%
for n = 1:length(startEventIdx) % step through each event time
    % find index of event n
    [~, eventTimeInd] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n))));
    [~, eventTimeEndInd] = min(abs(dataTime-eventTimes_risingEdge(endEventIdx(n))));
    eventData{n} = data(eventTimeInd:eventTimeEndInd);
    [~, segmentTimes{n}] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
    segmentTimes{n} = segmentTimes{n}-eventTimeInd;
end
%%
badTrial = zeros(size(eventData));
figure();
hold on;
colors = {'-r', '-g', '-b', '-y', '-k'};
for t=1:size(eventData,2)
    trialSig = eventData{t};
    rectified = movmean(abs(trialSig),50);
    baseline = nanmean(rectified(1:segmentTimes{t}(1)));
    Y = fft(trialSig);
    P2 = abs(Y/length(trialSig));
    P1 = P2(1:length(trialSig)/2 + 1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(length(trialSig)/2))/length(trialSig);
    if(any(find(f<15 & P1>4)))
        badTrial(t) = 1;
    elseif(range(abs(trialSig)) < 5)
        badTrial(t) = 2;
    elseif(baseline>80)
        badTrial(t) = 3;
    elseif(P1>7 & mod(f,notch)~=0 & mod(f,notch)<notch+.5 & mod(f,notch)>notch-.5)
        badTrial(t) = 4;
    end
    plot(f, P1, colors{badTrial(t)+1});
end

%% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE

% Initialize variables.
delimiter = '\t';

% Each line of text is read in as text (%s)
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

% Open the text file.
fileID = fopen(arduinoFilename,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
recordedTrial = [dataArray{1:end-1}];

% Clear temporary variables
clearvars delimiter formatSpec fileID dataArray;

% Clean up successfulTrial variable
recordedTrial = recordedTrial(ismember(recordedTrial(:,1), [graspConditions nonGraspConditions]),:); % remove lines that are not task related (pause, headers)
[row,~]=find(cellfun(@(x) removeShortSegs(x), recordedTrial(:,2:6)));

% Remove failed trials.
falseTrialIdx=cellfun(@str2num,recordedTrial(:,7),'UniformOutput',false);
correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
correctTrialIdx = setdiff(correctTrialIdx, row);
successfulTrial=recordedTrial(correctTrialIdx,1:7);

correctTrialCondition=successfulTrial(:,1);
badTrial = badTrial(correctTrialIdx);
eventData = eventData(correctTrialIdx);
segmentTimes = segmentTimes(correctTrialIdx);
%%
for i = 1:20
    if(~mod(i,5)==0)
        tempIdx=find(strcmp(correctTrialCondition,conditions(mod(i,5))));
        tempCondition=conditions{mod(i,5)};
        tempCondition=tempCondition(~isspace(tempCondition));
        sampleIdx = badTrial(tempIdx);
        sampleIdx = find(sampleIdx==0);
        plotIdx = randsample(length(sampleIdx), 1);
    else
        tempIdx = 1:length(badTrial);
        sampleIdx = find(badTrial==3);
        if(~any(sampleIdx))
            sampleIdx = find(badTrial==0);
            plotIdx = randsample(length(sampleIdx),1);
            tempCondition = correctTrialCondition(tempIdx(sampleIdx(plotIdx)));
        else
            plotIdx = randsample(length(sampleIdx),1);
            tempCondition = ['Bad ', correctTrialCondition{tempIdx(sampleIdx(plotIdx))}];
        end
    end
    s = subplot(4,5,i);
    
    sig = movmean(abs(eventData{tempIdx(sampleIdx(plotIdx))}),50);
    plot(s,sig);
    if(contains(tempCondition,'Photocell'))
        startWithdrawl = segmentTimes{tempIdx(sampleIdx(plotIdx))}(5);
    else
        startWithdrawl = segmentTimes{tempIdx(sampleIdx(plotIdx))}(6);
    end
    startReach = segmentTimes{tempIdx(sampleIdx(plotIdx))}(3);
    hold on;
    ymax = max(sig);
    plot([startReach, startReach], [0 ymax], '--r');
    plot([startWithdrawl, startWithdrawl], [0, ymax], '--r');
    xlabel('Frame');
    ylabel('uV');
    title(tempCondition);
    hold off;
end
figure('Name', sprintf('Session: %s, Channel: %d', sessionDate, dataChannel))
a = subplot(1,2,1);
xlabel('Frame');
ylabel('uV');
title(sprintf('Good Trials: %d out of %d', sum(~badTrial), length(badTrial)));
hold on;
b = subplot(1,2,2);
title(sprintf('Bad Trials: %d out of %d\n%d Bad Fourier (r), %d Dead (g),\n %d Baseline (b), %d Power (y)',...
    sum(badTrial>0), length(badTrial), sum(badTrial==1), sum(badTrial==2), sum(badTrial==3), sum(badTrial==4)));
xlabel('Frame');
ylabel('uV');
hold on;
for i = 1:size(badTrial,2)
    if(~badTrial(i))
        plot(a,movmean(abs(eventData{i}), 50));
    else
        plot(b,movmean(abs(eventData{i}), 50),colors{badTrial(i)});
    end
end
%savefig(['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',...
%    sessionDate,'\EMG\Results\Processed_', num2str(dataChannel)]);

function x = removeShortSegs(x)
if(isnumeric(str2double(x)))
    x = str2double(x)< 10;
else
    x = 0;
end
end