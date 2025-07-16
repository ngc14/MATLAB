clc
clear all
close all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
monkey = 'Gilligan';
date = '02_25_2021';
run = '0001';
stimChannel = 'SMA 2';
eventChannel = 'SMA 1';
time_before_event = -.1;
time_after_event = .1;
%%
dirName = ['\\pitt\sni\gharbawie\Lab\', monkey,'\All Data\', monkey,'_', ...
    date, '\EMG\Results_New\'];
directory = dir([dirName '*.mat']);
for f=1:length(directory)
    load([directory(f).folder, '\',directory(f).name])
    badTrials = EMG_getBadTrials(sortedEMGData.EMGData, sortedEMGData.SegTimes, sortedEMGData.SampleRate);
    if(any(badTrials))
        trialLengths = cellfun(@length, sortedEMGData.EMGData(find(badTrials)));
        nanTrials = arrayfun(@(a) nan(1,a), trialLengths, 'UniformOutput', false);
        sortedEMGData.EMGData(find(badTrials)) = nanTrials;
    end
    EMGSig{f} = sortedEMGData.EMGData;
    EMGTimes{f} = sortedEMGData.SegTimes;
    muscles{f} = sortedEMGData.Muscle;
end
disp('EMG data loaded');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['\\univ.pitt.edu\sni\Gharbawie\Lab\' monkey, '\All Data\', monkey,'_', date,'\Physiology\', monkey, '_', date, '_',run, '.nev'];
[~, hFile] = ns_OpenFile(filename);
eventEntityID = find(strcmp({hFile.Entity(:).Reason}, stimChannel));
% get event entity information
[~, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID);
% event data
numCount = eventEntityInfo.ItemCount;

[~, eventTimeStamps, eventData,~] = arrayfun(@(a) ns_GetEventData(hFile,eventEntityID,a), 1:numCount, 'UniformOutput', false);

% Only consider time stamps for rising edge (eventData = 0 on falling edge)
% This may not work with parallel input
tempValidTimeStampInds = find([eventData{:}] == 32767);
eventTimes = [eventTimeStamps{tempValidTimeStampInds}];

ns_CloseFile(filename);
disp('Stim data loaded');
%%
[ns_status, hFile] = ns_OpenFile(filename);
eventEntityID = find(strcmp({hFile.Entity(:).Reason}, eventChannel));

% get event entity information
[ns_RESULT, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID);

% event data
numCount = eventEntityInfo.ItemCount;
[~, eventTimeStamps, eventData,~] = arrayfun(@(a) ns_GetEventData(hFile,eventEntityID,a), 1:numCount, 'UniformOutput', false);
% Only consider time stamps for rising edge (eventData = 0 on falling edge)
% This may not work with parallel input
ns_CloseFile(filename);
disp('Event data loaded');
%%
filename = ['\\univ.pitt.edu\sni\Gharbawie\Lab\' monkey, '\All Data\', monkey,'_', date,'\EMG\', monkey, '_', date, '_',run, '.nf3'];

[ns_status, hFile] = ns_OpenFile(filename);

%% Determine correct entityID for desired datastream

% find index of specified electrode
EntityIndices = find(strncmp({hFile.Entity(:).Label} , 'hi-res', length('hi-res')));
fileTypeNums = [hFile.Entity(EntityIndices).FileType];
fileTypes = {hFile.FileInfo(fileTypeNums).Type};
entityID = EntityIndices(find(contains(fileTypes,filename(end-2:end))));
entityID = entityID(1);
[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID);
TimeStamps = hFile.FileInfo(hFile.Entity(entityID).FileType).TimeStamps;
numSamples = sum(TimeStamps(:,end));
analogInputDataTime_s = (0:numSamples-1)' ./ analogInfo.SampleRate;
if rem(length(analogInputDataTime_s),2)==1
    analogInputDataTime_s = analogInputDataTime_s(1:end-1);
end
ns_CloseFile(filename);
disp('EMG times loaded');
%%
%EMGSig = cellfun(@(a) [a{:}], EMGSig, 'UniformOutput', false);
eventTimeStamps = [eventTimeStamps{:}];
risingEdges = eventTimeStamps([eventData{:}] == 32767);
fallingEdges = eventTimeStamps([eventData{:}] == 0);
pulseLength=fallingEdges-risingEdges;
% Find start of correct trials
startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
tsInd = 1;
for n = 1:length(startEventIdx)
    if(risingEdges(startEventIdx(n))< eventTimes(tsInd) && risingEdges(endEventIdx(n))>eventTimes(tsInd))
        [~, eventTimeInd] = min(abs(analogInputDataTime_s-risingEdges(startEventIdx(n))));
        [~, stimTimes{n}] = min(abs(analogInputDataTime_s-eventTimes(tsInd)));
        stimTimes{n} = stimTimes{n}- eventTimeInd;
        tsInd = tsInd + 1;
    else
    end
end
failedTrials = find(pulseLength(endEventIdx)>.05);
startEventIdx(failedTrials) = [];
endEventIdx(failedTrials) = [];
stimTimes(failedTrials) = [];

%%
alignedTrial = cell(size(EMGSig));
tb = time_before_event * analogInfo.SampleRate;
ta = time_after_event * analogInfo.SampleRate;
for m = 1:length(muscles)
    s = subplot(1,length(muscles), m);
    title(s, muscles{m});
    hold on;
    alignedTrial{m} = cellfun(@(a,b) a(b+tb:b+ta), EMGSig{m}, stimTimes, 'UniformOutput', false);
    cellfun(@(a) plot(a), alignedTrial{m}, 'UniformOutput', false);
    matSig = reshape(cell2mat(alignedTrial{m})',[max(cellfun(@length,alignedTrial{m})),...
        length(alignedTrial{m})-sum(cellfun(@isempty,alignedTrial{m}))])';
    plot(nanmean(matSig), 'k', 'LineWidth', 2);
    xticks([0, round(max(cellfun(@length,alignedTrial{m}))/2), max(cellfun(@length,alignedTrial{m}))]);
    xticklabels({num2str(time_before_event), num2str(0), num2str(time_after_event)});
end