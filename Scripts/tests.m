clear sortedLFPData;

%sessionDate = '04_12_2019';
%sessionDate = '06_17_2019';
% sessionDate = '07_16_2019';
% sessionDate = '07_15_2019';
% sessionDate = '04_08_2019';
% sessionDate = '05_20_2019';
% sessionDate = '07_23_2019';
sessionDate = '04_22_2019';

delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
conds = {'Extra Small Sphere','Large Sphere','Photocell', 'Rest'};
relevantTrials = [];
events = repmat({{'TrialStart','GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}}, 3,1);
events(end+1) = {{'TrialStart','GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'}};
events = events';



% ns2 file
% [fname, pathname] = uigetfile('*.ns2', 'Select a ns2 file');
% filename = strcat(pathname, fname);
% results_folder = strcat(pathname, fname);
% fid = fopen(filename, 'r');
% if(fid == -1)
%     disp('cannot open file', 'r');
%     return
% end
%filepath='\\pitt\sni\gharbawie\Lab\Shan\Matlab Scriptes\New_data\extractEvent\Gilligan_04_12_2019_0001.ns2';
%filepath='G:\Lab\Shan\Matlab Scriptes\New_data\extractEvent\Gilligan_04_12_2019_0001.ns2';
%filepath='\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_04_12_2019\Physiology\Gilligan_04_12_2019.ns2';
%filepath='D:\Omar lab work folder\projects\new_structure\extractEvent\Physiology\Gilligan_04_12_2019_0001.ns2';
%filepath='\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_04_12_2019\Physiology\Gilligan_04_12_2019_0001.ns2';
%filepath='G:\Lab\Gilligan\All Data\Gilligan_04_12_2019\Physiology\Gilligan_04_12_2019_0001.ns2';
%filepath='G:\Lab\Gilligan\All Data\Gilligan_07_08_2019\Physiology\Gilligan_07_08_2019_0001.ns2';
% filepath='G:\Lab\Gilligan\All Data\Gilligan_07_16_2019\Physiology\Gilligan_07_16_2019_0001.ns2';
%filepath='G:\Lab\Gilligan\All Data\Gilligan_06_17_2019\Physiology\Gilligan_06_17_2019_0001.ns2';
% filepath='G:\Lab\Gilligan\All Data\Gilligan_07_15_2019\Physiology\Gilligan_07_15_2019_0001.ns2';
% filepath='G:\Lab\Gilligan\All Data\Gilligan_04_08_2019\Physiology\Gilligan_04_08_2019_0001.ns2';
% filepath='G:\Lab\Gilligan\All Data\Gilligan_05_20_2019\Physiology\Gilligan_05_20_2019_0001.ns2';
filepath='\\univ.pitt.edu\sni\Gharbawie\Lab\Gilligan\All Data\Gilligan_07_15_2019\Physiology\Gilligan_07_15_2019_0001.ns2';


eventChannel = 'SMA 1';
ChanMap=[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32];
for f = 28:32  % there are 32 data channel
    %load raw data and event times
    clear sortedLFPData;
    
    [data, dataTime, Fs, eventTimes_risingEdge, eventTimes_fallingEdge ] = ...
        loadRippleData(filepath,ChanMap(f),eventChannel);
    
    %% P1
    % data is the whole session
    % change the data into bandpass filtered data
    data_filt2300 = bandpass(data,[2,300],Fs);
%     data_filt24 = bandpass(data,[2,4],Fs);
%     data_filt58 = bandpass(data,[5,8],Fs);
%     data_filt914 = bandpass(data,[9,14],Fs);
%     data_filt215 = bandpass(data,[2,15],Fs);
%     data_filt315 = bandpass(data,[3,15],Fs);
%     data_filt15100 = bandpass(data,[15,100],Fs);
%     data_filt1530 = bandpass(data,[15,30],Fs);
%     data_filt3050 = bandpass(data,[30,50],Fs);
%     data_filt50100 = bandpass(data,[50,100],Fs);
%     data_filt100150 = bandpass(data,[100,150],Fs);

    %%

    % Find pulse lengths
    pulseLength=eventTimes_fallingEdge-eventTimes_risingEdge;
    % Find start and end of trials
%     pulseLength(2156) = 0.02;%7_15
    startEventIdx=find(pulseLength>0.0181 & pulseLength<0.024);
    endEventIdx = find(pulseLength>0.030 & pulseLength<0.070);
%     startEventIdx(142) = [];%7_15
%     endEventIdx(142) = [];%7_15
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
        for n = 1:length(startEventIdx)
            [~, eventTimeInd] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n))));
            [~, eventTimeEndInd] = min(abs(dataTime-eventTimes_risingEdge(endEventIdx(n))));
            
            % generate the event data
            %% p2
%             eventData{n} = data(eventTimeInd:eventTimeEndInd); 
            eventData_filt2300{n} = data_filt2300(eventTimeInd:eventTimeEndInd);%this will return the LFP trials 
%             eventData_filt24{n} = data_filt24(eventTimeInd:eventTimeEndInd);
%             eventData_filt58{n} = data_filt58(eventTimeInd:eventTimeEndInd);
%             eventData_filt914{n} = data_filt914(eventTimeInd:eventTimeEndInd);
%             eventData_filt215{n} = data_filt215(eventTimeInd:eventTimeEndInd);            
%             eventData_filt315{n} = data_filt315(eventTimeInd:eventTimeEndInd);
%             eventData_filt15100{n} = data_filt15100(eventTimeInd:eventTimeEndInd); 
%             eventData_filt1530{n} = data_filt1530(eventTimeInd:eventTimeEndInd);
%             eventData_filt3050{n} = data_filt3050(eventTimeInd:eventTimeEndInd);
%             eventData_filt50100{n} = data_filt50100(eventTimeInd:eventTimeEndInd);
%             eventData_filt100150{n} = data_filt100150(eventTimeInd:eventTimeEndInd);
            %%
            
            [~, segmentTimes{n}] = min(abs(dataTime-eventTimes_risingEdge(startEventIdx(n)+1:endEventIdx(n)-1)));
            segmentTimes{n} = segmentTimes{n}-eventTimeInd;
        end
        %% EXTRACT INFORMATION FROM ARDUINO OUTPUT FILE
        %arduinoFilename = ['\\pitt\sni\gharbawie\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',sessionDate,'\Gilligan_',sessionDate,'.txt'];
        arduinoFilename = ['\\univ.pitt.edu\sni\Gharbawie\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',sessionDate,'\Gilligan_',sessionDate,'.txt'];
        
        % Open the text file.
        fileID = fopen(arduinoFilename,'r');
        
        % Read columns of data according to the format.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
        
        % Close the text file.
        fclose(fileID);
        
        % Create output variable
        recordedTrial = [dataArray{1:end-1}];
        
        % remove lines that are not task related (pause, headers)
%         recordedTrial = recordedTrial(ismember(recordedTrial(:,1), conds),:);
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
        
        
        %% P3
%         eventData = eventData(correctTrialIdx);
        eventData_filt2300=eventData_filt2300(correctTrialIdx);
%         eventData_filt24=eventData_filt24(correctTrialIdx);
%         eventData_filt58=eventData_filt58(correctTrialIdx);
%         eventData_filt914=eventData_filt914(correctTrialIdx);
%         eventData_filt215=eventData_filt215(correctTrialIdx);
%         eventData_filt315=eventData_filt315(correctTrialIdx);
%         eventData_filt15100=eventData_filt15100(correctTrialIdx);
%         eventData_filt1530=eventData_filt1530(correctTrialIdx);
%         eventData_filt3050=eventData_filt3050(correctTrialIdx);
%         eventData_filt50100=eventData_filt50100(correctTrialIdx);
%         eventData_filt100150=eventData_filt100150(correctTrialIdx);
        %%
        
        
        segmentTimes = segmentTimes(correctTrialIdx);
        
        % the filtered data will be stored in the data structure
        %% P4
%         sortedLFPData.LFPData_org = eventData;
        sortedLFPData.LFPData_2300 = eventData_filt2300;
%         sortedLFPData.LFPData_24 = eventData_filt24;
%         sortedLFPData.LFPData_58 = eventData_filt58;
%         sortedLFPData.LFPData_914 = eventData_filt914;
%         sortedLFPData.LFPData_215 = eventData_filt215;
%         sortedLFPData.LFPData_315 = eventData_filt315;
%         sortedLFPData.LFPData_15100 = eventData_filt15100;
%         sortedLFPData.LFPData_1530 = eventData_filt1530;
%         sortedLFPData.LFPData_3050 = eventData_filt3050;
%         sortedLFPData.LFPData_50100 = eventData_filt50100;
%         sortedLFPData.LFPData_100150 = eventData_filt100150;
        %%
        
        
        sortedLFPData.SegTimes = segmentTimes;
        
        sortedLFPData.Date=sessionDate;
        sortedLFPData.DataChannel=f;
        sortedLFPData.EventChannel=eventChannel;
        sortedLFPData.SampleRate=Fs;
        
        sortedLFPData.ArduinoData = successfulTrial;
        sortedLFPData.Conditions = conds;
        sortedLFPData.ConditionSegments = events;
        % save(['D:\Omar lab work folder\projects\new_structure\extractEvent\output','sortedLFPData']); %% save the data, not rewrite
        %save(['G:\Lab\Shan\Matlab Scriptes\New_data\extractEvent\Method01\dataOutputs\E3_215_315\',strcat('sortedLFPData',sessionDate,'_channel_',num2str(f))],'sortedLFPData');
        %save(['G:\Lab\Shan\Matlab Scriptes\Ripple data\extractEvent\Method01_inUse\dataOutputs\41219_2300\',num2str(f)],'sortedLFPData');
%         save(['G:\Lab\Shan\Matlab Scriptes\Ripple data\extractEvent\Method01_inUse\dataOutputs\071519_2300\',num2str(f)],'sortedLFPData');
%         save(['G:\Lab\Shan\Matlab Scriptes\Ripple data\extractEvent\Method01_inUse\dataOutputs\040819_2300\',num2str(f)],'sortedLFPData');
%         save(['G:\Lab\Shan\Matlab Scriptes\Ripple data\extractEvent\Method01_inUse\dataOutputs\0715_19\',num2str(f)],'sortedLFPData');
        save(['G:\Lab\Shan\Matlab Scriptes\Ripple data\extractEvent\Method01_inUse\dataOutputs\92319_2300\',num2str(f)],'sortedLFPData');
    else
        disp('Mistmatch trial lengths');
    end
end