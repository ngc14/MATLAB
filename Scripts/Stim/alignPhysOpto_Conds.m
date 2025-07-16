clc
clear all
close all

monkey = 'Gilligan';
sessionDate = '06_22_2020';
run = '0003';
eventChannels = {'SMA 1', 'SMA 4'};
numColumns = 4;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds

chanMap = [1:2:32,2:2:32];
chanMap = [1:2:18,19:2:32,2:2:17];
chanMap = [32:-1:1];

sigma = .6; % smoothing window in ms

Kernel = (-3*sigma:3*sigma);
kernelBin = length(Kernel);
Kernel = (-kernelBin/2:kernelBin/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

buffer_time = 1;
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
eventTs = cell(1,length(eventChannels));
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filePath = ['S:\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\Gilligan_',sessionDate,'_',run,'.nev'];
[ns_status, hFile] = ns_OpenFile(filePath);

for l = 1:length(eventChannels)
    eventEntityID = find(cellfun(@strcmpi, {hFile.Entity.Reason}, repmat({eventChannels{l}}, size({hFile.Entity.Reason}))));
    
    [ns_RESULT, eventEntityInfo] = ns_GetEntityInfo(hFile, eventEntityID(end));
    
    % event data
    numCount = eventEntityInfo.ItemCount;
    eventData = NaN(1, numCount);
    eventTimeStamps = NaN(1, numCount);
    dataSize = NaN(1, numCount);
    
    % get actual event time stamps
    for i = 1:numCount
        [~, eventTimeStamps(i), eventData(i), ~] = ns_GetEventData(hFile, eventEntityID, i);
    end
    fallingEdges = eventTimeStamps(eventData==0);
    risingEdges = eventTimeStamps(eventData==32767);
    
    eventTimeStamps = risingEdges;
    fallingPulse = fallingEdges;
    eventTimeStamps([Inf,diff(eventTimeStamps)]<2) = [];
    
    fallingPulse = fallingPulse(end:-1:1);
    fallingPulse(abs([Inf,diff(fallingPulse)])<2) = [];
    fallingPulse = fallingPulse(end:-1:1);
    
    time_before_event{l} = -2;    % seconds before event
    eventTs{l} = eventTimeStamps;
    fallingPulses{l} = fallingPulse;
    
    disp('Loading Event Data... 100% Done.');
    if(strcmp(eventChannels{l}, 'SMA 1'))
        
        pulseLength=fallingEdges-risingEdges;
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
        startTrialTs = risingEdges(startEventIdx);
        endTrialTs = risingEdges(endEventIdx);
        
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
        recordedTrial = recordedTrial(ismember(recordedTrial(:,9), {'Close', 'Far'}),:);
        
        % Remove failed trials.
        falseTrialIdx=cellfun(@str2num,recordedTrial(:,8),'UniformOutput',false);
        correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
        
        successfulTrial=recordedTrial(correctTrialIdx,1:7);
    end
end
ns_CloseFile(filePath);

%%
for l = 2:length(eventChannels)
    trialPulses{l} = zeros(1,length(recordedTrial));
    trialPulsesEnd{l} = zeros(1,length(recordedTrial));
    for i = 1:length(eventTs{l})
        trialPulses{l}(find(eventTs{l}(i)>=startTrialTs & eventTs{l}(i)<=endTrialTs)) = eventTs{l}(i);
        trialPulsesEnd{l}(find(eventTs{l}(i)>=startTrialTs & eventTs{l}(i)<=endTrialTs)) = fallingPulses{l}(i);
    end
    segLengths = endEventIdx-startEventIdx;
    
    trialPulses{l}(segLengths==5) = startTrialTs(segLengths==5)+5;
    trialPulsesEnd{l}(segLengths==5) = startTrialTs(segLengths==5)+6;
    
    trialPulses{l} = trialPulses{l}(correctTrialIdx);
    trialPulsesEnd{l} = trialPulsesEnd{l}(correctTrialIdx);
    startTrialTs = startTrialTs(correctTrialIdx);
end
%%
for l = 1:length(chanMap)
    [dataTime, ids] = loadSpikeData(filePath,chanMap(l));
    sortedIDs = unique(ids(ids<255));
    for u = 1:length(sortedIDs)
        dataUnits = dataTime(ids==sortedIDs(u));
        spikeTimes{l,u} = dataUnits;
        spikeTimes{l,u} = dataTime;
    end
end

%%
numRows = ceil(length(chanMap)/numColumns);
conds = unique(recordedTrial(:,1));
condInds = cellfun(@(a) regexp(a, 'Stim_', 'end'), conds, 'UniformOutput', false);
noStimConds = cellfun(@(a,b) a(b+1:end), conds, condInds, 'UniformOutput', false);
numberOfConds = unique(noStimConds);
numberOfConds = numberOfConds(~cellfun(@isempty,numberOfConds));
colors = distinguishable_colors(5,{'k'});
for n = 1:length(conds)
    cond = conds{n};
    stimCond = 0;
    trialInds = find(strcmp(successfulTrial(:,1),cond));
    for e = 2:length(eventChannels)
        alignedTrial = cell(size(spikeTimes));
        stim_ts_cond = trialPulses{e}(trialInds);
        if(~any(stim_ts_cond))
            stim_ts_cond = startTrialTs(trialInds)+5;
            time_after_event{e} = 3;
        else
            time_after_event{e}  = round(mean(trialPulsesEnd{e}(trialInds)-stim_ts_cond)+buffer_time,3,'significant'); 
            
        end
        for s = 1:length(stim_ts_cond)
            centeredSpikes = cellfun(@(a) a-stim_ts_cond(s), spikeTimes, 'UniformOutput', false);
            centeredSpikes = cellfun(@(a) a(a>=time_before_event{e}-length(Kernel)/(2*1/bin) & a<=time_after_event{e}+length(Kernel)/(2*1/bin)), centeredSpikes, 'UniformOutput', false);
            [row, col] = find(~cellfun(@isempty, centeredSpikes));
            for i = 1:length(row)
                alignedTrial{row(i), col(i)}{end+1,:} = centeredSpikes{row(i), col(i)};
            end
        end
        if(contains(cond, numberOfConds) && contains(cond,'Stim') && ~contains(cond,'Rest'))
            figureInd = find(cellfun(@(a) contains(cond,a),numberOfConds));
            set(0, 'CurrentFigure', figures{figureInd});
            colorCond = colors(n-figureInd,:);
            stimCond = 1;
        else
            figures{n} = figure('units','normalized','outerposition',[0 0 1 1]);
            colorCond = 'k';
            title(cond);
        end
        plotMap = reshape(1:length(chanMap),[numRows,numColumns])';
        plotMap = reshape(plotMap,[1,length(chanMap)]);
        hold on;
        for r = 1:length(alignedTrial)
            s = subplot(numRows,numColumns,r);
            hold on;
            yPos = 0;
            %line(s,[time_before_event{e},time_after_event{e}],[yPos, yPos], 'Color','r', 'LineWidth', 2);
            
            %%PSTH
            for c = 1:size(alignedTrial,2)
                rasterUnit = alignedTrial{plotMap(r),c};
                if(~isempty(rasterUnit))
                    trialHists = cellfun(@(a) histcounts(a,...
                        time_before_event{e}-length(Kernel)/(2*1/bin):bin:...
                        time_after_event{e}+length(Kernel)/(2*1/bin)), rasterUnit,'UniformOutput', false);
                    trialHistsSmooth = cellfun(@(a) conv(a,Kernel), trialHists, 'UniformOutput', false);
                    trialHistsSmooth = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
                    totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
                    
                    axes(s);
                    shadedErrorBar(time_before_event{e}:bin:time_after_event{e},nanmean(totalHists,1)./bin, nanstd(totalHists,1)./bin,'lineProps', {'LineWidth', 2, 'Color', colorCond});
                    yPos = max(yPos,max(nanmean(totalHists,1))/bin);
                    
                end
                
            end
            
            %%RASTER
%                             for c = 1:size(alignedTrial,2)
%                                 rasterUnit = alignedTrial{plotMap(r),c};
%                                 if(~isempty(rasterUnit))
%                                     if(length(rasterUnit)>60)
%                                         rasterMap = 1:round(length(rasterUnit)/60):length(rasterUnit);
%                                         totalUnits = length(rasterMap);
%                                     else
%                                         rasterMap = 1:length(rasterUnit);
%                                         totalUnits = length(rasterUnit);
%                                     end
%                                     for t = 1:totalUnits
%                                         cellfun(@(a) arrayfun(@(a) line(s,[a,a], [yPos-1 yPos], 'Color', [0 0 0]), a), rasterUnit(rasterMap(t)));
%                                         yPos = yPos + 1;
%                                     end
%                                 end
%                             end
            %%
            if(stimCond || contains(cond, 'Rest'))
            ylim([0 yPos])
            xlim([time_before_event{e} time_after_event{e}]);
            set(s,'YTick',[0 yPos]);
            set(s,'YTickLabel',[0 round(yPos)]);
            set(s,'XTick',[time_before_event{e}, 0, time_after_event{e}-buffer_time,time_after_event{e}]);
            title(num2str(plotMap(r)));
            line(s,[0 0], [0, yPos], 'Color', colorCond, 'LineWidth', 1, 'LineStyle', '--');
            line(s,[time_after_event{e}-buffer_time time_after_event{e}-buffer_time ], [0, yPos], 'Color', colorCond, 'LineWidth', 1, 'LineStyle', '--');
            end
            end
        disp(cond);
        %close all;
    end
end