clc
clear all
close all

monkey = 'Gilligan';
sessionDate = '06_22_2020';
run = '0002';
eventChannels = {'SMA 1','SMA 2','SMA 3','SMA 4'};
numColumns = 4;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds

chanMap = [32:-1:1];
%chanMap = [18,15,17,16,22,11,21,12,31,2,29,9,32,1,20,10,30,4,19,7,27,3,25,8,28,13,23,5,26,6,24,14];

sigma = .2; % smoothing window in ms

Kernel = (-3*sigma:3*sigma);
kernelBin = length(Kernel);
Kernel = (-kernelBin/2:kernelBin/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

buffer_time = 1;
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
%     eventTimeStamps(1:2) = [];
%     eventData(1:2) = [];
    
    fallingEdges = eventTimeStamps(eventData==0);
    eventTimeStamps(eventData==0) = [];
    
    % get rid of pulsing edges
    timeDiff = diff(eventTimeStamps);
    fallingEdges([Inf,timeDiff]<2) = [];
    eventTimeStamps([Inf,timeDiff]<2) = [];
    
    %repeating conds
%     if(l==length(eventChannels))
%         timeDiff = diff(eventTimeStamps);
%         timeDiff = [Inf, timeDiff];
%         time_before_event{end+2} = -1;
%         time_before_event{end+1} = -1;
%         time_after_event{end+2} = 1+buffer_time;
%         time_after_event{end+1} = 1+buffer_time;
%         eventTs{end+1} = eventTimeStamps(timeDiff>7.5 & timeDiff<8.5);
%         nextInds = find(timeDiff>7.5 & timeDiff<8.5)+1;
%         if(nextInds(end)>length(eventTimeStamps))
%             nextInds = nextInds(1:end-1);
%         end
%         eventTs{end+1} = eventTimeStamps(nextInds);
%     end
%     
%     % condition repeat takes 8 seconds
%     if(strcmp(eventChannels{l}, 'SMA 3'))
%         fallingEdges([Inf, diff(fallingEdges)]<6) = [];
%         eventTimeStamps([Inf, diff(eventTimeStamps)]<6) = [];
%     else
%         nextInds = find(timeDiff>7.5 & timeDiff<8.5)-1;
%         if(nextInds(1)<1)
%             nextInds = nextInds(2:end);
%         end
%         fallingEdges = fallingEdges(nextInds);
%         eventTimeStamps = eventTimeStamps(nextInds);
%     end
    
    
    time_before_event{l} = -1;    % seconds before event
    time_after_event{l}  = round((fallingEdges(1)-eventTimeStamps(1))+buffer_time,3,'significant');     % seconds before event
    eventTs{l} = eventTimeStamps;
    
    disp('Loading Event Data... 100% Done.');
end
ns_CloseFile(filePath);
% if(strcmp(run,'0002'))
%     eventTs = [ {eventTs{1}-3}, eventTs];
%     time_before_event = [{-1}, time_before_event];
%     time_after_event = [{2}, time_after_event];
% end
%%
for l = 1:length(chanMap)
    [dataTime, ids] = loadSpikeData(filePath,chanMap(l));
    sortedIDs = unique(ids(ids<255));
    for u = 1:length(sortedIDs)
        dataUnits = dataTime(ids==sortedIDs(u));
        spikeTimes{l,u} = dataUnits;
        spikeTimes{l,1} = dataTime;
    end
end

%%
numRows = ceil(length(chanMap)/numColumns);
for e = 1:length(eventChannels)
    alignedTrial = cell(size(spikeTimes));
    stim_ts_cond = eventTs{e};
    for s = 1:length(stim_ts_cond)
        centeredSpikes = cellfun(@(a) a-stim_ts_cond(s), spikeTimes, 'UniformOutput', false);
        centeredSpikes = cellfun(@(a) a(a>=time_before_event{e}-length(Kernel)/(2*1/bin) & a<=time_after_event{e}+length(Kernel)/(2*1/bin)), centeredSpikes, 'UniformOutput', false);
        [row, col] = find(~cellfun(@isempty, centeredSpikes));
        for i = 1:length(row)
            alignedTrial{row(i), col(i)}{end+1,:} = centeredSpikes{row(i), col(i)};
        end
    end
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    plotMap = reshape(1:length(chanMap),[numRows,numColumns])';
    plotMap = reshape(plotMap,[1,length(chanMap)]);
    hold on;
    for r = 1:length(alignedTrial)
        s = subplot(numRows,numColumns,r);
        hold on;
        yPos = 0;
        line(s,[time_before_event{e},time_after_event{e}],[yPos, yPos], 'Color','r', 'LineWidth', 2);
        
        %%PSTH
                rasterUnit = alignedTrial{plotMap(r),1};
                trialHists = cellfun(@(a) histcounts(a,...
                    time_before_event{e}-length(Kernel)/(2*1/bin):bin:...
                    time_after_event{e}+length(Kernel)/(2*1/bin)), rasterUnit,'UniformOutput', false);
                trialHistsSmooth = cellfun(@(a) conv(a,Kernel), trialHists, 'UniformOutput', false);
                trialHistsSmooth = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
                totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
                plot(s,time_before_event{e}:bin:time_after_event{e},nanmean(totalHists,1)./bin, 'LineWidth', 2, 'Color', 'k');
                yPos = max(nanmean(totalHists,1))/bin;
        
        %%RASTER
%         for c = 1:size(alignedTrial,2)
%             rasterUnit = alignedTrial{plotMap(r),c};
%             if(~isempty(rasterUnit))
%                 if(length(rasterUnit)>60)
%                     rasterMap = 1:round(length(rasterUnit)/60):length(rasterUnit);
%                     totalUnits = length(rasterMap);
%                     
%                     for t = 1:totalUnits
%                         cellfun(@(a) arrayfun(@(a) line(s,[a,a], [yPos-1 yPos], 'Color', [0 0 0]), a), rasterUnit(rasterMap(t)));
%                         yPos = yPos + 1;
%                     end
%                 else
%                     totalUnits = length(rasterUnit);
%                     
%                     for t = 1:totalUnits
%                         cellfun(@(a) arrayfun(@(a) line(s,[a,a], [yPos-1 yPos], 'Color', [0 0 0]), a), rasterUnit(t));
%                         yPos = yPos + 1;
%                     end
%                 end
% 
%             end
%         end
        %%
        ylim([0 yPos])
        xlim([time_before_event{e} time_after_event{e}]);
        set(s,'YTick',[0 yPos]);
        set(s,'YTickLabel',[0 round(yPos)]);
        set(s,'XTick',[time_before_event{e}, 0, time_after_event{e}-buffer_time,time_after_event{e}]);
        title(num2str(plotMap(r)));
        line(s,[0 0], [0, yPos], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
        line(s,[time_after_event{e}-buffer_time time_after_event{e}-buffer_time ], [0, yPos], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
    end
    close all;
end
fclose('all');