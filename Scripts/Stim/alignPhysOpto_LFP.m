clc
clear all
close all

monkey = 'Gilligan';
sessionDate = '01_14_2020';
run = '0004';
eventChannels = {'SMA 1', 'SMA 2', 'SMA 3', 'SMA 4'};
numColumns = 4;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds

chanMap = [1:2:32,2:2:32];

buffer_time = 1;
low_freq = 30;
high_freq = 60;
step = 0.02;
sectionSize = .35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
eventTs = cell(1,length(eventChannels));
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filePath = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',sessionDate,'\Physiology\Gilligan_',sessionDate,'_',run,'.ns2'];
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
    eventTimeStamps(eventData==0) = [];
    timeDiff = diff(eventTimeStamps);
    fallingEdges(timeDiff<2) = [];
    timeDiff = [15 timeDiff];
    eventTimeStamps(timeDiff<2) = [];
    
    
    time_before_event{l} = -1;    % seconds before event
    time_after_event{l}  = round((fallingEdges(1)-eventTimeStamps(1))+buffer_time,3,'significant');     % seconds before event
    eventTs{l} = eventTimeStamps;
    
    disp('Loading Event Data... 100% Done.');
end
ns_CloseFile(filePath);
%%
for l = 1:length(chanMap)
    [LFP{l}, LFPTimes{l},Fs,~,~] = loadRippleData(filePath,chanMap(l));
end

%%
numRows = ceil(length(chanMap)/numColumns);
fRows = (high_freq-low_freq)*10+1;
fSteps = linspace(low_freq,high_freq,fRows);
for e = 1:length(eventChannels)
    alignedLFP = cell(size(LFPTimes));
    stim_ts_cond = eventTs{e};
    for s = 1:length(stim_ts_cond)
        centeredLFP = cellfun(@(a) a-stim_ts_cond(s), LFPTimes, 'UniformOutput', false);
        alignedLFP(s,:) = cellfun(@(a,b) b(a>=time_before_event{e} & a<=time_after_event{e}), centeredLFP,LFP, 'UniformOutput', false);
    end
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    plotMap = reshape(1:length(chanMap),[numRows,numColumns])';
    plotMap = reshape(plotMap,[1,length(chanMap)]);
    hold on;
    
    for r = 1:length(chanMap)
        s = subplot(numRows,numColumns,r);

        %LFP
        channelSig = alignedLFP(:,plotMap(r));
        cLength = mode(cellfun(@length, channelSig));
        channelSig = cellfun(@(a) a(1:cLength), channelSig, 'UniformOutput', false);
        
        for t=1:length(channelSig)
            start = 1;
            for k=1:step*Fs:length(channelSig{t})-sectionSize*Fs
                section = channelSig{t}(k:k+sectionSize*Fs);
                [pxx, ~] = pwelch(section,length(section), 0, fSteps, Fs);
                
                PS(:,start,t)=pxx;
                start = start + 1;
            end
        end
        plotPower = imgaussfilt(10*log10(median(PS,3)),3);
        imagesc(time_before_event{e}:1/Fs:time_after_event{e},[low_freq high_freq], plotPower)
        colormap('jet')
        set(gca,'YDir', 'normal');
        title(num2str(plotMap(r)));
        line(s,[0 0], [low_freq, high_freq], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
        line(s,[time_after_event{e}-buffer_time time_after_event{e}-buffer_time ], [low_freq, high_freq], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
        
        
        %         plot(s,time_before_event{e}:(1/Fs):time_after_event{e}-(1/Fs), ...
        %             mean(reshape(cell2mat(channelSig)', length(channelSig{1}),[]),2)','LineWidth', 2, 'Color', 'k');
        %         yPos = [min(mean(reshape(cell2mat(channelSig)', length(channelSig{1}),[]),2)), max(mean(reshape(cell2mat(channelSig)', length(channelSig{1}),[]),2))];
        %
        %
        %
        %         %%
        %         ylim([yPos(1) yPos(2)])
        %         xlim([time_before_event{e} time_after_event{e}]);
        %         set(s,'YTick',[yPos(1) yPos(2)]);
        %         set(s,'YTickLabel',[round(yPos(1)) round(yPos(2))]);
        %         set(s,'XTick',[time_before_event{e}, 0, time_after_event{e}-buffer_time,time_after_event{e}]);
        %         title(num2str(plotMap(r)));
        %         line(s,[0 0], [yPos(1), yPos(2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
        %         line(s,[time_after_event{e}-buffer_time time_after_event{e}-buffer_time ], [yPos(1), yPos(2)], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
    disp(r);
    minP(r) = min(min(min(plotPower)));
    maxP(r) = max(max(max(plotPower)));
    end
    for a = 1:length(chanMap)
        s = subplot(numRows,numColumns,a);
        set(gca, 'CLim',[min(minP) max(maxP)]);
    end
    close all;
end