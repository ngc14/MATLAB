clc
clear all
close all

monkey = 'Gilligan';
date = '04_18_2019';
run = '0001';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds

time_before_event = -0.15;    % seconds before event
time_after_event  = .5;     % seconds before event
%chanMap = [129:2:160, 130:2:160];
chanMap = 1:32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['\\pitt\sni\Gharbawie\Lab\' monkey, '\All Data\', monkey,'_', date,'\Physiology\', monkey, '_', date, '_',run];
NEV = openNEV([filename, '-sorted.nev']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodeIDs = unique([NEV.Data.Spikes.Electrode]);
for l = 1:length(electrodeIDs)
    elec_chan = electrodeIDs(l);
    elec_spikes_inds = find(NEV.Data.Spikes.Electrode==elec_chan);
    elec_units = NEV.Data.Spikes.Unit(elec_spikes_inds);
    elec_units = elec_units(elec_units~=255);
    for u = 1:max(elec_units)
        spikeTimes{l,u} = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Unit(elec_spikes_inds)==u));
    end
end
alignedTrial = cell(size(spikeTimes));
%%
risingEdges = NEV.Data.SerialDigitalIO.TimeStampSec(1:2:end);
fallingEdges = NEV.Data.SerialDigitalIO.TimeStampSec(2:2:end);
pulseLength = fallingEdges-risingEdges;
startEventIdx=find(pulseLength>0.019 & pulseLength<0.021);
endEventIdx = find(pulseLength>0.039 & pulseLength<0.061);
startEventIdx(152) = [];
%%
for s = 1:length(segTimes)
    centeredSpikes = cellfun(@(a) a-segTimes{s}(2), spikeTimes, 'UniformOutput', false);
    centeredSpikes = cellfun(@(a) a(a>=time_before_event & a<=time_after_event), centeredSpikes, 'UniformOutput', false);
    [row, col] = find(~cellfun(@isempty, centeredSpikes));
    for i = 1:length(row)
        alignedTrial{row(i), col(i)}{end+1,:} = centeredSpikes{row(i), col(i)};
    end
end

%%
yPos = 0;
for r = 1:size(alignedTrial,1)
    line([time_before_event,time_after_event],[yPos, yPos], 'Color',[.85 0 0], 'LineWidth', 1);
    yPos = yPos +1;
    textPos =  time_before_event - 0.15*diff([time_before_event,time_after_event]);
    %text(textPos,yPos,num2str(33-r));
    for d = 1:1
        [~,c] = max(cellfun(@(a) max(length(a)), alignedTrial(chanMap(33-r),:)));
        rasterUnit = alignedTrial{chanMap(33-r),c};
        if(~isempty(rasterUnit))
            %                 p = gca(figure());
            %                 trialHists = cellfun(@(a) histcounts(a,...
            %                     time_before_event-length(Kernel)/(2*1/bin):bin:...
            %                     time_after_event+length(Kernel)/(2*1/bin))./bin, rasterUnit,'UniformOutput', false);
            %                 trialHistsSmooth = cellfun(@(a) conv(a,Kernel), trialHists, 'UniformOutput', false);
            %                 trialHistsSmooth = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
            %                 totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
            %                 plot(p,time_before_event:bin:time_after_event,nanmean(totalHists,1), 'LineWidth', 2);
            for t = 1:length(rasterUnit)
                cellfun(@(a) arrayfun(@(a) line([a,a], [yPos-1 yPos+1], 'Color', [0 0 0], 'LineWidth', 1), a), rasterUnit(t));
                yPos = yPos + 2;
            end
            %line([time_before_event,time_after_event],[yPos, yPos], 'Color','m', 'LineWidth', 2);
        end
    end
    ylim([0 yPos]);
    xlim([time_before_event time_after_event]);
    line([reactT reactT], [0, yPos], 'Color', 'c', 'LineWidth', 2);
    line([reachT reachT], [0, yPos], 'Color', 'y', 'LineWidth', 2);
    line([graspT graspT], [0, yPos], 'Color', 'g', 'LineWidth', 2);
    line([liftT liftT], [0, yPos], 'Color', 'b', 'LineWidth', 2);
 end