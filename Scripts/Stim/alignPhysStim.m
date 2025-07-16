clc
clear all
close all

monkey = 'Gilligan';
date = '05_24_2019';
run = '0001';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds

time_before_event = -.05;    % seconds before event
time_after_event  = .05;     % seconds before event
chanMap = [129:2:160, 130:2:160];

sigma = 5; % smoothing window in ms

Kernel = (-3*sigma:3*sigma);
kernelBin = length(Kernel);
Kernel = (-kernelBin/2:kernelBin/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['\\pitt\sni\Gharbawie\Lab\' monkey, '\All Data\', monkey,'_', date,'\Physiology\', monkey, '_', date, '_stim_',run];
[n,names] = plx_chan_names([filename, '-sorted.plx']);
fileID = fopen([filename, '.txt']);
fileInfo = textscan(fileID, '%f%f%f%f%f%f%f%f','Delimiter', '\n', 'headerlines', 1);
fclose(fileID);
trialInfo = [fileInfo{1}, fileInfo{3}, fileInfo{4}, fileInfo{5}, fileInfo{6}, fileInfo{7}, fileInfo{8}];
condCombs = unique(trialInfo, 'rows');
for l = 1:length(chanMap)
    name = sprintf('elec%d', chanMap(l));
    elec_chan = find(contains(string(names), name));
    for u = 1:4
        [~,spikeTimes{l,u}] = plx_ts([filename, '-sorted.plx'],elec_chan,u);
    end
end
alignedTrial = cell(size(spikeTimes));
%%
stimElectrodes = unique(trialInfo(:,1));
stim_ts = cell(size(stimElectrodes));
for e = 1:length(stimElectrodes)
    stimE = stimElectrodes(e);
    name = sprintf('Nano2_stim.e%02.f',stimE);
    stim_chan = find(contains(string(names), name));
    [n, stim_ts_raw] = plx_ts([filename, '-sorted.plx'], stim_chan, 1);
    currTrialsInd = trialInfo(:,1)==stimE;
    pulses = trialInfo(currTrialsInd,2);
    tInd = 1;
    stim_ts{e}(end+1,:) = stim_ts_raw(1);
    for t = 1:length(pulses)-1
        stim_ts{e}(end+1,:) = stim_ts_raw(tInd+pulses(t));
        tInd = tInd + pulses(t);
    end
end
for e = 1:size(condCombs,1)
    cond = condCombs(e,:);
    name = sprintf('Nano2_stim.e%02.f',cond(1));
    stim_chan = find(contains(string(names), name));
    currTrials = trialInfo(trialInfo(:,1)==cond(1),:);
    currTrialInds = currTrials(:,2)==cond(2) &  currTrials(:,3)==cond(3) & ...
        currTrials(:,4)==cond(4) & currTrials(:,5)==cond(5) & currTrials(:,6)==cond(6) & ...
        currTrials(:,7)==cond(7);
    stim_ts_cond = stim_ts{find(cond(1)==stimElectrodes)}(currTrialInds);
    for s = 1:length(stim_ts_cond)
        centeredSpikes = cellfun(@(a) a-stim_ts_cond(s), spikeTimes, 'UniformOutput', false);
        centeredSpikes = cellfun(@(a) a(a>=time_before_event & a<=time_after_event), centeredSpikes, 'UniformOutput', false);
        [row, col] = find(~cellfun(@isempty, centeredSpikes));
        for i = 1:length(row)
            alignedTrial{row(i), col(i)}{end+1,:} = centeredSpikes{row(i), col(i)};
        end
    end
    s = subplot(1,size(condCombs,1),e);
    yPos = 0;
    for r = 1:size(alignedTrial,1)
        line(s,[time_before_event,time_after_event],[yPos, yPos], 'Color','r', 'LineWidth', 2);
        yPos = yPos +1;
        axisPos = get(s, 'Position');
        textPos =  time_before_event - 0.15*diff([time_before_event,time_after_event]);
        text(s,textPos,yPos,num2str(r));
        for c = 1:size(alignedTrial,2)
            rasterUnit = alignedTrial{r,c};
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
                    cellfun(@(a) arrayfun(@(a) line(s,[a,a], [yPos-1 yPos], 'Color', [0 0 0]), a), rasterUnit(t));
                    yPos = yPos + 1;
                end
                line(s,[time_before_event,time_after_event],[yPos, yPos], 'Color','m', 'LineWidth', 2);
            end
        end
    end
    ylim(s,[0 yPos]);
    set(s,'YTickLabel',[]);
    line(s,[0 0], [0, yPos], 'Color', 'b', 'LineWidth', 1);
    title(s,cond);
end