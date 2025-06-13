date = '06_03_2019';
monkey = 'Gilligan';
binSize=.01; % bin size in seconds
alignSegsConds = {["StartReach"],["StartReach"],["StartReach"]};
alignLims = {[-.5, 2.5]};
colors = [[224,144,38]./255; 0 1 0; 0 0 1; 1 0 1];
sigma = 5; % smoothing window in bin sizes
gapWind = 0.10;
MIN_BLOCKS_FOR_UNIT=15;
saveFolder = ['S:\Lab\ngc14\Working\EMG_UNITS\Stim_Triggered'];
conds = ["Extra Small Sphere", "Large Sphere", "Photocell"];
allUnits = cellfun(@(c) cell(1,length(alignLims)),conds,'UniformOutput',false);
sessionDir = dir(['S:\Lab\',monkey,'\All Data\',monkey,'_',date,'\Physiology\Results']);
sessionDir = sessionDir(~[sessionDir.isdir]);
[~,sortInd] = natsort({sessionDir.name});
sessionDir = sessionDir(sortInd);
avgSegs = cell(1,length(conds));
set(0, 'DefaultFigureRenderer', 'painters');
for f = 1:length(sessionDir)
    m = matfile_m([sessionDir(f).folder, '\', sessionDir(f).name]);
    m = m.sortedSpikeData;
    if(size(m.SegTimes,1)>0)
        load([sessionDir(f).folder,'\',sessionDir(f).name]);
        allGoodTrials = ~(cellfun(@(a,b) length(a) <=(b(end)-b(1)) ...
            | length(a)>200*(b(end)-b(1)) | any(isnan(a)), sortedSpikeData.SpikeTimes, sortedSpikeData.SegTimes));
        blockInds = cumsum(mod(1:length(sortedSpikeData.ArduinoData),length(conds))==1);
        unitTrials = num2cell(allGoodTrials.*blockInds,2);
        goodUnitsOnChannel = find(cellfun(@(tc) length(unique(tc))-1, unitTrials)>MIN_BLOCKS_FOR_UNIT);
        for u = 1:length(goodUnitsOnChannel)
            for c = 1:length(conds)
                alignSegs =  alignSegsConds{c};
                currTrialInds = find(cellfun(@(a) contains(a, conds{c}),sortedSpikeData.ArduinoData(:,1)));
                currTrialInds = currTrialInds(cellfun(@(a) ~any(isnan(a)), sortedSpikeData.SpikeTimes(goodUnitsOnChannel(u),currTrialInds)));
                currSpikes = sortedSpikeData.SpikeTimes(goodUnitsOnChannel(u),currTrialInds);
                currSegs = sortedSpikeData.SegTimes(goodUnitsOnChannel(u),currTrialInds);
                if(isempty(avgSegs{c}))
                    avgSegs{c} = currSegs;
                end
                for a = 1:length(alignSegs)
                    alignInd = find(strcmp(sortedSpikeData.ConditionSegments{c}, alignSegs(a)));
                    alignTimes = cellfun(@(s) getAlignedTimes(s,alignInd), currSegs, 'UniformOutput', false);
                    alignedSpikes = cellfun(@minus, currSpikes, alignTimes, 'UniformOutput', false);
                    trialHists = cellfun(@(s) histcounts(s,(alignLims{a}(1)):binSize:(alignLims{a}(end)))...
                        ./binSize, alignedSpikes,'UniformOutput', false);
                    trialHistsSmooth = cellfun(@(a) conv(a,gausswin(sigma)/sum(gausswin(sigma)),'same'), trialHists, 'UniformOutput', false);
                    allUnits{c}{a}{end+1} = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
                end
            end
        end
    end
end
%%
close all; set(0, 'DefaultFigureRenderer', 'painters');
figure();
normVals = max(cell2mat(cellfun(@(cp) cellfun(@(p) cellfun(@(n) max(mean(n,1,'omitnan')),p), cp(~isempty(cp)),...
    'UniformOutput',false), allUnits)'))';
maxFR=0; 
PSTHY = 1;
for c = 1:length(conds)
    alignSegs =  alignSegsConds{c};
    averageSegs = mean(cell2mat(transpose(avgSegs{c})),1,'omitnan');
    ax = subplot(3,2,(2*c)-1);
    hold on;
    ax2 = subplot(3,2,2*c);
    hold on;
    set(ax2,'ColorOrder',distinguishable_colors(length(normVals)));
    plotStart=0;
    for a = 1:length(alignSegs)
        secondsBeforeEvent = alignLims{a}(1);
        secondsAfterEvent = alignLims{a}(end);
        totalHists = cell2mat(cellfun(@(p) mean(p,1,'omitnan'),allUnits{c}{a},'UniformOutput',false)')./normVals;
        plot(ax2,plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-binSize),smoothdata(totalHists,2,'gaussian',7), ...
            'LineWidth', 1.5, 'LineStyle',':');
        plot(ax,plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-binSize),mean(totalHists,1,'omitnan'),...
            'Color',colors(c,:), 'LineWidth', 2);
        uE=mean(totalHists,1)+(std(totalHists,1)/sqrt(length(normVals)));
        lE=mean(totalHists,1)-(std(totalHists,1)/sqrt(length(normVals)));
        yP=[lE,fliplr(uE)];
        xP=[plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-...
            binSize),plotStart+fliplr(secondsBeforeEvent:binSize:secondsAfterEvent-binSize)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        d = patch(ax,xP,yP,colors(c,:));
        set(d,'facecolor',colors(c,:),'edgecolor','none','facealpha',.5)
        maxFR = max(maxFR,max(yP));
        alignInd = find(strcmp(sortedSpikeData.ConditionSegments{c}, alignSegs(a)));
        averageSegs = averageSegs - averageSegs(alignInd);
        beforeSegs = plotStart + cumsum(diff(averageSegs(alignInd:-1:1)));
        afterSegs = plotStart + cumsum(diff(averageSegs(alignInd:end)));
        aSegs = [beforeSegs plotStart afterSegs];
        plot(ax, [plotStart plotStart],[0 PSTHY],'Color', [0 0 0], 'LineStyle', '-');
        plot(ax2, [plotStart plotStart],[0 PSTHY],'Color', [0 0 0], 'LineStyle', '-');
        plotSegInds = (aSegs-plotStart<=secondsBeforeEvent | ...
            aSegs-plotStart>=secondsAfterEvent);
        plotSegInds(alignInd) = true;
        aSegs(plotSegInds) = [];
        for s = 1:length(aSegs)
            plot(ax, [aSegs(s) aSegs(s)],[0 PSTHY],'Color', [.15 .15 .15], 'LineStyle', '--');
            plot(ax2, [aSegs(s) aSegs(s)],[0 PSTHY],'Color', [.15 .15 .15], 'LineStyle', '--');
        end
        plotStart = plotStart + abs(secondsBeforeEvent) + gapWind + secondsAfterEvent;
    end
end


function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end