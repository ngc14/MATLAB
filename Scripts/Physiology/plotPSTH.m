function plotPSTH(savePath,conditionName,PSTHS,alignTimes,channels,bins)
conditionMappings = {'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'};
segLabs = repmat({{'RT', 'R', 'G', 'L','H', 'W', 'RH', 'SR', 'R'}},2,1);
segLabs(end+1) = {{'RT','R', 'G', 'H', 'W','RH', 'SR', 'R'}};
segLabs(end+1) = {{'H','RH', 'SR', 'R'}};
PSTHDisplayLimits = [-1, 2.5];

conditionInd = cellfun(@(a) strcmp(a, conditionName), conditionMappings);
figure('Name',conditionName,'units','normalized','outerposition',[0 0 1 1]);

subSZ = numSubplots(size(PSTHS,3));
prevChannel = NaN;
segInds = cellfun(@length, alignTimes)==length(segLabs{conditionInd});
segBins =  nanmean(reshape(cell2mat(alignTimes(segInds)'),sum(segInds),length(segLabs{conditionInd})),1);
segLabsLocs = segBins + [diff(segBins)/2, 0];
segLabsLocs(end) = segLabsLocs(end) + (PSTHDisplayLimits(end) + segLabsLocs(end))/2;
segsToPlot = segLabsLocs > PSTHDisplayLimits(1) & segLabsLocs<PSTHDisplayLimits(end);
segLabsLocs = segLabsLocs(segsToPlot);
segBins = segBins(segsToPlot);
segLabPlot = segLabs{conditionInd}(segsToPlot);
for u = 1:size(PSTHS,3)
    subplot(subSZ(1), subSZ(2),u);
    hold on;
    currPSTH = PSTHS(:,:,u);
    plot(bins,nanmean(currPSTH), 'LineWidth', 2);
    uE=nanmean(currPSTH,1)+nanstd(currPSTH,1);
    lE=nanmean(currPSTH,1)-nanstd(currPSTH,1);
    yP=[lE,fliplr(uE)];
    xP=[bins,fliplr(bins)];
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    d = patch(xP,yP,1);
    set(d,'edgecolor','none','facealpha',.5);
    if(prevChannel==channels(u))
        unitCount = unitCount + 1;
        prevChannel = channels(u);
    else
        unitCount = 1;
    end
    title(['Channel ', num2str(channels(u)), ' Unit ', ...
        num2str(unitCount)]);
    xlim(PSTHDisplayLimits);
    currYLim = get(gca,'ylim');
    ylim([0 currYLim(end)]);
    if(u<=subSZ(2))
        cellfun(@(a,b) text(a, currYLim(end), b, 'Color', 'r', 'FontSize', 12,...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center'),...
            num2cell(segLabsLocs), segLabPlot);
    end
    arrayfun(@(a) line([a a], [0 currYLim(end)], 'Color', 'r', 'LineStyle','--'),...
        segBins);
end
saveas(gcf,[savePath,'\',conditionName,'_Units.fig']);
saveas(gcf,[savePath,'\',conditionName,'_Units.png']);
pause(0.01);
close gcf;
end