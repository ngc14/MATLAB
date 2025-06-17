function f = plotPSTH(allUnits,alignSegsConds,avgSegs,conditionSegments,alignLims,binSize,colors,gapWind)
PSTHY=0;
for c = 1:length(alignSegsConds)
    alignSegs =  alignSegsConds{c};
    averageSegs = mean(cell2mat(transpose(avgSegs{c})),1,'omitnan');
    ax = nexttile();
    hold on;
    colororder(ax,colors{c});
    plotStart=0;
    for a = 1:length(alignSegs)
        secondsBeforeEvent = alignLims{a}(1);
        secondsAfterEvent = alignLims{a}(end);
        totalHists = cell2mat(cellfun(@(p) mean(p,1,'omitnan'),allUnits{c}{a},'UniformOutput',false)');
        if(~isempty(totalHists))
            plot(ax,plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-binSize),totalHists,...
                'LineWidth', 2);
            uE=cellfun(@(p) mean(p,1,'omitnan')+(std(p,1)/sqrt(size(p,1))),allUnits{c}{a},'UniformOutput',false)';
            lE=cellfun(@(p) mean(p,1,'omitnan')-(std(p,1)/sqrt(size(p,1))),allUnits{c}{a},'UniformOutput',false)';
            yP=cellfun(@(l,u) [l,fliplr(u)],lE,uE,'UniformOutput',false);
            xP=[plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-...
                binSize),plotStart+fliplr(secondsBeforeEvent:binSize:secondsAfterEvent-binSize)];
            ypInd = any(cell2mat(cellfun(@isnan,yP,'UniformOutput',false)),1);
            xP(ypInd)=[];
            yP = cell2mat(cellfun(@(y) y(~ypInd),yP,'UniformOutput',false));
            d = patch(ax,repmat(xP,size(yP,1),1)',yP',reshape(colors{c},size(colors{c},1),1,size(colors{c},2)));
            set(d,'edgecolor','none','facealpha',.5)
            PSTHY = max(PSTHY,max(totalHists,[],'all'));
        end
        alignInd = find(strcmp(conditionSegments{c}, alignSegs(a)));
        averageSegs = averageSegs - averageSegs(alignInd);
        beforeSegs = plotStart + cumsum(diff(averageSegs(alignInd:-1:1)));
        afterSegs = plotStart + cumsum(diff(averageSegs(alignInd:end)));
        aSegs = [beforeSegs plotStart afterSegs];
        plot(ax, [plotStart plotStart],[0 PSTHY],'Color', [0 0 0], 'LineStyle', '-');
        plotSegInds = (aSegs-plotStart<=secondsBeforeEvent | ...
            aSegs-plotStart>=secondsAfterEvent);
        plotSegInds(alignInd) = true;
        aSegs(plotSegInds) = [];
        for s = 1:length(aSegs)
            plot(ax, [aSegs(s) aSegs(s)],[0 PSTHY],'Color', [.15 .15 .15], 'LineStyle', '--');
        end
        plotStart = plotStart + abs(secondsBeforeEvent) + gapWind + secondsAfterEvent;
    end
    set(ax,'XLim',[secondsBeforeEvent,secondsAfterEvent]);
end
f=gcf();