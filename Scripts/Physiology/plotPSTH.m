function f = plotPSTH(allUnits,alignSegsConds,avgSegs,conditionSegments,alignLims,binSize,colors,gapWind)
PSTHY=0;
for c = 1:length(alignSegsConds)
    alignSegs =  alignSegsConds{c};
    ax{c} = nexttile();
    hold on;
    colororder(ax{c},colors{c});
    plotStart=0;
    for a = 1:length(alignSegs)
        condSegs = cellfun(@(i) i(1):binSize:(i(end)-binSize), (transpose(avgSegs{c})),'UniformOutput',false);
        secondsBeforeEvent =  mean(cellfun(@(i) i(1), condSegs),'omitnan');
        secondsAfterEvent = mean(cellfun(@(i) i(end),condSegs),'omitnan');
        nsteps = max(cellfun(@length,condSegs))+length(alignLims{a}(1):binSize:alignLims{a}(end)-binSize);
        averageSegs = mean(cell2mat(cellfun(@(bn) [bn(1),bn(end)],condSegs, 'UniformOutput',false)),1,'omitnan');
        totalHists = cell2mat(cellfun(@(p) mean(p,1,'omitnan'),allUnits{c}{a},'UniformOutput',false)')';
        if(~isempty(totalHists))
            plot(ax{c},plotStart+(linspace(secondsBeforeEvent,secondsAfterEvent,nsteps)),totalHists',...
                'LineWidth', 2);
            uE=cellfun(@(p) mean(p,1,'omitnan')+(std(p,1)/sqrt(size(p,1))),allUnits{c}{a},'UniformOutput',false)';
            lE=cellfun(@(p) mean(p,1,'omitnan')-(std(p,1)/sqrt(size(p,1))),allUnits{c}{a},'UniformOutput',false)';
            yP=cellfun(@(l,u) [l,fliplr(u)],lE,uE,'UniformOutput',false);
            xP=[plotStart+(linspace(secondsBeforeEvent,secondsAfterEvent,nsteps)),...
                fliplr(plotStart+(linspace(secondsBeforeEvent,secondsAfterEvent,nsteps)))];
            ypInd = any(cell2mat(cellfun(@isnan,yP,'UniformOutput',false)),1);
            xP(ypInd)=[];
            yP = cell2mat(cellfun(@(y) y(~ypInd),yP,'UniformOutput',false));
            d = patch(ax{c},repmat(xP,size(yP,1),1)',yP',reshape(colors{c},size(colors{c},1),1,size(colors{c},2)));
            set(d,'edgecolor','none','facealpha',.5)
            PSTHY = max(PSTHY,max(yP,[],'all'));
        end
        alignInd = find(strcmp(conditionSegments{c}, alignSegs(a)));
        averageSegs = averageSegs - averageSegs(alignInd);
        beforeSegs = plotStart + cumsum(diff(averageSegs(alignInd:-1:1)));
        afterSegs = plotStart + cumsum(diff(averageSegs(alignInd:end)));
        aSegs = [beforeSegs plotStart afterSegs];
        plot(ax{c}, [plotStart plotStart],[0 PSTHY],'Color', [0 0 0], 'LineStyle', '-');
        plotSegInds = (aSegs-plotStart<=secondsBeforeEvent | ...
            aSegs-plotStart>=secondsAfterEvent);
        plotSegInds(alignInd) = true;
        aSegs(plotSegInds) = [];
        for s = 1:length(aSegs)
            plot(ax{c}, [aSegs(s) aSegs(s)],[0 PSTHY],'Color', [.15 .15 .15], 'LineStyle', '--');
        end
        plotStart = plotStart + abs(secondsBeforeEvent) + gapWind + secondsAfterEvent;
    end
    set(ax{c},'XLim',[secondsBeforeEvent,secondsAfterEvent]);
end
cellfun(@(x) set(x,'YLim',[0 PSTHY]),ax);
f=gcf();