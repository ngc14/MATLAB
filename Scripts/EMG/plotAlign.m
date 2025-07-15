function alignedEMGSigs = plotAlign(trials, segs,  averageSegs, win, segInds,gap,Fs, plotColors)
plotStart = 0;
xTicks = [];
labels = [];
alignedEMGSigs = [];
for align = 1:length(segInds)
    if(align==1)
        plotStart = 0;
    else
        plotStart = plotStart + gap +win{align-1}(2) + abs(win{align}(1));
    end
    alginSegInd = segInds(align);
    alignedEMGSigs = [alignedEMGSigs, cellfun(@(a,b) a(b(alginSegInd)+win{align}(1):...
        b(alginSegInd)+win{align}(2)),trials, segs,'UniformOutput', false)'];
    
    cellfun(@(a,b) plot(plotStart+win{align}(1):plotStart+win{align}(2),...
        a, 'Color',b), alignedEMGSigs(:,align), plotColors);
    plot(plotStart+win{align}(1):plotStart+win{align}(2), ...
        mean(reshape([alignedEMGSigs{:,align}],size(alignedEMGSigs{1,align},2),size(trials,2))'),...
        'k', 'LineWidth', 2);
    beforeSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:-1:1)));
    afterSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:end)));
    avgSegs = [beforeSegs plotStart afterSegs];
    for s = 1:length(avgSegs)
        if(avgSegs(s)>=plotStart+win{align}(1) && avgSegs(s)<=plotStart+win{align}(2))
            plot([avgSegs(s) avgSegs(s)],[0 max([alignedEMGSigs{:,align}])],'r--')
        end
    end
    xTicks = [xTicks, plotStart+[win{align}(1) 0 win{align}(2)]];
    labels = [labels, [win{align}(1) 0 win{align}(2)]];
    if(align==1)
        plotLims(1) = avgSegs(alginSegInd)-.25*Fs;
    elseif(align==length(segInds))
        plotLims(2) = avgSegs(alginSegInd)+.5*Fs;
    end
end
end

