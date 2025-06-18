function figHandle = peakTimeHistograms(bins,utReps,sortValsNaN,segs,siteReps,saveDir,unitName)
figure('Units','normalized','Position',[0 0 1 1]); hold on;
rp = fieldnames(MotorMapping.repColors);
histNames = rp(ismember(string(rp),unique(utReps)));
histNames{end+1} = 'All';
ar = cell(1,length(histNames));
for r=1:length(histNames)
    ar{r} = subplot(length(histNames),1,r); hold on;
    if (r==length(histNames))
        jointInds = true(size(utReps));
        histogram(sortValsNaN(jointInds),bins,'FaceColor','k','EdgeColor','k');
        avgSegs = nanmedian(cell2mat(vertcat(segs(:))),1);
    else 
        jointInds = arrayfun(@(t) strcmp(t,histNames(r)), utReps);
        histogram(sortValsNaN(jointInds),bins,'FaceColor', MotorMapping.repColors.(histNames{r}),...
            'EdgeColor',MotorMapping.repColors.(histNames{r}));
        avgSegs = nanmedian(cell2mat(vertcat(segs(strcmp(utReps(siteReps),histNames{r})))),1);
    end
    title(strcat(unitName, ": n = ", num2str(nansum(~(isnan(sortValsNaN(jointInds))))),...
        " ", histNames{r}));
    scatter(nanmedian(sortValsNaN(jointInds)),0,'filled', 'CData', [1 0 0]);
    gscatter(quantile(sortValsNaN(jointInds),[0.25, 0.75])', [0 0]',[1,2]',...
        [1 0 0; 1 0 0],['<','>']',[5,5]');
    legend('off');
    currLim = ylim();
    for a=1:length(avgSegs)
        if(avgSegs(a)>=-1 & avgSegs(a)<=1)
            line([avgSegs(a),avgSegs(a)], currLim, 'Color', 'k', 'LineWidth', 2,'LineStyle', '--');
        end
    end
end
linkaxes([ar{:}]);
maxYlim = ylim();
cellfun(@(aa) set(findobj(aa.Parent, 'Type', 'Line','LineWidth', 2),'YData', maxYlim), ar, 'UniformOutput', false)
xlim([-1 1]);
saveFigures(gcf,saveDir,unitName,[]);
figHandle = ar;
end