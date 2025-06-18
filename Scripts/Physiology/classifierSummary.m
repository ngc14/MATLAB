mDirName = "S:\Lab\ngc14\Working\Both\Baseline_FR\Decoding\All\Bayes\";
load([mDirName+"accuracy.mat"]);
allSoma = somaUnits;
allSoma = cellfun(@(i) cell2mat(cellfun(@(d) reshape(mean(d,1,'omitnan'),...
    [ones(1,length(size(i))),size(d,2)]),i,'UniformOutput', false)), allSoma, 'UniformOutput',false);
allSoma = cat(5, allSoma{:});
%%
clr = {'r','b','k','g','m','k'};
symb = {'p', 'o', '*', '<'};
somatotopicLabsPlot = somatotopicLabs;
somatotopicLabsPlot = cellfun(@(s) s,somatotopicLabsPlot(1:end-2), 'UniformOutput',false);
fTypePlot = fTypes;
fTypePlot = arrayfun(@(s) s{1}(1)+"Units", fTypePlot);
subNames = transpose(somatotopicLabsPlot+fTypePlot);
subNames = arrayfun(@(p) (subNames +p), phaseNames(~strcmp(phaseNames, "Withdraw")), 'UniformOutput', false);
subNames = cat(3,subNames{:});
fGrouping = NaN(size(subNames));
for f = 1:length(fTypePlot)
    fGrouping(contains(subNames,fTypePlot(f))) = (f);
end
%%
close all;
for i = 1:size(subNames,3)
    xvals = NaN(size(subNames,1),size(subNames,2),1);
    for xa=1:length(somatotopicLabsPlot)
        xvals(:,xa) = repmat(xa,size(xvals(:,xa)));
    end
    % g = gscatter(categorical(reshape(somatotopicLabsPlot(xvals),1,[])),reshape(allThreshs(:,1:length(somatotopicLabsPlot),i),1,[]),...
    %     reshape(fGrouping(:,:,i),1,[]),[clr{1:length(fTypePlot)}],[],25,'off',[],'# of units to achieve 90%');
    %     {categorical(reshape(mGrouping(:,:,i,:),1,[])),categorical(max(mGrouping,[],'all')+reshape(fGrouping(:,:,i,:),1,[]))},...
    %     colorGroup,markerGroup,15,'off',[],'# of units to achieve 90%');

    f1 = figure(1);
    set(f1,'Units','normalized','Position',[0 0 1 1]);
    subplot(3,1,i);
    set(gca,'FontSize',16);
    hold on;
    title(phaseNames(i));
    threshCross = squeeze(allSoma(:,7,1:length(somatotopicLabsPlot),i,:));
    g = boxplot(reshape(threshCross,1,[]),{reshape(repmat(somatotopicLabsPlot(xvals),1,1,size(threshCross,3)),1,[]);...
        reshape(repmat(fTypePlot(fGrouping(:,:,i)),1,1,size(threshCross,3)),1,[])},'FactorGap',10,'Colors',[clr{1:length(fTypePlot)}],...
        'Notch','on','Symbol','','BoxStyle','filled','PlotStyle','traditional','Widths',1,'Labels',reshape(transpose(string(somatotopicLabsPlot)+fTypes),1,[]));
    xl = get(gca,'XTickLabel');
    xt = get(gca,'XTick');
    set(gca,'XTick', xt(floor(length(fTypes)/2):length(fTypes):end));
    set(gca,'XTickLabel',cellfun(@(s) s(1:end-4), xl(floor(length(fTypes)/2):length(fTypes):end),'Uniformoutput',false));
    % reshape(repmat(categorical(squeeze(mGrouping(:,:,i,:))).*categorical(max(mGrouping,[],'all')+squeeze(fGrouping(:,:,i,:))),1,1,1,size(threshCross,4)),1,[]),...
    % colorGroup,markerGroup,15,'off', [], 'Accuray evaluated at n = 10 units');
    if(i==size(subNames,3))
        p = [cellfun(@(s) scatter(NaN,NaN,55,s,'filled','o'), clr(1:length(fTypePlot)), 'UniformOutput',false)];
        legend([p{:}],[cellstr(fTypes)],'FontSize', 12,'location','southwest');
    else
        legend('off')
    end
    xticklabels(somatotopicLabsPlot);
    ylim([0.3 1]);
    set(gca,'view',[90 -90])
end
%%
saveFigures(gcf,mDirName, "Accuracy_Distribution",[]);