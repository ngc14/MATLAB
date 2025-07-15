function plotGroupedBars(condXrepXphase,saveIn,saveOn)
%%
saveDir = regexp(saveIn,'\','split');
saveName = saveDir(end);
saveDir = strjoin(saveDir(1:end-1),"\")+"\";
jNames = ["Arm", "Hand"];
jColor = flipud(cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)),jNames, 'UniformOutput', false)'));
condColor = flipud([.7 0 0; .8 .4 0; 0 0 .7;]);
%%
figure('Units','normalized','Position',[0 0 1 1]);
hold on
id = [];
condVals = {};
phaseGap = 1.75;
groupGap = 5;
for n = 1:length(condXrepXphase)
    condXrepPlot = cat(2,condXrepXphase{n});
    subInd = linspace(-.9,.9,length(condXrepPlot)+2);
    subInd = (subInd(2:end-1));
    condVals{n} = linspace((phaseGap*(n-1))*size(condXrepPlot{1},2),...
        (phaseGap*n*size(condXrepPlot{1},2))-1,size(condXrepPlot{1},2));
    condVals{n} = condVals{n}+(groupGap*(n-1));
    for r = 1:length(condXrepPlot)
        condXphasePlot = condXrepPlot{r};
        %         tb = notBoxPlot(condXphasePlot,condVals{n} + subInd(r),'style','patch','jitter', .5);
        %         arrayfun(@(a) set(a.data,'Marker', 'none'), tb, 'UniformOutput', false);
        pOffset = condVals{n}+subInd(r);
        %condXphasePlot = cell2mat(cellfun(@(m) mean(m,2,'omitnan'), condXphasePlot, 'UniformOutput',false));
        for p = 1:size(condXphasePlot,2)-1
            tb(p) = boxchart(condXphasePlot(:,p),'XData',repmat(pOffset(p),length(condXphasePlot(:,p)),1),...
                'Notch','on','BoxWidth',.5,'MarkerStyle','none');
        end
%         errorbar(condVals{n}+subInd(r),median(condXphasePlot,1,'omitnan'),std(condXphasePlot,0,1,'omitnan'),...
%             '.','Color',jColor(r,:));
        id(end+1:end+size(tb,2),:) = cell2mat(arrayfun(@(a) [n,r,a], 1:size(tb,2),'UniformOutput',false)');
    end
end
bx = (findall(gca,'Type','BoxChart'));
mx = findall(gca,'Tag','Median');
arrayfun(@(m) set(m,'Color',[0 0 0]), mx);
arrayfun(@(ur) set(bx(id(:,2)==ur),  'BoxEdgeColor',flipud(jColor(ur,:))), unique(id(:,2)));
arrayfun(@(uc) set(bx(id(:,1) ==uc),'BoxFaceColor',flipud(condColor(uc,:)), 'BoxFaceAlpha', 1),unique(id(:,1)));
%xlim([min(cellfun(@min,condVals))-1,max(cellfun(@max,condVals))+1]);
ylim([0 8])

ax = get(findall(gca,'Type', 'Axes'),'XAxis');
ax.Categories = categorical(min(cellfun(@min,condVals)):max(cellfun(@max,condVals)));
ax.TickValues = categorical(arrayfun(@num2str,cellfun(@(n) round(median(n)),condVals)-1, 'Uniformoutput', false));
xticklabels(["Precision", "Power", "Reach-only"]);
%title(jNames(n))
if(saveOn)
    saveFigures(gcf, saveDir,saveName,[]);
end
end