function fanoPlots(saveFig)
figure('Units', 'normalized', 'Position', [0 0 1 1]);
%%
colorPlots = [0 1 1;1 0 1;1 1 0;1 0 0;0 1 0];
axesHandles = {};
for a = 1:3
    if(a==1)
        activityCond = activityInd;
        activityType = 'Active';
    elseif(a==2)
        activityCond = ~activityInd;
        activityType = 'Inactive';
    else
        activityCond = true(size(activityInd));
        activityType = 'All';
    end
    for t = 1:length(masterTypes)
        plotNum = ((a-1) *length(masterTypes)) + t;
        disp(plotNum)
        axesHandles{plotNum} = subplot(3,length(masterTypes),plotNum);
        hold on;
        title([unitRow{2}, ' ', activityType]);
        for j = 1:length(jointName)
            jointInds = cellfun(@(a) strcmp(a, jointName{j}), masterJoints);
            allCondInd = activityCond & jointInds & masterTypes{t};
            condPSTF = cell2mat(masterFano(allCondInd)');
            plot(newBins,  nanmean(condPSTF,1),'Color', colorPlots(j,:), 'LineWidth', 1.5);
        end
        avgSegs = nanmean(reshape(cell2mat(masterSegs(jointInds & activityCond)'),...
            sum(jointInds & activityCond),length(events{condParamInd})),1);
        saveSegs{plotNum} = avgSegs;
        
    end
end
axesLinking = num2cell(reshape(1:length(axesHandles),3,length(masterTypes)),2);
yaxes = cellfun(@(d) [arrayfun(@(e) get(axesHandles{e},'YAxis'),d, 'UniformOutput', false)],axesLinking ,'UniformOutput', false);
xaxes = cellfun(@(d) [arrayfun(@(e) get(axesHandles{e},'XAxis'),d, 'UniformOutput', false)],axesLinking ,'UniformOutput', false);
cellfun(@(g) linkprop([cellfun(@(a) a(1), g)],'Limits'),xaxes, 'UniformOutput', false);
cellfun(@(g) linkprop([cellfun(@(a) a(1), g)],'Limits'),yaxes, 'UniformOutput', false);
arrayfun(@(a) set(axesHandles{a}, 'XLim', [-1 2.5]), cellfun(@(x) x(1), axesLinking));
maxLinked = cellfun(@(a) repmat(max(arrayfun(@(m) max(get(axesHandles{m},...
    'YLim')),a)),length(a),1), axesLinking, 'UniformOutput', false);
cellfun(@(x,s,m) arrayfun(@(t) plot(x, [t t], [0 m], 'k--'),s), axesHandles, ...
    saveSegs,num2cell(cell2mat(maxLinked))');
legend(axesHandles{1}, jointName, 'Location', 'best');

%
if(saveFigs)
    saveas(gcf,['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHs\',singleOrAll,...
        '\FRs\Quantification\Representation\',sessionOrUnitPSTHS,'_Discrete_XCorr'],'fig');
    saveas(gcf,['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHs\',singleOrAll,...
        '\FRs\Quantification\Representation\',sessionOrUnitPSTHS,'_Discrete_XCorr'],'png');
end