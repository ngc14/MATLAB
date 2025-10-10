varNames = ["Unit" "SiteNum" "Monkey" "Somatotopy" "Channel" "X" "Y" "Condition"...
     "TaskUnits" "DiffRest"];
rNames = ["RT_r", "RSpeed_r"];
conditions = params.condNames;
phaseNames = ["Go", "Reach", "Hold", "Withdraw","Reward"];
phaseAlignmentPoints = {["GoSignal","StartReach","StartHold","StartWithdraw","StartReward"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw","StartReward"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw","StartReward"],...
    ["GoSignal","StartReplaceHold","StartReward"]};
phaseWinSz = .2;
phaseWindows = repmat({{[-phaseWinSz*(3/4),phaseWinSz*(1/4)],[-phaseWinSz*(1/2),phaseWinSz*(1/2)],...
    [-phaseWinSz, 0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],[-phaseWinSz*(1/2),phaseWinSz*(1/2)]}},...
    1,length(conditions)-1);
phaseWindows(end+1) = {{[-phaseWinSz*(3/4),phaseWinSz*(1/4)],[-phaseWinSz*(1/4),phaseWinSz*(3/4)],...
    [-phaseWinSz*(1/4),phaseWinSz*(3/4)]}};
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\";
close all;
%%
mappedChannels = siteChannels;%cellfun(@(ch,l) ch(l(~isnan(l))), chMaps,siteChannels, 'Uniformoutput', false)';
avgSeg = cellfun(@(ct) cellfun(@(ca) cellfun(@(t) mean(t,1,'omitnan'), ca, 'UniformOutput',false),...
    ct, 'UniformOutput',false),siteSegs, 'UniformOutput',false);
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
sumSegs = cellfun(@(c) cellfun(@(n) NaN(size(n{1},1),length(maxSegL)), c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
for j = 1:size(sumSegs,2)
    for i = 1:length(mappedChannels)
        sumSegs{j}{i}(:,condSegMappedInds{j}) = siteSegs{j}{i}{1};
    end
end
[~,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,avgSeg,normPSTH,false,false);
avgPhase = cellfun(@(c) cellfun(@(a) median(cell2mat(reshape(cellfun(@cell2mat,a(1),'UniformOutput',false),1,1,[])),3,'omitnan'),...
    c, 'UniformOutput', false), avgPhase, 'UniformOutput',false);
taskUnits = cellfun(@(a,b) cell2mat(a) & repmat(sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2),1,size(b,2)), ...
    num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false);
%%
condPSTHS = cellfun(@(c) cellfun(@(cp) cell2mat(cp),c,'UniformOutput',false),normPSTH,'UniformOutput',false);
trialCondInfo = arrayfun(@(c) cellfun(@(s) s(strcmp(s(:,1),c),:), siteTrialInfo, 'UniformOutput',false)',conditions,'UniformOutput',false);
allPSTHS = cellfun(@(c,t) cellfun(@(r,n,i) permute(permute((any(i,2)./any(i,2)).*r,...
    [3 2 1]).*~isnan(cellfun(@str2double,n(:,end-1))),[3 2 1]),c,t,taskUnits,'UniformOutput',false),condPSTHS,trialCondInfo,'UniformOutput',false);
tBounds = cellfun(@(c) cellfun(@(s) s(:,1:3),c,'UniformOutput',false), sumSegs, 'UniformOutput',false);

RTMean = cellfun(@(c,cr) cellfun(@(r,tr) cellfun(@(s,p) mean(p(:,params.bins>s(:,1) & params.bins<s(:,2)),2,'omitnan'),...
    num2cell(tr,2),squeeze(num2cell(r,[1,2])),'UniformOutput',false),c,cr,'UniformOutput',false), allPSTHS, tBounds, 'UniformOutput',false);
RSpeedMean = cellfun(@(c,cr) cellfun(@(r,tr)cellfun(@(s,p) mean(p(:,params.bins>s(:,2) & params.bins<s(:,3)),2,'omitnan'),...
    num2cell(tr,2),squeeze(num2cell(r,[1,2])),'UniformOutput',false),c,cr,'UniformOutput',false), allPSTHS,tBounds,'UniformOutput',false);
RTr = cellfun(@(cS,cT) cellfun(@(ss,tE) cellfun(@(s) corr(s',diff(tE(:,1:end-1),1,2)./1000,'rows', 'complete'),...
    num2cell(cell2mat(ss'),2),'UniformOutput',false),cS,cT, 'UniformOutput',false),RTMean,tBounds,'UniformOutput',false);
RSpeedr =  cellfun(@(cS,cT) cellfun(@(ss,tE) cellfun(@(s) corr(s',diff(tE(:,2:end),1,2)./1000,'rows', 'complete'),...
    num2cell(cell2mat(ss'),2),'UniformOutput',false),cS,cT,'UniformOutput',false),RSpeedMean,tBounds,'UniformOutput',false);
Rs = {RTr,RSpeedr};
%%
tPhys = [];
condXphase = cellfun(@(pc) cell2mat(pc),cellfun(@(e) e(repmat(~isempty(e),size(e))),avgPhase,'UniformOutput', false),'UniformOutput',false);
condXphase = cellfun(@(s) [s,NaN(size(s,1),length(phaseNames)-size(s,2))], condXphase, 'UniformOutput',false);
%%
[~,sortCols]= cellfun(@(pa) sort([arrayfun(@(p) find(contains(phaseNames,p)), phaseNames(arrayfun(@(p) contains(char([pa{:}]),p),phaseNames))),...
setdiff(1:length(phaseNames),arrayfun(@(p) find(contains(phaseNames,p)), phaseNames(arrayfun(@(p) contains(char([pa{:}]),p),phaseNames))))]),phaseAlignmentPoints,'UniformOutput',false);
condXphase = cellfun(@(s,sc) s(:,sc), condXphase,sortCols,'UniformOutput',false);
for c = 1:length(conditions)
    condTable = table();
    AUCVals = condXphase{c};
    tUnits = cell2mat(cellfun(@(a) a(:,c), taskUnits, 'UniformOutput',false));
    tUnits(isnan(tUnits)) = 0;
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels)';    
    allReps =  mapSites2Units(condUnitMapping,siteRep');
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    mInds= contains(mLabs,"Skipper");
    unitLocation = mapSites2Units(condUnitMapping,num2cell([siteDateMap.x,siteDateMap.y],2));
    condTable.Unit = [1:length(mLabs)]';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), [1:length(condUnitMapping)]',condUnitMapping,'UniformOutput',false));
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(allReps);
    condTable.Channel =  [mappedChannels{:}]';
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    condTable.TaskUnits =  logical(tUnits);
    condTable.unitType = NaN(size(tUnits));
    for pn = 1:length(phaseNames)
        condTable.(phaseNames(pn)) = AUCVals(:,pn);
    end
    for r = 1:length(rNames)
        condTable.(rNames(r)) = cell2mat(vertcat(Rs{r}{c}{:}));
    end
    condTable.Properties.VariableNames = [varNames, phaseNames, rNames];
    tPhys = [tPhys;condTable];
end
plotNames = arrayfun(@(p) arrayfun(@(c) p+"_"+c{1}(1), conditions, 'UniformOutput', true), phaseNames, 'UniformOutput', false);
plotNames = [plotNames{:}];
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(...
    strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
%writetable(tPhys,savePath+'PMdTable.xlsx');
%%
plotGroupedBars(cellfun(@(n) num2cell(n,1),condXphase,'UniformOutput',false),savePath+"Units_Phys_FR_Box",false);
%%
cl = flipud(validatecolor(["#A2142F","#EDB120","#0072BD"],'multiple'));
lm = {};
for r = 1:length(rNames)
    figure(); hold on;
    title(rNames(r));
    for c = 1:length(conditions)-1
        xVals = (tPhys.X-min(tPhys.X)).*ImagingParameters.px2mm;
        yVals = tPhys.(strcat(rNames(r),"_",params.condAbbrev(conditions(c))));
        scatter(xVals,yVals,'MarkerFaceColor','flat','CData',cl(c,:),'Marker','o');
        lm{c} = fitlm(xVals,yVals);
        plot(0:1:max(xVals)+1,arrayfun(@(x) table2array(lm{c}.Coefficients(1,1)) + ...
            x*table2array(lm{c}.Coefficients(2,1)),0:1:max(xVals)+1),'Color',cl(c,:),'LineWidth',2);
        cc = coefCI(lm{c});
        plot(0:1:max(xVals)+1,arrayfun(@(x) cc(1,1) + x*cc(2,1),0:1:max(xVals)+1), 'Color',cl(c,:),'LineWidth',2,'LineStyle','--');
        plot(0:1:max(xVals)+1,arrayfun(@(x) cc(1,2) + x*cc(2,2),0:1:max(xVals)+1), 'Color',cl(c,:),'LineWidth',2,'LineStyle','--');
    end
    ylabel("r-value");
    xlabel("caudal to rostral (mm)");
    g = gca;
    g = g.Children;
    g = g(arrayfun(@(l) strcmp(l.Type,'line'),g));
    legend(g(arrayfun(@(l) strcmp(l.LineStyle,'-'),g)), cellfun(@(c,m) ...
        string(params.condAbbrev(c)+": R^2="+num2str(m.Rsquared.Ordinary,'%.4f')+", p="+ num2str(coefTest(m),'%.2f')), ...
        cellstr(conditions(1:end-1)),lm));
         ylim([-1 1]);
    saveFigures(gcf,savePath+"r-Plots\Max\",rNames(r),[]);
    f = gca(figure());
    boxchart(repmat(round(xVals),length(conditions)-1,1),cell2mat(arrayfun(@(c) ...
        tPhys.(strcat(rNames(r),"_",params.condAbbrev(c))),conditions(1:end-1),'Uniformoutput',false)'),...
        'GroupByColor',[ones(length(xVals),1);2*ones(length(xVals),1);3*ones(length(xVals),1)],...
        'Notch','on');
    hold on;
    plot([-1:1:max(round(xVals))+1], zeros(1,length(-1:1:max(round(xVals))+1)),'Color','k','LineStyle','--');
    f.ColorOrder = flipud(cl);
    legend(arrayfun(@(c) string(params.condAbbrev(c)), conditions(1:end-1)));
     ylabel("r-value");
     ylim([-1 1]);
    xlabel("caudal to rostral (mm)");
    xticks(0:1:max(round(xVals)));
    title(rNames(r));
    saveFigures(gcf,savePath+"r-Plots\Max\",rNames(r)+"_Box",[]);
    close all;
end
%%
close all;
allBounds = cellfun(@(cR) vertcat(cR{:}), num2cell(horzcat(tBounds{3}),2), 'UniformOutput', false);
rtTimes = cellfun(@(b)diff(b(:,1:end-1),1,2),allBounds, 'UniformOutput',false);
reachTimes =  cellfun(@(b) diff(b(:,2:end),1,2),allBounds, 'UniformOutput',false);
rtFRs = cellfun(@(v) cell2mat(cellfun(@transpose,vertcat(v{:}),'UniformOutput',false))',num2cell(horzcat(RTMean{3}),2),'UniformOutput',false);
reachFRs = cellfun(@(v) cell2mat(cellfun(@transpose,vertcat(v{:}),'UniformOutput',false))',num2cell(horzcat(RSpeedMean{3}),2),'UniformOutput',false);
[~,uI] = unique(tPhys.SiteNum);
siteColors = min(5,max(1,round(ImagingParameters.px2mm.*(tPhys.X(uI,:)-min(tPhys.X(uI,:))))));
cmap = distinguishable_colors(max(siteColors)+1);
for r = 1:2
    if(r==1)
        rPlotY = rtFRs;
        rPlotX = rtTimes;
    else
        rPlotY = reachFRs;
        rPlotX = reachTimes;
    end
    for c = 1:length(rPlotX)
        rYC = mean(rPlotY{c}(:,rPlotX{c}~=0).*1./~isinf(1./logical(rPlotY{c}~=0)),1,'omitnan');
        if(~all(isnan(rYC)))
            scatter(rPlotX{c}(rPlotX{c}~=0),rYC,25,cmap(siteColors(c),:),'filled'); hold on;
            lm = fitlm(rPlotX{c},rYC);
            cc = coefCI(lm);
            title([rNames(r)+" R^2="+num2str(lm.Rsquared.Ordinary,'%.4f')+", p="+num2str(coefTest(lm),'%.4f')]);
            ylabel("mean FR");
            xlabel("Duration (s)");
            ylim([0 65]);
            gc = gca();
            plot(gc.XTick,arrayfun(@(x) table2array(lm.Coefficients(1,1)) + ...
                x*table2array(lm.Coefficients(2,1)),gc.XTick),'LineWidth',2,'Color','k');
            plot(gc.XTick,arrayfun(@(x) cc(1,1) + x*cc(2,1),gc.XTick),'LineWidth',2,'LineStyle','--','Color','k');
            plot(gc.XTick,arrayfun(@(x) cc(1,2) + x*cc(2,2),gc.XTick),'LineWidth',2,'LineStyle','--','Color','k');
            xlim([min(rPlotX{c}(~isoutlier(rPlotX{c}))),max(rPlotX{c}(~isoutlier(rPlotX{c})))]);
            saveFigures(gcf,savePath+"r-Plots\Sessions\"+rNames(r)+"\","SiteNo"+num2str(siteDateMap{c,'Site'}),[]);
        end
        close all;
    end
    figure();hold on;
    colororder(cmap(1:end-1,:));
    s=cellfun(@(rx,ry,cl) scatter(rx(rx~=0),mean(ry(:,rx~=0).*1./~isinf(1./logical(ry~=0)),1,'omitnan'),...
        25,'filled','ColorVariable',cmap(cl,:)), rPlotX, rPlotY,num2cell(siteColors));
    allX = cell2mat(arrayfun(@(ss) ss.XData, s, 'UniformOutput',false)');
    allY = cell2mat(arrayfun(@(ss) ss.YData, s, 'UniformOutput',false)');
    lm = fitlm(allX,allY);
    cc = coefCI(lm);
    title([rNames(r)+" R^2="+num2str(lm.Rsquared.Ordinary,'%.4f')+", p="+num2str(coefTest(lm),'%.4f')]);
    ylabel("mean FR");
    xlabel("Duration (s)");
    ylim([0 65]);
    colormap(cmap(1:end-1,:));
    colorbar(gca,'eastoutside','Ticks',linspace(0,1,6),'TickLabels',[strcat('<',num2str(min(siteColors)));...
        string(arrayfun(@num2str, unique(siteColors(siteColors~=max(siteColors)))));strcat(num2str(max(siteColors)), '+')]);

    gc = gca();
    plot(gc.XTick,arrayfun(@(x) table2array(lm.Coefficients(1,1)) + ...
        x*table2array(lm.Coefficients(2,1)),gc.XTick),'LineWidth',2,'Color',cmap(end,:));
    plot(gc.XTick,arrayfun(@(x) cc(1,1) + x*cc(2,1),gc.XTick),'LineWidth',2,'LineStyle','--','Color',cmap(end,:));
    plot(gc.XTick,arrayfun(@(x) cc(1,2) + x*cc(2,2),gc.XTick),'LineWidth',2,'LineStyle','--','Color',cmap(end,:));
    xlim([min(allX(~isoutlier(allX))),max(allX(~isoutlier(allX)))]);
    saveFigures(gcf,savePath+"r-Plots\","Mean_"+rNames(r),[]);
end
%%
factorInd = find(contains(tPhys.Properties.VariableNames,phaseNames) & ~contains(tPhys.Properties.VariableNames,"_R"));
rModel = fitrm(tPhys(:,[find(strcmp(tPhys.Properties.VariableNames,"Somatotopy")),factorInd]),...
    char(join([join(tPhys.Properties.VariableNames(factorInd),','),...
    '~ 1 + Somatotopy'])),'WithinDesign',tPhys.Properties.VariableNames(:,factorInd));

stackVar = arrayfun(@(f) string(arrayfun(@(c) join([f,c],"_"),params.condAbbrev.values(...
    cellstr(conditions(~strcmp(conditions,"Rest")))))),["TaskUnits","DiffRest",phaseNames],'UniformOutput',false);
tPhysCases = stack(tPhys,stackVar,'NewDataVariableName',["TaskUnits","DiffRest",phaseNames],'IndexVariableName','Condition');
tPhysCases.Condition = categorical(params.condAbbrev.values(cellstr(conditions(tPhysCases.Condition-min(tPhysCases.Condition)+1)))');
tPhysCases = stack(tPhysCases,phaseNames,'NewDataVariableName','FR','IndexVariableName','Phase');
wd = cell2mat(arrayfun(@(c) strsplit(string(c)),reshape(unique(tPhysCases.Phase).*unique(tPhysCases.Condition)',[],1),...
    'UniformOutput',false));
wd(:,2) = strcat("_",wd(:,2));
rModel = fitrm(tPhys(:,[find(strcmp(tPhys.Properties.VariableNames,"Somatotopy")),factorInd]),...
    char(join([join(tPhys.Properties.VariableNames(factorInd),','),...
    '~ 1 + Somatotopy'])),'WithinDesign',array2table(wd,'VariableNames',{'Phase','Cond'}));
mc=multcompare(rModel,'Cond');
rn = ranova(rModel,'WithinModel','Phase*Cond');