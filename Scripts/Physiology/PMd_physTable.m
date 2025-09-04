pVal = 0.05;
varNames = ["Unit" "SiteNum" "Monkey" "Somatotopy" "Channel" "X" "Y" "Condition"...
     "TaskUnits" "DiffRest"];
repNames = ["Arm", "Hand", "Trunk"];
rNames = ["RT_r", "RSpeed_r"];
conditions = params.condNames;
taskAlign = params.PSTHAlignments;
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
allTaskInds = any(cell2mat(horzcat(tUnit{:})),2);
mappedChannels = cellfun(@(ch,l) ch{2}(l(~isnan(l))), chMaps,siteChannels, 'Uniformoutput', false)';
typeUnits = vertcat(restUnits{:});
avgSeg = cellfun(@(ct) cellfun(@(ca) cellfun(@(t) mean(t,1,'omitnan'), ca, 'UniformOutput',false),...
    ct, 'UniformOutput',false),siteSegs, 'UniformOutput',false);
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
[avgBaseline,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,avgSeg,normPSTH,false,false);
avgBase = cellfun(@cell2mat,cellfun(@(c) cellfun(@(a) median(cell2mat(a{1}),2,'omitnan'),...
    c, 'UniformOutput', false), avgBaseline, 'UniformOutput',false),'UniformOutput',false);
avgPhase = cellfun(@(c) cellfun(@(a) median(cell2mat(reshape(cellfun(@cell2mat,a(1),'UniformOutput',false),1,1,[])),3,'omitnan'),...
    c, 'UniformOutput', false), avgPhase, 'UniformOutput',false);
taskUnits = cell2mat(cellfun(@(a,b) cell2mat(a) & repmat(sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2),1,size(b,2)), ...
    num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false));
taskUnits(:,end+1) = any(taskUnits,2);
rawSp = cellfun(@(c) cellfun(@(a) a{1}, c,'UniformOutput', false),rawSpikes,'UniformOutput',false);
RTs = cellfun(@(c) cellfun(@(s) s{1}(:,2), c,'UniformOutput', false), siteSegs, 'UniformOutput',false);
Rspeeds = cellfun(@(c) cellfun(@(s) 1./(s{1}(:,3)), c,'UniformOutput', false), siteSegs, 'UniformOutput',false);
RTSpikes = cellfun(@(c,cr) cellfun(@(r,tr) num2cell(cell2mat(cellfun(@(a,t) sum(a<t), r,num2cell(repmat(tr',size(r,1),1)),...
    'UniformOutput',false)),2), c,cr,'UniformOutput',false), rawSp, RTs, 'UniformOutput',false);
RspeedSpikes = cellfun(@(c,cr,crs) cellfun(@(r,tr,trs) num2cell(cell2mat(cellfun(@(a,t,s) sum(a>t & a<s), r,...
    num2cell(repmat(tr',size(r,1),1)), num2cell(repmat(trs',size(r,1),1)),...
    'UniformOutput',false)),2), c,cr,crs,'UniformOutput',false), rawSp, RTs,Rspeeds, 'UniformOutput',false);
RTr = cellfun(@(cS,cT) cellfun(@(ss,t) cellfun(@(s) corr(s',t,'rows', 'complete'), ss,'UniformOutput',false),cS,cT, 'UniformOutput',false),...
    RTSpikes,RTs,'UniformOutput',false);
RSpeedr =  cellfun(@(cS,cT,cR) cellfun(@(ss,tE,tS) cellfun(@(s) corr(s',tE-tS,'rows', 'complete'), ss,'UniformOutput',false),cS,cT,cR, 'UniformOutput',false),...
    RspeedSpikes,RTs,Rspeeds,'UniformOutput',false);
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
    tUnits = taskUnits(:,c);
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
    condTable.unitType = typeUnits;
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
cl = validatecolor(["#A2142F","#0072BD","#EDB120","#77AC30"],'multiple');
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
    saveFigures(gcf,savePath+"r-Plots\",rNames(r),[]);
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