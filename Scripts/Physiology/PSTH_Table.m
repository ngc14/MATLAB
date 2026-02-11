conditions = ["Extra Small Sphere", "Large Sphere", "Photocell"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]}});
taskWindow =repmat({{[0.2, 0]}},1,length(conditions));
params = PhysRecording(string(conditions),.01,.15,-6,5,containers.Map(conditions,...
    {"StartReach","StartReach","StartReach"}));
phaseNames = categorical([ "Go", "Reach", "Hold","Withdraw"],'Ordinal',true);
phaseAlignmentPoints = {["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"]};
phaseWinSz = .2;
phaseWindows = repmat({{[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz*(5/4), -phaseWinSz*(1/4)],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions));
phaseWindows{end}{3} = [-.1 0];
pVal=0.05;
savePath = "S:\Lab\ngc14\Working\PSTHS\";
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","M1");
clear rawSpikes
%%
allCondCue = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,2)-4,'omitnan'),params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t) cellfun(@(a,n) max(1,mean(a(:,n(find(~isnan(n),1)):n(find(~isnan(n),1))+(3/params.binSize),:),2,'omitnan')),...
    p,t,'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),num2cell([allCondCue{:}],2),"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) cellfun(@(p,b)p./b,cp,nb,'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),...
    normBaseline,'Uniformoutput', false);
normPSTH = num2cell(vertcat(normPSTH{:}),1);
%%
sumSegs = cellfun(@(c) cellfun(@(n) [n{:}], c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
for c = 1:length(conditions)
    nanSegs = find(sum(isnan(cell2mat(sumSegs{c})),1)>size(cell2mat(sumSegs{c}),1)/2 & [true(1,length(maxSegL)-1) 0]);
    for a = 1:length(nanSegs)
        for n = 1:length(sumSegs{c})
            nextSeg = intersect(setdiff(1:size(sumSegs{c}{n},2),nanSegs),nanSegs(a)+1:size(sumSegs{c}{n},2));
            sumSegs{c}{n}(:,nanSegs(a)) = sumSegs{c}{n}(:,nextSeg(1));
        end
    end
end
sumSegs = cellfun(@(c)cellfun(@(s) {s},c,'UniformOutput',false),sumSegs,'UniformOutput',false);
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,sumSegs,siteTrialPSTHS,false,true);
[~,phaseFR] = calculatePhases(params,condPhaseAlign,phaseWindows,sumSegs,siteTrialPSTHS,false,true);
[~,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,sumSegs,normPSTH,false,false);
rAUC = cellfun(@(c) cellfun(@(a) a{1}{:,contains(string(phaseNames),"Reach")},c,'UniformOutput',false), phaseFR, 'UniformOutput',false);
gAUC = cellfun(@(c) cellfun(@(a) a{1}{:,contains(string(phaseNames),"Hold")},c,'UniformOutput',false), phaseFR, 'UniformOutput',false);
[gsubr,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials({{r}},{{g}},1,true,0.05),cr,cg,'UniformOutput',false),rAUC,gAUC, 'UniformOutput',false);
gsubr = cellfun(@(c) cellfun(@(v) mean([v{:}],2,'omitnan'), c, 'UniformOutput',false),gsubr,'UniformOutput',false);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
typeUnits = cellfun(@(rg,v,t) sum(cumprod([cell2mat(rg)>0 & cell2mat(v)<0, cell2mat(rg)>0 & cell2mat(v)>0, ...
    cell2mat(rg)==0 & cell2mat(t)==1] ==0,2),2)+1, rgInds,gsubr,tUnit, 'UniformOutput', false);
for t = 1:length(typeUnits)
    typeUnits{t}(typeUnits{t}>3) = 0;
end
condXphase = cellfun(@(pc) cell2mat(cellfun(@(v) permute(mean(cat(3,v{:}),2,'omitnan'),[1 3 2]),...
    cellfun(@(n) vertcat(n{:}),pc,'UniformOutput',false),'UniformOutput',false)),avgPhase,'UniformOutput',false);
condXphase = cellfun(@(s) [s,NaN(size(s,1),length(phaseNames)-size(s,2))], condXphase, 'UniformOutput',false);
chUnitMap = cellfun(@(m) m{end},chMaps','UniformOutput',false);
% chUnitMap(strcmp([siteDateMap.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};chUnitMap(strcmp([siteDateMap.Monkey],"Skipper")) = {[32:-1:1]};
%%
tPhys = table();
g = @(x,y,c)GetPointLineDistance(x,y,c(1),c(2),c(3),c(4));
for c = 1:length(conditions)
    condTable = table();
    AUCVals = condXphase{c};
    tUnits = cell2mat(tUnit{c});
    tUnits(isnan(tUnits)) = 0;
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels)';    
    allReps =  mapSites2Units(condUnitMapping,siteRep');
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    mInds= contains(mLabs,"Skipper");
    unitLocation = cell2mat(arrayfun(@(x,y,c) [g(x,y,OrthogonalLines(c).RCLine), g(x,y,OrthogonalLines(c).MLLine)],...
        siteDateMap.x,siteDateMap.y,siteDateMap.Monkey,'UniformOutput',false));
    unitLocation = mapSites2Units(condUnitMapping,num2cell([unitLocation(:,1),unitLocation(:,2)],2));
    condTable.Unit = [1:length(mLabs)]';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), [1:length(condUnitMapping)]',condUnitMapping,'UniformOutput',false));
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(allReps);
    condTable.Channel =  cell2mat(cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap,siteChannels, 'Uniformoutput', false))';
    condTable.X(mInds) = cellfun(@(x) x(1), unitLocation(mInds)) - min(cellfun(@(x) x(1), unitLocation(mInds)));
    condTable.X(~mInds) = cellfun(@(x) x(1), unitLocation(~mInds)) - min(cellfun(@(x) x(1), unitLocation(~mInds)));
    condTable.Y(mInds) =  cellfun(@(x) x(end), unitLocation(mInds)) - min(cellfun(@(x) x(end), unitLocation(mInds)));
    condTable.Y(~mInds) = cellfun(@(x) x(end), unitLocation(~mInds)) - min(cellfun(@(x) x(end), unitLocation(~mInds)));
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    condTable.TaskUnits =  logical(tUnits);
    condTable.TaskFR = cell2mat(cellfun(@(t) mean(t{1}{1},2,'omitnan'), taskFR{c},'UniformOutput',false));
    for pn = 1:length(phaseNames)
        condTable.(string(phaseNames(pn))) = AUCVals(:,pn);
    end
    condTable.rSI = condTable.Reach./(condTable.Reach+condTable.Hold);
    condTable.gSI = condTable.Hold./(condTable.Reach+condTable.Hold);
    condTable.rgSI = (condTable.Reach - condTable.Hold)./(condTable.Reach+condTable.Hold);
    condTable.unitType = typeUnits{c};
    tPhys = [tPhys;condTable];
end
plotNames = arrayfun(@(p) arrayfun(@(c) string(p)+"_"+c{1}(1), conditions, 'UniformOutput', true), phaseNames, 'UniformOutput', false);
plotNames = [plotNames{:}];
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
%%
varNames = tPhys.Properties.VariableNames(any(cell2mat(arrayfun(@(p) contains(tPhys.Properties.VariableNames,string(p)) & ~contains(tPhys.Properties.VariableNames,"_R"), phaseNames, 'UniformOutput', false)'),1));
phaseVals = tPhys{:,varNames};
condLabels = sum(cell2mat(cellfun(@(c,n) n.*contains(varNames,"_"+c),params.condAbbrev.values,num2cell(1:length(conditions)),'UniformOutput',false)'),1);
phaseLabels = sum(cell2mat(cellfun(@(c,n) n.*contains(varNames,string(c)),cellstr(phaseNames),num2cell(1:length(phaseNames)),'UniformOutput',false)'),1);
condIndex = any(tPhys{:,contains(tPhys.Properties.VariableNames,"Task") & ~contains(tPhys.Properties.VariableNames,"_R")},2) ...
    %& tPhys.Channel<=16;
condPhases = phaseVals(condIndex,:);
armPhases = condPhases(tPhys{condIndex,'Somatotopy'}=="Arm",:);
handPhases = condPhases(tPhys{condIndex,'Somatotopy'}=="Hand",:);
cl = [0 .5 0; .5 0 .5;0 0 0];
figure(); hold on;
s=swarmchart(reshape([repmat(condLabels+[-.3 -.3 -.3 -.1 -.1 -.1 .1 .1 .1 .3 .3 .3],size(armPhases,1),1);repmat(condLabels...
    +[-.3 -.3 -.3 -.1 -.1 -.1 .1 .1 .1 .3 .3 .3]+.05,size(handPhases,1),1)],[],1),reshape([armPhases;handPhases],[],1),...
    [],cl(reshape([ones(size(armPhases));2*ones(size(handPhases))],1,[]),:),'filled');
s.XJitter = 'rand'; s.XJitterWidth = .02;
[meanVals,gNames] = cellfun(@(c) groupsummary(cell2mat(arrayfun(@(v) tPhys{condIndex,string(v)},(varNames(contains(varNames,"_"+c))),'UniformOutput',false)),...
    {tPhys{condIndex,"Somatotopy"},any(tPhys{condIndex,contains(tPhys.Properties.VariableNames,"TaskUnit")},2)},"median"),params.condAbbrev.values,'UniformOutput',false);
%bx=violin(condPhases,'x',length(phaseNames)*(condLabels-1+repmat([-.3 -.1 .1 .3],1,length(conditions))),'facecolor',cl(phaseLabels,:));
ylim([0 8]);
xa = gca();
set(xa,'XTick',[1:1:length(conditions)]);
set(xa,'XTickLabel',params.condAbbrev.values);
gNames = cellfun(@(g) string(g{1})+"_"+string(g{2}), gNames, 'UniformOutput', false);
plotVals = cell2mat(cellfun(@(g,m) reshape(m((contains(g,"Arm") | contains(g,"Hand")) & contains(g,"true"),:),1,[]),gNames,meanVals,'UniformOutput',false));
scatter(unique([get(xa,"Children").XData],'sorted'),plotVals,200,'black',"_",'LineWidth',3);
saveFigures(gcf,savePath,"Normalized_Violins",[]);
%
allSegs = cellfun(@(c,n) cellfun(@(s,t) repmat(mean([s{:}],1,'omitnan'),size(t,1),1), c,n, 'UniformOutput',false), sumSegs,normPSTH,'UniformOutput',false);
psthLabs = arrayfun(@(c,s) repmat(c,size(cell2mat([s{:}]),1),1)+"_"+string(tPhys.Somatotopy)+"_"+tPhys.("unitType_"+c),string(params.condAbbrev.values),allSegs,'UniformOutput',false);%+"_"+string(tPhys.Channel>16)
allPSTHS = cellfun(@(p) num2cell(p,[2,3]),vertcat(normPSTH{:}),'UniformOutput',false);
for p = 1:length(psthLabs)
    psthLabs{p}(ismember(psthLabs{p},emptyConds{p}))="";
end
allLabs = unique([psthLabs{:}]); twoDim = length(conditions) * double(max(unique(arrayfun(@(r) length(regexp(r,"_",'match')),allLabs)))>1);
figure();plotJointPSTHS(params,{cell2mat(cellfun(@(m) mean(m,3,'omitnan'),vertcat(allPSTHS{:}),'UniformOutput',false)).*...
    double((0./(vertcat(tPhys.TaskUnits_ESS,tPhys.TaskUnits_LS,tPhys.TaskUnits_P)))+1)},{cell2mat(vertcat(allSegs{:}))},...
    vertcat(psthLabs{:}),(vertcat(tPhys.TaskUnits_ESS,tPhys.TaskUnits_LS,tPhys.TaskUnits_P)),[],[-.5 2.5],[.9 6],...
    cell2struct(reshape(repmat(repmat(num2cell(cl,2),1,length(conditions)),max(1,twoDim),1)',[],1),allLabs(allLabs~="")));
saveFigures(gcf,savePath,"Normalized_PSTHS",[]);