conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[0, 0]}},1,length(conditions));
params = PhysRecording(string(conditions),.01,.15,-6,5,containers.Map(conditions,...
    {"StartReach","StartReach","StartReach","GoSignal"}));
phaseNames = categorical([ "Go", "Reach", "Hold","Withdraw"],'Ordinal',true);
phaseAlignmentPoints = {["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReplaceHold","StartReward"]};
phaseWinSz = .2;
phaseWindows = repmat({{[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz, 0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions)-1);
phaseWindows(end+1) = {{[0, phaseWinSz], [-phaseWinSz,0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}};
pVal=0.05;
MIN_BLOCKS_FOR_UNIT = 13;
plotSessions = false;
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\";
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","M1");
clear rawSpikes
%%
allCondCue = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,2)+1,'omitnan'),params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t) mean(cell2mat(cellfun(@(a,n) max(1,mean(...
    a(:,n(~isnan(n)):n(~isnan(n))+(3/params.binSize),:),[2,3],'omitnan')),...
    p,t,'UniformOutput',false)),2,'omitnan'),num2cell([siteTrialPSTHS{:}],2),num2cell([allCondCue{:}],2),"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) cellfun(@(p)p./repmat(nb,1,1,size(p,3)),...
    cp,'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),normBaseline,'Uniformoutput', false);
%%
trialInfo = cellfun(@(c) cellfun(@(t) t(strcmp(t(:,1),c),:),siteTrialInfo,'UniformOutput',false)',conditions,'UniformOutput',false);
sumSegs = cellfun(@(c) cellfun(@(n) NaN(size(n,1),length(maxSegL)), c, 'UniformOutput',false), trialInfo,'UniformOutput',false);
trialInfo = cellfun(@(c) vertcat(c{:}),num2cell(horzcat(trialInfo{:}),2),'UniformOutput',false);
sumSegs = cellfun(@(c) cellfun(@(n) [n{:}], c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
for c = 1:length(conditions)-1
    nanSegs = find(sum(isnan(cell2mat(sumSegs{c}(~cellfun(@(a) all(isnan(a),'all'),sumSegs{c})))),1)>=...
        sum(cellfun(@(z) size(z,1),sumSegs{c}(~cellfun(@(a) all(isnan(a),'all'),sumSegs{c}))))/2);
    nanSegs= nanSegs(nanSegs~=length(maxSegL));
    for a = 1:length(nanSegs)
        for n = 1:length(sumSegs{c})
            nextSeg = intersect(setdiff(1:size(sumSegs{c}{n},2),nanSegs),nanSegs(a)+1:size(sumSegs{c}{n},2));
            sumSegs{c}{n}(:,nanSegs(a)) = sumSegs{c}{n}(:,nextSeg(1));
        end
    end
end
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,cellfun(@(c) cellfun(@(n) {n},c,'UniformOutput',false),...
    sumSegs,'Uniformoutput',false),siteTrialPSTHS,false,true);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),tb,tc,'UniformOutput',false),...
    taskBaseline,taskFR,'UniformOutput', false);
avgSeg = cellfun(@(ca) cellfun(@(t) mean(t,1,'omitnan'), ca, 'UniformOutput',false),sumSegs, 'UniformOutput',false);
[~,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,cellfun(@(c) cellfun(@(n) {n},c,'UniformOutput',false),...
    avgSeg,'UniformOutput',false),num2cell(vertcat(normPSTH{:}),1),false,true);
condXphase = cellfun(@(pc) cell2mat(cellfun(@(v) permute(mean(cat(3,v{:}),2,'omitnan'),[1 3 2]),cellfun(@(n) vertcat(n{:}),pc,'UniformOutput',false),'UniformOutput',false)),avgPhase,'UniformOutput',false);
condXphase = cellfun(@(s) [s,NaN(size(s,1),length(phaseNames)-size(s,2))], condXphase, 'UniformOutput',false);
%%
chUnitMap = cell(height(siteDateMap),1);
chUnitMap(strcmp([siteDateMap.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([siteDateMap.Monkey],"Skipper")) = {[32:-1:1]};
rAUC = cellfun(@(c) cellfun(@(a) a{1}{:,contains(string(phaseNames),"Reach")},c,'UniformOutput',false), avgPhase, 'UniformOutput',false);
gAUC = cellfun(@(c) cellfun(@(a) a{1}{:,contains(string(phaseNames),"Hold")},c,'UniformOutput',false), avgPhase, 'UniformOutput',false);
[gsubr,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials({{r}},{{g}},1,true,0.05),cr,cg,'UniformOutput',false),...
    rAUC,gAUC, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
gsubr = cellfun(@(c) cell2mat(cellfun(@(v) mean([v{:}],2,'omitnan'), c, 'UniformOutput',false)),gsubr,'UniformOutput',false);
rUnits = cellfun(@(rg,v) rg==1 & v<0, rgInds, gsubr, 'UniformOutput',false);
gUnits = cellfun(@(rg,v) rg==1 & v>0, rgInds, gsubr, 'UniformOutput',false);
bUnits = cellfun(@(r,g,t) ~r & ~g & cell2mat(t)==1, rUnits, gUnits, tUnit(1:length(rUnits)), 'UniformOutput',false);
typeUnits = cellfun(@(r,g,b) cellfun(@(a) find(a), num2cell([r,g,b],2),'Uniformoutput',false), rUnits,gUnits,bUnits, 'UniformOutput', false);
for t = 1:length(typeUnits)
    typeUnits{t}(cellfun(@isempty, typeUnits{t})) = {0};
    typeUnits{t} = cell2mat(typeUnits{t});
end
tPhys = table();
for c = 1:length(conditions)
    condTable = table();
    AUCVals = condXphase{c};
    tUnits = cell2mat(tUnit{c});
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
    condTable.Channel =  cell2mat(cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap',siteChannels, 'Uniformoutput', false))';
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    condTable.TaskUnits =  logical(tUnits);
    for pn = 1:length(phaseNames)
        condTable.(string(phaseNames(pn))) = AUCVals(:,pn);
    end
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
condIndex = any(tPhys{:,contains(tPhys.Properties.VariableNames,"Task") & ~contains(tPhys.Properties.VariableNames,"_R")},2); %& tPhys{:,'Channel'}>16;
condPhases = phaseVals(condIndex,:);
armPhases = condPhases(tPhys{condIndex,'Somatotopy'}=="Arm",:);
handPhases = condPhases(tPhys{condIndex,'Somatotopy'}=="Hand",:);
figure();
boxchart(reshape([repmat(condLabels+[-.3 -.3 -.3 -.1 -.1 -.1 .1 .1 .1 .3 .3 .3],size(armPhases,1),1);repmat(condLabels...
    +[-.3 -.3 -.3 -.1 -.1 -.1 .1 .1 .1 .3 .3 .3]+.1,size(handPhases,1),1)],1,[]),reshape([armPhases;handPhases],1,[]),'GroupByColor',...
    reshape(repmat(phaseLabels,size([armPhases;handPhases],1),1),1,[]),'Notch','on','MarkerStyle','none','BoxWidth',.4);
%bx=violin(condPhases,'x',4*(condLabels-1+repmat([-.3 -.1 .1 .3],1,3)),'facecolor',cm(phaseLabels,:));
ylim([0 5]);
set(gca,'XTick',[0:4:4*length(conditions)-1]);
set(gca,'XTickLabel',params.condAbbrev.values(conditions(1:length(conditions)-1)));
saveFigures(gcf,savePath,"Normalized_Violins",[]);
allSegs = cellfun(@(u,s) repmat(s,size(u,1),1), vertcat(condPSTHS{1:length(conditions)-1}),vertcat(segConds{1:length(conditions)-1}), 'UniformOutput',false);
psthLabs = arrayfun(@(c) cellfun(@(s) repmat(c,size(s,1),1),allSegs(1,:), 'UniformOutput',false), string(params.condAbbrev.values), 'UniformOutput',false);
plotJointPSTHS(params.bins,{cell2mat(horzcat(condPSTHS{1:length(conditions)-1})')},...
    {cell2mat(reshape(allSegs',1,[])')},cell2mat(horzcat(psthLabs{1:length(conditions)-1})'),...
    cell2mat(repmat(infoTable.TaskModulated,length(conditions)-1,1)),[], {[-1 3]},[1 2],cell2struct(num2cell(...
    distinguishable_colors(length(conditions),'r'),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath,"Normalized_PSTHS",[]);