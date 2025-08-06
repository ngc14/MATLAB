pVal = 0.05;
varNames = ["Unit" "SiteNum" "Monkey" "Somatotopy" "Channel"	"X"	"Y"	"Condition"...
     "TaskUnits" "Go" "Reach"	"Grasp"	"Withdraw", "rSI", "gSI", "rgSI","unitType"];
repNames = ["Arm", "Hand", "Trunk"];
% monkeyMapping = dictionary(true,"Gilligan",false,"Skipper");
conditions = [ "Extra Small Sphere","Large Sphere","Photocell"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartLift"]},{["GoSignal","StartLift"]},...
    {["GoSignal","StartHold"]}});
phaseNames = ["Go", "Reach", "Grasp", "Withdraw"];
savePath = "S:\Lab\ngc14\Working\Both\ChannelCheck\";
close all;
%%
m = matfile("S:\Lab\ngc14\Working\Both\Full_Baseline\phaseAnalysis.mat");
taskBaseline = m.taskBaseline;
taskFR = m.taskFR;
siteChannels = m.siteChannels;
siteChannels = siteChannels{2};
siteRep = m.siteRep;
siteDateMap = m.siteDateMap;
phaseFR = m.phaseFR;
avgPhase = m.avgPhase;
params = m.params;
avgSeg = m.avgSeg;
normPSTH = m.normPSTH;
clear m;
%%
[~,tUnitsCond] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,pVal),...
    pb,pc, 'UniformOutput', false),taskBaseline,taskFR,'UniformOutput',false);
tUnitsCond = cellfun(@cell2mat, tUnitsCond,'UniformOutput',false);
condRepInds = cellfun(@(t) arrayfun(@(s) t==1 & strcmp(mapSites2Units(...
    cellfun(@length ,siteChannels),siteRep'),s),repNames,'UniformOutput',false),tUnitsCond,'Uniformoutput',false);
tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(a) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
    a(cellfun(@(z) size(z,2), a)>=find(strcmp(phaseNames,pn))),'UniformOutput',false), s,'UniformOutput', false),...
    phaseFR,'UniformOutput', false), phaseNames, 'UniformOutput',false);
[~,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials(r,g,1,true,0.05),cr,cg,'UniformOutput',false),...
    tPhase{strcmp(phaseNames,"Reach")},tPhase{strcmp(phaseNames,"Grasp")}, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
%rgInds{cellfun(@isempty,rgInds)} = NaN(size(rgInds{find(~cellfun(@isempty,rgInds),1)}));
AUCVals = cellfun(@(c) [cell2mat(c),NaN(size(cell2mat(c),1),length(phaseNames)-size(cell2mat(c),2))], avgPhase,'UniformOutput',false);
rAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Reach")), AUCVals, 'UniformOutput',false);
gAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Grasp")), AUCVals, 'UniformOutput',false);
rUnits = cellfun(@(rg,r,g) rg==1 & r > g, rgInds, rAUC, gAUC, 'UniformOutput',false);
gUnits = cellfun(@(rg,r,g) rg==1 & r < g, rgInds, rAUC, gAUC, 'UniformOutput',false);
bUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, tUnitsCond(1:length(rUnits)), 'UniformOutput',false);
chUnitMap = cell(height(siteDateMap),1);
chUnitMap(strcmp([siteDateMap.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([siteDateMap.Monkey],"Skipper")) = {[32:-1:1]};
mappedChannels = cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap,siteChannels, 'Uniformoutput', false);
typeUnits = cellfun(@(r,g,b) cellfun(@(a) find(a), num2cell([r,g,b],2),'Uniformoutput',false), rUnits,gUnits,bUnits, 'UniformOutput', false);
for t = 1:length(typeUnits)
    typeUnits{t}(cellfun(@isempty, typeUnits{t})) = {0};
    typeUnits{t} = cell2mat(typeUnits{t});
end
%%
tPhys = [];
g = @(x,y,c)GetPointLineDistance(x,y,c(1),c(2),c(3),c(4));
condXphase = cellfun(@(pc) cell2mat(pc),cellfun(@(e) e(repmat(~isempty(e),size(e))),avgPhase,'UniformOutput', false),'UniformOutput',false);
condXphase = cellfun(@(s) [s,NaN(size(s,1),length(phaseNames)-size(s,2))], condXphase, 'UniformOutput',false);
condXphase = cellfun(@(cp,t) [cp,t], condXphase,typeUnits, 'UniformOutput',false);
for c = 1:length(conditions)
    condTable = table();
    AUCVals = condXphase{c};
    tUnits = tUnitsCond{c};
    tUnits(isnan(tUnits)) = 0;
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels);    
    allReps =  mapSites2Units(condUnitMapping,siteRep');
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    mInds= contains(mLabs,"Skipper");
    unitLocation = cell2mat(arrayfun(@(x,y,c) [g(x,y,OrthogonalLines(c).RCLine), g(x,y,OrthogonalLines(c).MLLine)],...
        siteDateMap.x,siteDateMap.y,siteDateMap.Monkey,'UniformOutput',false));
    unitLocation = mapSites2Units(condUnitMapping,num2cell([unitLocation(:,1),unitLocation(:,2)],2));
    % unitLocation(mInds,:) = round(transformPointsForward(mtform.tform,unitLocation(mInds,:)));

    condTable.Unit = [1:length(mLabs)]';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), [1:length(condUnitMapping)]',condUnitMapping,'UniformOutput',false));
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(allReps);
    condTable.Channel =  [mappedChannels{:}]';
    condTable.X = unitLocation(:,1);
    condTable.X(~mInds) = condTable.X(~mInds) - min(condTable.X(~mInds));
    condTable.X(mInds) = condTable.X(mInds) - min(condTable.X(mInds));
    condTable.Y = unitLocation(:,2);
    condTable.Y(~mInds) = condTable.Y(~mInds) - min(condTable.Y(~mInds));
    condTable.Y(mInds) = condTable.Y(mInds) - min(condTable.Y(mInds));
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    condTable.Condition = categorical(cellstr(repmat(conditions{c}(1),length(mLabs),1)));
    condTable.TaskUnits =  logical(tUnits);
    condTable.Go = AUCVals(:,strcmp(phaseNames,"Go"));
    condTable.Reach = AUCVals(:,strcmp(phaseNames,"Reach"));
    condTable.Grasp = AUCVals(:,strcmp(phaseNames,"Grasp"));
    condTable.Withdraw = AUCVals(:,strcmp(phaseNames,"Withdraw"));
    condTable.rSI = condTable.Reach./(condTable.Reach + condTable.Grasp);
    condTable.gSI = condTable.Grasp./(condTable.Reach + condTable.Grasp);
    condTable.rgSI = (condTable.Reach - condTable.Grasp)./(condTable.Reach + condTable.Grasp);
    condTable.unitType = typeUnits{c};
    % condTable.ISOI = logical(activityInds);
    condTable.Properties.VariableNames = varNames;
    tPhys = [tPhys;condTable];
end
taskAlignments = cellfun(@(tc,pc) cellfun(@(pa) arrayfun(@(s) find(strcmp(string(tc), s)),pa),pc,...
    'UniformOutput', false),params.condSegMap.values,taskAlign.values,'UniformOutput',false);
avgTask = cellfun(@(cn,ct,cp) cellfun(@(tp,tt)cellfun(@(ap,at) mean(ap(:,max(1,...
    findBins(at(cp{1}(1)),params.bins)):max(length(ap),findBins(at(cp{1}(end)),params.bins))),2,'omitnan'),...
    tp,tt,'UniformOutput', false),cn, ct,'UniformOutput', false),normPSTH(1:length(taskAlignments)),...
    avgSeg(1:length(taskAlignments)),taskAlignments,"UniformOutput",false);
meanTask = cellfun(@(c) cell2mat(cellfun(@(s) median(cat(3,s{:}),3,'omitnan'),c,'UniformOutput', false)), avgTask, 'UniformOutput',false);
plotNames = arrayfun(@(p) arrayfun(@(c) p+"_"+c{1}(1), conditions, 'UniformOutput', true), [phaseNames,"rSI","gSI","Task"], 'UniformOutput', false);
plotNames = [plotNames{:}];
tPhys = unstack(tPhys,varNames(9:end),"Condition");
tPhys = addvars(tPhys,meanTask{1},meanTask{2},meanTask{3}, 'NewVariableNames',plotNames(end-(length(conditions)-1):end));
%tPhys = addvars(tPhys, NaN(height(tPhys),1),'After','Y','NewVariableNames',"BIN");
%tPhys = addvars(tPhys, NaN(height(tPhys),1),'After','XB','NewVariableNames',"YB");
%quads = locationBins(tPhys,repNames);
% tPhys(~(tPhys.TaskUnits_E | tPhys.TaskUnits_L | tPhys.TaskUnits_P),:) = [];
tPhys = addvars(tPhys,mode(table2array(tPhys(:,{'unitType_E','unitType_L','unitType_P'})),2),'NewVariableNames','Type');
percs = table2array(rowfun(@(a,s,m) sum(a==tPhys.SiteNum & s==tPhys.Type)/sum(s==tPhys.Type), tPhys,...
    'InputVariables',{'SiteNum','Type','Monkey'}, 'OutputVariableNames', "Percentage"));
tPhys= addvars(tPhys,100.*percs,'After','Type','NewVariableNames','Percentage');

% conds x rep x (phase)
allChannels = cell2mat(mappedChannels')';
condXrepXphase = cellfun(@(p,ic) cellfun(@(i) [p(i==1 & allChannels<=16,:);NaN(max(cellfun(@(s)...
    sum(s,'omitnan'),ic))-sum(i==1),size(p,2))],ic(1:2),'UniformOutput',false),...
    condXphase,condRepInds,'UniformOutput',false)';
rgPlot = cellfun(@(r) cellfun(@(m) [m(:,strcmp(phaseNames,"Reach"))./...
    sum(m(:,contains(phaseNames,["Reach","Grasp"])),2),m(:,strcmp(phaseNames,"Grasp"))./...
    sum(m(:,contains(phaseNames,["Reach","Grasp"])),2)], r,'UniformOutput',false),condXrepXphase,'UniformOutput',false);
% writetable(tPhys,savePath+'Task_Units_FR_AVGPSTH.xlsx');
%%
plotGroupedBars(condXrepXphase,savePath+"UnmappedSuperficial_Units_Phys_FR_Box",true);
%%
trialFRS = cellfun(@(c) cellfun(@(s) cellfun(@(a) cat(3,a{:}), s, ...
    'UniformOutput',false),c,'UniformOutput',false), phaseFR(1:3),'UniformOutput',false);
trialFRS = cellfun(@(c) cellfun(@(s) median(cat(4,s{:}),4,'omitnan'), c,'UniformOutput', false),trialFRS,'UniformOutput',false);
trialFRS = cellfun(@(n) cellfun(@(t) num2cell(t,[2,3]),n,'UniformOutput',false), trialFRS, 'UniformOutput', false);
trialFRS = cellfun(@(c) vertcat(c{:}), trialFRS, 'UniformOutput',false);
trialFRS = cellfun(@(p) permute(p,[2 3 1]),[trialFRS{:}],'UniformOutput',false);
unitNum = transpose(1:size(trialFRS,1));
maxTrialNum = max(cellfun(@(s) size(s,1), trialFRS),[],2);
id = cell2mat(arrayfun(@(u,t) repmat(u,t,1), unitNum, maxTrialNum,'UniformOutput',false));
trialNum = num2cell(repmat(maxTrialNum,1,size(trialFRS,2)));
missingTrialInds = cellfun(@(s,m) size(s,1) ~= m, trialFRS,trialNum);
trialFRS(missingTrialInds) = cellfun(@(t,m) [t;NaN(m-size(t,1),size(t,2))], trialFRS(missingTrialInds), trialNum(missingTrialInds), 'UniformOutput',false);
trialFRS = cellfun(@(a) vertcat(a{:}),num2cell(trialFRS,1),'UniformOutput',false);

condPhaseNames = arrayfun(@(c) [arrayfun(@(p) strjoin([c,p],""), phaseNames)],conditions,'UniformOutput',false);
condPhaseNames = [condPhaseNames{:}];
trialFRSTable = addvars(array2table([trialFRS{:}],'VariableNames',condPhaseNames),id,'Before',1,'NewVariableName','UnitTrials');
varTableNames = trialFRSTable.Properties.VariableNames;
trialFRSTable = stack(trialFRSTable,arrayfun(@(p) find(contains(varTableNames,p)), phaseNames,'Uniformoutput', false),...
    'ConstantVariables',1,'NewDataVariableName',phaseNames, 'IndexVariableName',"Condition");
trialFRSTable.Condition = arrayfun(@(p) varTableNames{p}(1), trialFRSTable.Condition);
%%
rModel = fitrm(tPhys(all(isfinite(table2array(tPhys(:,11:end))),2),[3,11:end]),['Reach_E ,' ...
    'Reach_L , Reach_P , Grasp_E , Grasp_L , Grasp_P ,Withdraw_E , Withdraw_L ,'...
    'Withdraw_P ~ 1 + Somatotopy'],'WithinDesign',w);
mn = manova(tPhys,'Go,Reach,Grasp,Withdraw ~ 1 + Rep + Condition + Rep*Condition');
tb = groupmeans(mn,"Rep",MarginalMean=false);
ts = grpstats(tPhys,"Rep",["mean","sem"],"DataVars",2:5);
mc = multcompare(mn,["Rep"],"CriticalValueType","bonferroni");
mc = multcompare(mn,["Condition"],"CriticalValueType","bonferroni");

figure();
hold on;
aa = gca();
ls=plotprofile(mn,mn.FactorNames);
aa.ColorOrder = [cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)),repNames, 'UniformOutput', false)')];
aa.LineStyleOrder = [ "-",":","--"];
aa.LineStyleCyclingMethod = 'beforecolor';
set(gcf,'Renderer','painters');
orderCats = renamecats(categorical(string(aa.XAxis.Categories')),mn.ResponseNames);
aa.XAxis.Categories = orderCats;
aa.XAxis.TickValues = orderCats;
arrayfun(@(l) set(l,"LineStyleMode","auto"),ls);
arrayfun(@(l) set(l,"XData",orderCats),ls);
arrayfun(@(l) set(l,"YData",l.YData([1 3 2 4])),ls);
drawnow
legend
%%
monkey = "Gilligan";
monkeyCompDir = "S:\Lab\ngc14\Working\Both\Baseline_FR\Stats\" + monkey;
factorN = dir(monkeyCompDir+"\*.xls");
spOpts = detectImportOptions(monkeyCompDir + "\" + factorN(1).name);
spOpts.VariableNamesRange = '2:2';
spOpts.DataRange = 'A4';
for t = 1:length(factorN)
    sp = readtable(string([factorN(t).folder,'\',factorN(t).name]),spOpts);
    emptyRowVals =cell2mat(cellfun(@str2num,unique(table2array(sp(:,1)),'stable'),'UniformOutput',false));
    emptyRowInds = find(~cellfun(@isempty,table2array(sp(:,1))));
    emptyRowInds = emptyRowInds(diff(emptyRowInds)==mode( diff(emptyRowInds)))
    for p = 1:length(emptyRowInds)-1
        sp(emptyRowInds(p):emptyRowInds(p+1)-1,1) = array2table(repmat({emptyRowVals(p)},emptyRowInds(p+1)-emptyRowInds(p),1));
    end
    allVars = sp.Properties.VariableNames(1:find(contains(sp.Properties.VariableNames, 'MeanDifference_I_J_'))-1);
    sigInd = find(contains(sp.Properties.VariableNames, "Sig"));
    interactionName = cellfun(@(n) string(n), allVars(~cellfun(@(a) contains(a,'Var'), allVars)));
    siteInd = arrayfun(@(n) contains(n,"Site"), interactionName);
    interactionName(siteInd) = "Site";

    interactionName = strjoin(interactionName,'x');
end
%%
condAbbrevs = string(cellfun(@(r) erase(strjoin(r)," "),regexp(conditions,"[A-Z]",'match','matchcase'),'UniformOutput',false));
for p = 1:length(plotNames)
    figure('Units','normalized','Position',[0 0 1 1]);
    for r = 1:length(repNames)
        if(r<length(repNames))
            currRep = repNames(r);
            repInds = strcmp(string(tPhys.Somatotopy),repNames(r));
        else
            currRep = "Both";
            repInds = contains(string(tPhys.Somatotopy),["Arm","Hand"]);
        end
        currPlot = tPhys.(plotNames(p));
        currPlot = currPlot(repInds);
        XBP = tPhys.XB;
        XBP = XBP(repInds);
        YBP = tPhys.YB;
        YBP = YBP(repInds);
        subplot(1,2,1);
        boxchart(XBP,currPlot ,'Notch', 'on');
        ylabel(plotNames(p));
        title("X" + currRep);
        xlabel("Caudal to Rostral");
        ylim([min(currPlot(~isoutlier(currPlot))), max(currPlot(~isoutlier(currPlot)))])
        subplot(1,2,2);
        boxchart(YBP,currPlot, 'Notch', 'on');
        ylabel(plotNames(p));
        xlabel("Medial to Lateral");
        title("Y" + currRep);
        ylim([min(currPlot(~isoutlier(currPlot))), max(currPlot(~isoutlier(currPlot)))])

        condCase = regexp(plotNames(p),'_','split');
        currPath = savePath + condAbbrevs(startsWith(conditions,condCase)) + "\XY\" + currRep+"\";
        if(~exist(currPath,'dir'))
            mkdir(currPath);
        end
        saveFigures(gcf, currPath, plotNames(p), []);
    end
    close all;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phaseSubVals, sigs] = ttestTrials(dist1,dist2,taskPhase,paired,pVal)
if(paired)
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)>=floor(length(dist1)/2);
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}./...
        d1{taskPhase}),dist1,dist2,'UniformOutput',false)',1);
else
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest2(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)'>=floor(length(dist1)/2);
    smallestTrialCount = min(size(dist1{1}{taskPhase},2),size(dist2{1}{taskPhase},2));
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}(:,1:smallestTrialCount)-...
        d1{taskPhase}(:,1:smallestTrialCount)),dist1,dist2,'UniformOutput',false)',1);
end
phaseSubVals = cellfun(@(m) nanmean(cat(3,m{:}),3), phaseSubVals, 'UniformOutput', false);
sigs = double(sigs);
for s = 1:length(phaseSubVals)
    if(isnan(phaseSubVals{s}))
        sigs = NaN(size(sigs));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quadR = locationBins(tbl,repNames)
XVALS = tbl.X;
YVALS = tbl.Y;
XHISTBINS = [min(tbl.X),linspace(110,260,5),max(tbl.X)];
YHISTBINS = [min(tbl.Y),linspace(125,425,5),max(tbl.Y)];
for r = 1:length(repNames)-1
    currRep = strcmp(string(tbl.Somatotopy),repNames(r));
    XHISTBINS = [min(XVALS(currRep)),quantile(XVALS(currRep), [.55]), max(XVALS(currRep))];
    YHISTBINS = [min(YVALS(currRep)),quantile(YVALS(currRep), [.55]), max(YVALS(currRep))];

    [quads,cn,wsum,bsum] = kmeans(unique([XVALS(currRep),YVALS(currRep)],'rows'),4);
    quads = cell2mat(cellfun(@(c,q) repmat(q,1,sum(all([XVALS(currRep),YVALS(currRep)]==c,2)))', ...
        num2cell(unique([XVALS(currRep),YVALS(currRep)],'rows'),2),num2cell(quads), 'UniformOutput',false));
    [G1,G2] = meshgrid(1:768,1:768);
    grid = [G1(:),G2(:)];
    idx2Region = kmeans(grid,4,'MaxIter',1,'Start',cn);
    figure;
    gscatter(grid(:,1),grid(:,2),idx2Region,colormap(hsv(4)),'..');
    hold on;
    plot(XVALS(currRep),YVALS(currRep),'k*','MarkerSize',5);
    axis ij;
    tbl.BIN(currRep) = quads;
    title(repNames(r));
    legend(arrayfun(@(s) strjoin([num2str(s),": ",num2str(sum(quads==s)),",",...
        num2str(100*(sum(quads==s)/length(quads)),2),"%"]),unique(quads),'UniformOutput',false));
    quadR{r} = quads;
    %saveFigures(gcf,savePath,repNames(r),[]);
end
end