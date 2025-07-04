conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions(1:end-1),{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]}});
taskWindow = {[0.2, 0]};
alignLimits = {[-.75, 1]};
pVal=0.05;
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\";
monkey = "Gilligan";
MIN_BLOCKS_FOR_UNIT = 15;
params = PhysRecording(string(conditions),.01,.15,-6,5,containers.Map(conditions,...
    {["StartReach"],["StartReach"],["StartReach"],"GoSignal"}));
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions] = getAllSessions(params,"Single","PMd");
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs(1:end-1),siteTrialPSTHS(1:end-1),false,true);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(params.bins,mean(t(:,1),'omitnan')-3),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) ...
    permute(mean(c(:,s:s+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),a(~isnan(n)),...
    num2cell(n(~isnan(n))),'UniformOutput',false),[1,1,sum(~isnan(n))])),3,'omitnan'));NaN(all(isnan(n)).*size(a{1},1),1)],p,t,...
    'UniformOutput',false),siteTrialPSTHS,allCondSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp{:}),repmat(nb,1,size(vertcat(cp{:}),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
%%
allGoodTrials = cellfun(@(ct,cs) cellfun(@(a,b) ~(cellfun(@length,a{1}) <=(repmat(b{1}(:,end)-b{1}(:,1),1,size(a{1},1))')...
    | cellfun(@length,a{1})>200*(repmat(b{1}(:,end)-b{1}(:,1),1,size(a{1},1))') | ...
    cellfun(@(e) any(isnan(e)),a{1})), ct,cs,'UniformOutput',false), rawSpikes, siteSegs,'UniformOutput',false);
goodUnits = cellfun(@(tn) cell2mat(cellfun(@(s)sum(s,2), tn,'UniformOutput',false)),...
    num2cell(cat(2,allGoodTrials{:}),2),'UniformOutput',false);
[tVals,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
tUnits = cellfun(@(a,b) any(cell2mat(a),2) & sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2), ...
    num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false);
condUnitMapping = cellfun(@(si) size(si,2),siteChannels{2})';
%%
close all;
[allTrials,allPSTHS]= deal(cell(1,1));
[condInds,allSiteInds] = deal([]);
for c =1:length(conditions)
    condSegs = siteSegs{c};
    trialSegs = vertcat(condSegs{:});
    mMapping = cell2mat(arrayfun(@(m,n) ones(1,m)*n,condUnitMapping,...
        1:length(condUnitMapping),'UniformOutput',false));
    sessionInds = ~isnan(mMapping) & vertcat(tUnits{:})';
    siteUnitMods = mMapping;
    siteUnitMods(~sessionInds) = 0;
    [siteIndsN,siteInds] = unique(siteUnitMods);
    siteInds = siteInds(siteIndsN>0);
    siteUnits = mapSites2Units(condUnitMapping,cellfun(@(c) string(datetime(c,'Format','MMMM_dd')),siteDateMap.Date','UniformOutput',false));
    siteUnits(siteUnitMods==0) = "";

    plotColors =cell2struct(num2cell(distinguishable_colors(length(siteDateMap.Date)),2),...
        arrayfun(@(s) char(datetime(s,'Format','MMMM_dd')),string(siteDateMap.Date),'UniformOutput',false));
    condPSTHS = num2cell(cellfun(@(m) mean(m,3,'omitnan'),vertcat(normPSTH{c}{:}),'UniformOutput',false),1);
    %unitPSTHS = cellfun(@(m,p) (m./m).*vertcat(p{:}),tUnits,condPSTHS,'UniformOutput',false);
    unitPSTHS = cellfun(@(n) cellfun(@(t) squeeze(t)',num2cell(n,[2,3]),'UniformOutput',false), ...
        vertcat(normPSTH{c}{:}),'Uniformoutput',false);
    % for s = 1:length(unitPSTHS)
    %     currsiteUnits =  find(siteUnitMods==s);
    %     if(~isempty(currsiteUnits))            
    %         condtrialPSTHS = {cell2mat(unitPSTHS{s}(tUnits{s}))};
    %         allTrialSegs = trialSegs{s};
    %         plotColors =cell2struct(num2cell(distinguishable_colors(length(currsiteUnits)),2),...
    %             string(arrayfun(@(t) char(datetime(siteDateMap.Date{s},'Format','MMMM_dd'))+"_"+t,currsiteUnits,'UniformOutput',false)));
    %         trialSiteUnits = cellfun(@(n) cellstr(repmat(string(n),size(allTrialSegs,1),1)),fieldnames(plotColors),'UniformOutput',false);
    %         allTrialLabels = cellfun(@string,vertcat(trialSiteUnits{:}));            
    %         plotJointPSTHS(params.bins,condtrialPSTHS,{repmat(allTrialSegs,length(currsiteUnits),1)},....
    %             allTrialLabels,true(size(allTrialLabels)),[],  alignLimits,[0 25],plotColors);
    %         saveFigures(gcf,savePath+"\"+string(datetime(siteDateMap.Date{s},'Format','MMMM_dd'))+"\",...
    %             params.condAbbrev(params.condNames(c))+"_PSTH",[]);
    %     end
    % end
    siteUnits(siteUnitMods>0) = "All";
    allPSTHS = cellfun(@(a,b) vertcat(a,b), allPSTHS, {cellfun(@(u,t) u(t,:), condPSTHS{1},tUnits,'UniformOutput',false)}, 'UniformOutput', false);
    condInd = repmat(string(params.condAbbrev(params.condNames(c))),1,length(siteUnits));
    condInd(siteUnitMods==0) = "";
    condInds = [condInds,condInd(condInd~="")];
    allSiteInds = [allSiteInds;cellfun(@sum,tUnits)];
    allTrials = cellfun(@(a,d) vertcat(a,cellfun(@(m,u)repmat(mean(m,1,'omitnan'),sum(u),1),...
        d,tUnits,'UniformOutput',false)), allTrials,{trialSegs}, 'UniformOutput',false);
end 
allTrials = cellfun(@(c) cellfun(@(t) [t(:,1:3), NaN(size(t,1),size(t,2)==8),t(:,4:end)],c,'UniformOutput',false),allTrials, 'UniformOutput',false);
allTrials = cellfun(@(c) cellfun(@(t) [t, NaN(size(t,1),5*double(size(t,2)==4))],c,'UniformOutput',false),allTrials, 'UniformOutput',false);
plotJointPSTHS(params.bins,{vertcat(allPSTHS{:}{:})},{vertcat(allTrials{:}{:})},condInds,cumsum(allSiteInds),[],  alignLimits,[0 15],...
    cell2struct(num2cell(distinguishable_colors(length(conditions)),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath,"All_PSTH",[]);
%% task phase units
maxCondUnits = max(cellfun(@length, taskFR));
reachVgrasp = zeros(1,maxCondUnits);
taskU = false(1,maxCondUnits);
[SIAvg,SI] = deal(zeros(1,maxCondUnits));
for u = 1:maxCondUnits
    baselinePhase = cell2mat(cellfun(@(bf) bf{u}, baselineFR', 'UniformOutput', false));
    taskPhase = cell2mat(cellfun(@(tf) tf{u}, taskFR', 'Uniformoutput', false));
    reachPhase = cell2mat(cellfun(@(tr) tr{u}, phaseFR(1:2,3), 'UniformOutput', false));
    graspPhase = cell2mat(cellfun(@(tg) tg{u}, phaseFR(1:2,4), 'UniformOutput', false));
    [~,pT] = ttest2(reachPhase(~isnan(reachPhase)), graspPhase(~isnan(graspPhase)));
    reachVgrasp(u) = (pT<pVal) + (nanmean(reachPhase) > nanmean(graspPhase));
    [~,pTU] = ttest2(baselinePhase(~isnan(baselinePhase)), taskPhase(~isnan(taskPhase)));
    taskU(u) = pTU < pVal;
    SIAvg(u) = (nanmean(reachPhase)-nanmean(graspPhase))/(nanmean(reachPhase) + nanmean(graspPhase));
    SI(u) = nanmean((reachPhase-graspPhase)./(reachPhase+graspPhase));
    %     SIspheres = cellfun(@(rc,gc) cellfun(@(r,g) (nanmean(r)-nanmean(g))/...
    %         (nanmean(r)+nanmean(g)), rc, gc,'UniformOutput',true),phaseFR(:,3),phaseFR(:,4),'UniformOutput', false);
    %     SIspheres = cell2mat(SIspheres(1:2));
    %     [~,largestInd] = max(abs(SIspheres),[],1);
    %     SIspheres = SIspheres(sub2ind(size(SIspheres),largestInd,1:length(SIspheres)));
    %     SIspheres(~taskUnitInds) = NaN;
end
% task modulated unit distribution
figTaskD = modulatedUnitsPerRep(repSave,{taskU}," of units","Task modulated (all conditions)");
saveFigures(figTaskD, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Task Modulated\'],[fileSaveName],[])
figCondMods = modulatedUnitsPerRep(repSave,{logical(condModUnits)},' of units',"Units modulated differently across condition (task phase)");
saveFigures(figCondMods,[saveDir,'Representation Distribution\Min Representation\Condition Effect\'],[fileSaveName],[]);


%% compare condition FR Changes between representations and phases
[figJoint,figPhaseCond] = plotConditionTimings(conditions,masterPSTH,...
    repSave,phaseFR,baselineFR,masterActivity);
if(saveFigs)
    saveName = [saveDir,'Conds\',monkey,'_',sessionOrUnitPSTHS,'_',singleOrAll];
    saveFigures(figJoint,[saveDir,'Conds\Min Representations\'],['CondsByJoint_',fileSaveName],[]);
    saveFigures(figPhaseCond,[saveDir,'Conds\Min Representations\'],['CondsByPhase_',fileSaveName],[]);
end
%% plot SI representation distribution
reachPhaseConds = cellfun(@(cIn) cellfun(@nanmean, cIn), phaseFR(1:2,3), 'UniformOutput', false)
graspPhaseConds = cellfun(@(cIn) cellfun(@nanmean, cIn), phaseFR(1:2,4), 'UniformOutput', false)
rgCondsAvg = cellfun(@(r,g) (r-g)./(r+g), reachPhaseConds, graspPhaseConds, 'UniformOutput', false);
rgConds = cellfun(@(rin,gin) cellfun(@(r,g) nanmean((r-g)./(r+g)), rin, gin, 'UniformOutput', true),  phaseFR(1:2,3),  phaseFR(1:2,4), 'UniformOutput', false);
reachVgraspConds = cellfun(@(rin,gin) cellfun(@(r,g) ttest2(r(~isnan(r)),...
    g(~isnan(g))),rin,gin,'UniformOutput', false), phaseFR(1:2,3), phaseFR(1:2,4), 'UniformOutput', false);
reachVgraspConds = padMissingCells(reachVgraspConds,  cellfun(@(cin) cellfun(@isempty, cin),reachVgraspConds, 'UniformOutput', false));
reachVgraspConds = cellfun(@(c) cell2mat(c), reachVgraspConds, 'UniformOutput', false);
reachVgraspConds = cellfun(@(c,r,g) c + (nanmean(r) > nanmean(g)), reachVgraspConds, reachPhaseConds, graspPhaseConds, 'UniformOutput', false);

figSIT = modulatedUnitsPerRep(repSave,{SI},...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"SI trial by trial then averaged from spheres");
ylim([-.1 .1]);
saveFigures(figSIT, [saveDir, 'Representation Distribution\Min Representation\SI\'],...
    ['Combined_',fileSaveName,'_Trial_Averaged'],[])
figSIU = modulatedUnitsPerRep(repSave,{SIAvg},...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"SI averaged per unit from spheres");
ylim([-.1 .1]);
saveFigures(figSIU, [saveDir, 'Representation Distribution\Min Representation\SI\'],...
    ['Combined_',fileSaveName,'Averaged_Per_Unit'],[])


figSIT2 = modulatedUnitsPerRep(repSave,rgCondsAvg',...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"SI averaged per unit per sphere");
ylim([-.1 .1]);
saveFigures(figSIT2, [saveDir, 'Representation Distribution\Min Representation\SI\'],...
    ['Spheres_',fileSaveName,'_Trial_Averaged'],[])
figSIU2 = modulatedUnitsPerRep(repSave,rgConds',...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"SI trial by trial then averaged per sphere");
ylim([-.1 .1]);
saveFigures(figSIU2, [saveDir, 'Representation Distribution\Min Representation\SI\'],...
    ['Spheres_',fileSaveName,'Averaged_Per_Unit'],[])
allRepNames = fieldnames(repColors);

f1 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
f2 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
f3 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
f4 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;

for r = 1:length(allRepNames)
    currInds = cellfun(@(s) strcmp(s,allRepNames{r}), repSave);
    for f = 1:4
        if(f==1)
            figure(f1);
            plotVals = {SIAvg(currInds & any(cell2mat(taskUnitInds(1:2)')))};
        elseif(f==2)
            figure(f2);
            plotVals = {SI(currInds & any(cell2mat(taskUnitInds(1:2)')))};
        elseif(f==3)
            figure(f3);
            SI1 = cellfun(@(r,g) (r-g)./(r+g), reachPhaseConds, graspPhaseConds, 'UniformOutput', false);
            plotVals = cellfun(@(s,t) s(currInds & t), SI1,taskUnitInds(1:2)', 'UniformOutput', false);
        elseif(f==4)
            figure(f4);
            SI2 = cellfun(@(rIn, gIn) cellfun(@(r,g) nanmean((r-g)./(r+g)), ...
                rIn, gIn), phaseFR(1:2,3), phaseFR(1:2,4), 'UniformOutput', false);
            plotVals = cellfun(@(s,t) s(currInds & t), SI2, taskUnitInds(1:2)', 'UniformOutput', false);
        end
        for v = 1:size(plotVals,1)
            subplot(2,4,r);hold on;
            h = histogram(plotVals{v}(~isnan(plotVals{v})),linspace(-.5,.5,10),'FaceColor',repColors.(allRepNames{r}));
            title([allRepNames{r}, ', n=', num2str(sum(currInds))]);
            xlim([-.5 .5]);
            xticks([-.5 -.25 0 .25 .5]); ylim([0 125]);
            xlabel(sprintf(['Grasp ', blanks(10), 'SI', blanks(10),'Reach']));
            if(r==1 | r==5)
                ylabel('# of units');
            end
        end
        axx=gca;
        histData = get(axx.Children(strcmp(get(axx.Children,'Type'),'histogram')),'Data');
        if(f>2)
            arrayfun(@(a,c) set(a,'EdgeColor', c{1}), axx.Children(strcmp(get(axx.Children,'Type'),'histogram')),{[1 0 0]; [1 .6 0]})
            histData = cellfun(@nanmean,histData);%sum((h.BinEdges(1:end-1)+h.BinWidth) .* h.Values)/sum(h.Values)
            if(histData(1)>histData(2))
                cellfun(@(h,n,a) text(h,100,n,'HorizontalAlignment', a), num2cell(histData'),{'S', 'L'}, {'left','right'})
            else
                cellfun(@(h,n,a) text(h,100,n,'HorizontalAlignment', a), num2cell(histData'),{'S', 'L'}, {'right','left'})
            end
            arrayfun(@(hx,c) plot(repmat(hx,1,2),[0 125], 'LineStyle', '--',...
                'Color',c{1},'LineWidth', 1),histData,{[1 0 0]; [1 .6 0]});
        else
            histData = cellfun(@nanmean,{histData});%sum((h.BinEdges(1:end-1)+h.BinWidth) .* h.Values)/sum(h.Values)
            arrayfun(@(hx) plot(repmat(hx,1,2),[0 125], 'LineStyle', '--',...
                'Color',[1 1 0],'LineWidth', 1),histData);
        end
        if(f==1)
            saveFigures(gcf, [saveDir, 'Representation Distribution\Min Representation\SI\ByRep\'],['Combined_',fileSaveName,'_Trial_By_Trial_Average'],[])
        elseif(f==2)
            saveFigures(gcf, [saveDir, 'Representation Distribution\Min Representation\SI\ByRep\'],['Combined_',fileSaveName,'_Averaged_Per_Unit'],[])
        elseif(f==3)
            saveFigures(gcf, [saveDir, 'Representation Distribution\Min Representation\SI\ByRep\'],['Spheres_',fileSaveName,'_Trial_By_Trial_Phase_Average'],[])
        elseif(f==4)
            saveFigures(gcf, [saveDir, 'Representation Distribution\Min Representation\SI\ByRep\'],['Spheres_',fileSaveName,'_Averaged_Per_Unit_Phase'],[])
        end
    end
end


figReachDA = modulatedUnitsPerRep(repSave,cellfun(@(r) r==2, reachVgraspConds(1), 'UniformOutput', false),...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"Reach > grasp phase (small sphere)");
saveFigures(figReachDA, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Spheres\Reach Modulated\'],['Small_Reach_Exclusive_',fileSaveName],[])
figReachD2 = modulatedUnitsPerRep(repSave,cellfun(@(r) r==1, reachVgraspConds(1), 'UniformOutput', false),...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"Grasp > reach phase (small sphere)");
saveFigures(figReachD2, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Spheres\Grasp Modulated\'],['Small_Grasp_Exclusive_',fileSaveName],[])

figReachDA = modulatedUnitsPerRep(repSave,cellfun(@(r) r==2, reachVgraspConds(2), 'UniformOutput', false),...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"Reach > grasp phase (large sphere)");
saveFigures(figReachDA, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Spheres\Reach Modulated\'],['Large_Reach_Exclusive_',fileSaveName],[])
figReachD2 = modulatedUnitsPerRep(repSave,cellfun(@(r) r==1, reachVgraspConds(2), 'UniformOutput', false),...
    sprintf("Reach \n \n \n SI \n \n \n Grasp"),"Grasp > reach phase (large sphere)");
saveFigures(figReachD2, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Spheres\Grasp Modulated\'],['Large_Grasp_Exclusive_',fileSaveName],[])

allRepNames = fieldnames(repColors);

figReachD = modulatedUnitsPerRep(repSave,{reachVgrasp==2}," of units","Reach > grasp phase (combined spheres)");
saveFigures(figReachD, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\Reach Modulated\'],['Combined_Reach_Exclusive_',fileSaveName],[])
figGraspD = modulatedUnitsPerRep(repSave,{reachVgrasp==1}," of units","Grasp > reach phase (combined spheres)");
saveFigures(figGraspD, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\Grasp Modulated\'],['Combined_Grasp_Exclusive_',fileSaveName],[])

reachCounts = arrayfun(@(ind) sum(reachVgrasp'==2 & unitSessionLabels==ind),unitSessionLabels);
graspCounts = arrayfun(@(ind) sum(reachVgrasp'==1 & unitSessionLabels==ind),unitSessionLabels);
reachPercs = 100*(reachCounts./arrayfun(@(ind) sum(unitSessionLabels==ind),unitSessionLabels));
graspPercs = 100*(graspCounts./arrayfun(@(ind) sum(unitSessionLabels==ind),unitSessionLabels));

figReachDP = modulatedUnitsPerRep(repSave,{reachPercs},...
    "% of units","Average % of units with reach > grasp phase (combined spheres)");
ylim([0 60]);
saveFigures(figReachDP, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\Reach Modulated\'],...
    ['_Reach_Exclusive_',fileSaveName,'_AvgPerc'],[])
figGraspDP = modulatedUnitsPerRep(repSave,{graspPercs},...
    "% of units","Average % of units with grasp > reach phase (combined spheres)");
ylim([0 60])
saveFigures(figGraspDP, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\Grasp Modulated\'],...
    ['_Grasp_Exclusive_',fileSaveName,'_AvgPerc'],[])

figTaskPercAvg = modulatedUnitsPerRep(repSave,{reachCounts-graspCounts},{'Grasp', 'Reach'}, "# of reach units - # of grasp units per site");
ylim([-10 0])
saveFigures(figTaskPercAvg, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\R_G\'],[fileSaveName,'_PercSite'],[])
figTaskPercAvg = modulatedUnitsPerRep(repSave,{reachPercs-graspPercs},{'Grasp', 'Reach'}, "% of reach units - % of grasp units per site");
ylim([-50 0]);
saveFigures(figTaskPercAvg, [saveDir, 'Representation Distribution\Min Representation\Phase Effect\Combined\R_G\'],[fileSaveName,'_PercSite'],[]);

percRFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,reachPercs);
title([fileSaveName,'_Reach_Perc'])
saveFigures(percRFig, [saveDir, 'Count Maps\Phase Effect\Reach Mod\'],['Reach_Exclusive_',fileSaveName,'_Perc'],[])
percGFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,graspPercs);
title([fileSaveName,'_Grasp_Perc'])
saveFigures(percGFig, [saveDir, 'Count Maps\Phase Effect\Grasp Mod\'],['Grasp_Exclusive_',fileSaveName,'_Perc'],[])

countRFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,reachCounts);
title([fileSaveName,'_Reach_Count'])
countGFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,graspCounts);
title([fileSaveName,'_Grasp_Count'])
saveFigures(countRFig, [saveDir, 'Count Maps\Phase Effect\Reach Mod\'],['Reach_Exclusive_',fileSaveName,'_Count'],[]);
saveFigures(countGFig, [saveDir, 'Count Maps\Phase Effect\Grasp Mod\'],['Grasp_Exclusive_',fileSaveName,'_Count'],[])

reachVgrasp(reachVgrasp==2) = -1;
reachVgrasp(reachVgrasp==0) = NaN;
subMap = arrayfun(@(ind) nansum(-reachVgrasp(unitSessionLabels==ind)),unitSessionLabels);
subMapFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,subMap);
title([fileSaveName,'_Reach-Grasp'])
saveFigures(subMapFig, [saveDir, 'Count Maps\Phase Effect\R_G\'],['ReachEx-GraspEx_',fileSaveName],[])
%%


for c = 1:length(movementConds)
    figPSTHS = plotJointPSTHS(bins,masterPSTH{c},masterSegs{c},...
        repSave,masterActivity{c},sessionInds);
    taskSessionInds = logical([1; any(diff(cell2mat(siteLocation(taskUnitInds{c}==1)')),2)]);
    figPSTHSTask = plotJointPSTHS(bins,masterPSTH{c}(taskUnitInds{c}==1),...
        masterSegs{c}(taskUnitInds{c}==1),repSave(taskUnitInds{c}==1),...
        masterActivity{c}(taskUnitInds{c}==1),taskSessionInds);
    saveFigures(figPSTHS, [saveDir, 'PSTHS\Min Representation\', condAbbrev{c},'\PSTHS\'],[fileSaveName],[])
    saveFigures(figPSTHSTask, [saveDir, 'PSTHS\Min representation\', condAbbrev{c},'\PSTHS\'],[fileSaveName],[])
end
cCount = arrayfun(@(ind) sum(condModUnits' & unitSessionLabels==ind),unitSessionLabels);
cPercs = 100*(cCount./arrayfun(@(ind) sum(unitSessionLabels==ind), unitSessionLabels));
countCFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,cCount);
title([fileSaveName,'_CondEffect_Count'])
saveFigures(countCFig, [saveDir, 'Count Maps\Condition Effect\'],[fileSaveName,'_Count'],[])
percCFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,cPercs);
title([fileSaveName,'_CondEffect_Perc'])
saveFigures(percCFig, [saveDir, 'Count Maps\Condition Effect\'],[fileSaveName,'_Perc'],[])

tCount = arrayfun(@(ind) sum(taskU' & unitSessionLabels==ind), unitSessionLabels);
tPercs = 100*(tCount./arrayfun(@(ind) sum(unitSessionLabels==ind), unitSessionLabels));
countTFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,tCount);
title([fileSaveName,'_Task_Count'])
saveFigures(countTFig, [saveDir, 'Count Maps\Phase Effect\Task Mod\'],['TaskMod_',fileSaveName,'_Count'],[])
percTFig = mapUnitVals(verticies, vCells, vMask, siteMasks,sessionInds,tPercs);
title([fileSaveName,'_Task_Perc'])
saveFigures(percTFig, [saveDir, 'Count Maps\Phase Effect\Task Mod\'],['TaskMod_',fileSaveName,'_Perc'],[])
%% PSTH FR phase changes, peak rise and fall time, AUC
for c = 1:length(movementConds)
    [figPhase,figTimes,figAUC] = plotPSTHTimings(bins,masterPSTH{c},...
        repSave,baselineFR{c},masterPhaseBins(c,:) ,phaseFR(c,:),masterActivity{c},alignmentPoint);
    [figPhaseTask,figTimesTask,figAUCTask] = plotPSTHTimings(bins,...
        masterPSTH{c}(taskUnitInds{c}==1),repSave(taskUnitInds{c}==1),...
        baselineFR{c}(taskUnitInds{c}==1),cellfun(@(b) b(taskUnitInds{c}==1), masterPhaseBins(c,:), 'UniformOutput', false),...
        cellfun(@(p) p(taskUnitInds{c}==1), phaseFR(c,:), 'UniformOutput', false),...
        masterActivity{c}(taskUnitInds{c}==1),alignmentPoint);
    if(saveFigs)
        if(strcmp(sessionOrUnitPSTHS, 'Session'))
            saveFolder = [saveDir,'\PSTHS\Min Representation\',condAbbrev{c},'\Session\'];
        else
            saveFolder = [saveDir,'\PSTHS\Min Representation\',condAbbrev{c},'\'];
        end
        saveFigures(figPSTHS,[saveFolder,'PSTHS\'],['PSTHS_',fileSaveName],[]);
        saveFigures(figPhase,[saveFolder,'Phase\'],['Phase_',fileSaveName],[]);
        saveFigures(figTimes,[saveFolder,'Times\'],['PSTHTimes_',fileSaveName],[]);
        saveFigures(figAUC,[saveFolder,'AUC\'],['AUC_',fileSaveName],[]);
        saveFigures(figPSTHSTask,[saveFolder,'PSTHS\'],['PSTHS_',fileSaveName,'_Task'],[]);
        saveFigures(figPhaseTask,[saveFolder,'Phase\'],['Phase_',fileSaveName,'_Task'],[]);
        saveFigures(figTimesTask,[saveFolder,'Times\'],['PSTHTimes_',fileSaveName,'_Task'],[]);
        saveFigures(figAUCTask,[saveFolder,'AUC\'],['AUC_',fileSaveName,'_Task'],[]);
    end
end









%% task phase units
taskFRBaseline = cellfun(@(tbs) tbs(1,:), taskFRChanges, 'UniformOutput', false);
taskFRPhase = cellfun(@(tbs) tbs(2,:), taskFRChanges, 'UniformOutput', false);
maxCondUnits = max(cellfun(@length, taskFRChanges));
taskUnits = false(1,maxCondUnits);
reachUnits = false(1,maxCondUnits);
graspUnits = false(1,maxCondUnits);
for u = 1:maxCondUnits
    currBase = cell2mat(cellfun(@(tb) tb{u}, taskFRBaseline, 'UniformOutput', false)');
    currPhase = cell2mat(cellfun(@(tp) tp{u}, taskFRPhase, 'UniformOutput', false)');
    reachPhase = cell2mat(cellfun(@(tr) tr{u}, masterFRChanges(:,3), 'UniformOutput', false));
    graspPhase = cell2mat(cellfun(@(tg) tg{u}, masterFRChanges(:,4), 'UniformOutput', false));
    [~,pT] = ttest2(currBase(~isnan(currBase)), currPhase(~isnan(currPhase)));
    taskUnits(u) = pT<pVal;
    [~,pR] = ttest2(currBase(~isnan(currBase)), reachPhase(~isnan(reachPhase)));
    reachUnits(u) = pR<pVal;
    [~,pG] = ttest2(currBase(~isnan(currBase)), graspPhase(~isnan(graspPhase)));
    graspUnits(u) = pG<pVal;
end

figTaskUnits = modulatedUnitsPerRep({masterJoints},{taskUnits},...
    repNames,"","Task modulated (combined conditions)");
saveFigures(figTaskUnits, [saveDir, 'Representation Distribution\Task Modulated\'],[fileSaveName],[])
figReachUnits = modulatedUnitsPerRep({masterJoints},{reachUnits},...
    repNames,"","Reach modulated (combined conditions)");
saveFigures(figReachUnits, [saveDir, 'Representation Distribution\Reach Modulated\'],[fileSaveName],[])
figGraspUnits = modulatedUnitsPerRep({masterJoints},{graspUnits},...
    repNames,"","Grasp modulated (combined conditions)");
saveFigures(figGraspUnits, [saveDir, 'Representation Distribution\Grasp Modulated\'],[fileSaveName],[])


taskCounts = arrayfun(@(ind) sum(taskUnits' & unitSessionLabels==ind),...
    unitSessionLabels);
countFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,taskCounts,cm);
saveFigures(countFig, [saveDir, 'Count Maps\Task Mod\'],[fileSaveName,'_Count'],[])

reachCounts = arrayfun(@(ind) sum(reachUnits' & unitSessionLabels==ind),...
    unitSessionLabels);
countRFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,reachCounts,cm);
saveFigures(countRFig, [saveDir, 'Count Maps\Reach Mod\'],[fileSaveName,'_Count'],[])

graspCounts = arrayfun(@(ind) sum(taskUnits' & unitSessionLabels==ind),...
    unitSessionLabels);
countGFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,graspCounts,cm);
saveFigures(countGFig, [saveDir, 'Count Maps\Grasp Mod\'],[fileSaveName,'_Count'],[])

taskPercs = arrayfun(@(ind) round(100*sum(taskUnits' & ...
    unitSessionLabels==ind)/sum(unitSessionLabels==ind)), unitSessionLabels);
percFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,taskPercs,cm);
saveFigures(percFig, [saveDir, 'Count Maps\Task Mod\'],[fileSaveName,'_Perc'],[])
reachPercs = arrayfun(@(ind) round(100*sum(reachUnits' & ...
    unitSessionLabels==ind)/sum(unitSessionLabels==ind)), unitSessionLabels);
percRFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,reachPercs,cm);
saveFigures(percRFig, [saveDir, 'Count Maps\Reach Mod\'],[fileSaveName,'_Perc'],[])
graspPercs = arrayfun(@(ind) round(100*sum(graspUnits' & ...
    unitSessionLabels==ind)/sum(unitSessionLabels==ind)), unitSessionLabels);
percGFig = mapUnitVals(verticies, vCells, vMask, siteMasks,...
    sessionInds,graspPercs,cm);
saveFigures(percGFig, [saveDir, 'Count Maps\Grasp Mod\'],[fileSaveName,'_Perc'],[])

reachPercsT = arrayfun(@(ind) round(100*sum(reachUnits' & ...
    (taskUnits' & unitSessionLabels==ind))/sum(unitSessionLabels==ind & taskUnits')), unitSessionLabels);
percRTFig= modulatedUnitsPerRep({masterJoints},{reachPercsT},repNames,"% of task units","Reach modulated (combined conditions)");
saveFigures(percRTFig, [saveDir, 'Representation Distribution\Reach Percentage\'],[fileSaveName,'_PercOfTask'],[])
graspTPercs = arrayfun(@(ind) round(100*sum(graspUnits' & ...
    (taskUnits' & unitSessionLabels==ind))/sum(unitSessionLabels==ind & taskUnits')), unitSessionLabels);
percGTFig = modulatedUnitsPerRep({masterJoints},{graspTPercs},repNames,"% of task units","Grasp modulated (combined conditions)");
saveFigures(percGTFig, [saveDir, 'Representation Distribution\Grasp Percentage\'],[fileSaveName,'_PercOfTask'],[])



figTaskPercAvg = modulatedUnitsPerRep({masterJoints},{taskPercs},repNames,"", "Avg % of task units in a session per representation");
saveFigures(figTaskPercAvg, [saveDir, 'Representation Distribution\Task Percentage\'],[fileSaveName],[])
figTaskPercAvg = modulatedUnitsPerRep({masterJoints},{reachPercs},repNames,"", "Avg % of reach units in a session per representation");
saveFigures(figTaskPercAvg, [saveDir, 'Representation Distribution\Reach Percentage\'],[fileSaveName],[])
figTaskPercAvg = modulatedUnitsPerRep({masterJoints},{graspPercs},repNames,"", "Avg % of grasp units in a session per representation");
saveFigures(figTaskPercAvg, [saveDir, 'Representation Distribution\Grasp Percentage\'],[fileSaveName],[])

%% map units that modulated more for their respective condition
for s = 1:size(condCombs,1)
    % calculate percentage of total units that demonstrate significant
    % modulation between conditions
    condModCountsVS = arrayfun(@(ind) sum(condUnits{s}' & ...
        unitSessionLabels==ind), unitSessionLabels);
    condModPercsVS= arrayfun(@(ind) 100*(sum(condUnits{s}' & ...
        unitSessionLabels==ind))/sum(unitSessionLabels==ind),unitSessionLabels);
    figureCountVS =  mapUnitVals(verticies, vCells, vMask, siteMasks,...
        sessionInds{maxCondInd},condModCountsVS,cm);
    figurePercVS =  mapUnitVals(verticies, vCells, vMask, siteMasks,...
        sessionInds{maxCondInd},condModPercsVS,cm);
    % overlay imaging results onto figure
    activityIm = im2double(imresize(imread([sessionDirPrefix,...
        '\Mapping\tTests\Cond_Subtracted\',condAbbrevVS,'.png']),...
        [size(MMColorReference,1),size(MMColorReference,2)]));
    for cl = 1:length(condColor)
        activityIm = sqrt(sum((activityIm-permute(repmat([condColor{cl}].',1,...
            size(activityIm,2),size(activityIm,1)),[3 2 1])).^2,3));
        activityIm = (activityIm - min(activityIm(:)))/...
            (max(activityIm(:))-min(activityIm(:)))<= tol;
        activityIm = imdilate(activityIm, strel('square', lineWidth));
        boundaries = bwboundaries(activityIm);
        condAbbrevVS = strjoin(cellfun(@(cn) cn([0,regexp(cn, '\s')]+1),...
            condCombs(s,:), 'UniformOutput', false),'_');

        countIm = copyobj(gca(figureCountVS),figure()); hold on;
        cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth),...
            boundaries);
        percIM = copyobj(gca(figurePercVS),figure()); hold on;
        cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth),...
            boundaries);
        if(saveFigs)
            if(cl==1)
                saveFigures(figureCountVS,[saveDir  '\Count Maps\'...
                    condAbbrevVS,'\'], [fileSaveName,'_Count']);
                saveFigures(figurePercVS,[saveDir '\Count Maps\'...
                    condAbbrevVS,'\'], [fileSaveName,'_Perc']);
            end
            saveFigures(countIm,[saveDir, '\Count Maps\' ...
                condAbbrevVS,'\Imaging\'], [fileSaveName,'_Count']);
            saveFigures(percIM,[saveDir,'\Count Maps\' condAbbrevVS,...
                '\Imaging\'], [fileSaveName,'_Perc']);
        end
    end
end
for p = 1:length(phaseNames)

    jP =  plotJointPSTHS(bins,masterPSTH{d},masterSegs{d},simpRep,masterActivity{d},sessionInds);
    saveFigures(jP,['S:\Lab\ngc14\Figures\Physiology\Results\PSTH\',condAbbrev{d},'\'],saveName, []);
    masterPhaseIndices = cellfun(@(pIn) cellfun(@(p) round(nanmean(p,1)), pIn, 'UniformOutput', false), masterPhaseBins(d,p), 'UniformOutput',false);
    trialAUC = cellfun(@(rin, gin) cellfun(@(r,g) nanmean((r-g)),rin,gin), phaseFR(d,p), baselineFR(d),'UniformOutput', false);
    baseP = cellfun(@(psth) trapz(binSize,psth(baselineStartBin:baselineEndBin)./(baselineEndBin-baselineStartBin)), masterPSTH{d}, 'UniformOutput',  true);
    phaseP = cellfun(@(psth,r) trapz(binSize,psth(round(max(1,nanmean(r(:,1)))):round(min(length(bins),nanmean(r(:,end))))))./(diff(nanmean(r,1))),...
        masterPSTH{d}, masterPhaseBins{d,p},'UniformOutput', true);
    unitAUC = arrayfun(@(r,g) (r-g), phaseP, baseP, 'UniformOutput', true);

    trialAUC = cellfun(@(t) t.*tkInds, trialAUC, 'UniformOutput', false);
    unitAUC = unitAUC.*tkInds;

    repT = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
    repU = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
    s1 = []; s2 = [];
    for s = 1:length(repNames)
        figure(repT);
        s1(s) =subplot(1,4,s);
        hold on;
        title(repNames{s});
        bar(1:length(phaseNames), cellfun(@(p) nanmean(p(jointInds{s})), trialAUCS), 'FaceColor', repColors.(repNames{s}));
        errorbar(1:length(phaseNames), cellfun(@(p) nanmean(p(jointInds{s})), trialAUCS),...
            cellfun(@(p) nanstd(p(jointInds{s}))/sqrt(sum(~isnan(p(jointInds{s})))), trialAUCS), 'k','linestyle', 'none');
        xticks(1:length(phaseNames));
        xticklabels(phaseNames);

        figure(repU)
        s2(s) = subplot(1,4,s);
        hold on;
        title(repNames{s});
        bar(1:length(phaseNames), cellfun(@(p) nanmean(p(jointInds{s})), unitAUCS), 'FaceColor', repColors.(repNames{s}));
        errorbar(1:length(phaseNames), cellfun(@(p) nanmean(p(jointInds{s})), unitAUCS),...
            cellfun(@(p) nanstd(p(jointInds{s}))/sqrt(sum(~isnan(p(jointInds{s})))), unitAUCS), 'k','linestyle', 'none');
        xticks(1:length(phaseNames));
        xticklabels(phaseNames);
    end
    linkaxes(s1);
    linkaxes(s2);
    saveFigures(repT,['S:\Lab\ngc14\Figures\Physiology\Results\SI\ByRep\',condAbbrev{d},'\Trial\'],[saveName,'_Task_Trial_',condAbbrev{d}] ,[]);
    saveFigures(repU,['S:\Lab\ngc14\Figures\Physiology\Results\SI\ByRep\',condAbbrev{d},'\Unit\'],[saveName, '_Task_Unit_',condAbbrev{d}] ,[]);

    trialAUCFig = modulatedUnitsPerRep(repSave,trialAUC,string(['AUC']),string(sprintf(phaseNames{p})));
    unitAUCFig = modulatedUnitsPerRep(repSave, mat2cell(unitAUC,1),string(['AUC']),string(sprintf(phaseNames{p})));

    saveFigures(trialAUCFig,['S:\Lab\ngc14\Figures\Physiology\Results\AUC\',condAbbrev{d},'\Trial\'],[saveName,'_Task_Trial'] ,[]);
    saveFigures(unitAUCFig,['S:\Lab\ngc14\Figures\Physiology\Results\AUC\',condAbbrev{d},'\Unit\'],[saveName,'_Task_Unit'] ,[]);

    figMapUnit = mapUnitVals(verticies,vCells,vMask,siteMasks,sessionInds,unitAUC);
    hold on;
    weighting = cellfun(@(c) normalize(c(sessionInds), 'range'), unitAUCS, 'UniformOutput', false)
    centr =nansum((cell2mat(siteLocation(sessionInds)').*weighting{p}')/nansum(weighting{p}'),1);
    scatter(centr(1), centr(2),250,'b','filled','square');
    figMapTrial = mapUnitVals(verticies,vCells,vMask,siteMasks,sessionInds,trialAUC{1});
    saveFigures(figMapTrial,['S:\Lab\ngc14\Figures\Physiology\Results\Maps\AUC\',condAbbrev{d},'Trial\'],[saveName,'_', phaseNames{p},'_Task_Trial'] ,[]);
    saveFigures(figMapUnit,['S:\Lab\ngc14\Figures\Physiology\Results\Maps\AUC\',condAbbrev{d},'Unit\'],[saveName,'_',phaseNames{p} '_Task_Unit'] ,[]);
end
condModUnits = [];
condCompUnits = repmat({false(1,length(sessionInds))}, 1,size(condCombs,1));
for u = 1:length(sessionInds)
    currUnitFR = cellfun(@(p) p{u} ,taskFR, 'UniformOutput', false);
    maxTrials = max(cellfun(@length,currUnitFR));
    paddedUnitFR = cell2mat(cellfun(@(uFR) [uFR; ...
        NaN(maxTrials-length(uFR),1)],currUnitFR, 'UniformOutput', false));
    condModUnits(u) = anova1(paddedUnitFR,[],'off')<pVal;
    if(condModUnits(u))
        for s = 1:size(condCombs,1)
            dist1 = paddedUnitFR(:,find(strcmp(movementConds,condCombs{s,1})));
            dist2 = paddedUnitFR(:,find(strcmp(movementConds,condCombs{s,2})));
            condCompUnits{s}(u) = anova1([dist1,dist2],[],'off')<pVal/size(condCombs,1);
        end
    end
end
%% clustering
cluster_PSTHS(saveFigs, monkey, singleOrAll, sessionOrUnitPSTHS, ...
    jointName, jointClass, masterPSTH, activation)
unitPhaseVals = cellfun(@(a,b,c) nanmean(nanmean(a(:,c(1):c(2)),2)./...
    nanmean(a(:,b(1):b(2)),2)),condPSTHS,baselineSegs(condInd),currSegs);

%
cMapJointEnd = cellfun(@(jc) rgb2hsv(jointColors.(jc)),repNames,'UniformOutput', false);
cMapJoint = cellfun(@(jc) [repmat(jc(1),1,nBins); ...
    linspace(jc(2)/nBins,jc(2),nBins);...
    linspace(jc(end)-(jc(end)*(1/nBins)),jc(end),nBins)]', cMapJointEnd, 'UniformOutput', false);
cMapJoint(contains(repNames, {'Arm', 'Hand'})) = cellfun(@flipud,...
    cMapJoint(contains(repNames, {'Arm', 'Hand'})),'UniformOutput', false);
cMapJoint = cell2mat(cMapJoint);
cMapJoint = hsv2rgb(cMapJoint);
mappedInd = cellfun(@(jn) find(strcmp(repNames,jn))-1, masterJoints{maxCondInd});
condModCounts = condModCounts +(mappedInd*length(countRange));

%% emg/phys
%%
events = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}});
units = {};
muscleCorrs = {};
binSize=.01; % bin size in seconds
secondsBeforePSTHAlignmentPoint = -6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 7; % time (in seconds) after alignement point to create PSTH
inclusiveBins = secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values
bins = inclusiveBins(1:end-1);
muscles = dir('*.mat');
for c = 1:32
    l = load(['S:\Lab\Gilligan\All Data\Gilligan_04_22_2019\Physiology\Results\Gilligan_04_22_2019_',num2str(c),'.mat']);
    if(~isfield(l.sortedSpikeData, 'label'))
        labs = l.label;
    else
        labs = l.sortedSpikeData.label;
    end
    allGoodTrials = ~(cellfun(@(a,b) length(a) <=(b(end)-b(1)) ...
        | length(a)>200*(b(end)-b(1)) | any(isnan(a)), l.sortedSpikeData.SpikeTimes,l.sortedSpikeData.SegTimes));
    goodUnits = find(sum(allGoodTrials,2)>length(allGoodTrials)/4);
    units(end+1:end+length(find(strcmp(labs(goodUnits),"s"))),:) = l.sortedSpikeData.SpikeTimes(find(strcmp(labs(goodUnits),"s")),:);
end
%%
ard = l.sortedSpikeData.ArduinoData;
units = cellfun(@(u,t) u - t(2), units, repmat( l.sortedSpikeData.SegTimes,size(units,1),1),'UniformOutput', false);
units = cellfun(@(a) histcounts(a,inclusiveBins)./binSize, units, 'UniformOutput', false);
smoothHists = cellfun(@(a) conv(a, gausswin(.15/.01)/sum(gausswin(.15/.01)), 'same'), units, 'UniformOutput', false);
CSegs = cellfun(@(ar) values(events,{ar}), ard(:,1), 'UniformOutput', false);
restInds = cellfun(@(cc) strcmp(cc, 'Rest'), ard(:,1));
smoothHists = smoothHists(:,~restInds);
CSegs = CSegs(~restInds);
l.sortedSpikeData.SegTimes(~restInds) = cellfun(@(a) a-a(2), l.sortedSpikeData.SegTimes(~restInds)','UniformOutput', false);
tAlign = cellfun(@(ar,s) [findBins(bins,(s(strcmp('StartReach',ar{1})))-.5),...
    findBins(bins,s(strcmp('StartReplaceHold',ar{1})))], ...
    CSegs,l.sortedSpikeData.SegTimes(~restInds)', 'UniformOutput', false);
muscCorrs = {};
muscWin = {};
for m = 1:length(muscles)
    subplot(4,3,m); hold on; title(muscles(m).name(1:end-4));
    currMuscle = load(muscles(m).name);
    currMuscle.sortedEMGData.EMGData(212) = [];
    currMuscle.sortedEMGData.SegTimes(212) = [];
    mAlign = cellfun(@(em,s) [s(strcmp('StartReach', em{1})),...
        s(strcmp('StartReplaceHold',em{1}))], CSegs,currMuscle.sortedEMGData.SegTimes(~restInds)', 'UniformOutput', false);
    muscWin{m} = cellfun(@(e,w) e(w(1):w(end)), currMuscle.sortedEMGData.EMGData(~restInds),mAlign', 'UniformOutput', false);
    unitCorrs = {};
    for u = 1:size(units,1)
        unitCorrs{u} = cellfun(@(sh,ta,mu) corrcoef(interp1(1:length(mu),mu,...
            1:(ta(end)-ta(1)+1)),sh(ta(1):ta(end))),...
            smoothHists(u,:),tAlign',muscWin{m},'UniformOutput', false);
    end
    for c = 1 :length(l.sortedSpikeData.Conditions)-1
        condInds = cellfun(@(cc) strcmp(cc,l.sortedSpikeData.Conditions{c}), ard(~restInds,1));muscWin
        cellfun(@(sh, ta) plot(1:diff(ta)+1,sh(ta(1):ta(end)),'r'), smoothHists(u,condInds),tAlign(condInds)')
        cellfun(@(ta,mu) plot(interp1(1:length(mu),mu,1:(ta(end)-ta(1)+1)),'b'),tAlign(condInds)',muscWin{m}(condInds),'UniformOutput', false);
    end
    muscCorrs(m,:) = unitCorrs;
end



%% correlate units to EMG averages
uPhases = {};
unitCorrs = {};

for c = 1:length(movementConds)
    currCondSegs = values(events,movementConds(c));
    currCondSegs = currCondSegs{1};
    condMuscles = cellfun(@(m) m{c}, mPhases,'UniformOutput',false);
    avgUnitSegs = cellfun(@(s) nanmean(s,1), masterSegs{c}, 'UniformOutput',false);
    %     if(strcmp(movementConds{c}, 'Photocell'))
    %         phaseAlignments{end} = 'StartHold';
    %     else
    %         phaseAlignments{end} = 'StartLift';
    %     end
    uwin = cellfun(@(u) cellfun(@(s,w) findBins(bins,u(strcmp(s,currCondSegs)))...
        + w/binSize, phaseAlignments, phaseWindows, 'UniformOutput', false),...
        avgUnitSegs, 'UniformOutput',false);
    %uPhases{c} = cellfun(@(u, uwi) cellfun(@(p) u(p(1):p(end)) , uwi, ...
    %'UniformOutput', false), masterPSTH{c}, uwin, 'UniformOutput', false);
    uPhases{c} = cellfun(@(u,uwi) u(uwi{1}(1):uwi{end}(end)), masterPSTH{c}, uwin, 'Uniformoutput', false);
    minInterp = min(cellfun(@length,condMuscles));
    %     condMuscles = cellfun(@(cm) interpft(cm,minInterp), condMuscles, 'UniformOutput', false);
    %     unitCorrs{c} = cellfun(@(u) cellfun(@(m) corrcoef(interpft(m,length(u)),u),...
    %         condMuscles, 'UniformOutput', false), uPhases{c}, 'UniformOutput', false);
    %     unitCorrs{c} = cellfun(@(ur) cellfun(@(r) r(2), ur, 'UniformOutput', true),unitCorrs{c},'UniformOutput', false);
    %     for m = 1:length(modelNames)
    %         currMuscleCorrs = cellfun(@(u) u(m), unitCorrs{c}, 'UniformOutput', true);
    %         figMapMCorr = mapUnitVals(vXY,vMask,siteMasks,sessionInds,currMuscleCorrs,true,10,[-1,1]);
    %         saveFigures(figMapMCorr,[saveDirPath,'Maps\Aligned\'],[saveName,'_',modelNames{m}],[]);
    %         corrDistsFig = modulatedUnitsPerRep(repSave,currMuscleCorrs,{'Task'},string(modelNames(m)),repColors);
    %         saveFigures(corrDistsFig,[saveDirPath,'Aligned\'],[saveName,'_',modelNames{m}],[]);
    %     end
    %% EMG PSTH.corr
    allAvgM = {};
    for f = 1:numUnits
        for m = 1:7
            figure('Units','normalized','Position',[0 0 1 1])
            hold on
            for c = 1:3
                subplot(2,2,c)
                hold on;
                interpTo =taskFR{c}{f}(condInds{c});
                interpFrom =  mTrials{m}{c};
                largestSeg = max(cellfun(@length,interpTo));
                interpTrials = cellfun(@(t,f) interpft(f,length(t)),...
                    interpTo,interpFrom,'UniformOutput', false);
                un = nanmean(cell2mat(cellfun(@(b,ii) nanPadArrays(...
                    interpft(b,length(ii)),largestSeg), taskFR{c}{:,wM(c)}(condInds{c}),interpTrials,'UniformOutput', false)'),1);
                allAvgM{m} = nanmean(cell2mat(cellfun(@(b,ii) nanPadArrays(...
                    interpft(b,length(ii)),largestSeg), mTrials{m}{c},interpTrials,'UniformOutput', false)'),1);


                yyaxis right
                cellfun(@(ts) plot(FRBins(1:length(ts)),ts, 'Color', [randi([50 100],[1,1]),...
                    randi([0 50],[1,2])]./100,'LineWidth', .5,'LineStyle',...
                    '--','Marker','none'),interpTo);
                avgS = cell2mat(cellfun(@(t) nanPadArrays(t,largestSeg),interpTo,'UniformOutput', false)');
                plot(FRBins(1:length(avgS)),nanmean(avgS,1),'Color',[.5 0 0],'LineWidth',2,'LineStyle', '-','Marker','none');

                yyaxis left
                cellfun(@(t) plot(FRBins(1:length(t)),t,'Color',[randi([0 50],[1,2]),...
                    randi([50 100],[1,1])]./100,'LineWidth',.5,'LineStyle', '--',...
                    'Marker','none'),interpTrials);
                avgI = cell2mat(cellfun(@(t) nanPadArrays(t,largestSeg),interpTrials,'UniformOutput',false)');
                plot(FRBins(1:length(avgI)),nanmean(avgI,1),'Color',[0 0 .5],'LineStyle','-','LineWidth',2,'Marker','none');

                corrs = cellfun(@(t,u) corrcoef(t,u),interpTrials,...
                    interpTo,'UniformOutput',false);
                allCorrs{c}{f,m} = cellfun(@(r) r(2), corrs);
                avgCorrHold = corrcoef(nanmean(avgS,1),nanmean(avgI,1));
                avgCorr{c}{f,m} = avgCorrHold(2);
                title([string(movementConds{c}); strcat("Avg trial correlation = ", num2str(nanmean(allCorrs{c}{f,m}),2));...
                    strcat("Avg correlation = ",num2str(avgCorr{c}{f,m},2))]);
            end
            if(~exist(['Trials\',modelNames{m},'\'],'dir'))
                mkdir(['Trials\',modelNames{m},'\'])
            end
            saveas(gcf,['Trials\',modelNames{m},'\',num2str(f),'.png'])
            close all;
        end
        figure('Units','normalized','Position',[0 0 1 1])
        hold on
        [wmV,wM] = cellfun(@(ac) max(cellfun(@(m) max(m),ac(f,:),'UniformOutput',true)), allCorrs,'UniformOutput',true);
        [wmV2,wM2] = cellfun(@(av) max(abs([av{f,:}])), avgCorr,'UniformOutput',true);
        for c = 1:3
            subplot(2,2,c); hold on;

            title([strcat("Best muscle (trial corr = ", num2str(wmV(c),2),"): ",...
                string(modelNames{wM(c)})); strcat("Best muscle (avg corr = ", ...
                num2str(wmV2(c),2), "): ",string(modelNames{wM2(c)}))]);
            yyaxis left
            plot(FRBins(1:length(un)),un, 'LineStyle','-','LineWidth',2,'Color','k','Marker','none');
            yyaxis right;
            for m = 1:7
                plot(FRBins(1:length(allAvgM{m})),allAvgM{m},'LineStyle','--','LineWidth',1,'Color',musclesColors(m,:),'Marker','none')
            end
            plot(FRBins(1:length(allAvgM{wM(c)})),allAvgM{wM(c)},'LineStyle','-','LineWidth',2,'Color',musclesColors(wM(c),:),'Marker','none')
            plot(FRBins(1:length(allAvgM{wM2(c)})),allAvgM{wM2(c)},'LineStyle','-','LineWidth',2,'Color',musclesColors(wM2(c),:),'Marker','none')

        end
        if(~exist(['Trials\Winning Muscle\'],'dir'))
            mkdir(['Trials\Winning Muscle\'])
        end

    end
    uPhases{c} = cellfun(@(u, uwi) cellfun(@(p) u(p(1):p(end)) , uwi, ...
        'UniformOutput', false), masterPSTH{c}, uwin, 'UniformOutput', false);
end

%% plot unit psths
%% plot joint PSTHS

close all
events = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}});

alignmentPoint = {{'StartReach', 'StartLift', 'StartWithdraw'}, {'StartReach', 'StartLift', 'StartWithdraw'},...
    {'StartReach', 'StartHold', 'StartWithdraw'}, {'GoSignal','StartReplaceSuccess', 'StartReward'}};
saveDirPath = 'S:\Lab\ngc14\Figures\Physiology\Results\Phase\';
saveName = [monkey,'_','Unit','_','Single','_Aligned'];
alignLimits = {[-.6,.15],[-.6, .25], [-.25, .6],[-.15,.6]};
alignLimits = repmat({[-3,3]},1,4)
zeroBinInd = find(bins==0);
alignmentGap = 25

condAlignment =  cellfun(@(a) find(strcmp(events(conditions{c}),a)), alignmentPoint{c});
PSTH = cellfun(@(a,w) a(:,(fix(w(1)/binSize)+zeroBinInd):...
    fix((w(end)/binSize)+zeroBinInd),:),unitPSTHS, alignLimits, 'UniformOutput', false);
if(strcmp(monkey,'Gilligan'))
    colorRange = [1 30];
    FRLim = [5 20]
else
    colorRange = [1 30];
    FRLim = [15 35];
end
allSegs = masterSegs{c};
for u = 1:size(PSTH{1},1)
    figure(); hold on;
    plotStart = 0;
    plotted = false(1,9);
    xAlignTicks = {};
    for a = 1:length(PSTH)
        xAlignTicks{a} = plotStart+[1:size(PSTH{a},2)];
        plot(xAlignTicks{a},nanmean(squeeze(PSTH{a}(u,:,:)),2), 'LineWidth',2,'Color', 'b')
        ue = nanmean(squeeze(PSTH{a}(u,:,:)),2) + (nanstd(squeeze(PSTH{a}(u,:,:)),0,2)/sqrt(size(PSTH{a},3)));
        le = nanmean(squeeze(PSTH{a}(u,:,:)),2) - (nanstd(squeeze(PSTH{a}(u,:,:)),0,2)/sqrt(size(PSTH{a},3)));
        yP = [le',fliplr(ue')];
        xP = [xAlignTicks{a},fliplr(xAlignTicks{a})];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        d = patch(xP,yP,1);
        set(d,'edgecolor','none','facealpha',.5,'facecolor','r');
        avgSegs = nanmean(allSegs{a}{u},1);
        for s = 1:length(avgSegs)
            if(avgSegs(s)>=alignLimits{a}(1) && ...
                    avgSegs(s)<=alignLimits{a}(end) && ...
                    (~plotted(s) ||avgSegs(s)==plotStart))
                if(avgSegs(s)==0)
                    plotColor = 'k';
                else
                    plotColor = 'r';
                end
                %plotted(s) = true;
                pSeg = find(isalmost(alignLimits{a}(1):binSize:...
                    alignLimits{a}(end),avgSegs(s),binSize/1.99),1);
                plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 FRLim(end)],'Color',plotColor,'LineStyle','--')
            end
        end
        if(a==length(PSTH))
            allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs([pd(1):binSize:pd(end)])...
                ==min(abs([pd(1):binSize:pd(end)]))), ta(end)], xAlignTicks, ...
                alignLimits, 'UniformOutput', false);
            xticks([allXTicks{:}]);
            allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                "0"; num2str(pd(end),'%.2f')],alignLimits, 'UniformOutput', false);
            xticklabels([allLabels{:}]);
        end
        plotStart = plotStart + size(PSTH{a},2) + alignmentGap;
    end
    saveFigures(gcf,['\\univ.pitt.edu\nb\Gharbawie\Lab\ngc14\Figures\Physiology\Results\Single PSTHS\05_17_2019\',condAbbrev{c},'\'], num2str(u),[]);
    if(u==size(PSTH{1},1))

        figure(); hold on;
        plotStart = 0;
        plotted = false(1,9);
        xAlignTicks = {};
        for a = 1:length(PSTH)
            xAlignTicks{a} = plotStart+[1:size(PSTH{a},2)];
            plot(xAlignTicks{a},nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1), 'LineWidth',2,'Color', 'b')
            ue = squeeze(nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1))+  squeeze(nanmean(nanstd(squeeze(PSTH{a}(:,:,:)),0,3)/sqrt(size(PSTH{a},3)),1));
            le = squeeze(nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1)) - squeeze(nanmean(nanstd(squeeze(PSTH{a}(:,:,:)),0,3)/sqrt(size(PSTH{a},3)),1));
            yP = [le,fliplr(ue)];
            xP = [xAlignTicks{a},fliplr(xAlignTicks{a})];
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            d = patch(xP,yP,1);
            set(d,'edgecolor','none','facealpha',.5,'facecolor','r');
            avgSegs = nanmean(allSegs{a}{u},1);
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=alignLimits{a}(1) && ...
                        avgSegs(s)<=alignLimits{a}(end) && ...
                        (~plotted(s) ||avgSegs(s)==plotStart))
                    if(avgSegs(s)==0)
                        plotColor = 'k';
                    else
                        plotColor = 'r';
                    end
                    %plotted(s) = true;
                    pSeg = find(isalmost(alignLimits{a}(1):binSize:...
                        alignLimits{a}(end),avgSegs(s),binSize/1.99),1);
                    plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 FRLim(end)],'Color',plotColor,'LineStyle','--')
                end
            end
            if(a==length(PSTH))
                allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs([pd(1):binSize:pd(end)])...
                    ==min(abs([pd(1):binSize:pd(end)]))), ta(end)], xAlignTicks, ...
                    alignLimits, 'UniformOutput', false);
                xticks([allXTicks{:}]);
                allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                    "0"; num2str(pd(end),'%.2f')],alignLimits, 'UniformOutput', false);
                xticklabels([allLabels{:}]);
            end
            plotStart = plotStart + size(PSTH{a},2) + alignmentGap;
        end
    end
    saveFigures(gcf,['\\univ.pitt.edu\nb\Gharbawie\Lab\ngc14\Figures\Physiology\Results\Single PSTHS\05_17_2019\',condAbbrev{c},'\'], 'Session',[]);

end

% figMapG = mapUnitVals(vXY,vMask,vMask,sessionInds,int32(unitPhaseInds{2}.*tks),true,5,[1 10]);
% saveFigures(figMapG,[saveDirPath,'Maps\Count\',condAbbrev{c},'\'],[saveName,'_Grasp_units'],[]);
% figMapRGx = mapUnitVals(vXY,vMask,vMask,sessionInds,int32(unitPhaseInds{3}.*tks),true,5,[-5 5]);
% saveFigures(figMapRGx,[saveDirPath,'Maps\Count\',condAbbrev{c},'\'],[saveName,'_Go_units'],[]);
%saveFigures(gcf,[saveDirPath,condAbbrev{c},'\'],saveName,[]);
close all;
%encodingPhase{c} = phaseAssign;
%
% %
% % figMapGUnit = mapUnitVals(vXY,vMask,siteMasks,sessionInds,(encodingPhase{c}==2),true,5,[1,7]);
% % unitGPhase = modulatedUnitsPerRep(repSave,{(encodingPhase{c}==2)},["G units"],"G");









%                         [~, phaseFRTask, baseline] = ...
%                             calcFRChanges(PSTHCalc,[repmat(baselineStartBin,...
%                             sum(condInds),1),repmat(baselineEndBin,...
%                             sum(condInds),1)],[cellfun(@(a) ...
%                             findBins(bins,a(taskStartInd)+taskTimeRange(1)),...
%                             alignTimes(condInds)); cellfun(@(a) ...
%                             findBins(bins,a(taskEndInd)+taskTimeRange(end)),...
%                             alignTimes(condInds))]');

%
%                         baseline = cellfun(@(pc,bc,be) cell2mat(cellfun(@(p,s,e) ...
%                             (trapz(p(:,s:e),2)),pc,bc,be,...
%                             'UniformOutput', false)),PSTHCalc,baselineStart,...
%                             baselineEnd,'UniformOutput', false);

%                         phaseFRTask = mat2cell(mat2cell(...
%                             [phaseFRTask{:}],ones(1,numUnits),...
%                             cellfun(@(e,s) e+1-s, pE,pS)),ones(1,numUnits));



%
%     if(c<4)
%         rSpeed =cellfun(@(p) cellfun(@(a) 1./(abs(diff([a(:,reachAlignment)],1,2))), ...
%             p, 'UniformOutput', false), masterSegs{c}, 'UniformOutput', false);
%         rFR = cellfun(@(t) t, reachFR{c},'UniformOutput', false);
%         allTrialReps = vertcat(cellfun(@(r,s) repmat(string(r),length(s),1), simpRep,rSpeed{2},'UniformOutput',false));
%         allTrialReps = vertcat(allTrialReps{:});
%         [unitB,~,unitR,~,unitStats] = cellfun(@(us,uf) regress(uf', [ones(length(us),1), vertcat(us)]), rSpeed{2},reachFR{c}{2}', 'UniformOutput', false);
%         unitCorrs = cellfun(@(us,uf) corrcoef(us, uf'), rSpeed{2}, rFR{2}', 'UniformOutput', false);
%         unitCorrs = cellfun(@(r) r(2), unitCorrs, 'UniformOutput', false);
%         unitStats =cellfun(@(s) s(1), unitStats, 'UniformOutput', true);
%
%         scatter(vertcat(rSpeed{2}{:}), [rFR{2}{:}]',36,cell2mat(cellfun(@(s) repColors.(s), allTrialReps, 'Uniformoutput', false)),'filled');
%         hold on;
%         ylabel('FR');
%         xlabel('Speed');
%         allCorr = corrcoef(vertcat(rSpeed{2}{:}), [rFR{2}{:}]');
%         title(['All units: r^2 = ', num2str(allCorr(2),3), '; mean r^2 = ',...
%             num2str(nanmean(unitStats),3)]);
%
%         pAll = polyfit(vertcat(rSpeed{2}{:}),[(rFR{2}{:})]',1);
%         fAll = polyval(pAll,vertcat(rSpeed{2}{:}));
%         plot(vertcat(rSpeed{2}{:}),fAll,'k-')
%         saveFigures(gcf,[saveDirPath,'Reach_Speed\',condAbbrev{c},'\'],'All_Units',[]);
%         close all;
%         for j = 1:length(jointName)
%             figure(); hold on;
%             currJ = cellfun(@(s) strcmp(s, jointName{j}), simpRep, 'UniformOutput', true);
%             scatter(vertcat(rSpeed{2}{currJ}), [rFR{2}{currJ}]',36,repColors.(jointName{j}),'filled');
%             [b,~,r,~,stats] = regress([rFR{2}{currJ}]', ...
%                 [ones(sum(cellfun(@length,rSpeed{2}(currJ))),1), vertcat(rSpeed{2}{currJ})])
%             ylabel('FR');
%             xlabel('Speed');
%             title([jointName{j}, ' units: r^2 = ', num2str(stats(1),3),...
%                 '; mean r^2 = ',num2str(nanmean(unitStats(currJ)),3), '; corrcoef: = ', ...
%                 num2str(nanmean(unitStats(currJ)),2)]);
%             p = polyfit(vertcat(rSpeed{2}{currJ}),[(rFR{2}{currJ})]',1);
%             f = polyval(p,vertcat(rSpeed{2}{currJ}));
%             plot(vertcat(rSpeed{2}{currJ}),f,'k-')
%             saveFigures(gcf,[saveDirPath,'Reach_Speed\',condAbbrev{c},'\'],[jointName{j},'_Units'],[]);
%             close all
%         end
%         unitStats(~tks) = NaN;
%         rFig = modulatedUnitsPerRep(currRep,unitStats,...
%             "r-squared", "Task Units", repColors);
%         rMapFig = mapUnitVals(vXY,vMask,siteMasks,sessionInds,...
%             unitStats,false,5,[]);
%         saveFigures(rFig, [saveDirPath,'\Reach_Speed\', condAbbrev{c},'\'],[saveName],[])
%         saveFigures(rMapFig,[saveDirPath,'Maps\Reach_Speed\',condAbbrev{c},'\'],...
%             [saveName],[]);
%     end
%
%
%     figure('Units','normalized', 'Position',[0 0 1 1]);
%     hold on;
%     raSpeed = (vertcat(rSpeed{2}{:}));
%     raFR = horzcat(rFR{2}{:})';
%     for j = 1:length(jointName)
%             jointInds = cell2mat(cellfun(@(a,s) repmat(strcmp(a,jointName{j}), size(s,1),1), simpRep,masterSegs{1}{2},'UniformOutput', false)');
%             subplot(2,2,j);hold on;
%             scatter(raSpeed(jointInds), raFR(jointInds),36,'MarkerEdgeColor',...
%                 repColors.(jointName{j}));
%             rVal = corr(raSpeed(jointInds), raFR(jointInds));
%
%     p = polyfit(raSpeed,raFR,1);
%     f = polyval(p,raSpeed);
%      plot(raSpeed,f,'k-')
%
%             title([jointName{j}, ': r = ', num2str(rVal,3)]);
%     end
%
%
%




%     condEvents = eventsAll(conditions{c});
%     if(strcmp(conditions{c}, 'Photocell'))
%         condPhaseAlignments = {{'GoSignal', 'StartHold'},{'GoSignal'},{'StartReach'},{'StartHold'}};
%         condPhaseWindows = phaseWindows;
%     elseif(strcmp(conditions{c}, 'Rest'))
%         condPhaseAlignments = {{'GoSignal', 'StartReplaceHold'}, {'GoSignal'}};
%         condPhaseWindows = {[0 0], [-1 1]};
%     else
%         condPhaseAlignments = phaseAlignments;
%         condPhaseWindows = phaseWindows;
%     end
%
%     baselineAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
%         baselineAlignments,'UniformOutput', true);
%     psthPhaseAlignment = cellfun(@(pa) cellfun(@(a) find(strcmp(condEvents,a)),...
%         pa,'UniformOutput',true),condPhaseAlignments, 'UniformOutput', false);
%     % find start and end index in (1) each PSTH for (2) each phase
%     psthPhaseEnds = cellfun(@(at) cellfun(@(a) cellfun(@(pa,pw) num2cell(...
%         findBins(a(:,pa)+pw,bins),2),psthPhaseAlignment,condPhaseWindows,'UniformOutput',false),at,'UniformOutput', false),...
%         masterSegs{c},'UniformOutput', false);
%     % create mask that NaNs out all values outside of interrogated phase
%     % for each (1) PSTH for each (2) phase
%     psthPhases = cellfun(@(pr) cellfun(@(ps) cellfun(@(p) ...
%         NaNMaskBins(p,length(bins)),ps, 'UniformOutput',false),...
%         pr,'UniformOutput', false), psthPhaseEnds, 'UniformOutput', false);
%     % create equivalent sized start and end points in the baseline window
%     % for each (1) PSTH for each (2) phase
%     psthBaselineEnds = cellfun(@(at,pp) cellfun(@(a,ps) cellfun(@(p) num2cell(...
%         max(ones(size(p,1),2), findBins(a(:,baselineAlignment)+baselineWindow,bins)...
%         -[zeros(size(p,1),1),nansum(p,2)]),2),ps,'UniformOutput',false),at,pp,'UniformOutput', false),...
%         masterSegs{c},psthPhases, 'UniformOutput', false);
%
%     baselinePSTHSampAUC{c} = cellfun(@(tb,ph,pc) cellfun(@(pb,hh,pt) cellfun(@(p,w)...
%         AUCBaselineBootstrap(p,pt,w),pb,hh, 'UniformOutput', false),tb,ph,pc,...
%         'UniformOutput', false),psthBaselineEnds,psthPhases,masterTrialPSTH{c},'UniformOutput',false);
%
%     psthPhasesAUC{c} = cellfun(@(up,ua) cellfun(@(ps,p) cellfun(@(i) permute(trapz(p.*...
%         ~isnan(permute(repmat(i,[1,1,size(p,1)]),[3 2 1])),2),[1 3 2]),ps,'UniformOutput', false),up,ua,'UniformOutput',...
%         false),psthPhases,masterTrialPSTH{c},'UniformOutput', false);
%%
function bin = findBins(allBins, timePoint)
binSize = mode(diff(allBins));
if(isempty(timePoint))
    bin = NaN;
else
    bin = find(isalmost(allBins,timePoint,binSize/1.99),1);
end
if(isempty(bin))
    bin = NaN;
end
end