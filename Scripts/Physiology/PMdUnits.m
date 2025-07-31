conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions(1:end-1),{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]}});
taskWindow = {[-0.3, 0]};
alignLimits = {[-.75, 2.5]};
pVal=0.05;
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\";
monkey = "Gilligan";
MIN_BLOCKS_FOR_UNIT = 13;
params = PhysRecording(string(conditions),.01,.15,-5,7,containers.Map(conditions,...
    {["GoSignal"],["GoSignal"],["GoSignal"],"GoSignal"}));
plotUnits = false;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps] = getAllSessions(params,"Single","PMd");
trialFR = cellfun(@(ct,cs,ta) cellfun(@(a,b) cell2mat(cellfun(@(m,tt) ...
    squeeze(mean(m(:,max(1,tt(1)):max(1,tt(end)),:),2,'omitnan').*(1*~all(isnan(tt)))),...
    squeeze(num2cell(a{1},[1,2])),cellfun(@(e) [cell2mat(e),repmat([NaN,NaN],isempty(cell2mat(e)),1)],num2cell(...
    [arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(1)),'UniformOutput',false),...
    arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(end)),'UniformOutput',false)],...
    2),'UniformOutput',false),'UniformOutput',false)'),ct,cs,'UniformOutput',false), siteTrialPSTHS(1:end-1), siteSegs(1:end-1),...
    cellfun(@(c,t) arrayfun(@(e) find(strcmp(c,e)),t{1}), params.condSegMap.values(taskAlign.keys),taskAlign.values,...
    'UniformOutput',false),'UniformOutput',false);
goodFR = cellfun(@(c) cellfun(@(s) s>2 & s<200,c,'UniformOutput',false),trialFR,'UniformOutput',false);
goodUnits = cellfun(@(tn) cell2mat(cellfun(@(s)sum(s,2), tn,'UniformOutput',false)),...
    num2cell(cat(2,goodFR{:}),2),'UniformOutput',false);
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs(1:end-1),siteTrialPSTHS(1:end-1),false,true);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,1),'omitnan')-3,params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) ...
    permute(mean(c(:,s:s+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),a(~isnan(n)),...
    num2cell(n(~isnan(n))),'UniformOutput',false),[1,1,sum(~isnan(n))])),3,'omitnan'));NaN(all(isnan(n)).*size(a{1},1),1)],p,t,...
    'UniformOutput',false),siteTrialPSTHS,allCondSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp{:}),repmat(nb,1,size(vertcat(cp{:}),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
[tVals,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
[~,restFR] = calculatePhases(params,containers.Map(["Rest"],{["GoSignal","StartReplaceHold"]}),taskWindow,siteSegs(end),siteTrialPSTHS(end),false,true);
[rVals,rUnit] = cellfun(@(tc) cellfun(@(r,c) ttestTrials(r,c,1,false,pVal),...
    restFR{1},tc,'UniformOutput',false),taskFR,'UniformOutput',false);
restUnits = cellfun(@(d) any(vertcat(d{:}),1)', num2cell([rUnit{:}],2), 'UniformOutput', false);
allMaps = cellfun(@(d) d{end}, chMaps, 'UniformOutput',false);
taskUnits = cellfun(@(a,b) any(cell2mat(a),2) & sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2), ...
    num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false);
unitChannels = cellfun(@(sc,am) arrayfun(@(s) find(s==am),sc),siteChannels{1},allMaps','UniformOutput',false);
unit2SiteMap=cell2mat(arrayfun(@(m,n) ones(1,m)*n,cellfun(@(s)size(s,2),unitChannels)',1:length(unitChannels),'UniformOutput',false));
%%
maxClusters = 10;
for c =1:length(conditions)
    tUnits = cell2mat(taskUnits);
    taskSiteInds = find(cellfun(@any,arrayfun(@(a) tUnits(unit2SiteMap==a),...
        min(unit2SiteMap):max(unit2SiteMap),'UniformOutput', false))'); %[19,20,21];
    tUnits = taskUnits(taskSiteInds);
    maxUnitFR = cell2mat(cellfun(@(m) cell2mat(m(taskSiteInds)), maxCondsFR, 'UniformOutput',false));
    unitPSTHS = cell2mat(cellfun(@(m,i) (i./i).*mean(m,3,'omitnan'),vertcat(normPSTH{c}{taskSiteInds}),tUnits,'UniformOutput',false));
    unitPSTHS =  unitPSTHS./maxUnitFR(:,c);
    clusterP = evalclusters(unitPSTHS,'kmeans','CalinskiHarabasz','KList',1:maxClusters);
    close all;
    plot(1:maxClusters,clusterP.CriterionValues);
    saveFigures(gcf,savePath+params.condAbbrev(params.condNames(c))+"\","EvaluationValues",[]);
    for s = 1:length(clusterP.InspectedK)
        close all;
        sCluster{s} = kmeans(unitPSTHS,s);
        clusterGroups = arrayfun(@(f) "A_"+num2str(f),sCluster{s});
        siteUnitSegs = cellfun(@(si,tu) repmat(mean(si{1},1,'omitnan'),length(tu),1),...
            siteSegs{c}(taskSiteInds),tUnits,'UniformOutput',false);
        clusterColors = cell2struct(num2cell(distinguishable_colors(s),2),...
            unique(clusterGroups(~cellfun(@(cs) contains(cs,"NaN"),clusterGroups))));
        plotJointPSTHS(params.bins,{unitPSTHS},{cell2mat(siteUnitSegs)},...
            clusterGroups,cell2mat(tUnits)',[],alignLimits,[0,1],clusterColors);
        if(clusterP.OptimalK==s)
            saveFigures(gcf,savePath+params.condAbbrev(params.condNames(c))+"\","Optimal"+num2str(s)+"_PSTH",[]);
        else
            saveFigures(gcf,savePath+params.condAbbrev(params.condNames(c))+"\",num2str(s)+"_PSTH",[]);
        end
    end
end
%%
[allTrials,allPSTHS]= deal(cell(1,1));
[condInds,allSiteInds,allRestInds] = deal([]);
maxCondsFR = cellfun(@(c) cellfun(@(d) max(mean(d{1},3,'omitnan'),[],2,'omitnan'),...
    c,'UniformOutput',false),siteTrialPSTHS, 'UniformOutput',false);
for c =1:length(conditions)
    close all;
    tUnits = cell2mat(taskUnits);
    [~,siteUnitMods] = unique(unit2SiteMap);
    taskSiteInds = find(cellfun(@any,arrayfun(@(a) tUnits(unit2SiteMap==a),...
        min(unit2SiteMap):max(unit2SiteMap),'UniformOutput', false))'); %[19,20,21];
    tUnits = taskUnits(taskSiteInds);
    rUnits = restUnits(taskSiteInds);
    siteDates = siteDateMap(taskSiteInds,:);
    siteUnitMods = siteUnitMods(taskSiteInds);
    unitChannelMaps = unitChannels(taskSiteInds);
    siteCondSegs = vertcat(siteSegs{c}{taskSiteInds});
    unitLocation = mapSites2Units(cellfun(@length,unitChannelMaps),siteLocation(taskSiteInds));
    siteUnitNo = mapSites2Units(cellfun(@length,unitChannelMaps),[siteDateMap{taskSiteInds,'Site'}]);
    siteUnitNames = "SiteNo"+arrayfun(@num2str,siteUnitNo,'UniformOutput',false);
    maxUnitFR = cell2mat(cellfun(@(m) cell2mat(m(taskSiteInds)), maxCondsFR, 'UniformOutput',false));
    if(contains(conditions(c),'sphere','IgnoreCase',true))
        trialSegs = cellfun(@(t) [t(:,1:2), NaN(size(t,1),double(size(t,2)==8)), t(:,3:end)], siteCondSegs,'UniformOutput',false);
    elseif(strcmp(conditions(c),"Photocell"))
        trialSegs = cellfun(@(t) [t(:,1:3), NaN(size(t,1),double(size(t,2)==8)), t(:,4:end)], siteCondSegs,'UniformOutput',false);
    else
        trialSegs = cellfun(@(t) [t(:,1:4), NaN(size(t,1),5*double(size(t,2)==4)), t(:,5:end)],siteCondSegs,'UniformOutput',false);
    end
    siteUnitSegs = cellfun(@(si,tu) repmat(mean(si,1,'omitnan'),length(tu),1),siteCondSegs,tUnits,'UniformOutput',false);
    unitLocation = unitLocation-min(unitLocation(:,:)).*ImagingParameters.px2mm;

    unitPSTHS = cellfun(@(m,i) (i./i).*mean(m,3,'omitnan'),vertcat(normPSTH{c}{taskSiteInds}),tUnits,'UniformOutput',false);
    sitePSTHS = cellfun(@(n,ti) cellfun(@(t) squeeze(t)',num2cell(n(ti,:,:),[2,3]),'UniformOutput',false),...
        vertcat(normPSTH{c}{taskSiteInds}),tUnits,'Uniformoutput',false);
    if(plotUnits)
        plotJointPSTHS(params.bins,{vertcat(unitPSTHS{:})},{cell2mat(siteUnitSegs)},...
            siteUnitNames,cell2mat(tUnits)', [],alignLimits,[0 25],cell2struct(num2cell(...
            distinguishable_colors(length(siteUnitSegs)),2),arrayfun(@(t) "SiteNo"+num2str(siteDateMap{t,'Site'}),taskSiteInds)));
        saveFigures(gcf,strcat(savePath,"PSTHS\Session_PSTHS\"),strcat("Task_Units_",params.condAbbrev(char(conditions(c))),"_"),[]);
        for s = 1:length(sitePSTHS)
            ucs = unitChannelMaps{s}(tUnits{s});
            unitIndex = histcounts(ucs,[unique(ucs),max(ucs)+1]);
            utc = arrayfun(@(aa,bb) [string(aa),string(arrayfun(@(bi) string([num2str(aa),'_',num2str(bi)]),2:bb,'UniformOutput',false))],...
                unique(ucs),unitIndex,'UniformOutput',false);
            unitLabels = string(arrayfun(@(t) "SiteNo"+num2str(siteDates{s,'Site'})+"_"+t,cell2mat(utc),'UniformOutput',false))';
            plotColors =cell2struct(num2cell(distinguishable_colors(length(ucs)),2),unitLabels);
            unitTrialLabels = cell2mat(arrayfun(@(n) repmat(n,size(siteCondSegs{s},1),1),unitLabels,'UniformOutput',false));
            siteTrialSegs = repmat(siteCondSegs(s),length(unitLabels),1);
            plotJointPSTHS(params.bins,{cell2mat(sitePSTHS{s})},{cell2mat(siteTrialSegs)},...
                unitTrialLabels,true(fliplr(size(unitTrialLabels))), [],alignLimits,[0 25],plotColors);
            saveFigures(gcf,savePath+"PSTHS\Units\"+string(datetime(siteDates.Date{s},'Format','MMMM_dd_yyyy'))+"\",...
                params.condAbbrev(params.condNames(c))+"_PSTH",[]);
            close all;
        end
    end
    taskPSTHS = cellfun(@(c,t) c(t,:), unitPSTHS,tUnits, 'UniformOutput',false);
    for t = 1:2
        if(t==1)
            typeName = "EqRest";
        else
            typeName = "DiffRest";
        end
        condGroupUnits = cellfun(@(u,d) d(u)==(t-1), tUnits, rUnits,'UniformOutput',false);
        unitAvgPSTHS = cellfun(@(m,p) (m./m).*p,condGroupUnits,taskPSTHS,'UniformOutput',false);
        if(plotUnits)
            plotJointPSTHS(params.bins,{vertcat(unitAvgPSTHS{:})},{cell2mat(cellfun(@(s,u)repmat(mean(s,1,'omitnan'),length(u),1),...
                siteCondSegs,condGroupUnits,'UniformOutput',false))},siteUnitNames(vertcat(tUnits{:})),...
                cell2mat(condGroupUnits)', [],alignLimits,[0 15],cell2struct(num2cell(...
                distinguishable_colors(length(condGroupUnits)),2),"SiteNo"+num2str(siteDates{:,'Site'})));
            saveFigures(gcf,strcat(savePath,"PSTHS\Session_PSTHS\"),strcat(params.condAbbrev(params.condNames(c)),"_",typeName),[]);
        end
    end
    allSiteInds = [allSiteInds;tUnits];
    allRestInds = [allRestInds;cellfun(@(u,d) d(u)==0, tUnits, rUnits,'UniformOutput',false)];
    condInds = [condInds,repmat(string(params.condAbbrev(params.condNames(c))),1,length(siteUnitNames))];
    allPSTHS = cellfun(@(a,b) vertcat(a,b), allPSTHS, {unitPSTHS}, 'UniformOutput', false);
    allTrials = cellfun(@(a,d) vertcat(a,cellfun(@(m) vertcat(m{:}),d,'UniformOutput',false)),...
        allTrials,{cellfun(@(su,tu) repmat({mean(su,1,'omitnan')},length(tu),1),trialSegs,tUnits,...
        'UniformOutput',false)},'UniformOutput',false);
    
    siteUnitInds = siteUnitNo;
    siteUnitInds(~ismember(1:length(siteUnitNo),siteUnitMods)) = NaN;
    siteUnitNames(~cell2mat(tUnits)) = "";
    allUnitSegs = vertcat(siteUnitSegs{:});

    heatmap_distribution_plots(siteUnitNames,cell2mat(tUnits)',(siteUnitInds'),...
    {allUnitSegs},params.bins,{cell2mat(unitPSTHS)},maxUnitFR(:,c),unitLocation,...
    alignLimits,strcat(savePath,"Heatmaps\FR\"),conditions{c});
    close all;
end
%%
allPSTHSCond = vertcat(allPSTHS{:}{:});
allTrialsCond = vertcat(allTrials{:}{:});
allTaskInds = vertcat(allSiteInds{:});
%allPSTHSCond = [abs(diff(allPSTHSCond,1,2)),zeros(size(allPSTHSCond,1),1)];

plotJointPSTHS(params.bins,{allPSTHSCond},{allTrialsCond},condInds,allTaskInds,[], alignLimits,[0 10],...
    cell2struct(num2cell(distinguishable_colors(length(conditions)),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath+"PSTHS\","All_PSTH",[]);

allRestTaskInds = vertcat(allRestInds{:});
restPSTHSCond = (allRestTaskInds./allRestTaskInds) .* allPSTHSCond(allTaskInds,:);
restTrialsCond = (allRestTaskInds./allRestTaskInds) .*allTrialsCond(allTaskInds,:);
plotJointPSTHS(params.bins,{restPSTHSCond},{restTrialsCond},condInds(allTaskInds),allRestTaskInds,...
    [],  alignLimits,[0 1],cell2struct(num2cell(distinguishable_colors(length(conditions)),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath+"PSTHS\","All_PSTH_EqRest",[]);

allRestTaskInds = ~vertcat(allRestInds{:});
restPSTHSCond = (allRestTaskInds./allRestTaskInds) .* allPSTHSCond(allTaskInds,:);
restTrialsCond = (allRestTaskInds./allRestTaskInds) .*allTrialsCond(allTaskInds,:);
plotJointPSTHS(params.bins,{restPSTHSCond},{restTrialsCond},condInds(allTaskInds),allRestTaskInds,...
    [], alignLimits,[0 1],cell2struct(num2cell(distinguishable_colors(length(conditions)),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath+"PSTHS\","All_PSTH_DiffRest",[]);