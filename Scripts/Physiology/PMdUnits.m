conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[0, 0]}},1,length(conditions));
alignLimits = {[-1.5, 2]};
pVal=0.05;
saveDir = "S:\Lab\ngc14\Working\PMd\Task_Units\";
params = PhysRecording(string(conditions),.01,.15,-6,15,containers.Map(conditions,...
    {"StartReach","StartReach","StartReach","GoSignal"}));
plotUnits = false;
MIN_BLOCKS_FOR_UNIT = 13;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","PMd");
clear rawSpikes;
%%
trialFR = cellfun(@(ct,cs,ta) cellfun(@(a,b) cell2mat(cellfun(@(m,tt) ...
    squeeze(mean(m(:,max(1,tt(1)):max(1,tt(end)),:),2,'omitnan').*(1*~all(isnan(tt)))),...
    squeeze(num2cell(a,[1,2])),cellfun(@(e) [cell2mat(e),repmat([NaN,NaN],isempty(cell2mat(e)),1)],num2cell(...
    [arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(1)),'UniformOutput',false),...
    arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(end)),'UniformOutput',false)],...
    2),'UniformOutput',false),'UniformOutput',false)'),ct,cs,'UniformOutput',false), siteTrialPSTHS, siteSegs,...
    cellfun(@(c,t) arrayfun(@(e) find(strcmp(c,e)),t{1}), params.condSegMap.values(taskAlign.keys),taskAlign.values,...
    'UniformOutput',false),'UniformOutput',false);
goodFR = cellfun(@(c) cellfun(@(s) s>2 & s<200,c,'UniformOutput',false),trialFR,'UniformOutput',false);
goodUnits = cellfun(@(tn) cell2mat(cellfun(@(s)sum(s,2), tn,'UniformOutput',false)),...
    num2cell(cat(2,goodFR{:}),2),'UniformOutput',false);
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,1),'omitnan')+2,params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t) mean(cell2mat(cellfun(@(a,n) max(1,mean(...
    a(:,n(~isnan(n)):n(~isnan(n))+(1/params.binSize),:),[2,3],'omitnan')),p,t,...
    'UniformOutput',false)),2,'omitnan'),num2cell([siteTrialPSTHS{:}],2),num2cell([allCondSegs{:}],2),"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) cellfun(@(p)p./repmat(nb,1,1,size(p,3)),cp,...
    'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),normBaseline,'Uniformoutput', false);
[tVals,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
taskUnits = cellfun(@(a,b) any(cell2mat(a),2) & sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2), ...
    num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false);
allMaps = cellfun(@(d) d{end}, chMaps, 'UniformOutput',false);
unitChannels = siteChannels;
unit2SiteMap=cell2mat(arrayfun(@(m,n) ones(1,size(m{1},2))*n,unitChannels,1:length(unitChannels),'UniformOutput',false));
maxCondsFR = cellfun(@(c) cellfun(@(d) max(mean(d,3,'omitnan'),[],2,'omitnan'),...
    c,'UniformOutput',false),normPSTH, 'UniformOutput',false);
% normPSTH = num2cell([siteTrialPSTHS{:}],2);
%
[~,restFR] = calculatePhases(params,containers.Map(["Rest"],{["GoSignal","StartReplaceHold"]}),taskWindow,siteSegs(end),siteTrialPSTHS(end),false,true);
[rVals,rUnit] = cellfun(@(tc) cellfun(@(r,c) ttestTrials(r,c,1,false,pVal),...
    restFR{1},tc,'UniformOutput',false),taskFR,'UniformOutput',false);
restUnits = cellfun(@(d) any(vertcat(d{:}),1)', num2cell([rUnit{:}],2), 'UniformOutput', false);
%%
maxClusters = 10;
for c =1:length(conditions)-1
    tUnits = cell2mat(taskUnits);
    taskSiteInds = find(cellfun(@any,arrayfun(@(a) tUnits(unit2SiteMap==a),...
        min(unit2SiteMap):max(unit2SiteMap),'UniformOutput', false))');
    tUnits = taskUnits(taskSiteInds);
    maxUnitFR = cell2mat(cellfun(@(m) cell2mat(m(taskSiteInds)), maxCondsFR, 'UniformOutput',false));
    unitPSTHS = cell2mat(cellfun(@(m,i) (i./i).*sqrt(mean(m,3,'omitnan')),vertcat(normPSTH{c}{taskSiteInds}),tUnits,'UniformOutput',false));
    for m = 1:maxClusters
        [solutions(:,m),~,sumd{m},~] = kmeans(unitPSTHS,m,'Distance','correlation','Replicates',5,'Options',statset('UseParallel',1));
    end
    clusterP = evalclusters(solutions,'kmeans','CalinskiHarabasz','KList',1:maxClusters);
    close all;
    plot(1:maxClusters,clusterP.CriterionValues);
    hold on;
    xlim([1,maxClusters]);
    scatter(clusterP.OptimalK,clusterP.CriterionValues(clusterP.OptimalK),'r','filled','o');
    saveFigures(gcf,saveDir+"Clustering\Correlation\"+params.condAbbrev(params.condNames(c))+"\","EvaluationValues",[]);
    wrapPlots = 3;
    for s = 1:length(clusterP.InspectedK)
        close all;
        figure('Units','normalized','Position',[0 0 1 1]);
        sCluster{s} = kmeans(unitPSTHS,s);
        clusterGroups = arrayfun(@(f) "Cluster_"+num2str(f),sCluster{s});
        siteUnitSegs = cellfun(@(si,tu) repmat(mean(si{1},1,'omitnan'),length(tu),1),...
            siteSegs{c}(taskSiteInds),tUnits,'UniformOutput',false);
        clusterColors = cell2struct(num2cell(distinguishable_colors(s),2),...
            unique(clusterGroups(~cellfun(@(cs) contains(cs,"NaN"),clusterGroups))));
        for i = 1:s
            subplot(round(s/wrapPlots)+1,min(s,wrapPlots),i);
        end
        plotJointPSTHS(params.bins,{unitPSTHS},{cell2mat(siteUnitSegs)},...
            clusterGroups,cell2mat(tUnits)',[],alignLimits,[0,1],clusterColors);
        hold on;
        clusterColors = cell2struct(repmat({[0.8 0.8 0.8]},length(fieldnames(clusterColors)),1),fieldnames(clusterColors));
        restPSTHS = cell2mat(cellfun(@(m,i) (i./i).*sqrt(mean(m,3,'omitnan')),...
            vertcat(normPSTH{end}{taskSiteInds}),tUnits,'UniformOutput',false));
        plotJointPSTHS(params.bins,{restPSTHS},{cell2mat(siteUnitSegs)},...
            clusterGroups,cell2mat(tUnits)',[],alignLimits,[0,1],clusterColors);
        if(clusterP.OptimalK==s)
            saveFigures(gcf,saveDir+"Clustering\Correlation\"+params.condAbbrev(params.condNames(c))+"\",num2str(s)+"_Optimal_PSTH",[]);
        else
            saveFigures(gcf,saveDir+"Clustering\Correlation\"+params.condAbbrev(params.condNames(c))+"\",num2str(s)+"_PSTH",[]);
        end
    end
end
%%
savePath = saveDir + "PSTHS\";
[allTrials,allPSTHS]= deal(cell(1,1));
[condInds,allSiteInds,allRestInds] = deal([]);
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
for c =1:length(conditions)
    close all;
    tUnits = cell2mat(taskUnits);
    [~,siteUnitMods] = unique(unit2SiteMap);
    taskSiteInds = find(cellfun(@any,arrayfun(@(a) tUnits(unit2SiteMap==a),...
        min(unit2SiteMap(unit2SiteMap~=0)):max(unit2SiteMap),'UniformOutput', false))');
    % taskSiteInds = intersect(taskSiteInds,find(cellfun(@(s) s(1), siteLocation)<515 & cellfun(@(s) s(1), siteLocation)<600));
    tUnits = taskUnits(taskSiteInds);
    rUnits = restUnits(taskSiteInds);
    siteDates = siteDateMap(taskSiteInds,:);
    siteUnitMods = siteUnitMods(taskSiteInds);
    unitChannelMaps = unitChannels(taskSiteInds);
    siteCondSegs = cellfun(@(a)  NaN(sum(strcmp(a(:,1),conditions(c)) & ...
        (~isnan(str2double(a(:,end-1))) | cellfun(@isempty,a(:,end-1)))),length(maxSegL)),siteTrialInfo(taskSiteInds)','UniformOutput',false);
    unitLocation = mapSites2Units(cellfun(@length,unitChannelMaps)',siteLocation(taskSiteInds));
    siteUnitNo = mapSites2Units(cellfun(@length,unitChannelMaps)',[siteDateMap{taskSiteInds,'Site'}]);
    siteUnitNames = "SiteNo"+arrayfun(@num2str,siteUnitNo,'UniformOutput',false);
    maxUnitFR = cell2mat(cellfun(@(m) cell2mat(m), maxCondsFR(taskSiteInds), 'UniformOutput',false));
    for j = 1:size(siteCondSegs,1)
        siteCondSegs{j}(:,condSegMappedInds{c}) = siteSegs{c}{taskSiteInds(j)}{1};
    end
    goodTrials=cellfun(@(s) sum(isnan(s),2)==mode(sum(isnan(s),2)) |...
        sum(isnan(s),2)==mode(sum(isnan(s),2))-1,siteCondSegs,'UniformOutput',false);
    siteUnitSegs = cellfun(@(si,tu,g) repmat(mean(si(g,:),1,'omitnan'),length(tu),1),siteCondSegs,tUnits,goodTrials,'UniformOutput',false);
    unitLocation = unitLocation-min(unitLocation(:,:)).*ImagingParameters.px2mm;

    unitPSTHS = cellfun(@(m,i,g) (i./i).*mean(m{c}(:,:,g),3,'omitnan'),num2cell(vertcat(normPSTH{taskSiteInds}),2),tUnits,goodTrials,'UniformOutput',false);
    sitePSTHS = cellfun(@(n,ti,g) cellfun(@(t) squeeze(t)',num2cell(n{c}(ti,:,g),[2,3]),'UniformOutput',false),...
        num2cell(vertcat(normPSTH{taskSiteInds}),2),tUnits,goodTrials,'Uniformoutput',false);
    if(plotUnits)
        plotJointPSTHS(params.bins,{vertcat(unitPSTHS{:})},{cell2mat(siteUnitSegs)},...
            siteUnitNames,cell2mat(tUnits), [],alignLimits,[0 5],cell2struct(num2cell(...
            distinguishable_colors(length(siteUnitSegs)),2),arrayfun(@(t) "SiteNo"+num2str(siteDateMap{t,'Site'}),taskSiteInds)));
        saveFigures(gcf,strcat(savePath,"Session_PSTHS\"),strcat("Task_Units_",params.condAbbrev(char(conditions(c))),"_"),[]);
        close all;
        for s = 1:length(sitePSTHS)
            ucs = unitChannelMaps{s}(tUnits{s});
            unitIndex = histcounts(ucs,[unique(ucs),max(ucs)+1]);
            utc = arrayfun(@(aa,bb) [string(aa),string(arrayfun(@(bi) string([num2str(aa),'_',num2str(bi)]),2:bb,'UniformOutput',false))],...
                unique(ucs),unitIndex,'UniformOutput',false);
            unitLabels = string(arrayfun(@(t) "SiteNo"+num2str(siteDates{s,'Site'})+"_"+t,cell2mat(utc),'UniformOutput',false))';
            plotColors =cell2struct(num2cell(distinguishable_colors(length(ucs)),2),unitLabels);
            unitTrialLabels = cell2mat(arrayfun(@(n) repmat(n,sum(goodTrials{s}),1),unitLabels,'UniformOutput',false));
            siteTrialSegs = repmat(siteCondSegs(s),length(unitLabels),1);
            plotJointPSTHS(params.bins,{cell2mat(sitePSTHS{s})},{cell2mat(siteTrialSegs)},...
                unitTrialLabels,true(size(unitTrialLabels)), [],alignLimits,[0 5],plotColors);
            saveFigures(gcf,savePath+"Units\"+string(datetime(siteDates.Date{s},'Format','MMMM_dd_yyyy'))+"\",...
                params.condAbbrev(params.condNames(c))+"_PSTH",[]);
            close all;
        end
    end
    allSiteInds = [allSiteInds;tUnits];
    allRestInds = [allRestInds;cellfun(@(u,d) d(u)==0, tUnits, rUnits,'UniformOutput',false)];
    condInds = [condInds,repmat(string(params.condAbbrev(params.condNames(c))),1,length(siteUnitNames))];
    allPSTHS = cellfun(@(a,b) vertcat(a,b), allPSTHS, {unitPSTHS}, 'UniformOutput', false);
    allTrials = cellfun(@(a,d) vertcat(a,cellfun(@(m) vertcat(m{:}),d,'UniformOutput',false)),...
        allTrials,{cellfun(@(su,tu) repmat({mean(su,1,'omitnan')},length(tu),1),siteCondSegs,tUnits,...
        'UniformOutput',false)},'UniformOutput',false);
    
    siteUnitInds = siteUnitNo;
    siteUnitInds(~ismember(1:length(siteUnitNo),siteUnitMods)) = NaN;
    siteUnitNames(~cell2mat(tUnits)) = "";
    allUnitSegs = vertcat(siteUnitSegs{:});

    % heatmap_distribution_plots(siteUnitNames,cell2mat(tUnits)',(siteUnitInds'),...
    % {allUnitSegs},params.bins,{cell2mat(unitPSTHS)},maxUnitFR(:,c),unitLocation,...
    % alignLimits,strcat(savePath,"Heatmaps\FR\"),conditions{c});
    close all;
end
%allPSTHSCond = [abs(diff(allPSTHSCond,1,2)),zeros(size(allPSTHSCond,1),1)];
allPSTHSCond = vertcat(allPSTHS{:}{:});
allTrialsCond = vertcat(allTrials{:}{:});
allTaskInds = vertcat(allSiteInds{:});
restUnitInds = strcmp(condInds,"R")';
restPSTHS = cell2mat(vertcat(repmat(allPSTHS{1}(end-(length(taskSiteInds)-1):end),length(conditions),1)));
%allPSTHSCond = allPSTHSCond-restPSTHS;

figure('Units','normalized','Position',[0 0 1 1]);
for a = 1:length(conditions)-1; subplot(1,length(conditions)-1,a); end
plotJointPSTHS(params.bins,{repmat(allPSTHSCond(restUnitInds,:),length(conditions)-1,1)},...
    {repmat(allTrialsCond(restUnitInds,:),length(conditions)-1,1)},cell2mat(arrayfun(@(s) ...
    repmat("R"+num2str(s),sum(restUnitInds),1),1:length(conditions)-1,'UniformOutput',false)'),...
    repmat(ones(sum(restUnitInds),1),length(conditions)-1,1),[],alignLimits,[-.5 5], ...
    cell2struct(repmat({[.6 .6 .6]},length(conditions)-1,1), ["R1","R2","R3"]));
plotJointPSTHS(params.bins,{allPSTHSCond(~restUnitInds,:)},{allTrialsCond(~restUnitInds,:)},condInds(~restUnitInds)',...
    allTaskInds(~restUnitInds),[], alignLimits,[-.5 5],cell2struct(num2cell(...
    distinguishable_colors(length(conditions),'r'),2),string(params.condAbbrev.values)));
saveFigures(gcf,savePath,"All_PSTH_RestSub",[]);

if(plotUnits)
    figure('Units','normalized','Position',[0 0 1 1]);
    for a = 1:length(conditions)-1; subplot(1,length(conditions)-1,a); end
    restUnitInds = restUnitInds(allTaskInds);
    condInds = condInds(allTaskInds);
    allRestTaskInds = vertcat(allRestInds{:});
    restPSTHSCond = (allRestTaskInds./allRestTaskInds) .* allPSTHSCond(allTaskInds,:);
    restTrialsCond = (allRestTaskInds./allRestTaskInds) .*allTrialsCond(allTaskInds,:);
    plotJointPSTHS(params.bins,{repmat(restPSTHSCond(restUnitInds,:),length(conditions)-1,1)},...
        {repmat(restTrialsCond(restUnitInds,:),length(conditions)-1,1)},cell2mat(arrayfun(@(s) ...
        repmat("R"+num2str(s),sum(restUnitInds),1),1:length(conditions)-1,'UniformOutput',false)'),...
        repmat(ones(sum(restUnitInds),1),length(conditions)-1,1),[],alignLimits,[0 5], ...
        cell2struct(repmat({[.4 .4 .4]},length(conditions)-1,1), ["R1","R2","R3"]));
    plotJointPSTHS(params.bins,{restPSTHSCond(~restUnitInds,:)},{restTrialsCond(~restUnitInds,:)},...
        condInds(~restUnitInds)',allRestTaskInds(~restUnitInds),[],alignLimits,[0 5],...
        cell2struct(num2cell(distinguishable_colors(length(conditions),'r'),2),string(params.condAbbrev.values)));
    saveFigures(gcf,savePath,"All_PSTH_EqRest",[]);

    figure('Units','normalized','Position',[0 0 1 1]);
    for a = 1:length(conditions)-1; subplot(length(conditions)-1,1,a); end
    allRestTaskInds = ~vertcat(allRestInds{:});
    restPSTHSCond = (allRestTaskInds./allRestTaskInds) .* allPSTHSCond(allTaskInds,:);
    restTrialsCond = (allRestTaskInds./allRestTaskInds) .*allTrialsCond(allTaskInds,:);
    plotJointPSTHS(params.bins,{repmat(restPSTHSCond(restUnitInds,:),length(conditions)-1,1)},...
        {repmat(restTrialsCond(restUnitInds,:),length(conditions)-1,1)},cell2mat(arrayfun(@(s) ...
        repmat("R"+num2str(s),sum(restUnitInds),1),1:length(conditions)-1,'UniformOutput',false)'),...
        repmat(ones(sum(restUnitInds),1),length(conditions)-1,1),[],alignLimits,[0 5], ...
        cell2struct(repmat({[.4 .4 .4]},length(conditions)-1,1), ["R1","R2","R3"]));
    plotJointPSTHS(params.bins,{restPSTHSCond(~restUnitInds,:)},{restTrialsCond(~restUnitInds,:)},...
        condInds(~restUnitInds)',allRestTaskInds(~restUnitInds),[],alignLimits,[0 5],...
        cell2struct(num2cell(distinguishable_colors(length(conditions),'r'),2),string(params.condAbbrev.values)));
    saveFigures(gcf,savePath,"All_PSTH_DiffRest",[]);
end