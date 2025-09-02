conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]}, {["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[-0.3, 0]}},1,length(conditions));
pVal=0.05;
monkey = "Gilligan";
params = PhysRecording(string(conditions),.01,.15,-5,11,containers.Map(conditions,...
    {["GoSignal"],["GoSignal"],["GoSignal"],"GoSignal"}));
alignLimits = {[-.75, 10.5]};
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
savePath = "S:\Lab\ngc14\Working\PMd\Errors\";
plotUnits = false;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,trialInfo] = getAllSessions(params,"Single","PMd");
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
siteTrialSegs = cellfun(@cell2mat,horzcat(siteSegs{:}),'UniformOutput',false);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,1),'omitnan')-3,params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) ...
    permute(mean(c(:,s:s+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),a(~isnan(n)),...
    num2cell(n(~isnan(n))),'UniformOutput',false),[1,1,sum(~isnan(n))])),3,'omitnan'));NaN(all(isnan(n)).*size(a{1},1),1)],p,t,...
    'UniformOutput',false),siteTrialPSTHS,allCondSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp{:}),repmat(nb,1,size(vertcat(cp{:}),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
normPSTH = horzcat(normPSTH{:});
tUnit = cat(2,tUnit{:});
%%
unit2SiteMap=cell2mat(cellfun(@(m,n) ones(1,size(m,2))*n,siteChannels,num2cell(1:length(siteChannels)),'UniformOutput',false));
[~,siteUnitMods] = unique(unit2SiteMap);
taskUnits = cell2mat(cellfun(@(a) any(cell2mat(a),2), num2cell(tUnit,2),'Uniformoutput',false));
taskSiteInds = find(cellfun(@any,arrayfun(@(a) taskUnits(unit2SiteMap==a),...
    min(unit2SiteMap(unit2SiteMap~=0)):max(unit2SiteMap),'UniformOutput', false))');
siteTrialSegs = cellfun(@(n) NaN(size(n,1),length(maxSegL)), siteTrialSegs(taskSiteInds,:),'UniformOutput',false);
for i = 1:size(siteTrialSegs,1)
    for j = 1:size(siteTrialSegs,2)
        siteTrialSegs{i,j}(:,condSegMappedInds{j}) = siteSegs{j}{taskSiteInds(i)}{1}; 
    end
end
siteDates = siteDateMap(taskSiteInds,:);
siteUnitMods = siteUnitMods(taskSiteInds);
unitChannelMaps = siteChannels(taskSiteInds)';
unitLocation = mapSites2Units(cellfun(@length,unitChannelMaps),siteLocation(taskSiteInds));
unitLocation = unitLocation-min(unitLocation(:,:)).*ImagingParameters.px2mm;
siteUnitNo = mapSites2Units(cellfun(@length,unitChannelMaps),[siteDateMap{taskSiteInds,'Site'}]);
siteUnitNo = "SiteNo"+arrayfun(@num2str,siteUnitNo,'UniformOutput',false);
taskUnitsCond = cellfun(@(a) logical(cell2mat(a)),num2cell(tUnit(taskSiteInds,:),2), 'UniformOutput',false);
trialInfo = cellfun(@(c) cellfun(@(t) t(strcmp(t(:,1),c),:), trialInfo(taskSiteInds), 'UniformOutput', false),conditions,'UniformOutput',false);
normPSTH  = normPSTH(taskSiteInds,:);
%%
taskUnitsAll = cellfun(@(a) any(a,2), taskUnitsCond,'UniformOutput',false);
allTrialInfo = cellfun(@(r) vertcat(r{:}),num2cell(cat(1,trialInfo{:})',2),'UniformOutput',false);
allTrials = cellfun(@(r) cell2mat(reshape(r,[],1)),num2cell(siteTrialSegs,2),'UniformOutput',false);
allPSTHS =cellfun(@(r) cell2mat(reshape(r,1,1,[])),num2cell(...
    cellfun(@cell2mat,normPSTH,'UniformOutput',false),2),'UniformOutput',false);
failTypes = cellfun(@(t) string(unique(t(cellfun(@(n) isnan(str2double(n)) & ~isempty(n),...
    t(:,7)),7))), allTrialInfo,'UniformOutput',false);
failTypes = unique(lower(vertcat(failTypes{:})));
trialInds = cellfun(@(a) ismember(lower(a(:,end-1)),failTypes),allTrialInfo, 'UniformOutput',false);
allPSTHS = cellfun(@(m,i,a) (i./i).*m(:,:,a),allPSTHS,taskUnitsAll,trialInds,'UniformOutput',false);
allTrialInfo = cellfun(@(t,a) t(a,:),allTrialInfo,trialInds,'UniformOutput',false);
allTrials = cellfun(@(t,a) t(a,:),allTrials,trialInds,'UniformOutput',false);
siteLabels = siteUnitNo(ismember(cellfun(@(s)string(s(7:end)),...
    cellstr(siteUnitNo)),string(table2array(siteDateMap(:,'Site')))));
for c =1:length(failTypes)
    condInds = reshape(cellfun(@(t) strcmpi(t(:,end-1),failTypes(c)),allTrialInfo,'UniformOutput',false),[],1);
    unitPSTHS = cellfun(@(m,i,t) (i./i).*mean(m(:,:,t),3,'omitnan'),allPSTHS,...
        taskUnitsAll,condInds,'UniformOutput',false);
    trialPSTHS = cellfun(@(n,ti,r) squeeze(num2cell((r./r).*permute(n(ti,:,:),[3 2 1]),[1,2])),...
       allPSTHS,taskUnitsAll,condInds,'Uniformoutput',false);
    if(plotUnits)
        plotPSTHS(params.bins,{vertcat(unitPSTHS{:})},{cell2mat(cellfun(@(m,r) repmat(mean(m,1,'omitnan'),...
            size(r,1),1),allTrials,taskUnitsCond,'UniformOutput',false))},...
            siteLabels,any(vertcat(taskUnitsCond{:}),2)',alignLimits,[0 10],...
            cell2struct(num2cell(distinguishable_colors(height(siteDates)),2),...
            arrayfun(@(t)"SiteNo"+num2str(t),[siteDates{:,'Site'}])));
        saveFigures(gcf,strcat(savePath,"PSTHS\Session_PSTHS\"),strcat("Task_Units_",failTypes(c)),[]);
        close all;
        for s = 1:size(trialPSTHS,1)
            sPSTHS = condInds{s};
            if(any(sPSTHS))
            ucs = siteChannels{s}(taskUnitsAll{s});
            unitIndex = histcounts(ucs,[unique(ucs),max(ucs)+1]);
            utc = arrayfun(@(aa,bb) [string(aa),string(arrayfun(@(bi) string([num2str(aa),'_',num2str(bi)]),...
                2:bb,'UniformOutput',false))],unique(ucs),unitIndex,'UniformOutput',false);
            unitLabels = string(arrayfun(@(t) "SiteNo"+num2str(siteDates{s,'Site'})+"_"+t,cell2mat(utc),'UniformOutput',false))';
            plotColors =cell2struct(num2cell(distinguishable_colors(length(ucs)),2),unitLabels);
            unitTrialLabels = cell2mat(arrayfun(@(n) repmat(n,size(allTrials{s},1),1),unitLabels,'UniformOutput',false));
            plotPSTHS(params.bins,{cell2mat(trialPSTHS{s})},{repmat(allTrials{s},length(ucs),1)},...
                unitTrialLabels,true(fliplr(size(unitTrialLabels))),alignLimits,[0 10],plotColors);
            saveFigures(gcf,savePath+"PSTHS\Units\"+string(datetime(siteDates.Date{s},'Format','MMMM_dd_yyyy'))+"\",...
                strcat("Task_Units_",failTypes(c))+"_PSTH",[]);
            end
            close all;
        end
    end
    close all;
end

function figHandle = plotPSTHS(bins,PSTHIn,allSegsIn,allRepsIn,siteInds,...
    PSTHDisplayLimits,FRLimIn,plotColors)
if(~exist('FRLimIn', 'var'))
    FRLim = [5 30];
else
    FRLim = FRLimIn;
end
alignmentGap = .1;
segColors = ['k','r'];
PSTH = PSTHIn;
allSegs = allSegsIn;
allGroups = allRepsIn;
plotNames = fieldnames(plotColors);
groupName = natsort(plotNames(matches(unique(allGroups(allGroups~="")),plotNames)));

zeroBinInd = find(bins==0);
binSize = mode(diff(bins));
alignmentGap = alignmentGap/binSize;

PSTH =  cellfun(@(t,w) t(:,(fix(w(1)/binSize)+zeroBinInd):...
    fix((w(end)/binSize)+zeroBinInd)),PSTH,PSTHDisplayLimits,'UniformOutput',false);
g = groot;
figHandle = g.CurrentFigure;
if(~isempty(figHandle))
    pos = cell2mat(arrayfun(@(n) get(n,'Position'), get(figHandle,'Children'), 'UniformOutput', false));
    wrapPlots = numel(unique(pos(:,1)));
    nrows = numel(unique(pos(:,2)));
else
    figure();
    figHandle = gcf;
    wrapPlots = min(length(groupName),3);
    nrows = ceil(length(groupName)/wrapPlots);
end
xAlignTicks = {};
maxPlot = 0;
for j = 1:length(groupName)
    groupInds = arrayfun(@(s) cellfun(@(a) strcmp(a,groupName{j}),s), allGroups,'UniformOutput',true);
    groupPSTH = cellfun(@(t) t(groupInds,:), PSTH, 'UniformOutput', false);
    groupSegs = cellfun(@(t) t(strcmp(allGroups,groupName{j}) & siteInds',:), allSegs,'UniformOutput',false);
    subplot(nrows,wrapPlots,j);hold on;
    titleName = replace(groupName{j},"_", " ");
    titleName = strcat(titleName, " (n= ", num2str(size(groupSegs{1},1)),")");
    if(sum(siteInds)~=length(allGroups))
        titleName = strcat(titleName{:}(1:end-1), " of ", num2str(sum(groupInds)),")");
    end
    title(titleName);
    if(sum(groupInds)>0)
        plotStart = 0;
        for a = 1:size(groupPSTH,2)
            currSegs = groupSegs{a};
            currJointAlign = groupPSTH{a};
            xAlignTicks{a} = plotStart+(1:size(currJointAlign,2));
            plot(xAlignTicks{a},nanmean(currJointAlign,1), 'LineWidth',2,'Color', plotColors.(groupName{j}));
            uE=nanmean(currJointAlign,1)+(nanstd(currJointAlign,0,1)/sqrt(sum(groupInds)));
            lE=nanmean(currJointAlign,1)-(nanstd(currJointAlign,0,1)/sqrt(sum(groupInds)));
            yP=[lE,fliplr(uE)];
            xP=[xAlignTicks{a},fliplr(xAlignTicks{a})];
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            d = patch(xP,yP,1);
            set(d,'edgecolor','none','facealpha',.5,'facecolor',plotColors.(groupName{j}));
            if(~isempty(yP))
                maxPlot = max(maxPlot,FRLim(end)*max(1,ceil(max(yP)/FRLim(end))));
            end
            avgSegs = nanmean(currSegs,1);
            if(a==1)
                plotted = false(1,size(currSegs,2));
            end
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=PSTHDisplayLimits{a}(1) && ...
                        avgSegs(s)<=PSTHDisplayLimits{a}(end) && ...
                        (~plotted(s) ||avgSegs(s)==plotStart))
                    if(avgSegs(s)==0)
                        plotColor = segColors(1);
                    else
                        plotColor = segColors(end);
                    end
                    %plotted(s) = true;
                    pSeg = find(isalmost(PSTHDisplayLimits{a}(1):binSize:...
                        PSTHDisplayLimits{a}(end),avgSegs(s),binSize/1.99),1);
                    plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 maxPlot],...
                        'Color',plotColor,'LineStyle','--');
                end
            end
            if(a==size(groupPSTH,2))
                allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs(...
                    pd(1):binSize:pd(end))==min(abs(pd(1):binSize:pd(end))))-1,...
                    ta(end)], xAlignTicks,PSTHDisplayLimits, 'UniformOutput', false);
                xticks([allXTicks{:}]);
                allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                    "0"; num2str(pd(end),'%.2f')],PSTHDisplayLimits, 'UniformOutput', false);
                xticklabels([allLabels{:}]);
            end
            plotStart = plotStart + size(currJointAlign,2) + alignmentGap;
        end
    end
end
set(figHandle.Children,'YLim',[FRLim(1),maxPlot]);
end