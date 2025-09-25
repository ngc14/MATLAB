conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[-0.3, 0]}},1,length(conditions));
pVal=0.05;
alignLimits = [-3, 10.5];
params = PhysRecording(string(conditions),.01,.15,-5,15,containers.Map(conditions,...
    {"GoSignal","GoSignal","GoSignal","GoSignal"}));
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\Errors\";
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","PMd");
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
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
%%
trialInfo = cellfun(@(c) cellfun(@(t) t(strcmp(t(:,1),c),:),siteTrialInfo','UniformOutput',false),conditions,'UniformOutput',false);
siteTrialSegs = cellfun(@(c) cellfun(@(n) NaN(size(n,1),length(maxSegL)), c, 'UniformOutput',false), trialInfo,'UniformOutput',false);
trialInfo = cellfun(@(c) vertcat(c{:}),num2cell(horzcat(trialInfo{:}),2),'UniformOutput',false);
for j = 1:size(siteTrialSegs,2)
    for i = 1:length(siteTrialInfo)
        siteTrialSegs{j}{i}(:,condSegMappedInds{j}) = siteSegs{j}{i}{1};
    end
end
siteTrialSegs=cellfun(@(v) vertcat(v{:}),num2cell(cat(2,siteTrialSegs{:}),2),'UniformOutput',false);
failedTrials = cell(height(siteDateMap),1);
for s = 1:length(trialInfo)
    numNans = cumsum(isnan(siteTrialSegs{s}),2);
    cInds = contains(trialInfo{s}(:,1),'Sphere');
    tInds = isnan(str2double(trialInfo{s}(:,end-1))) & ~isempty(trialInfo{s}(:,end-1));
    failedTrials{s} = strings(length(tInds),1);
    failedTrials{s}(cInds & tInds) = replace(lower(string(trialInfo{s}(cInds & tInds,end-1))),{' ', '-'},'_');
    failedTrials{s}(cInds & ~tInds & (numNans(:,find(contains(maxSegL,'Replace'),1))<=1 & ...
        numNans(:,end)>1)) = 'failed_to_replace';
    failedTrials{s}(cInds & failedTrials{s}=='' & numNans(:,end)<=1) = "success";
end
failTypes = unique(cell2mat(failedTrials));
failTypes = failTypes(arrayfun(@(a) find(contains(failTypes,a)),["start","reach","contact","lift","hold","replace","success"]));
infoTable = table();
infoTable.SiteNo = "SiteNo" + table2array(siteDateMap(:,'Site'));
infoTable.TaskModulated = cellfun(@(c) any(horzcat(c{:}),2), num2cell(cat(2,tUnit{:}),2), 'UniformOutput',false);
infoTable.Condition = cellfun(@(t) string(t(:,1)), trialInfo,'UniformOutput',false);
infoTable.SegTimes = siteTrialSegs;
infoTable.Outcomes = failedTrials;
infoTable.PSTHS = cellfun(@(r,i)(i./i).*cell2mat(reshape(r,1,1,[])),num2cell(...
    cellfun(@cell2mat,normPSTH,'UniformOutput',false),2),infoTable.TaskModulated,'UniformOutput',false);
%%
for f = 1:length(failTypes)
    failTrials = cellfun(@(s) strcmp(s,failTypes(f)), infoTable.Outcomes,'UniformOutput',false);
    if(f==length(failTrials))
        otherTrials = cellfun(@(s) sum(contains(s,failTypes(1:end-1)))/length(failTypes)-1, infoTable.Outcomes,'UniformOutput',false);
        failTrials = cellfun(@(f) f.*ismember(find(f),sample(find(f),o)), failTrials,otherTrials,'UniformOutput',false);
    end
    sessionPSTHS = cellfun(@(p,f) p(:,:,f), infoTable.PSTHS,failTrials,'UniformOutput',false);
    sessionSegs = cellfun(@(s,f) s(f,:), infoTable.SegTimes,failTrials, 'UniformOutput',false);
    sessionNo = cellfun(@(s,f) string(repmat(s,sum(f),1)), infoTable.SiteNo, failTrials, 'UniformOutput',false);
    close all;
    plotPSTHS(params.bins,sessionPSTHS,sessionSegs,sessionNo,[-3,1.5],[0 4],...
        cell2struct(num2cell(distinguishable_colors(height(infoTable)),2),infoTable.SiteNo));
    saveFigures(gcf,savePath+"Session_PSTHS\","Traces_"+failTypes(f),[]);
end
%%
close all;
failColors = cell2struct(num2cell(distinguishable_colors(length(failTypes)),2),string(failTypes));
totals = cellfun(@(h) histcounts(categorical(h)), infoTable.Outcomes, 'UniformOutput', false);
%sampleNum = cellfun(@(m) ceil(mean(m(1:end-1))), totals);
sampleNum(isnan(sampleNum)) = 0;
plotOutcomes = infoTable.Outcomes;
for p = 1:length(plotOutcomes)
     f = contains(infoTable.Condition{p},"Sphere") & strcmp(plotOutcomes{p},"success") & ~all(isnan(infoTable.SegTimes{p}),2);
     plotOutcomes{p}(strcmp(plotOutcomes{p},"success")) = "";
     plotOutcomes{p}(randsample(find(f),sampleNum(p))) = "success";
end
plotPSTHS(params.bins,infoTable.PSTHS,infoTable.SegTimes,plotOutcomes,alignLimits,[0 10],failColors);
lastAx = gca;
lastAx = gcf;
lastAx = lastAx.Children;
lastAx = lastAx(end);
xPos = get(lastAx.Children,'XData');
xPos = unique(cell2mat(xPos(cellfun(@length,xPos)==2)));
newAx = axes('Position',get(lastAx,'Position'),'XAxisLocation','top','Color','none','FontWeight','bold');
newAx.XTick = xPos./(max(lastAx.XLim)-min(lastAx.XLim));
newAx.XTickLabel = arrayfun(@(s) string(erase(s,'Start')), maxSegL);
newAx.XTickLabelRotation = 25;
newAx.YTick = [];
saveFigures(gcf,savePath+"Errors\","Averages",[]);


function figHandle = plotPSTHS(bins,PSTHIn,allSegsIn,groupInds,...
    PSTHDisplayLimits,FRLimIn,plotColors)
if(~exist('FRLimIn', 'var'))
    FRLim = [5 30];
else
    FRLim = FRLimIn;
end
segColors = ['k','r'];
PSTH = PSTHIn;
allSegs = allSegsIn;
allGroups = groupInds;
plotNames = fieldnames(plotColors);
groupName = plotNames(cellfun(@(c) contains(string(c),cell2mat(allGroups)),plotNames));
zeroBinInd = find(bins==0);
binSize = mode(diff(bins));
plotBins = (fix(PSTHDisplayLimits(1)/binSize)+zeroBinInd):...
    (fix(PSTHDisplayLimits(end)/binSize)+zeroBinInd);
PSTH =  cellfun(@(t) t(:,plotBins,:),PSTH,'UniformOutput',false);
g = groot;
figHandle = g.CurrentFigure;
if(~isempty(get(figHandle,'Children')))
    pos = cell2mat(arrayfun(@(n) get(n,'Position'), get(figHandle,'Children'), 'UniformOutput', false));
    wrapPlots = numel(unique(pos(:,1)));
    nrows = numel(unique(pos(:,2)));
else
    figure();
    figHandle = gcf;
    wrapPlots = min(length(groupName),3);
    nrows = ceil(length(groupName)/wrapPlots);
    wrapPlots = nrows; nrows = 3;
   %wrapPlots=1;nrows=length(groupName);
end
maxPlot = 0;
for j = 1:length(groupName)
    condInds = cellfun(@(a) strcmp(a,groupName{j}), allGroups,'UniformOutput',false);
    groupPSTH = cellfun(@(t,g) mean(t(:,:,g),3,'omitnan'), PSTH,condInds, 'UniformOutput', false);
    tUnits = cellfun(@(g) ~all(isnan(g),[2 3]), groupPSTH, 'UniformOutput',false);
    subplot(nrows,wrapPlots,j);hold on;
    title(replace(strcat(groupName{j},": ", num2str(sum(cell2mat(tUnits))), " units, ",... %num2str(length(allGroups))," sites, "
        num2str(sum(cell2mat(condInds))), " trials, ",num2str(sum(cellfun(@sum,condInds).*cellfun(@sum,tUnits))),...
        " instances"),"_", " "));
    tUnits = cell2mat(tUnits);
    if(any(cellfun(@any,condInds)))
        xAlignTicks = 0+(1:length(plotBins));
        plot(xAlignTicks,nanmean(cell2mat(groupPSTH),1), 'LineWidth',2,'Color', plotColors.(groupName{j}));
        uE=nanmean(cell2mat(groupPSTH),1)+sum(cell2mat(cellfun(@(g) ...
            nanstd(g,0,1),groupPSTH,'UniformOutput',false))./sum(tUnits),1,'omitnan');
        lE=nanmean(cell2mat(groupPSTH),1)-sum(cell2mat(cellfun(@(g) ...
            nanstd(g,0,1),groupPSTH,'UniformOutput',false))./sum(tUnits),1,'omitnan');
        yP=[lE,fliplr(uE)];
        xP=[xAlignTicks,fliplr(xAlignTicks)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        %d = patch(xP,yP,1);
        %set(d,'edgecolor','none','facealpha',.5,'facecolor',plotColors.(groupName{j}));
        plot(xAlignTicks,cell2mat(groupPSTH),'LineWidth',1)%,'Color', plotColors.(groupName{j}));
        if(~isempty(yP))
            maxPlot = max(maxPlot,FRLim(end)*max(1,ceil(mean(prctile(cell2mat(groupPSTH),95,2),'omitnan')/FRLim(end))));
        end
    end
end
for j = 1:length(groupName)
    groupSegs = cell2mat(cellfun(@(t,g) t(g,:), allSegs,cellfun(@(a) ...
        strcmp(a,groupName{j}), allGroups,'UniformOutput',false),'UniformOutput',false));
    subplot(nrows,wrapPlots,j);hold on;
    c = sum(~isnan(groupSegs),2);
    groupSegs(c~=mode(c),mode(c)) = groupSegs(sub2ind(size(groupSegs),find(c~=mode(c)),c(c~=mode(c))));
    avgSegs = nanmean(groupSegs,1);
    lastInd = sub2ind(size(groupSegs),1:size(groupSegs,1),c');
    for s = 1:mode(c)
        if(avgSegs(s)>=PSTHDisplayLimits(1) && ...
                avgSegs(s)<=PSTHDisplayLimits(end))
            if(avgSegs(s)==0 || s>=mode(c))
                plotColor = segColors(1);
            else
                plotColor = segColors(end);
            end
            pSeg = find(isalmost(PSTHDisplayLimits(1):binSize:...
                PSTHDisplayLimits(end),avgSegs(s),binSize/1.99),1);
            plot([xAlignTicks(pSeg) xAlignTicks(pSeg)],[0 maxPlot],...
                'Color',plotColor,'LineStyle','--');
        end
    end
    lowCI= find(isalmost(PSTHDisplayLimits(1):binSize:PSTHDisplayLimits(end),...
        max(PSTHDisplayLimits(1),prctile(groupSegs(lastInd),5)),binSize/1.99),1);
    highCI= find(isalmost(PSTHDisplayLimits(1):binSize:PSTHDisplayLimits(end),...
        min(PSTHDisplayLimits(end),prctile(groupSegs(lastInd),95)),binSize/1.99),1);
    if(~isempty(lowCI) && ~isempty(highCI))
    patch([lowCI, highCI, highCI, lowCI],[0 0 maxPlot maxPlot],plotColors.(groupName{j}),'FaceAlpha',0.25,'LineStyle','none');
    end
    allXTicks = PSTHDisplayLimits(1):binSize:PSTHDisplayLimits(end);
    xticks([0,find(isalmost(mod(allXTicks,1),1,binSize)),xAlignTicks(end)]);
    allLabels = [num2str(PSTHDisplayLimits(1),'%.2f'), arrayfun(@num2str,...
        round(allXTicks(isalmost(mod(allXTicks,1),1,binSize))),'UniformOutput',false),...
        num2str(PSTHDisplayLimits(end),'%.2f')];
    xticklabels(allLabels);
    set(gca,'XMinorTick','on');
    set(get(gca,'XAxis'),'MinorTickValues',find(isalmost(mod(allXTicks,1),.5,binSize/1.99)));
end
figHandle = figHandle.Children;
set(figHandle(strcmp(get(figHandle,'Type'),'axes')),'YLim',[FRLim(1),maxPlot]);
end