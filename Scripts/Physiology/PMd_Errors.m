conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[-0.3, 0]}},1,length(conditions));
alignLimits = [-5.5, 14];
pVal=0.05;
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\Errors\";
params = PhysRecording(string(conditions),.01,.15,-6,15,containers.Map(conditions,...
    {"GoSignal","GoSignal","GoSignal","GoSignal"}));
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","PMd");
clear rawSpikes
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,1)+2,'omitnan'),params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t) mean(cell2mat(cellfun(@(a,n) max(1,median(...
    a(:,n(~isnan(n)):n(~isnan(n))+(1/params.binSize),:),[2,3],'omitnan')),...
    p,t,'UniformOutput',false)),2,'omitnan'),num2cell([siteTrialPSTHS{:}],2),num2cell([allCondSegs{:}],2),"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) cellfun(@(p)p./repmat(nb,1,1,size(p,3)),...
    cp,'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),normBaseline,'Uniformoutput', false);
normPSTH = vertcat(normPSTH{:});
normPSTH = horzcat(siteTrialPSTHS{:});
%%
trialInfo = cellfun(@(c) cellfun(@(t) t(strcmp(t(:,1),c),:),siteTrialInfo,'UniformOutput',false)',conditions,'UniformOutput',false);
siteTrialSegs = cellfun(@(c) cellfun(@(n) NaN(size(n,1),length(maxSegL)), c, 'UniformOutput',false), trialInfo,'UniformOutput',false);
trialInfo = cellfun(@(c) vertcat(c{:}),num2cell(horzcat(trialInfo{:}),2),'UniformOutput',false);
for j = 1:size(siteTrialSegs,2)
    for i = 1:length(siteTrialInfo)
        siteTrialSegs{j}{i}(:,condSegMappedInds{j}) = siteSegs{j}{i}{1};
    end
end
siteTrialSegs=cellfun(@(v) vertcat(v{:}),num2cell(cat(2,siteTrialSegs{:}),2),'UniformOutput',false);
failedTrials = cell(height(siteDateMap),1);
allPSTHS = cellfun(@(r) cell2mat(reshape(r,1,1,[])),num2cell(normPSTH,2),'UniformOutput',false);
for s = 1:length(trialInfo)
    badTrials = all(isnan(siteTrialSegs{s}),2);
    siteTrialSegs{s} = siteTrialSegs{s}(~badTrials,:);
    trialInfo{s} = trialInfo{s}(~badTrials,:);
    allPSTHS{s} = allPSTHS{s}(:,:,~badTrials);
    numNans = cumsum(isnan(siteTrialSegs{s}),2);
    tInds = isnan(str2double(trialInfo{s}(:,end-1))) & ~cellfun(@isempty,trialInfo{s}(:,end-1));
    failedTrials{s} = strings(length(tInds),1);
    failedTrials{s}(tInds) = replace(lower(string(trialInfo{s}(tInds,end-1))),{' ', '-'},'_');
    failedTrials{s}(~tInds & (numNans(:,find(contains(maxSegL,'Replace'),1))<=1 & ...
        numNans(:,end)>1)) = 'failed_to_replace';
    failedTrials{s}(failedTrials{s}=='' & (numNans(:,end)<=1 |...
        cellfun(@(t) strcmp(t,"Rest"),trialInfo{s}(:,1)) & numNans(:,end)==5)) = "success";
end
failTypes = unique(cell2mat(failedTrials));
failTypes = failTypes(arrayfun(@(a) find(endsWith(failTypes,a),1),["success",...
    "start","reach","contact","lift","hold","replace"]));
infoTable = table();
infoTable.SiteNo = "SiteNo" + table2array(siteDateMap(:,'Site'));
infoTable.TaskModulated = cellfun(@(c) any(horzcat(c{:}),2), num2cell(cat(2,tUnit{:}),2), 'UniformOutput',false);
infoTable.Condition = cellfun(@(t) string(t(:,1)), trialInfo,'UniformOutput',false);
infoTable.SegTimes = siteTrialSegs;
infoTable.Outcomes = failedTrials;
infoTable.PSTHS = cellfun(@(r,i)(i./i).*r,allPSTHS,infoTable.TaskModulated,'UniformOutput',false);
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
    plotPSTHS(params.bins,sessionPSTHS,sessionSegs,sessionNo,alignLimits,[0 4],...
        cell2struct(num2cell(distinguishable_colors(height(infoTable)),2),infoTable.SiteNo));
    saveFigures(gcf,savePath+"Session_PSTHS\","Averages_"+failTypes(f),[]);
end
%%
close all;
failColors = cell2struct(num2cell(distinguishable_colors(length(failTypes),[0 0 0; 1 1 1]),2),string(failTypes)); 
% cell2struct(num2cell(distinguishable_colors(length(conditions),[1 0 0; 0 0 0]),2),replace(conditions," ","_"))
totals = cellfun(@(h) histcounts(categorical(h)), infoTable.Outcomes, 'UniformOutput', false);
sampleNum = cellfun(@(m) ceil(mean(m(1:end-1))), totals);
sampleNum(isnan(sampleNum)) = 0;
plotOutcomes = infoTable.Outcomes;
for p = 1:length(plotOutcomes)
     f = contains(infoTable.Condition{p},"Large Sphere") & strcmp(plotOutcomes{p},"success") & ~all(isnan(infoTable.SegTimes{p}),2);
     plotOutcomes{p}(strcmp(plotOutcomes{p},"success")) = "";%replace(infoTable{p,'Condition'}{1}(strcmp(plotOutcomes{p},'success'))," ","_");
     plotOutcomes{p}(randsample(find(f),sampleNum(p))) = "success";
end
plotPSTHS(params.bins,infoTable.PSTHS,infoTable.SegTimes,plotOutcomes,alignLimits,[4 12],failColors);
lastAx = gca;
lastAx = gcf;
lastAx = lastAx.Children;
lastAx = lastAx(end);
xPos = get(lastAx.Children,'XData');
xPos = unique(cell2mat(xPos(cellfun(@length,xPos)==2)));
newAx = axes('Position',get(lastAx,'Position'),'XAxisLocation','top','Color','none','FontWeight','bold');
newAx.Position(2) = newAx.Position(2)-.002;
newAx.XTick = (xPos-min(lastAx.XLim))./(max(lastAx.XLim)-min(lastAx.XLim));
newAx.XTickLabel = arrayfun(@(s) string(erase(s,'Start')), maxSegL);
newAx.XTickLabelRotation = 45;
newAx.YTick = [];
saveFigures(gcf,savePath+"Errors\","Averages",[]);


function figHandle = plotPSTHS(bins,PSTHIn,allSegsIn,groupInds,...
    PSTHDisplayLimits,FRLimIn,plotColors)
if(~exist('FRLimIn', 'var'))
    FRLim = [5 30];
else
    FRLim = FRLimIn;
end
PSTH = PSTHIn;
allSegs = allSegsIn;
allGroups = groupInds;
plotNames = fieldnames(plotColors);
groupName = plotNames(cellfun(@(c) contains(string(c),cell2mat(allGroups)),plotNames));
binSize = mode(diff(bins));
plotBins = findBins(PSTHDisplayLimits(1),bins):1:findBins(PSTHDisplayLimits(end),bins);
bins=bins(plotBins);
PSTH =  cellfun(@(t) t(:,plotBins,:),PSTH,'UniformOutput',false);
gr = groot;
figHandle = gr.CurrentFigure;
if(~isempty(get(figHandle,'Children')))
    pos = cell2mat(arrayfun(@(n) get(n,'Position'), get(figHandle,'Children'), 'UniformOutput', false));
    wrapPlots = numel(unique(pos(:,1)));
    nrows = numel(unique(pos(:,2)));
else
    figure();
    figHandle = gcf;
    nrows = min(length(groupName),3);
    wrapPlots = ceil(length(groupName)/nrows);
    wrapPlots=1;nrows=length(groupName);
end
maxPlot = 0;
for j = 1:length(groupName)
    condInds = cellfun(@(a) strcmp(a,groupName{j}), allGroups,'UniformOutput',false);
    groupPSTH = cellfun(@(t,g) mean(t(~all(isnan(t),[2 3]),:,g),3,'omitnan'),PSTH,condInds, 'UniformOutput', false);
    tUnits = cell2mat(cellfun(@(t) sum(~all(isnan(t),[2,3]).*size(t,3)>=1), groupPSTH, 'UniformOutput',false));
    subplot(nrows,wrapPlots,j);hold on;
    title(replace(strcat(groupName{j},": ", num2str(sum(tUnits)), " units, ",... %num2str(length(allGroups))," sites, "
        num2str(sum(cellfun(@sum,condInds))), " trials, ",num2str(sum(cellfun(@sum,condInds).*tUnits)), ...
        " instances"),"_", " "));
    groupPSTH = cellfun(@(g) reshape(g,[],size(g,2)),groupPSTH(tUnits~=0),'UniformOutput',false);
    if(any(cellfun(@any,condInds)))
        xAlignTicks = bins;
        avgPSTH = mean(cell2mat(groupPSTH),1,'omitnan');
        sessionSTD = cell2mat(cellfun(@(g,s)nanstd(g,0,1)./s,groupPSTH,num2cell(tUnits(tUnits~=0)),'UniformOutput',false));
        uE=avgPSTH+median(sessionSTD,1,'omitnan');
        lE=avgPSTH-median(sessionSTD,1,'omitnan');
        yP=[lE,fliplr(uE)];
        xP=[xAlignTicks,fliplr(xAlignTicks)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        d = patch(xP,yP,1);
        set(d,'edgecolor','none','facealpha',.5,'facecolor',plotColors.(groupName{j}));
        %plot(xAlignTicks,cell2mat(groupPSTH),'LineWidth',.5)%,'Color', plotColors.(groupName{j}));
        if(~isempty(yP))
            maxPlot = max(maxPlot,FRLim(end)*max(1,ceil(mean(prctile(cell2mat(groupPSTH),75,2),'omitnan')/FRLim(end))));
        end
        plot(xAlignTicks,avgPSTH, 'LineWidth',2,'Color', plotColors.(groupName{j}));
    end
end
for j = 1:length(groupName)
    groupSegs = cell2mat(cellfun(@(t,g) t(g,:), allSegs,cellfun(@(a) ...
        strcmp(a,groupName{j}), allGroups,'UniformOutput',false),'UniformOutput',false));
    avgSegs = mean(groupSegs,1,'omitnan');
    subplot(nrows,wrapPlots,j);hold on;
    c = sum(~isnan(groupSegs),2)-1;
    groupSegs(c<mode(c)-1,:) = NaN;
    c = c-(avgSegs(mode(c))==0);
    % mapC = find(~all(isnan(groupSegs),1)); groupSegs(c~=mode(c),mode(c)) = groupSegs(sub2ind(size(groupSegs),find(c~=mode(c)),c(c~=mode(c))));
    lastInd = sub2ind(size(groupSegs),1:size(groupSegs,1),c');
    for s = 1:mode(c)
        lineWidth = 1;
        lineStyle = '--';
        plotColor = 'k';
        if(avgSegs(s)>=(PSTHDisplayLimits(1) - .1) && ...
                avgSegs(s)<=PSTHDisplayLimits(end))
            if(avgSegs(s)==0)
                lineWidth = 2;
                lineStyle=':';
            elseif(s==mode(c))
                plotColor = rgb2hsv(plotColors.(groupName{j}));
                plotColor(2) = .7;
                plotColor(end) = .6;
                plotColor = hsv2rgb(plotColor);
                lineWidth = 3;
                lineStyle=':';
            end
            plot([avgSegs(s) avgSegs(s)],[0 maxPlot],...
                'Color',plotColor,'LineStyle',lineStyle,'LineWidth',lineWidth);
        end
    end
    lowCI= max(PSTHDisplayLimits(1),prctile(groupSegs(lastInd),5));
    highCI= min(PSTHDisplayLimits(end),prctile(groupSegs(lastInd),95));
    if(~isempty(lowCI) && ~isempty(highCI))
        patch([lowCI, highCI, highCI, lowCI],[0 0 maxPlot maxPlot],...
            plotColors.(groupName{j}),'FaceAlpha',0.15,'LineStyle','none');
    end
    allXTicks = bins;
    xticks([bins(1),bins(isalmost(mod(allXTicks,1),1,binSize))]);
    allLabels = [num2str(bins(1),'%.2f'), arrayfun(@num2str,...
        round(allXTicks(isalmost(mod(allXTicks,1),1,binSize))),'UniformOutput',false),...
        num2str(PSTHDisplayLimits(end),'%.2f')];
    xticklabels(allLabels);
    set(gca,'XMinorTick','on');
    set(get(gca,'XAxis'),'MinorTickValues',bins(isalmost(mod(allXTicks,1),.5,binSize/1.99)));
end
figHandle = figHandle.Children;
set(figHandle(strcmp(get(figHandle,'Type'),'axes')),'YLim',[FRLim(1),maxPlot]);
set(figHandle(strcmp(get(figHandle,'Type'),'axes')),'XLim',[bins(1)-binSize,bins(end)+binSize]);
end