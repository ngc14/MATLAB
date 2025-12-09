function figHandle = plotJointPSTHS(params,PSTHIn,allSegsIn,allRepsIn,siteInds,...
    activityIn,PSTHDisplayLimits,FRLimIn,plotColors)
if(~exist('FRLimIn', 'var'))
    FRLim = [5 30];
else
    FRLim = FRLimIn;
end
alignmentGap = .1;
segColors = {[0 0 0],[.3 .3 .3]};
phaseWinSz = .2;
pw = {[0,phaseWinSz],[-phaseWinSz, 0]};
pa = cellstr(["GoSignal","StartHold"]);
[~,maxSegL]= max(cellfun(@length,params.condSegMap.values));
maxSegL= cell2mat(params.condSegMap.values({params.condNames(maxSegL)}));
activityInd = activityIn;
PSTH = PSTHIn;
allSegs = allSegsIn;
allReps = allRepsIn;

if(~exist('plotColors','var'))
   if(any(contains(allReps,["Reach","Grasp","Both"])))
       plotColors  = struct("Both", [1 .9 0], "Reach", [1 0 0],...
           "GraspEx",[.75 .895 .5], "Grasp", [0 0 1]);
    else
       plotColors = MotorMapping.repColors;
    end
end
plotNames = fieldnames(plotColors);
jointName = natsort(plotNames(matches(unique(allReps(allReps~="")),plotNames)));

zeroBinInd = find(params.bins==0);
binSize = params.binSize;
alignmentGap = alignmentGap/binSize;

PSTH =  cellfun(@(t,w) t(:,(fix(w(1)/binSize)+zeroBinInd):...
    fix((w(end)/binSize)+zeroBinInd)),PSTH,PSTHDisplayLimits,'UniformOutput',false);
g = groot;
figHandle = g.CurrentFigure;
reuseAxes = ~isempty(figHandle);
if(reuseAxes)
    pos = cell2mat(arrayfun(@(n) get(n,'Position'), get(figHandle,'Children'), 'UniformOutput', false));
    wrapPlots = numel(unique(pos(:,1)));
    nrows = numel(unique(pos(:,2)));
else
    figHandle = gcf;
    wrapPlots = min(length(jointName),5);
    nrows = ceil(length(jointName)/wrapPlots);
end
xAlignTicks = {};
maxPlot = 0;
trialSamples = sum(siteInds)~=length(allReps);
for j = 1:length(jointName)
    jointInds = arrayfun(@(s) cellfun(@(a) strcmp(a,jointName{j}),s), allReps,'UniformOutput',true);
    if(sum(jointInds)>0)
        jointPSTH = cellfun(@(t) t(jointInds,:), PSTH, 'UniformOutput', false);
        jointSegs = cellfun(@(t) t(strcmp(allReps,jointName{j}) & siteInds,:), allSegs,'UniformOutput',false);
        %     [~, ~, activeJointInds] = intersect(find(activityInd),jointInds(activityInd));
        %     [~,~,inactiveJoints] = intersect(find(~activityInd),find(jointInds));
        %% plot PSTHS
        subplot(nrows,wrapPlots,j);hold on;
        titleName = replace(jointName{j},"_", " ");
        if(reuseAxes)
            oldName = get(gca,'Title');
            oldName = oldName.String;
            titleName = strcat(titleName," (n= ",num2str(min(cellfun(@(s) str2num(s{1}),...
                regexp(oldName,'[=]\s(\d*)','tokens')),sum(~all(isnan(jointPSTH{1}),2)))),")");
        else
            titleName = strcat(titleName, " (n= ", num2str(sum(~all(isnan(jointPSTH{1}),2))),")");
        end
        if(trialSamples)
            titleName = strcat(titleName{:}(1:end-1), " of ", num2str(sum(jointInds)),")");
        end
        title(titleName);
        plotStart = 0;
        for a = 1:size(jointPSTH,2)
            currSegs = jointSegs{a};
            currJointAlign = jointPSTH{a};
            xAlignTicks{a} = plotStart+(1:size(currJointAlign,2));
            meanTrace = mean(currJointAlign,1,'omitnan');
            plot(xAlignTicks{a},meanTrace, 'LineWidth',2,'Color', plotColors.(jointName{j}));
            SEM = nanstd(currJointAlign,0,1);
            if(~trialSamples)
                SEM = SEM/sqrt(sum(~all(isnan(currJointAlign),2)));
            end
            uE=meanTrace+SEM;
            lE=meanTrace-SEM;
            yP=[lE,fliplr(uE)];
            xP=[xAlignTicks{a},fliplr(xAlignTicks{a})];
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            d = patch(xP,yP,1);
            set(d,'edgecolor','none','facealpha',.15,'facecolor',plotColors.(jointName{j}));
            if(~isempty(yP))
                groupMax = max(FRLim(end),FRLim(end)*ceil(quantile(yP,.8)/FRLim(end)));
            end
            avgSegs = nanmean(currSegs,1);
            if(sum(~isnan(avgSegs))==6)
                plotted = true(1,size(currSegs,2));
                pa = {};
                pw = {};
            end
            if(a==1)
                plotted = false(1,size(currSegs,2));
                maxSegNames=maxSegL;
                patches = cellfun(@(i,w) findBins(avgSegs(find(contains(maxSegNames,i),1))+w,...
                     PSTHDisplayLimits{a}(1):binSize:PSTHDisplayLimits{a}(end)), pa,pw,'UniformOutput',false);
            end
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=PSTHDisplayLimits{a}(1) && ...
                        avgSegs(s)<=PSTHDisplayLimits{a}(end) && ...
                        (~plotted(s) ||avgSegs(s)==plotStart))
                    if(avgSegs(s)==0)
                        plotColor = segColors{1};
                    else
                        plotColor = segColors{end};
                    end
                    %plotted(s) = true;
                    pSeg = find(isalmost(PSTHDisplayLimits{a}(1):binSize:...
                        PSTHDisplayLimits{a}(end),avgSegs(s),binSize/1.99),1);
                    plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[FRLim(1) groupMax],...
                        'Color',plotColor,'LineStyle','--');
                end
            end
            if(a==size(jointPSTH,2))
                allXTicks = cellfun(@(ta,pd) [find(mod(pd(1):.01:pd(end),1)==0),...
                    length(ta)],xAlignTicks,PSTHDisplayLimits,'UniformOutput',false);
                allXTicks = unique(cell2mat(allXTicks),'stable');
                xticks(allXTicks(1:1:end));
                allLabels = arrayfun(@(pd)num2str(pd,'%.2f'),...
                    unique([PSTHDisplayLimits{a}(1),ceil(PSTHDisplayLimits{a}(1)):1:...
                    floor(PSTHDisplayLimits{a}(end)),PSTHDisplayLimits{a}(end)]), 'UniformOutput', false);
                xticklabels(allLabels);
            end
            plotStart = plotStart + size(currJointAlign,2) + alignmentGap;
        end
    end
    set(gca,'XLim',[allXTicks(1), allXTicks(end)]);
    set(gca,'YLim',[FRLim(1),groupMax]);
    maxPlot = max(maxPlot,groupMax);
    cellfun(@(cr) patch([cr,fliplr(cr)],[FRLim(1) FRLim(1) groupMax groupMax],[.5 .5 .5],'FaceAlpha',.25,'EdgeColor','none'),patches);
end
%set(figHandle.Children,'YLim',[FRLim(1),maxPlot]);
end