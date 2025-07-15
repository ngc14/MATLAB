function figHandle = plotJointPSTHS(bins,PSTHIn,allSegsIn,allRepsIn,siteInds,...
    activityIn,PSTHDisplayLimits,FRLimIn,plotColors)
if(~exist('FRLimIn', 'var'))
    FRLim = [5 30];
else
    FRLim = FRLimIn;
end
wrapPlots = 3;
alignmentGap = .1;
segColors = ['k','r'];

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
jointName =plotNames(matches(unique(allReps),plotNames));

zeroBinInd = find(bins==0);
binSize = mode(diff(bins));
alignmentGap = alignmentGap/binSize;

PSTH =  cellfun(@(t,w) t(:,(fix(w(1)/binSize)+zeroBinInd):...
    fix((w(end)/binSize)+zeroBinInd)),PSTH,PSTHDisplayLimits,'UniformOutput',false);
figHandle=figure('Units', 'normalized', 'Position', [0 0 1 1]); hold on;
xAlignTicks = {};
for j = 1:length(jointName)
    jointInds = arrayfun(@(s) cellfun(@(a) strcmp(a,jointName{j}),s), allReps,'UniformOutput',true);
    jointPSTH = cellfun(@(t) t(jointInds,:), PSTH, 'UniformOutput', false);
    jointSegs = cellfun(@(t) t(strcmp(allReps(siteInds),jointName{j}),:), allSegs,'UniformOutput',false);
    %     [~, ~, activeJointInds] = intersect(find(activityInd),jointInds(activityInd));
    %     [~,~,inactiveJoints] = intersect(find(~activityInd),find(jointInds));
    %% plot PSTHS
    subplot(ceil(length(jointName)/wrapPlots),wrapPlots,j);hold on;
    titleName = replace(jointName{j},"_", " ");
    titleName = strcat(titleName, " (n= ", num2str(size(jointSegs{1},1)),")");
    if(sum(siteInds)~=length(allReps))
        titleName = strcat(titleName{:}(1:end-1), " of ", num2str(sum(jointInds)),")");
    end
    title(titleName);
    if(sum(jointInds)>0)
        plotStart = 0;
        for a = 1:size(jointPSTH,2)
            currSegs = jointSegs{a};
            currJointAlign = jointPSTH{a};
            xAlignTicks{a} = plotStart+(1:size(currJointAlign,2));
            plot(xAlignTicks{a},nanmean(currJointAlign,1), 'LineWidth',2,'Color', plotColors.(jointName{j}));
            %             if(~isempty(activeJointInds))
            %                 plot(xAlignTicks{a},nanmean([jointPSTH{a}(activeJointInds,:)],1),...
            %                     'LineWidth',2,'Color', plotColors.(jointName{j}), 'LineStyle','--');
            %             end
            %             if(~isempty(inactiveJoints))
            %                 plot(xAlignTicks{a},nanmean([jointPSTH{a}(inactiveJoints,:)],1),...
            %                     'LineWidth',2,'Color', plotColors.(jointName{j}), 'LineStyle','--');
            %             end
            uE=nanmean(currJointAlign,1)+(nanstd(currJointAlign,1)/sqrt(sum(jointInds)));
            lE=nanmean(currJointAlign,1)-(nanstd(currJointAlign,1)/sqrt(sum(jointInds)));
            yP=[lE,fliplr(uE)];
            xP=[xAlignTicks{a},fliplr(xAlignTicks{a})];
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            d = patch(xP,yP,1);
            set(d,'edgecolor','none','facealpha',.5,'facecolor',plotColors.(jointName{j}));
            
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
                    plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 2*FRLim(end)],...
                        'Color',plotColor,'LineStyle','--');
                end
            end
            if(a==size(jointPSTH,2))
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
    if(max(yP)>FRLim(end))
        ylim([FRLim(1),FRLim(end)*2]);
    else
        ylim(FRLim);
    end
end
end