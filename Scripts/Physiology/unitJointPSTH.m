function figHandle = unitJointPSTH(bins,PSTHIn,allRepsIn,sortVals,allSegsIn,...
    sessionInds,activityIn,PSTHDisplayLimits)
alignmentGap = .1;
FRLim = [0.2 .8];

plotColors = fieldnames(MotorMapping.repColors);
zeroBinInd = find(bins==0);
binSize = mode(diff(bins));
alignmentGap = alignmentGap/binSize;

activityInd = activityIn;
PSTH = PSTHIn;
allSegs = allSegsIn;
allReps = allRepsIn;
plotJoints =plotColors(ismember(plotColors,unique(allReps)));

PSTH =  cellfun(@(t,w) t(:,(fix(w(1)/binSize)+zeroBinInd):...
    fix((w(end)/binSize)+zeroBinInd)),PSTH,PSTHDisplayLimits,'UniformOutput',false);

xAlignTicks = {};
for j = 1:length(plotJoints)+1
    if(j>length(plotJoints))
        jointInds = true(size(allReps));
        jointSegs = allSegs;
        plotTitle = "All reps";
    else
        jointInds = arrayfun(@(s) cellfun(@(a) strcmp(a,plotJoints{j}),s), allReps,'UniformOutput',true);
        jointSegs = cellfun(@(t) t(strcmp(allReps(sessionInds),plotJoints{j})), allSegs,'UniformOutput',false);
        plotTitle = plotJoints{j};
    end
    jointPSTH = cellfun(@(t) t(jointInds,:), PSTH, 'UniformOutput', false);
    sortJoint = sortVals(jointInds);
    numUnits = nansum(cell2mat(cellfun(@(j)~any(isnan(j),2), jointPSTH,...
        'UniformOutput', false)),2)>size(jointPSTH,2)/2;
    timespan = round(nansum((cellfun(@range,PSTHDisplayLimits)./binSize)+1)...
        + (alignmentGap * (size(PSTHIn,2)-2)));
    jointIm = Inf(sum(numUnits),timespan);
    figHandle{j} = figure('Units','normalized','Position',[0 0 1 1]);
    hold on;
    title(strcat(plotTitle, ': n = ', num2str(sum(numUnits))))%, ' units, ', num2str(length(activeJointInds)), ' Active']);
    %% plot PSTHS
    if(sum(numUnits)>1)
        plotStart = 0;
        for a = 1:size(jointPSTH,2)
            currSegs = cell2mat(jointSegs{a});
            currJointAlign = jointPSTH{a}(numUnits,:);
            if(a==1)
                [sortedVals,sortInds] = sort(sortJoint(numUnits),'ascend');
            end
            xAlignTicks{a} = plotStart+[1:size(currJointAlign,2)];
            jointIm(:,round(xAlignTicks{a})) = currJointAlign(sortInds,:);
            plotStart = plotStart + size(currJointAlign,2) + alignmentGap;
        end
        imagesc(jointIm,FRLim);
        axis tight
        cm = colormap('jet');
        colormap(cm);
        plotStart = 0;
        for a = 1:size(jointPSTH,2)
            avgSegs = nanmean(cell2mat(jointSegs{a}),1);
            if(a==1)
                plotted = false(1,size(currSegs,2));
            end
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=PSTHDisplayLimits{a}(1) && ...
                        avgSegs(s)<=PSTHDisplayLimits{a}(end) && ...
                        (~plotted(s) ||avgSegs(s)==plotStart))
                    if(avgSegs(s)==0)
                        plotColor = 'k';
                    else
                        plotColor = [.75 .75 .75];
                    end
                    plotted(s) = true;
                    pSeg = find(isalmost(PSTHDisplayLimits{a}(1):binSize:...
                        PSTHDisplayLimits{a}(end),avgSegs(s),binSize/1.99),1);
                    line([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],...
                        [0 sum(numUnits)],'Color',plotColor,'LineStyle','--',...
                        'LineWidth',2);
                end
            end
            if(a==size(jointPSTH,2))
                if(sum(diff(sortVals)~=0) < length(sessionInds))
                    gaps = isinf(jointIm(1,:));
                    startgaps = find(diff(gaps)==1);
                    endgaps = find(diff(gaps)==-1)+1;
                    arrayfun(@(l) arrayfun(@(s,e) line([s e], [l l], 'Color',[.4 .4 .4],...
                        'LineStyle','-','LineWidth',1),startgaps,endgaps),sessionInds);
                end
                allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs([pd(1):binSize:pd(end)])...
                    ==min(abs([pd(1):binSize:pd(end)])))-1, ta(end)], xAlignTicks, ...
                    PSTHDisplayLimits, 'UniformOutput', false);
                set(gca,'xtick',[allXTicks{:}]);
                allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                    ""; num2str(pd(end),'%.2f')],PSTHDisplayLimits, 'UniformOutput', false);
                set(gca,'xticklabels',[allLabels{:}]);
                xtickangle(-25);
                
                ylim(gca,[1,sum(numUnits)]);
                nSteps = 4-(3*(sum(numUnits)<3)-((sum(numUnits)<3)*sum(numUnits)));
                set(gca,'ytick',linspace(1,sum(numUnits),nSteps));
                set(gca,'yticklabels',round(sortedVals(fix(linspace(1,length(sortedVals),nSteps))),4));
                set(gca,'YDir','reverse');
                cb = colorbar;
                cb.Location = 'westoutside';
                imChild = allchild(gca);
                set(imChild(strcmp(get(imChild,'Type'),'image')),'AlphaData', ~isinf(jointIm))
            end
            plotStart = plotStart + size(currJointAlign,2) + alignmentGap;
        end
    end
end
end