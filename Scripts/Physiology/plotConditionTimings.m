function [figHandle1,figHandle2] = ...
    plotConditionTimings(conditions,masterPSTH,masterJoints,masterFRChanges,masterBaseline,activityInd)
alphaVal = 0.05;
bufferWindow = 0.5;
FRLim = [0 9];
phaseNames = {'Go', 'Task', 'Reach', 'Grasp'};
titleNames = containers.Map([1 2 3], {'All', 'Active', 'Inactive'});
repColors = struct('Arm', [.75 .75 .75], 'Hand', [.35 .35 .35], 'Forelimb',...
    [.5 .5 .5], 'Mixed', [121 76 92]./255,'Trunk', [.5 .25 .1],...
    'Face',[155 9 144]./255,'Axial', [102 0 36]./255);

repColors = struct('Arm', [.75 .75 .75], 'Hand', [.35 .35 .35], 'Trunk', [.5 .25 .1],'Face',[155 9 144]./255);
condColors = struct('ExtraSmallSphere', [1 0 0], 'LargeSphere', [1 1/1.5 0], 'Photocell', [0 0 1]);

currentJoints = unique(masterJoints);
plotJoints = fieldnames(repColors)';
jointNames = plotJoints(contains(plotJoints,unique(currentJoints)));
jointPhaseVals = cell(length(jointNames),length(phaseNames),3);
sAx = [];
figHandle1 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
figure(figHandle1);
for j = 1:length(jointNames)
    jointCondInds = cellfun(@(m) strcmp(m,jointNames{j}), masterJoints);
    jointCondPSTH = cellfun(@(a,b) a(jointCondInds), masterPSTH, 'UniformOutput', false);
    jointCondPSTH = cellfun(@(a) cell2mat(a'), jointCondPSTH, 'UniformOutput', false);
    activeJointCond = cellfun(@(a) find(a(jointCondInds)), activityInd, 'UniformOutput', false);
    inactiveJointCond = cellfun(@(a) find(~a(jointCondInds)), activityInd, 'UniformOutput', false);
    activityCondInds= {cellfun(@(j) 1:size(j,1), jointCondPSTH, 'UniformOutput', false),activeJointCond,inactiveJointCond};
    for a = 1:size(activityCondInds,2)
        sigPhase = false(1,size(masterFRChanges,1));
        currActive = activityCondInds{a}';
        for p = 1:length(phaseNames)
            FRChanges = cellfun(@(pC,bC) cellfun(@(p,b) nanmean(p./b), pC, bC, 'UniformOutput', true), masterFRChanges(:,p),masterBaseline', 'UniformOutput', false);
            FRChanges = cellfun(@(fr,ain) fr(ain), FRChanges,currActive,'UniformOutput', false);
            maxArrayLength = max(cellfun(@length,FRChanges));
            paddedFRChanges = cell2mat(cellfun(@(a) [a,NaN(1,maxArrayLength-length(a))], FRChanges, 'UniformOutput', false))';
            sigPhase(p) = anova1(paddedFRChanges,[],'off')<alphaVal;
            jointPhaseVals{j,p,a} = paddedFRChanges;
        end
        subInd = ((a-1)*length(jointNames))+j;
        sAx(subInd) = subplot(size(activityCondInds,2),length(jointNames),subInd); hold on;
        title([jointNames{j}, ',' titleNames(a)], 'Color', repColors.(jointNames{j}));
        bR = bar(cell2mat(cellfun(@nanmean, jointPhaseVals(j,:,a)','UniformOutput', false)),'FaceColor', 'flat');
        pause(0.01);
        xVals = cell2mat(arrayfun(@(bc) get(bc,'XData')' + [bc.XOffset],bR, 'UniformOutput', false));
        errorbar(reshape(xVals,1,numel(xVals)),[bR.YData],reshape(cell2mat(cellfun(@(f) ...
            nanmean(f)/(sqrt(length(f))), jointPhaseVals(j,:,a),'UniformOutput', false)),1,numel(xVals)),'k', 'linestyle', 'none');
        for c = 1:length(conditions)
            bR(c).CData = repmat(condColors.(strrep(conditions{c}, ' ','')),size(bR(c).CData,1),1);
        end
        arrayfun(@(a) sigline([xVals(a,1),xVals(a,3)], ['p < ', num2str(alphaVal)]), find(sigPhase));
    end
end

linkaxes(sAx);
xlim([bufferWindow,length(phaseNames)+bufferWindow]);
arrayfun(@(ax) xticks(ax,1:length(phaseNames)),sAx);
arrayfun(@(ax) xticklabels(ax,phaseNames),sAx);
ylim(FRLim);


sAx2 = [];
figHandle2 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
figure(figHandle2);
for p = 1:length(phaseNames)
    for a = 1:size(activityCondInds,2)
        sigJoint = false(length(jointNames),1);
        subInd = ((a-1)*length(phaseNames))+p;
        sAx2(subInd) = subplot(size(jointPhaseVals,3),length(phaseNames),subInd); hold on;
        title([phaseNames{p}, ', ' titleNames(a)]);
        jointFRChanges = [jointPhaseVals(:,p,a)'].';
        currMeans = cell2mat(cellfun(@nanmean, jointFRChanges, 'UniformOutput', false));
        bR = bar(currMeans,'FaceColor',  'flat');
        pause(0.01);
        xVals = cell2mat(get(bR,'XData'))' + [bR.XOffset];
        currErrs = cell2mat(cellfun(@(a) nanstd(a)/sqrt(length(a)),jointFRChanges, 'UniformOutput', false));
        errorbar(xVals(:),[bR.YData].',currErrs(:), 'k', 'linestyle', 'none');
        
        sigPhase = cellfun(@(an) (anova1(an,[],'off'))<alphaVal,...
            jointFRChanges,'UniformOutput', true)';
        for c = 1:length(conditions)
            bR(c).CData = repmat(condColors.(strrep(conditions{c}, ' ','')),...
                size(bR(c).CData,1),1);
        end
        arrayfun(@(a) sigline([xVals(a,1),xVals(a,3)], ['p < ', num2str(alphaVal)]), find(sigPhase));
    end
end
linkaxes(sAx2);
xlim([bufferWindow,length(jointNames)+bufferWindow])
arrayfun(@(ax) xticks(ax,1:length(jointNames)),sAx2);
arrayfun(@(ax) xticklabels(ax,jointNames),sAx2);
ylim(FRLim);