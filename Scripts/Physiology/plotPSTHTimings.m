function [figHandle1,figHandle2,figHandle3] = ...
    plotPSTHTimings(bins,masterPSTH,masterJoints,baselineFR,phaseBins,phaseFR,activityInd,alignmentPointName)
startTaskBin = find(bins==-1)-1; % 1000 ms before alignment
taskBinEndBuffer = 50; %500 ms after start of hold
alphaVal = 0.05;
wrapPlots = 4;
phaseCompare = {{'Reach', 'Grasp'}};
jointCompare = {{'Arm', 'Hand'}};
phaseNames = {'Go', 'Task', 'Reach', 'Grasp'};
repColors = struct('Arm', [.75 .75 .75], 'Hand', [.35 .35 .35], 'Forelimb',...
    [.5 .5 .5], 'Mixed', [121 76 92]./255,'Trunk', [.5 .25 .1],...
    'Face',[155 9 144]./255,'Axial', [102 0 36]./255);
repColors = struct('Arm', [.75 .75 .75], 'Hand', [.35 .35 .35], 'Trunk', [.5 .25 .1],'Face',[155 9 144]./255);

plotReps = fieldnames(repColors)';

AUCLim = 2200;
FRLim = 10;
FRLim = FRLim/1.2;
timeLim = .45;
timeLim = timeLim/1.2;

figHandle1 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
figHandle2 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
figHandle3 = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
[maxTP, maxTPErrs, maxD1, maxD1Errs, minD1, minD1Errs, AUCs, AUCErrs] = deal([]);
allActiveJointInds = containers.Map();
jointVals = containers.Map();

activityInd = activityInd(~cellfun(@isempty,masterPSTH));
masterPSTH = masterPSTH(~cellfun(@isempty,masterPSTH));
masterJoints = masterJoints(~cellfun(@isempty,masterJoints));
phaseFR = phaseFR(~cellfun(@isempty,phaseFR));
currentJoints = contains(plotReps,unique(masterJoints));
jointName = plotReps(currentJoints);
[sigPeakActivity,sigMinActivity,sigMaxActivity,sigAUC] = deal(false(1,length(jointName)));
for j = 1:length(jointName)
    jointInds = cellfun(@(a) strcmp(a,jointName{j}), masterJoints);
    jointPSTH = masterPSTH(jointInds);
    jointPSTH = cell2mat(jointPSTH');
    activeJointInds = find(activityInd(jointInds));
    inactiveJointInds = find(~activityInd(jointInds));
    
    [allFR, allErr]= deal([]);
    phaseVals = containers.Map();
    activePhaseInds = containers.Map();
    sigActivity = false(1,length(phaseNames));
    for p = 1:length(phaseNames)
        FRChange = cellfun(@(p,b) nanmean(p./b), phaseFR{p}(jointInds), baselineFR(jointInds));%cell2mat(phaseFR{p})./cell2mat(baselineFR');
        if(~isempty(activeJointInds) & ~isempty(inactiveJointInds))
            [~, pVal] = ttest2(FRChange(inactiveJointInds(~isnan(FRChange(inactiveJointInds)))), ...
                FRChange(activeJointInds(~isnan(FRChange(activeJointInds)))));
            sigActivity(p) = pVal < alphaVal;
        end
        if(contains(phaseNames{p}, [phaseCompare{:}]))
            phaseVals(phaseNames{p}) = FRChange;
            activePhaseInds(phaseNames{p}) = activeJointInds;
        end
        allFR(p,:) = [nanmean(FRChange),nanmean(FRChange(activeJointInds)),nanmean(FRChange(inactiveJointInds))];
        allErr(p,:) = [nanstd(FRChange)/sqrt(sum(~isnan(FRChange))), ...
            nanstd(FRChange(activeJointInds))/sqrt(sum(~isnan(FRChange(activeJointInds)))),...
            nanstd(FRChange(inactiveJointInds))/sqrt(sum(~isnan(FRChange(inactiveJointInds))))];
    end
    %% plot phase changes
    figure(figHandle1)
    subplot(ceil(length(jointName)/wrapPlots),wrapPlots,j);hold on;
    title([jointName{j}, ': n = ', num2str(sum(jointInds)), ' sites, ', num2str(length(activeJointInds)), ' Active']);
    b = bar(allFR, 1, 'FaceColor', repColors.(jointName{j}));
    nbars = size(allFR,2);
    pause(0.01);
    xVals = arrayfun(@(bx) bx.XData+bx.XOffset,b, 'UniformOutput', false);
    for i = 1:nbars
        errorbar(xVals{i}', allFR(:,i), allErr(:,i), 'k', 'linestyle', 'none');
    end
    if(sum(strcmp({b.Type}, 'bar'))==3)
        b(2).LineStyle = '--';
        b(2).LineWidth = 3;
        b(3).LineStyle = ':';
        b(3).LineWidth = 3;
    end
    for s = 1:length(phaseNames)
        if(sigActivity(s))
            sigline([xVals{2}(s) xVals{3}(s)], sprintf('AC/IN \n p < %0.2f', alphaVal), -1)
        end
    end
    for ps = 1:length(phaseCompare)
        phaseG1 = phaseCompare{ps}(1);
        phaseG2 = phaseCompare{ps}(end);
        allVals1 = values(phaseVals,phaseG1);
        allVals2 = values(phaseVals,phaseG2);
        allVals1 = allVals1{:};
        allVals2 = allVals2{:};
        for a = 1:3
            if(a==1)
                currInds1 = ~isnan(allVals1);
                currInds2 = ~isnan(allVals2);
                sigName = 'ALL';
            elseif(a==2)
                currInds1 = values(activePhaseInds,phaseG1);
                currInds2 = values(activePhaseInds,phaseG2);
                activeInds1 = currInds1{:};
                activeInds2 = currInds2{:};
                currInds1 = activeInds1(~isnan(allVals1(activeInds1)));
                currInds2 = activeInds2(~isnan(allVals2(activeInds2)));
                sigName = 'ACTIVE';
            else
                currInds1 = setdiff(1:length(allVals1),activeInds1);
                currInds2 = setdiff(1:length(allVals2),activeInds2);
                currInds1 = currInds1(~isnan(allVals1(currInds1)));
                currInds2 = currInds2(~isnan(allVals2(currInds2)));
                sigName = 'INACTIVE';
            end
            if(~isempty(currInds1) & ~isempty(currInds2))
                [~, pVal] = ttest2(allVals1(currInds1), allVals2(currInds2));
                if(pVal < alphaVal)
                    sigline([xVals{a}(find(strcmp(phaseNames,phaseG1))),...
                        xVals{a}(find(strcmp(phaseNames,phaseG2)))],...
                        sprintf('%s \n  p < %0.2f', sigName,alphaVal), ...
                        str2double(sprintf('1.%d',(3*(a-1))))*min(FRLim-1,max(...
                        allFR(find(contains(phaseNames,[phaseG1,phaseG2])),:),[],'all')));
                end
            end
        end
    end
    xticks(1:length(phaseNames));
    xticklabels(phaseNames);
    yticks(0:1:10);
    if(j>1)
        yticklabels([]);
    end
    ylim([0 FRLim+1]);
    %% store info for PSTH timings
    PSTHTimes = jointPSTH;
    currBins = cellfun(@(b) round(nanmean(b,1)), phaseBins{2}(jointInds), 'UniformOutput', false);
    AUCJoints = cellfun(@(a,b) trapz(abs(a(max(1,b(1)):min(b(end),length(a))))), num2cell(PSTHTimes,2)', currBins, 'UniformOutput', true);
    AUCs(j,:) = [nanmean(AUCJoints), nanmean(AUCJoints(activeJointInds)), nanmean(AUCJoints(inactiveJointInds))];
    AUCErrs(j,:) = [nanstd(AUCJoints)/sqrt(sum(~isnan(AUCJoints))),...
        nanstd(AUCJoints(activeJointInds))/sqrt(sum(~isnan(AUCJoints(activeJointInds)))),...
        nanstd(AUCJoints(inactiveJointInds))/sqrt(sum(~isnan(AUCJoints(inactiveJointInds))))];
    % restrict to task window
    PSTHTimes(:,1:startTaskBin) = NaN;
    for u = 1:size(PSTHTimes,1)
        PSTHTimes(u,(max(1,currBins{u}(end))+taskBinEndBuffer):end) = NaN;
    end
    badPSTHS = ~any(~isnan(PSTHTimes),2);
    
    [~,maxInds] = max(abs(PSTHTimes),[],2);
    maxInds = bins(maxInds);
    maxInds(badPSTHS) = NaN;
    maxTP(j,:) = [nanmean(maxInds),nanmean(maxInds(activeJointInds)),nanmean(maxInds(inactiveJointInds))];
    maxTPErrs(j,:) = [nanstd(maxInds)/sqrt(sum(~isnan(maxInds))),...
        nanstd(maxInds(activeJointInds))/sqrt(sum(~isnan(maxInds(activeJointInds)))),...
        nanstd(maxInds(inactiveJointInds))/sqrt(sum(~isnan(maxInds(inactiveJointInds))))];
    
    d1PSTHS = diff(PSTHTimes,[],2);
    [~,d1MaxInds] = max(d1PSTHS,[],2);
    d1MaxInds = bins(d1MaxInds);
    d1MaxInds(badPSTHS) = NaN;
    maxD1(j,:) = [nanmean(d1MaxInds),nanmean(d1MaxInds(activeJointInds)),nanmean(d1MaxInds(inactiveJointInds))];
    maxD1Errs(j,:) = [nanstd(d1MaxInds)/sqrt(sum(~isnan(d1MaxInds))),...
        nanstd(d1MaxInds(activeJointInds))/sqrt(sum(~isnan(d1MaxInds(activeJointInds)))),...
        nanstd(d1MaxInds(inactiveJointInds))/sqrt(sum(~isnan(d1MaxInds(inactiveJointInds))))];
    for u = 1:size(d1PSTHS,1)
        if(~badPSTHS(u))
            d1PSTHS(u,startTaskBin:d1MaxInds(u)) = NaN;
        end
    end
    [~,d1MinInds] = min(d1PSTHS,[],2);
    d1MinInds = bins(d1MinInds);
    d1MinInds(~any(~isnan(d1PSTHS),2)) = NaN;
    minD1(j,:) = [nanmean(d1MinInds),nanmean(d1MinInds(activeJointInds)),...
        nanmean(d1MinInds(inactiveJointInds))];
    minD1Errs(j,:) = [nanstd(d1MinInds)/sqrt(sum(~isnan(d1MinInds))),...
        nanstd(d1MinInds(activeJointInds))/sqrt(sum(~isnan(d1MinInds(activeJointInds)))),...
        nanstd(d1MinInds(inactiveJointInds))/sqrt(sum(~isnan(d1MinInds(inactiveJointInds))))];
    
    if(~isempty(activeJointInds) & ~isempty(inactiveJointInds))
        [~, pVal] = ttest2(maxInds(inactiveJointInds(~isnan(maxInds(inactiveJointInds)))),...
            maxInds(activeJointInds(~isnan(maxInds(activeJointInds)))));
        sigPeakActivity(j) = pVal < alphaVal;
        [~, pVal] = ttest2(d1MaxInds(inactiveJointInds(~isnan(d1MaxInds(inactiveJointInds)))),...
            d1MaxInds(activeJointInds(~isnan(d1MaxInds(activeJointInds)))));
        sigMaxActivity(j) = pVal < alphaVal;
        [~, pVal] = ttest2(d1MinInds(inactiveJointInds(~isnan(d1MinInds(inactiveJointInds)))),...
            d1MinInds(activeJointInds(~isnan(d1MinInds(activeJointInds)))));
        sigMinActivity(j) = pVal < alphaVal;
        [~, pVal] = ttest2(AUCJoints(inactiveJointInds(~isnan(AUCJoints(inactiveJointInds)))),...
            AUCJoints(activeJointInds(~isnan(AUCJoints(activeJointInds)))));
        sigAUC(j) = pVal < alphaVal;
    end
    if(contains(jointName{j}, [jointCompare{:}]))
        jointVals(jointName{j}) = [{maxInds},{d1MaxInds},{d1MinInds},{AUCJoints}];
        allActiveJointInds(jointName{j}) = activeJointInds;
    end
end

%% plot PSTH timings
for f = 1:4
    if(f==1)
        graphName = 'PSTH Peaks';
        plotVals = maxTP;
        errVals = maxTPErrs;
        activityComp = sigPeakActivity;
    elseif(f==2)
        graphName = 'PSTH Max Rise Time';
        plotVals = maxD1;
        errVals = maxD1Errs;
        activityComp = sigMaxActivity;
    elseif(f==3)
        graphName = 'PSTH Max Fall Time';
        plotVals = minD1;
        errVals = minD1Errs;
        activityComp = sigMinActivity;
    elseif(f==4)
        graphName = 'AUC of Task Window';
        plotVals = AUCs;
        errVals = AUCErrs;
        activityComp = sigAUC;
    end
    if(f<=3)
        figure(figHandle2);
        subplot(1,3,f);
        labelVal = min(-0.01,min(plotVals(j,:),[],'all'));
        yAxLim = timeLim;
    else
        figure(figHandle3);
        labelVal = -250;
        yAxLim = AUCLim;
    end
    title(graphName);
    hold on;
    b = bar(plotVals,1, 'FaceColor', 'flat');
    nbars = size(plotVals,2);
    pause(0.01);
    xVals = arrayfun(@(bx) bx.XData+bx.XOffset,b, 'UniformOutput', false);
    for i = 1:nbars
        errorbar(xVals{i}', plotVals(:,i), errVals(:,i), 'k', 'linestyle', 'none');
        b(i).CData = cell2mat(cellfun(@(a) repColors.(a), jointName, 'UniformOutput', false)');
    end
    if(sum(strcmp({b.Type}, 'bar'))==3)
        b(2).LineStyle = '--';
        b(2).LineWidth = 3;
        b(3).LineStyle = ':';
        b(3).LineWidth = 3;
    end
    for j = 1:length(jointName)
        if(activityComp(j))
            sigline([xVals{2}(j), xVals{3}(j)], sprintf('AC/IN \n p < %0.2f', alphaVal), labelVal);
        end
    end
    for jc = 1:length(jointCompare)
        jointG1 = jointCompare{jc}(1);
        jointG2 = jointCompare{jc}(end);
        allVals1 = values(jointVals,jointG1);
        allVals2 = values(jointVals,jointG2);
        allVals1 = allVals1{1}{f};
        allVals2 = allVals2{1}{f};
        for a = 1:3
            if(a==1)
                currInds1 = ~isnan(allVals1);
                currInds2 = ~isnan(allVals2);
                sigName = 'ALL';
            elseif(a==2)
                currInds1 = values(allActiveJointInds,jointG1);
                currInds2 = values(allActiveJointInds,jointG2);
                activeInds1 = currInds1{:};
                activeInds2 = currInds2{:};
                currInds1 = activeInds1(~isnan(allVals1(activeInds1)));
                currInds2 = activeInds2(~isnan(allVals2(activeInds2)));
                sigName = 'ACTIVE';
            else
                currInds1 = setdiff(1:length(allVals1),activeInds1);
                currInds2 = setdiff(1:length(allVals2),activeInds2);
                currInds1 = currInds1(~isnan(allVals1(currInds1)));
                currInds2 = currInds2(~isnan(allVals2(currInds2)));
                sigName = 'INACTIVE';
            end
            if(~isempty(currInds1) & ~isempty(currInds2))
                [~, pVal] = ttest2(allVals1(currInds1), allVals2(currInds2));
                if(pVal < alphaVal)
                    sigline([xVals{a}(find(strcmp(jointName,jointG1))),...
                        xVals{a}(find(strcmp(jointName,jointG2)))],...
                        sprintf('%s \n p < %0.2f', sigName,alphaVal), ...
                        str2double(sprintf('1.%d',2*(a-1)))*min(yAxLim,max(...
                        plotVals(find(contains(jointName,[jointG1,jointG2])),:),[],'all')));
                end
            end
        end
        xticks(1:length(jointName));
        xticklabels(jointName);
        if(f<=3)
            ylabel(['Seconds from ', alignmentPointName]);
            ylim([-.5, .5]);
            yticks(-.5:.25:.5);
        else
            ylim([0 AUCLim]);
            yticks(0:500:AUCLim);
        end
    end
end