close all;
clear all;
conds = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};
graspEvents = {'TrialStart','GoSignal','StartReach','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward',};
nonGraspEvents = {'TrialStart','GoSignal','StartReach','StartHold','StartWithdraw',...
    'StartReplaceHold','StartReplaceSuccess','StartReward',};
condSegs = {graspEvents, graspEvents, nonGraspEvents};
muscles = {'Deltoid','Bicep', 'Tricep', 'Extensor', 'ProximalRadialFlexor',...
    'D4D5Extensor','DistalFlexor'};
printNames = { 'Deltoid', 'Biceps', 'Triceps','WristEx', 'WristFlex', 'D4D5Ex', 'DigitFlex'};
directory = dir('S:\Lab\Gilligan\Sorted EMG Data\Results\Models\Averages_Raw\');
alignSegs = {'GoSignal'};
xLims = {[-.15 1.5]};
gap = .25;
Fs = 2000;
xLims = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);

% figure('units','normalized','outerposition',[0 0 1 1]);
% ax1 = axes('Position',[0 0 1 1],'Visible','off');
% axes(ax1)
% tx = text(.03,0.5,'Voltage (uV)', 'FontSize', 26, 'FontWeight', 'bold', 'HorizontalAlignment', 'Center');
% set(tx, 'Rotation', 90);
% text(.5,.03,'Time (ms)', 'FontSize', 26, 'FontWeight', 'bold','HorizontalAlignment', 'Center');
% axes(ax1);
% hf = tight_subplot(length(muscles),1,[.03, 0.03],[.09 .07],[.06 .04]);
figure('units','normalized','outerposition',[0 0 1 1]);
ax2 = axes('Position',[0 0 1 1],'Visible','off');

names = {directory.name};
firstInd = cellfun(@(a) regexp(a, '_', 'split'), names, 'UniformOutput', false);
firstInd = cellfun(@(a) a{1}, firstInd,'UniformOutput', false);
numMuscles = length(unique(firstInd(3:end)));

hf = cellfun(@(a) figure(), muscles, 'UniformOutput', false);
hf = cellfun(@gca, hf);
arrayfun(@(ax) hold(ax, 'on'), hf, 'UniformOutput', false);
axes(ax2)
hm = tight_subplot(2,2,[.1 .06],[.1 .1],[.07 .01]);
arrayfun(@(ax) hold(ax, 'on'), hm, 'UniformOutput', false);
segNames = {'RT', 'R', 'L', 'H', 'W'};
times = [];
timeName = {};

for c = 1:length(conds)
    cond = conds{c};
    timeFiles = find(cellfun(@(a) contains(a, 'Times') & contains(a, cond), names));
    for s = 1:length(timeFiles)
        time = load([directory(timeFiles(s)).folder, '\', directory(timeFiles(s)).name]);
        timeName{end+1} = directory(timeFiles(s)).name;
        if(strcmp(cond, 'Photocell'))
            time = nanmean(time.allSegs,1);
            time = [time(1:3),NaN, time(4:end)];
            times(end+1,:) = time;
        else
            times(end+1,:) = nanmean(time.allSegs,1);
        end
    end
end
[timeName, ind, ~] = unique(timeName);
times = times(ind,:);
colors = {'r', 'g', 'b', 'y'};
colorM = colormap('colorcube');
colorM = [colorM(1,:);colorM(4,:);colorM(7,:);colorM(10,:);colorM(13,:);...
    colorM(21,:);colorM(27,:);colorM(44,:);colorM(50,:);];


for d=1:length(muscles)
    for c = 1:length(conds)
        cond = conds{c};
        EMGfiles = find(cellfun(@(a) endsWith(a, ['_', cond '.mat']), names));
        di = find(startsWith({directory(EMGfiles).name}, muscles{d}));
        ti = find(cellfun(@(a) strcmp(a,[muscles{d},'_',cond,'_Times.mat']), names));
        tMusc = find(contains(timeName, ['_',cond,'_Times.mat']));
        tCond = find(contains(timeName, muscles{d}));
        sig = load([directory(EMGfiles(di)).folder, '\', directory(EMGfiles(di)).name]);
        sig = sig.rawTrials;
        currTimes = load([directory(ti).folder, '\', directory(ti).name]);
        currTimes = currTimes.allSegs;
        muscTimes = times(tMusc,:);
        condTimes = times(tCond,:);
        %condTimes = condTimes([1 4 2 3],:);
        plotStart = 0;
        xTicks = [];
        labels = [];
        for align = 1:length(alignSegs)
            if(align==1)
                plotStart = 0;
            else
                plotStart = plotStart + (gap*Fs) +xLims{align-1}(2) + abs(xLims{align}(1));
            end
            alignSegInd = find(strcmp(condSegs{c},alignSegs{align}))-1;
            alignedSigs = sig(:,currTimes(:,alignSegInd) + xLims{align}(1):...
                currTimes(:,alignSegInd)+xLims{align}(2));
            avgSig = nanmean(alignedSigs,1);
            xVals = plotStart+xLims{align}(1):plotStart+xLims{align}(2);
            plot(hf(d),xVals,avgSig,'Color',colors{c}, 'LineWidth', 2);
            plot(hf(d),xVals,avgSig-nanstd(alignedSigs),'Color',colors{c}, 'LineWidth', .3,'LineStyle', '--');
            plot(hf(d),xVals,avgSig+nanstd(alignedSigs),'Color',colors{c}, 'LineWidth', .3,'LineStyle', '--');
            plot(hm(c),xVals,avgSig,'Color',colorM(d,:), 'LineWidth', 2);
            
            [maxVal, maxInd] = max(avgSig);
            maxInd = xVals(maxInd);
            %             scatter(hf(d),maxInd, maxVal, 120, colors{c}, 'filled');
            %             scatter(hm(c),maxInd, maxVal, 75, colorM(d,:), 'filled');
            
            beforeMSegs = plotStart + cumsum(diff(muscTimes(:,alignSegInd:-1:1),1,2),2 ,'omitnan', 'reverse');
            afterMSegs = plotStart + cumsum(diff(muscTimes(:,alignSegInd:end),1,2),2, 'omitnan');
            avgMSegs = median([beforeMSegs repmat(plotStart,size(beforeMSegs,1),1) afterMSegs]);
            
            beforeCSegs = plotStart + cumsum(diff(condTimes(1:2,alignSegInd:-1:1),1,2),2, 'omitnan', 'reverse');
            afterCSegs = plotStart + cumsum(diff(condTimes(1:2,alignSegInd:end),1,2),2, 'omitnan');
            avgCSegs = median([beforeCSegs repmat(plotStart,size(beforeCSegs,1),1) afterCSegs]);
            for s = 1:length(avgMSegs)
                beforePSegs = plotStart + cumsum(diff(condTimes(3,alignSegInd+1:-1:1),1,2),2, 'omitnan', 'reverse');
                afterPSegs = plotStart + cumsum(diff(condTimes(3,alignSegInd+1:end),1,2),2, 'omitnan');
                avgPSegs = [beforePSegs repmat(plotStart,size(beforePSegs,1),1) afterPSegs];
                
                if(avgCSegs(s)>=plotStart && avgCSegs(s)<=plotStart+xLims{align}(2))
                    plot(hf(d),[avgCSegs(s) avgCSegs(s)],[0 150],'k--')
                end
                if(avgMSegs(s)>=plotStart && avgMSegs(s)<=plotStart+xLims{align}(2) && d==1)
                    plot(hm(c),[avgMSegs(s) avgMSegs(s)],[0 150],'k--')
                end
            end
            xTicks = [xTicks, plotStart+[xLims{align}(1) 0 xLims{align}(2)]];
            labels = [labels, [xLims{align}(1) 0 xLims{align}(2)]];
        end
        
        if(c==1)
            title(hf(d), printNames{d});
            set(hf(d), 'YLim', [0 150]);
            set(hf(d), 'FontSize', 24);
            set(hf(d), 'XTick', xTicks);
            set(hf(d), 'XTickLabel', arrayfun(@num2str, (labels./Fs*1000), 'UniformOutput', false));
            for fL = 1:length(conds)
                hfL(fL) = plot(hf(d),NaN,NaN,'Color', colors{fL}, 'LineWidth', 2);
            end
            legend(hf(d),hfL, conds, 'FontSize', 14, 'Location', 'north', 'Orientation', 'horizontal', 'AutoUpdate', 'off');
        end
        
        if(d==1)
            title(hm(c), cond);
            set(hm(c), 'FontSize', 24);
            set(hm(c), 'YLim', [0 150]);
            set(hm(c), 'XTick', xTicks);
            set(hm(c), 'XTickLabel', arrayfun(@num2str, (labels./Fs*1000), 'UniformOutput', false));
            if(c==1)
                for mL = 1:size(printNames,2)
                    hmL(mL) = plot(hm(c),NaN,NaN,'Color', colorM(mL,:), 'LineWidth', 2);
                end
                legend(hm(c),hmL, printNames, 'FontSize', 14, 'Location', 'north', 'Orientation', 'horizontal','AutoUpdate', 'off');
            end
        end
        
    end
end