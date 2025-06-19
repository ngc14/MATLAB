clear all;
close all;
monkey = 'Gilligan';
conds = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};
events = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess'}});
modelNames = {'Deltoid','Bicep', 'Tricep', 'Extensor', 'ProximalRadialFlexor',...
    'D4D5Extensor','DistalFlexor'};
modelNames = {'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor',...
    'Digit Extensor','Digit Flexor'};
groupings = {{'Deltoid'},{'Biceps'},{'Triceps'},{'Wrist Extensor'},{'Wrist Flexor'},{'Digit Extensor'},{'Digit Flexor'}};
  groupings = {{'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor','Digit Extensor','Digit Flexor'}};
printNames = { 'Deltoid', 'Biceps', 'Triceps','WristEx', 'WristFlex', 'DigitEx', 'DigitFlex'};
Fs = 2000;
pVal = 0.05;
phaseWindow = .5;

phaseWindow = phaseWindow*Fs;
alignPoint = ["GoSignal","StartReplaceHold"];
alignWindow = [0 0];
comps = combnk(1:length(conds),2);
comps(:, [1 2]) = comps([1 3 2], [2 1]);
directory = dir(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\']);

[times,allTrialDates,allAUC] = deal(repmat({[]},length(conds),length(modelNames)));
[traceAvg,traceErr] = deal(NaN(length(conds),length(modelNames)));
for d=1:length(modelNames)
    EMGfiles = find(cellfun(@(a) endsWith(a, '.mat'), {directory.name}));
    di = find(startsWith({directory(EMGfiles).name}, modelNames{d}));
    s = matfile([directory(EMGfiles(di)).folder, '\', directory(EMGfiles(di)).name]);
    conditions = s.Conditions;
    for c = 1:length(conds)
        %         trials = load([trialDir, modelNames{d}, '_', conds{c},'.mat']);
        %         trials = trials.rawTrials;
        %         trials = cumsum([1 diff(cellfun(@length,trials))~=0]);
        trialArea =[];
        cond = find(strcmp(conditions, conds{c}));
        datesR = s.dates(cond,1);
        datesR = datesR{:};
        currDates = cellfun(@(d) datetime(d,'InputFormat','MM_dd_uuuu'), datesR);
        currSegEvents = values(events,conds(c));
        currSegEvents = currSegEvents{1};

        muscTimes = s.segs(c,1);
        muscTimes = cellfun(@(m) m(1:min(length(m),length(currSegEvents))), muscTimes{1}, 'UniformOutput',false);
        muscTimes = cellfun(@(m) [m(1:2), NaN(1,double(length(m)<length(currSegEvents))), m(3:end)], muscTimes,'UniformOutput',false);

        dateTrials = find(currDates>=datetime(2018,12,01) & currDates<datetime(2021,01,01));
        muscTimes = vertcat(muscTimes{dateTrials});
        meanMuscTimes = nanmean(muscTimes);
        sig =  s.rawTrials(c,1);
        sig = sig{1}(dateTrials);

        %         sig = abs(sig);
        %         sig = filter(1/100*ones(1,100),1,sig, [],2);
        if(iscell(sig))
%             maxL = max(cellfun(@length,sig));
%             sig = cellfun(@(s) nanPad(s,maxL),sig,'UniformOutput', false);
%             sig = cell2mat(sig');
        end

        alignRange = horzcat(muscTimes(:,find(strcmp(currSegEvents,alignPoint(1)))) +alignWindow(1),...
            muscTimes(:,find(strcmp(currSegEvents,alignPoint(end))))+ alignWindow(end));

        badTimes = any(isnan(alignRange),2) | all(isnan(muscTimes),2);
        currDates = currDates(~badTimes);
        sig = sig(~badTimes);
        alignRange = alignRange(~badTimes,:);
        timesC{c} = muscTimes(~badTimes,:);
        muscTimes = [];
        for t = 1:size(sig,2)
            trialArea(t) = trapz(sig{t}([alignRange(t,1):alignRange(t,2)]))/Fs;

            %/sum(~isnan(sig(s,(startTime(s):stopTime(s)))));
        end
        %
        traceAvgC(c) =  nanmean(trialArea);
        traceErrC(c) = nanstd(trialArea,[],2)./sqrt(size(trialArea,2));
        allTrialDatesC{c} = currDates; 
        allAUCc{c} = trialArea;
        disp([c,d]);
    end
times(:,d) = timesC;
traceAvg(:,d) = traceAvgC';
traceErr(:,d) = traceErrC';
allTrialDates(:,d) = allTrialDatesC;
allAUC(:,d) = allAUCc;
end
%%
AUCs = [];
AUCerr=[];
for d = 1:size(allTrialDates,2)
    [trialIDs, utInd]= unique([allTrialDates{:,d}]);
    matchSegs = cellfun(@(cc) arrayfun(@(s) cc==s,trialIDs,'UniformOutput',false), ...
    allTrialDates(:,d),'UniformOutput',false);
    matchSegs = vertcat(matchSegs{:});
    normAUC = {};
    for s = 1:length(trialIDs)
        matchDates = nanmean(cell2mat(cellfun(@(a,b) a(b), allAUC(:,d), matchSegs(:,s),'UniformOutput', false)'));
        normAUC{s} = cellfun(@(a,b) a(b)/matchDates, allAUC(:,d), matchSegs(:,s),'UniformOutput', false);
    end
    allNAUC(:,d) = cellfun(@cell2mat, num2cell([normAUC{:}],2), 'UniformOutput', false);
end
%%
AUCS = {};
allMWeights = cellfun(@length,allNAUC)./sum(cellfun(@length,allNAUC),2);
allMs = mean(sum(cellfun(@(m) mean(m,'omitnan'),allNAUC).*allMWeights,2),'omitnan');
for c = 1:length(conds)
    for g = 1:length(groupings)
        groupAssignment = find(cellfun(@(a) ismember(a, groupings{g}), modelNames));
        maxGLength = max(cellfun(@length, allNAUC(c,groupAssignment)));
        %groupAUCS = nansum(cell2mat(cellfun(@(ga) [ga NaN(1,maxGLength-length(ga))], allNAUC(c,groupAssignment),'UniformOutput',false)'),1);
        %AUCs{c} = groupAUCS(~isnan(groupAUCS) & groupAUCS>0);
        %AUCs(c,g) = mean(groupAUCS(~isnan(groupAUCS)));
        %AUCerr(c,g) = nanstd(groupAUCS(~isnan(groupAUCS)));
        groupAUCS = allNAUC(c,groupAssignment);
        mWeights(c,g) = sum(cellfun(@length, groupAUCS));
        mu(c,g) = sum(cellfun(@(g) mean(g,'omitnan'), groupAUCS).*(mWeights(c,g)./sum(mWeights(c,g))));
        sumSqr(c,g) = sum(cellfun(@(a) sumsqr(a), groupAUCS) .* (mWeights(c,g)./sum(mWeights(c,g))));
    end
end
Y = sum(sumSqr);
T = sum((mu.^2).*mWeights);
A = sum(mWeights)*sum(mu.*(mWeights./sum(mWeights)))^2;
SSW = Y-A;
SSB = A-T;
SST = Y-T;
MSB = sum((AUCs-allMs).^2)/(length(conds)-1);
MSW = sum(AUCerr)/(mean(sum(cellfun(@length,allNAUC),1))-2);
F = MSW/MSB;

% %allAreas = cellfun(@(a,b) a./b', allAreas, allMeans, 'UniformOutput', false);
% AUCs = cellfun(@nanmean, trialAreas);
% AUCerr = cellfun(@(a) nanstd(a)/sqrt(length(a)), trialAreas);

% AUCs = traceAvg;
% AUCerr = traceErr;

figure('Units','normalized' ,'Position',[0 0 1 1])
hbar = bar(AUCs,'FaceColor','flat');
hbar(1).CData = [1 .5 0];
ctr = [];
for m = 1:length(hbar)
    ctr(m,:) = bsxfun(@plus, hbar(1).XData, hbar(m).XOffset);
end
hold on
% hbar(1).FaceColor = [1 0 0];
% hbar(2).FaceColor = [.75 0 0];
% hbar(3).FaceColor = [.5 0 0];
% hbar(4).FaceColor = [0 1 0];
% hbar(5).FaceColor = [0 .5 0];
% hbar(6).FaceColor = [0 0 1];
% hbar(7).FaceColor = [0 0 .5];
legend(["Total"],'AutoUpdate','off');

errorbar(ctr, AUCs', AUCerr', '.k')
xticklabels(conds);
ylabel('Area Under Curve');
for cm = 1:length(comps)
%      aDiff = abs(diff(AUCs(comps(cm,:),1)))/sqrt(sum(AUCerr(comps(cm,:),1).^2)/2);
%      hDiff = abs(diff(AUCs(comps(cm,:),2)))/sqrt(sum(AUCerr(comps(cm,:),2).^2)/2);
    tDiff = abs(diff(sum(AUCs(comps(cm,:),:),2)))/sqrt(sum(AUCerr(comps(cm,:),:).^2,'all')/2);
%     xArm = ctr(1,comps(cm,:));
%     xHand = ctr(2,comps(cm,:));
    xTotal = comps(cm,:);
    yPos = [20,20];
    sT = cm*1.5;
%     line(xArm,yPos+sT,'Color', [1 0 0], 'LineWidth', 2, 'LineStyle','--');
%     line(xHand,yPos+sT,'Color', [0 .5 .5], 'LineWidth', 2,'LineStyle','--');
     line(xTotal,yPos+sT+1, 'Color', 'k', 'LineWidth',2);
     text(sum(xTotal)/2, yPos(1)+.51+sT+1,['d: ',num2str(tDiff,2)],'HorizontalAlignment','center','FontWeight','bold','FontSize',15);
%     text(min(xHand),yPos(1)+.51+sT, ['d: ',num2str(aDiff,2)],'HorizontalAlignment','left','FontWeight', 'bold', 'FontSize', 15);
%     text(max(xArm),yPos(1)+.51+sT, ['d: ',num2str(hDiff,2)],'HorizontalAlignment','right','FontWeight', 'bold', 'FontSize', 15);
end
% areas1 = areas;
% errors1 = errors;
% % clear errors areas;
% areas(1,:) = areas1(1,:)-areas1(3,:);
% areas(2,:) = areas1(2,:)-areas1(3,:);
% areas(3,:) = areas1(1,:)-areas1(2,:);
% errors(1,:) = errors1(1,:)-errors1(3,:);
% errors(2,:) = errors1(2,:)-errors1(3,:);
% errors(3,:) = errors1(1,:)-errors1(2,:);
saveFigures(gcf,['S:\Lab\ngc14\Working\', monkey,'\'],'EMG_AUC', []);

%
% allComps = [];
% for d = 1:length(modelNames)
%     if(anova1(cell2mat(trialAreas(:,d)'),[repmat([1],1,length(trialAreas{1,d})), ...
%             repmat([2],1,length(trialAreas{2,d})), repmat([3],1,length(trialAreas{3,d}))],'off')<pVal)
%         for c = 1:length(comps)
%             if(anova1(cell2mat(trialAreas(comps(c,:),d)'), [repmat([1],1,length(trialAreas{comps(c,1),d})),...
%                     repmat([2],1,length(trialAreas{comps(c,2),d}))],'off')<pVal/length(comps))
%                 allComps(c,d) = anova1(cell2mat(trialAreas(comps(c,:),d)'), [repmat([1],1,length(trialAreas{comps(c,1),d})),...
%                     repmat([2],1,length(trialAreas{comps(c,2),d}))], 'off');
%                 if(nanmean(trialAreas{comps(c,2),d}) > nanmean(trialAreas{comps(c,1),d}))
%                     allComps(c,d) = Inf;
%                 end
%             else
%                 allComps(c,d)= Inf;
%             end
%         end
%     end
% end
% 
% CIs = [];
% figure();
% for d = 1:length(modelNames)
%     for c = 1:length(comps)
%         areasComps = [];
%         labels = [];
%         if(nanmean(trialAreas{comps(c,2),d}) > nanmean(trialAreas{comps(c,1),d}))
%             areasComps = cell2mat(trialAreas(comps(c,:),d));
%             labels =  [repmat('2',length(trialAreas{comps(c,1),d}),1);repmat('1',length(trialAreas{comps(c,2),d}),1)];
%             tabRes = dabest(table(labels,areasComps));
%             CIs(c,d) = tabRes.Value(3);
%         else
%             CIs(c,d) = Inf;
%         end
%     end
% end
% maxTrials = max(max(cellfun(@length,trialAreas)));
% %
% customCMap = [1 0 0;.7 .4 0;1 .64 0;0 1 0;0 .5 0;0 0 1;0 0 .5];
% nanPadAA=cellfun(@(a) nanPad(a,maxTrials), trialAreas,'UniformOutput', false);
% nanPadAA= cell2mat(nanPadAA);
% bGroups = zeros(size(nanPadAA));
% bp = gca(figure('Units', 'normalized', 'Position', [0 0 1 1]));
% 
% xPosB = .5:1/(d-1):1.5;
% hold on;
% for c = 0:length(conds)-1
%     cF = figure('Units', 'normalized', 'Position',  [0 0 1 1]);
%     avgSegs = nanmean(times{c+1});
%     alignTime = cellfun(@(s) avgSegs(strcmp(s,alignPoint)), values(events,conds(c+1)));
%     hold on;
%     title(conds{c+1});
%     ylabel('microVolts');
%     xlabel('Time from reach onset (s)')
%     cellfun(@(s,n) plot(s,'Color',n, 'LineWidth',2), traceAvg(c+1,:)', num2cell(customCMap,2));
%     lg = legend(modelNames);
%     lg.AutoUpdate = 'Off';
%     cellfun(@(s,e,n) plot(s+e,'Color', n, 'LineWidth', 1, 'LineStyle', ':'), traceAvg(c+1,:)',traceErr(c+1,:)', num2cell(customCMap,2));
%     cellfun(@(s,e,n) plot(s-e,'Color', n, 'LineWidth', 1,'LineStyle', ':'), traceAvg(c+1,:)', traceErr(c+1,:)',num2cell(customCMap,2));
%     arrayfun(@(t) line([t,t],[0 250],'Color','r','LineStyle','--'), avgSegs);
%     set(gca, 'XTick', [xPlotLims(1):Fs:xPlotLims(end)] + alignTime);
%     set(gca, 'XTickLabel', arrayfun(@num2str, [xPlotLims(1):Fs:xPlotLims(end)]./Fs, 'UniformOutput', false));
%     set(gca, 'XLim', xPlotLims+alignTime);
%     set(gca, 'YLim', [0 100]);
% 
% 
%     start = 1+(c*length(modelNames)) + (c*2);
%     boxplot(bp,nanPadAA((c*maxTrials)+1:((c*maxTrials)+1)+maxTrials-1,:),...
%         'BoxStyle', 'filled','Colors', customCMap, 'Symbol', '', 'Width', 2,...
%         'Positions', start:start+d-1,'Notch', 'marker');
%     if(c==0)
%         lg2 = legend(printNames);
%         lg2.AutoUpdate = 'Off';
%     end
% end
% xlim(bp,[0 26]);
% xticks(bp,[4 13 22]);
% xticklabels(bp,[{'Precision'}, {'Power'}, {'Reach-only'}]);
% ylabel('SI (r - g)');
% ylim(bp,[-1 1]);
% line(bp,[0 26], [0 0], 'Color', 'k')
% gs = gscatter(NaN(7,1),NaN(7,1),1:7',customCMap);
% legend(gs,printNames,'Location', 'best');
% 
% saveFigures(gcf,['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\SI\'],'MuscleSI', []);
% figure('Units', 'normalized', 'Position', [0 0 1 1]);
% hbarG = bar(areasG);
% for m = 1:size(areasG',1)
%     ctrG(m,:) = bsxfun(@plus, hbarG(1).XData, hbarG(m).XOffset');
% end
% hold on
% hbarG(1).FaceColor = [1 0 0];
% %hbarG(2).FaceColor = [.7 .4 0];
% %hbarG(3).FaceColor = [1 .64 0];
% %hbarG(4).FaceColor = [0 1 0];
% %hbarG(5).FaceColor = [0 .5 0];
% hbarG(2).FaceColor = [0 0 1];
% %hbarG(7).FaceColor = [0 0 .5];
% errorbar(ctrG', areasG, errorsG, '.k')
% legend([{'Proximal'}, {'Distal'}],'Location','best');
% xticklabels([{'Precision'}, {'Power'}, {'Reach-only'}]);
% ylabel('SI (r-g)');
% hold off;
% saveFigures(gcf,['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\SI\'],'GroupedSI', []);

function saveFigures(fHandle,saveDir,saveName,imgTrace)
if(~exist(saveDir, 'dir'))
    mkdir(saveDir)
end
if(strcmp(get(fHandle,'type'),'figure'))
    set(fHandle,'Renderer','painters');
end
saveas(fHandle,[saveDir,saveName,'.epsc']);
saveas(fHandle,[saveDir,saveName,'.png']);
saveas(fHandle,[saveDir,saveName,'.fig']);
if(~isempty(imgTrace))
    boundaries = bwboundaries(imgTrace);
    fHandle2 = copyobj(gca(fHandle),figure()); hold on;
    cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth), boundaries);
    saveFigures(fHandle2,[saveDir,'Imaging\'],saveName,[]);
end
end

function padded = nanPad(x,m)
padded = x;
padded(end+1:end+(m-length(x))) = nan(m-length(x),1);
end