clear all;
close all;
monkey = 'Gilligan';
conds = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};

events = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}});

modelNames = {'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor',...
    'Digit Extensor','Digit Flexor'};
printNames = { 'Deltoid', 'Biceps', 'Triceps','WristEx', 'WristFlex', 'DigitEx', 'DigitFlex'};
directory = dir(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\']);

Fs = 2000;
phaseWindow = .2;
phaseWindow = phaseWindow*Fs;
EMGfiles = find(cellfun(@(a) endsWith(a, '.mat'), {directory.name}));
loadVars = {'dates','rawTrials','segs'};
phaseWindows = {[0, phaseWindow],[-phaseWindow/2,phaseWindow/2],[-phaseWindow, 0]};
condPhaseAlign = containers.Map(conds,{["GoSignal","StartReach","StartLift"],...
    ["GoSignal","StartReach","StartLift"],["GoSignal","StartReach","StartHold"]});

condTraces = {};
condSegs = {};
for c = 1:length(conds)
    currCond = {};
    currSegs = {};
    muscDates = {};
    mDates = {};
    mTrials = {};
    mSegs = {};
    phases = condPhaseAlign(conds{c});
    parfor d = 1:length(modelNames)
        di = find(startsWith({directory(EMGfiles).name}, modelNames{d}));
        s = load([directory(EMGfiles(di)).folder, '\', directory(EMGfiles(di)).name],loadVars{:});
        mDates(d) = s.dates(c,1);
        mTrials(d) = s.rawTrials(c,1);
        mSegs(d) = s.segs(c,1);
    end
    clear s;
    %%
    allTrace = {}; 
    allSegs = {};
    allDates = unique(string([mDates{:}]));
    validDates = allDates(cellfun(@(d) all(cellfun(@(m) any(find(datenum(string(d))==datenum(string(m)))), mDates)), allDates));
    for a = 1:length(validDates)
        dateInds = cellfun(@(d,t) datenum(string(d))==datenum(validDates(a)) ...
            & cellfun(@(tt) ~all(isnan(tt)),t)', mDates,mTrials, 'UniformOutput',false);
        currSegs = cellfun(@(s,d) cell2mat(s(d)'), mSegs,dateInds,'UniformOutput', false);
        dateInds = cellfun(@(d) find(d), dateInds, 'UniformOutput',false);
        trialIDs = unique(cell2mat(currSegs'),'row');
        trialSessionSegs = num2cell(trialIDs,2);
        matchSegs = cellfun(@(cc) cellfun(@(s) find(all(cc==s,2)), ...
            trialSessionSegs,'UniformOutput',false), currSegs,'UniformOutput',false);
        matchSegs = [matchSegs{:}];
        allTrialInds = all(~cellfun(@isempty,matchSegs),2);
        matchSegs = num2cell(cell2mat(matchSegs(allTrialInds,:)),1);
        if(~isempty(matchSegs))
            allMuscTrials = cellfun(@(t,d,ai) cell2mat(t(d(ai))'), mTrials,dateInds,matchSegs,'UniformOutput',false);
            allTrace{a} = nanmean(cat(3,allMuscTrials{:}),3);
            allSegs{a} = cell2mat(trialSessionSegs(allTrialInds));
            allMuscTrials = [];
        end
        disp(a)
    end
    condTraces{c}= cell2mat(cellfun(@(p) nanPad(p,max(cellfun(@length,allTrace))),allTrace,'UniformOutput', false)');
    condSegs{c} = cell2mat(allSegs');
    %%
    AUC = {};
    for p = 1:length(phases)
        alignPoint = strcmp(events(conds{c}),phases{p});
        alignWindow = phaseWindows{p};
        tTimes = condSegs{c}(:,alignPoint)+alignWindow;
        AUC{p} = cellfun(@(t,r) trapz(1/Fs,t(r(1):r(2))), num2cell(condTraces{c},2),num2cell(tTimes,2), 'UniformOutput', true);
    end
    hx(c) = subplot(2,2,c);
    AUCs = cellfun(@nanmean, AUC);
    AUCerr = cellfun(@(a) nanstd(a), AUC);
    hbar = bar(AUCs','FaceColor','flat');
    hold on
    errorbar(hbar.XData', AUCs', AUCerr', '.k');
    xticklabels(["Go","Reach","Grasp"]);
    title([conds{c}," (", num2str(length(AUC{1}))," trials)"]);
    xlabel('Avg Muscle Activity');
    ylabel('Area Under Curve');
end
linkaxes(hx);
saveFigures(gcf,['S:\Lab\ngc14\Figures\Physiology\Results\EMG\AUC\'],[monkey,'_Phases'], []);

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
padded(:,end+1:end+(m-length(x))) = nan(size(padded,1),m-length(x));
end