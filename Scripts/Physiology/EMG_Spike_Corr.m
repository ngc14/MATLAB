function EMG_Spike_Corr(monkey,sessionDate)%% EMG XCorr analysis
if(~exist('monkey', 'var'))
    monkey ='Gilligan';
end
if(~exist('sessionDate', 'var'))
    sessionDate ='11_19_2019';
end
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'mm_dd_yyyy';
else
    dateFormat = 'yyyy_mm_dd';
end
params = PhysRecording(["Extra Small Sphere", "Large Sphere", "Photocell", "Rest"]);
musclesColors = containers.Map({'Deltoid','Biceps', 'Triceps', 'WristExtensor', 'WristFlexor',...
    'DigitExtensor','DigitFlexor'}, {[1 0 0],[.7 .4 0],[1 .64 0],[0 1 0],[0 .5 0],[0 0 1],[0 0 .5]});
musclesColors = containers.Map({'Arm','Hand'}, {[1 0 0], [0 0 1]});
trialCorr = true;
singleOrAll = 'Single';
sessionDateF = datestr(sessionDate,dateFormat);
startAlign = {'StartReach'};
phaseAlignments = {{'GoSignal', 'StartReach', 'StartLift', 'StartWithdraw'},...
    {'GoSignal','StartReach', 'StartLift', 'StartWithdraw'},...
    {'GoSignal','StartReach', 'StartHold', 'StartWithdraw'},...
    {'GoSignal','StartReplaceHold','StartReplaceSuccess', 'StartReward'}};
phaseAlignments = {{'GoSignal' , 'StartLift'}, {'GoSignal' , 'StartLift'},...
    {'GoSignal' , 'StartHold'},{'GoSignal' , 'StartReplaceHold'}};
phaseWindows = {[0,1],[-.55, .45], [-.90, .10],[-.10,.90]};
phaseWindows = {[-0.35 0]};
evokedResponseCategorization = containers.Map(...
    {'Deltoid','Biceps', 'Triceps', 'WristExtensor', 'WristFlexor','DigitExtensor','DigitFlexor'},...
    {'Arm', 'Arm','Arm', 'Hand','Hand','Hand', 'Hand'});
EMGdir = dir(['S:\Lab\', monkey, '\All Data\', monkey,'_', sessionDateF,'\EMG\Results_New\*.mat']);
saveDir = 'S:\Lab\ngc14\Meeting\EMG_UNITS\';
if(trialCorr)
   saveDir = [saveDir,'Trials\']; 
else
    saveDir = [saveDir,'Session\'];
end
%%
[spikes,times,~,condTrials,sConds,~,~,~] =...
    getSessionInfo2(['S:\Lab\',monkey,'\All Data\', monkey,'_',...
    sessionDateF],singleOrAll);

zeroBin = find(params.bins(1:end-1)==0);
numUnits = size(spikes,1);
trialCondInds =  cellfun(@(cn) strcmp(condTrials(:,1),cn)', sConds,'UniformOutput', false); 

%%
modelNames = unique(evokedResponseCategorization.values(evokedResponseCategorization.keys()));
muscles = evokedResponseCategorization.keys();
alignedMuscles = {};
allSegs = {};
for m = 1:length(muscles)
    loadFile = matfile([EMGdir(m).folder,'\',EMGdir(m).name]);
        allTrials = loadFile.rawTrials';
    if(m==1)
        Fs = loadFile.Fs;
        phaseSegInd = cellfun(@(sg,ph) cellfun(@(p) find(strcmp(sg,p)),ph,'UniformOutput',true),...
            loadFile.ConditionSegs', phaseAlignments,'UniformOutput', false);
        alignInd = cellfun(@(sg)  find(strcmp(sg,startAlign)),...
            loadFile.ConditionSegs','UniformOutput', false);
                alignInd{end} = 1;

    end
    segTimes = loadFile.segs';
    nanTrials = cellfun(@(tc) cellfun(@(t) all(isnan(t)),tc,'UniformOutput',false) ,segTimes,'UniformOutput',false);
    segTimes = cellfun(@(sc,nc) cellfun(@(s,n) nansum(vertcat(s,repmat(n,1,length(s)).*Fs)) ,sc,nc,'UniformOutput',false), ...
        segTimes,nanTrials,'UniformOutput',false);
    alignedEMG =  cellfun(@(tc,ac,sc) cellfun(@(t,a) cellfun(@(pw) t(max(1,a(sc(1))+(Fs*pw(1))):...
        min(length(t),(a(sc(end))+(Fs*pw(end))))),phaseWindows,'UniformOutput', false), ...
        tc,ac,'UniformOutput', false),allTrials,segTimes,phaseSegInd, 'UniformOutput', false);
    alignedEMG = cellfun(@(tc) cellfun(@(t) cellfun(@(a,pw) [a,NaN(1,...
        (range(pw)*Fs)-(length(a)-1))],t,phaseWindows,'UniformOutput',false),...
        tc, 'UniformOutput',false),alignedEMG, 'UniformOutput', false);
    alignedEMG = cellfun(@(ac) vertcat(ac{:}), alignedEMG, 'UniformOutput', false);
    NaNPaddedEMG = cellfun(@(ac) cellfun(@(av) cell2mat(cellfun(@(a) nanPadArrays(a,...
        max(cellfun(@length,av))), av, 'UniformOutput', false)),...
        num2cell(ac,1),'UniformOutput',false), alignedEMG,'UniformOutput',false);
    
    alignedMuscles{m} = NaNPaddedEMG;
    allSegs{m} = segTimes;
end
groupedMuscles = {};
alignTimes = cellfun(@(ci,ac) cellfun(@(t) t(ac),...
    times(ci), 'UniformOutput', false), trialCondInds,alignInd,'UniformOutput', false);
for c = 1:length(params.condNames)
   [ trialIDs utInd,~]= unique(cell2mat(cellfun(@(a) vertcat(a{c}{:}), allSegs, 'UniformOutput', false)'),'row');
    trialSessionSegs = num2cell(trialIDs,2);
    matchSegs = cellfun(@(cc) cellfun(@(s) find(all(cell2mat(cc{c}')==s,2),1), ...
        trialSessionSegs,'UniformOutput',false), allSegs,'UniformOutput',false);
    matchSegs = [matchSegs{:}];
    allTrialInds = all(~cellfun(@isempty,matchSegs),2);
     uTrialMatch{c} = utInd(allTrialInds);
    matchSegs = num2cell(cell2mat(matchSegs(allTrialInds,:)),1);
    if(~isempty(matchSegs))
        allMuscTrials{c} = cellfun(@(t,m) cellfun(@(ta) ta(m,:),t{c},'UniformOutput',false), alignedMuscles,matchSegs,'UniformOutput',false);
    end
    for m = 1:length(modelNames)
        mGroupInd = find(cellfun(@(s) strcmp(modelNames{m}, evokedResponseCategorization(s)), muscles));
        
        groupedMuscles{m}(c) = cellfun(@(a) num2cell(nanmean(cat(3,a{:}),3),2),...
            squeeze(num2cell(permute(cat(3,allMuscTrials{c}{mGroupInd}),...
            [1,3,2]),[1,2])),'UniformOutput', false)';
        if(~trialCorr)
            groupedMuscles{m}(c) = cellfun(@(g) nanmean(g,1), groupedMuscles{m,c},'UniformOutput',false);
        end
    end
end
%%
allUnitWinds = {};
uCorrs ={};
uLag = {};
phaseAligned = cellfun(@(t,a,ps) cellfun(@(tt,aa) cellfun(@(pw) (tt(ps)+pw)-aa,...
    phaseWindows,'Uniformoutput',false), times(t),a, 'UniformOutput',false), trialCondInds,alignTimes,phaseSegInd,'UniformOutput',false);
for u = 1:numUnits
    %%
    alignSpikes = cellfun(@(at,ci) cellfun(@(s, a)...
        histcounts(s-a,params.bins)./params.binSize,...
        spikes(u,ci),at,'UniformOutput', false),alignTimes,trialCondInds, 'UniformOutput', false);
    alignedHists = cellfun(@(hc) vertcat(hc{:}), alignSpikes, 'UniformOutput', false);
    smoothAligned = cellfun(@(hc,u) cellfun(@(a) conv(a,gausswin(params.sigma)/...
        sum(gausswin(params.sigma)),'same'), num2cell(hc,2),...
        'UniformOutput', false),alignedHists, uTrialMatch,'UniformOutput', false);
    segmentedPSTH = cellfun(@(s,p,u) cellfun(@(pw,ss) cell2mat(cellfun(@(ps)  ss(:,find(isalmost(...
        params.bins,ps(1),params.binSize/1.99)):find(isalmost(...
        params.bins,ps(end),params.binSize/1.99))),pw,'UniformOutput',false)),...
        p(u(u<=size(p,2)))',s(u(u<=size(p,2))), 'UniformOutput',false),...
        smoothAligned,phaseAligned,uTrialMatch, 'UniformOutput',false);
     NaNPaddedPSTH = cellfun(@(ac) cellfun(@(av) nanPadArrays(av,...
        max(cellfun(@length,ac))), ac, 'UniformOutput', false),...
        smoothAligned,'UniformOutput',false);
%     if(trialCorr)
%         alignedPSTH = cellfun(@(rc)  cellfun(@(r) reshape(cell2mat(r),size(r,1),[]),...
%             num2cell(rc,1),'UniformOutput',false),smoothAligned, 'UniformOutput',false);
%     else
%          alignedPSTH = cellfun(@(rc) num2cell(squeeze(nanmean(cell2mat(...
%             cellfun(@(r) reshape(r,1,1,[]),rc,'UniformOutput', false)),1)),2)',...
%             smoothAligned,'UniformOutput',false);
%     end
    allUnitWinds(u,:) = cellfun(@(ac,ic) ac(ic(ic<=size(ac,1))), NaNPaddedPSTH, uTrialMatch, 'UniformOutput',false);
    
    mCorr = {};
    mLag = {};
    for g = 1:length(modelNames)
        %interp
        [corrValues,tLags] = cellfun(@(mc) cellfun(@(mp,cp) cellfun(@(m,c) xcorr(c,...
            interpft(g,length(c)),'unbiased'),num2cell(cell2mat(mp),2),...
            num2cell(cell2mat(cp),2),'UniformOutput', false),mc,allUnitWinds(u,:),'UniformOutput',false),...
            groupedMuscles(g),'UniformOutput', false);%cellfun(@(m,u)xcorr(u,interpft(m,length(u)),'biased'), mc,uc, 'UniformOutput', false),...
        [muCorrs,mulag] = cellfun(@(cs) cellfun(@(cp) cellfun(@(c) ...
            max(c(true(1,round(length(c)/2)))),cp,'UniformOutput', true),...
            cs,'UniformOutput', false),corrValues,'UniformOutput', false);
        [mCorr(g,:), maxInd] = cellfun(@(cc) nanmax(cc,[],2), [muCorrs{:}], 'UniformOutput', false);
        mulag = [mulag{:}];
        maxMuLag = cellfun(@(mx,mu) mu(sub2ind(size(mu),1:size(mu,1),mx'))', maxInd, mulag, 'UniformOutput',false);
        mLag(g,:) = cellfun(@(t,m) t{1}(m),[tLags{:}],maxMuLag,'UniformOutput',false);
    end
    uCorrs{u} = cellfun(@nanmedian,mCorr,'UniformOutput',true);
    uLag{u} = cellfun(@nanmedian,mLag,'UniformOutput', true);
end
groupedMuscles = cellfun(@(g) cellfun(@(ga,b) cellfun(@(a) [interpft(a,length(cell2mat(b))), ...
    NaN(size(a,1),50)],ga,'UniformOutput',false),g,allUnitWinds(1,:),...
    'UniformOutput', false),groupedMuscles,'UniformOutput', false);
%%
[maxCorr,maxInd] = cellfun(@(u) max(u,[],1),uCorrs,'UniformOutput',false);
xE = [0,cellfun(@(r) range(r)/params.binSize, phaseWindows)];
xVals{1} = xE(1)+1:xE(2);
for x = 2:length(xE)-1
    xVals{x} = (50+xVals{x-1}(end)+1):(xVals{x-1}(end)+xE(x+1)+50);
end
%%
for u = 1:numUnits
    figure('Units', 'normalized', 'Position', [0 0 1 1]); hold on;
    ax= {};
    for c = 1:length(sConds)
        ax{c} = subplot(2,2,c); hold on;
        title({string(sConds{c}); strcat("Best correlation: ", string(modelNames{maxInd{u}(c)})," ",...
            num2str(maxCorr{u}(c),'% .02f')); strcat("Lag: ", num2str(uLag{u}(maxInd{u}(c),c)*10), " ms")});
        psthNaNGap = cellfun(@(a) [a, NaN(size(a,1),50)],[allUnitWinds{u,c}(:)],'UniformOutput',false);
        unitCondPSTH = nanmean(cell2mat(psthNaNGap)',2)';

        yyaxis right
        plot(unitCondPSTH,'Color','k','LineWidth',2,'LineStyle','-','Marker','none');
        
        yyaxis left
        cellfun(@(m,mc) plot(nanmean(cell2mat(m{c}),1),'Color',mc,'LineStyle','-',...
            'LineWidth',2,'Marker','none'),groupedMuscles,values(musclesColors,modelNames),'UniformOutput', false);
        
        g =cellfun(@(cl) plot(NaN,NaN,'Color',cl,'LineWidth', 2,...
            'LineStyle', '-','Marker','none'),values(musclesColors,modelNames));
        legend(g,cellfun(@(m,c) [m,' r= ', num2str(c,'%.2f')],modelNames,...
            num2cell(uCorrs{u}(:,c))','UniformOutput', false),'AutoUpdate','off');
        xticks(cell2mat(cellfun(@(p) [p(1), p(floor(length(p)/2)), p(end)], xVals,'UniformOutput', false)));
        xticklabels([cell2mat(cellfun(@(p) [p(1), NaN, p(end)], phaseWindows,'UniformOutput',false))]);
    end
    if(~exist([saveDir,sessionDateF,'\'],'dir'))
        mkdir([saveDir,sessionDateF,'\']);
    end
    linkprop(cellfun(@(aa) aa.YAxis(1), ax),'Limits');
    linkprop(cellfun(@(aa) aa.YAxis(2), ax),'Limits');
    saveas(gcf,[saveDir,sessionDateF,'\',num2str(u),'.png']);
    close all;
end
end
%%

function paddedArr = nanPadArrays(arr,m)
paddedArr = arr;
paddedArr(end+1:m) = NaN;
end

function bin = findBins(allBins, timePoint)
binSize = mode(diff(allBins));
if(isempty(timePoint))
    bin = NaN;
else
    bin = find(isalmost(allBins,timePoint,binSize/1.99),1);
end
if(isempty(bin))
    bin = NaN;
end
end
