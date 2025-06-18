%% load data
loadPath = "S:\Lab\ngc14\Working\Both\Baseline_FR\";%"E:\Data\"
phaseNames = ["Go", "Reach", "Grasp", "Withdraw"];
fTypes = ["Reach","Grasp","Both","Task"];
condLabels = {'E','L','P'};
countsOrFR = "counts";
cMatFigs = false;
numRuns = 500;
m = matfile(loadPath+"phaseAnalysis_Face.mat");
mtform = load('S:\Lab\ngc14\Working\S2G_Reference_tform.mat');
avgPhase = m.avgPhase;
siteSegs = m.siteSegs;
siteChannels = m.siteChannels;
siteChannels = siteChannels{2};
phaseFR = m.phaseFR;
phaseBaseline = m.phaseBaseline;
taskBaseline = m.taskBaseline;
sdm = m.siteDateMap;
taskFR = m.taskFR;
siteTrialPSTHS = m.normPSTH;
siteImaging = m.siteActiveInd;
chUnitMap = cell(height(sdm),1);
chUnitMap(strcmp([sdm.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([sdm.Monkey],"Skipper")) = {[32:-1:1]};
mappedChannels = cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap,siteChannels, 'Uniformoutput', false);
mappedChannels = [mappedChannels{:}];
params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5);
add_ndt_paths_and_init_rand_generator
phaseWinSz = 0.2;
phaseWindows = {[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz, 0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]};
phaseAlignmentPoints ={["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"]};
phaseInds = cellfun(@(p,cc) arrayfun(@(b)find(strcmp(string(cc),b)),...
    p),phaseAlignmentPoints,values(params.condSegMap,...
    params.condSegMap.keys()),'UniformOutput',false);
psthPhaseEnds = cellfun(@(a,pa) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,pw) ...
    findBins(ap(:,p)+pw,params.bins),num2cell(pa),phaseWindows,'UniformOutput', false),ps,'UniformOutput',false),...
    a, 'UniformOutput', false),siteSegs,phaseInds,'UniformOutput', false);
psthPhaseEnds = cellfun(@(n) vertcat(n{:}), psthPhaseEnds, 'UniformOutput',false);
maxTrialsSize = num2cell(max(cell2mat(cellfun(@(c) cellfun(@(s) size(s{1},3),c), siteTrialPSTHS(1:length(condLabels)), 'UniformOutput',false)),[],2));
%% trial by trial condition labels organized by unit
FRWind = cellfun(@(sc,pi) cellfun(@(s,p,m) cellfun(@(tp,pe) cellfun(@(pp) [cellfun(@(t,i) t(:, max(1,i(1)):min(size(t,2),i(end))), ...
    squeeze(num2cell(tp,[1 2])),num2cell(pp,2), 'UniformOutput',false);repmat({NaN(size(tp,1),1)},m-size(pp,1),1)], pe, 'UniformOutput',false),s,p,'UniformOutput',false),...
    sc,num2cell(pi,2),maxTrialsSize, 'UniformOutput',false), siteTrialPSTHS(1:length(condLabels)),psthPhaseEnds, 'UniformOutput',false);
[maxVal,maxInd] = cellfun(@(c) cellfun(@(a) cellfun(@(p) cellfun(@(s) cellfun(@(st) max(st,[],2), s, 'UniformOutput',false),...
    p, 'UniformOutput',false), a, 'UniformOutput',false), c, 'UniformOutput',false), FRWind, 'UniformOutput',false);
maxVal = cellfun(@(c) cellfun(@(a) cellfun(@(p) cellfun(@(t)cell2mat(t'),p,'UniformOutput',false),a, 'UniformOutput',false),c,'UniformOutput',false), maxVal, 'UniformOutput',false);
maxVal = cellfun(@(c) cellfun(@(s) cellfun(@(m) median(m,3,'omitnan'), cellfun(@cell2mat,num2cell(cat(3,s{:}),3),'UniformOutput',false), 'UniformOutput',false),...
    c, 'UniformOutput',false), maxVal, 'UniformOutput',false);
maxVal = cellfun(@(v) cellfun(@(f) num2cell([f{:}],2), num2cell(vertcat(v{:}),1), 'UniformOutput',false), num2cell([maxVal{:}],2), 'UniformOutput',false);
maxVal = cellfun(@(m) vertcat(m{:}), num2cell(vertcat(maxVal{:}),1), 'UniformOutput',false);
maxVal = [maxVal{:}];
maxInd = cellfun(@(c,pc) cellfun(@(a,e) cellfun(@(p,i) cell2mat(p'),...%+(zeros(size(p{1},1),size(p,1))./double(cell2mat(p')~=1))+[i(:,1)',NaN(1,size(p,1)-size(i,1))], ...
    a,e,'UniformOutput',false),vertcat(c{:}),pc,'UniformOutput',false),maxInd, psthPhaseEnds,'UniformOutput',false);
maxInd = cellfun(@(c) cellfun(@(s) cellfun(@(m) fix(median(cat(3,m{:}),3,'omitnan')), ...
    cellfun(@squeeze,num2cell(cat(3,s{:}),3),'UniformOutput',false), 'UniformOutput',false),...
    num2cell(c,2), 'UniformOutput',false), maxInd, 'UniformOutput',false);
maxInd = cellfun(@(v) cellfun(@(f) num2cell([f{:}],2), num2cell(vertcat(v{:}),1), 'UniformOutput',false), num2cell([maxInd{:}],2), 'UniformOutput',false);
maxInd = cellfun(@(m) vertcat(m{:}), num2cell(vertcat(maxInd{:}),1), 'UniformOutput',false);
maxInd = [maxInd{:}];
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), m.siteRep'));
monkeyUnitInd = cellstr(mapSites2Units(cellfun(@length, siteChannels), sdm.Monkey));
siteImaging = cellfun(@(si) mapSites2Units(cellfun(@length, siteChannels), si), siteImaging, 'UniformOutput',false);
xyLoc = num2cell(mapSites2Units(cellfun(@length,siteChannels),num2cell([sdm.x,sdm.y],2)),2);
xyLoc(strcmp(monkeyUnitInd,"Skipper")) = cellfun(@(s) round(transformPointsForward(mtform.tform,s)), xyLoc(strcmp(monkeyUnitInd,"Skipper")),'UniformOutput',false);
trialCondTable = getTrialPhaseTable(phaseFR(1:length(condLabels)),phaseNames,condLabels,sdm);
for a = length(phaseNames):-1:1
    trialCondTable = addvars(trialCondTable, cell2mat(maxVal(:,a)')', 'NewVariableNames',char("MaxFR_"+phaseNames(a)), 'Before',1);
    trialCondTable  = addvars(trialCondTable, cell2mat(maxInd(:,a)')', 'NewVariableNames',char("MaxInd_"+phaseNames(a)), 'Before',1);
end
trialConds = trialCondTable(:,strcmp(trialCondTable.Properties.VariableNames,'Condition'));
trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t,:)),unique(trialCondTable.Unit),'UniformOutput',false);
clear siteTrialPSTHS FRWind
%% only use units that are task modulated according to the task phase window
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(condLabels)),taskFR(1:length(condLabels)),'UniformOutput',false);
taskUnits = cell2mat(cellfun(@cell2mat, taskUnits,'UniformOutput',false));
tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(a) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
    a,'UniformOutput',false), s,'UniformOutput', false),phaseFR(1:length(condLabels)),'UniformOutput', false), phaseNames, 'UniformOutput',false);
[~,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials(r,g,1,true,0.05),cr,cg,'UniformOutput',false),...
    tPhase{strcmp(phaseNames,"Reach")},tPhase{strcmp(phaseNames,"Grasp")}, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
%%
AUCVals = cellfun(@(c) cell2mat(cellfun(@(r) median(cat(3,r{:}),3), c, 'UniformOutput', false)), avgPhase(1:length(condLabels)),'UniformOutput',false);
rAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Reach")), AUCVals, 'UniformOutput',false);
gAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Grasp")), AUCVals, 'UniformOutput',false);
rUnits = cellfun(@(rg,r,g) rg==1 & r > g, rgInds, rAUC, gAUC, 'UniformOutput',false);
gUnits = cellfun(@(rg,r,g) rg==1 & r < g, rgInds, rAUC, gAUC, 'UniformOutput',false);
bothUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, num2cell(taskUnits,1), 'UniformOutput',false);
fUnits = {cell2mat(rUnits), cell2mat(gUnits), cell2mat(bothUnits)};
fUnits{end+1} = taskUnits;
lUnits = {repmat(mappedChannels'<16,1,3)};
lUnits{end+1} = repmat(mappedChannels'>=16,1,3);
tUnitsInd = sum(taskUnits,2)>0;

trialConds = trialConds(tUnitsInd);
unitSomatotopy = unitSomatotopy(tUnitsInd);
monkeyUnitInd = string(monkeyUnitInd(tUnitsInd));
xyLoc  = xyLoc(tUnitsInd);
maxVal = maxVal(tUnitsInd,:);
maxInd = maxInd(tUnitsInd,:);
fUnits = cellfun(@(f) f(tUnitsInd,:), fUnits, 'UniformOutput',false);
%fUnits = cellfun(@(f) cellfun(@(l) f & l, lUnits, 'UniformOutput', false), fUnits, 'UniformOutput',false);
%% get trial FR or trial spike counts
if(strcmpi(countsOrFR,'FR'))
    normFR = trialCondTable(:, contains(trialCondTable.Properties.VariableNames,phaseNames));
    normFR = arrayfun(@(t) table2array(normFR(trialCondTable.Unit==t,:)),unique(trialCondTable.Unit),'UniformOutput',false);
    trialTaskFR = cellfun(@(c) cellfun(@(s) median(cell2mat(cat(3,s{:})),3,'omitnan'), c,'UniformOutput',false), taskFR(1:3), 'UniformOutput',false);
    trialTaskFR = cellfun(@(s)num2cell(horzcat(s{:}),2)',num2cell([trialTaskFR{:}],2),'UniformOutput',false);
    trialTaskFR = cellfun(@transpose,[trialTaskFR{:}]','UniformOutput',false);
    trialTaskFR = trialTaskFR(tUnitsInd);
    allUnitsTrials = normFR;
elseif(strcmpi(countsOrFR,'counts'))
    % counts
    clear trialCondTable
    rawSpikes = m.rawSpikes;
    baselineWindow = [-5,-1];
    baselineAlignPhases = ["GoSignal", "GoSignal"];
    baselineInds = cellfun(@(cc) arrayfun(@(b) find(strcmp(string(cc),b)),...
        baselineAlignPhases,'UniformOutput', true),values(params.condSegMap,...
        params.condSegMap.keys()),'UniformOutput',false);
    baselineEnds = cellfun(@(a,pa) cellfun(@(ps) ...
        findBins(ps{1}(:,pa)+baselineWindow,params.bins), a, 'UniformOutput', false),...
        siteSegs,baselineInds,'UniformOutput', false);
    taskWindow = {[phaseWinSz, 0]};
    taskAlignmentPoints = {["GoSignal" "StartLift"],["GoSignal","StartLift"],...
        ["GoSignal","StartHold"]};
    %   baselineTable = getTrialPhaseTable(phaseBaseline(1:length(condLabels)), "Baseline", condLabels,sdm);
    %   taskCondTable = getTrialPhaseTable(taskFR(1:length(condLabels)),"Task",condLabels,sdm);
    %   rawBCounts = cellfun(@(p,s) cellfun(@(t,a) cellfun(@(h) bootstrapBaseline(h,t), a,'UniformOutput',false),...
    %     p(~cellfun(@isempty,p)),num2cell(s(~cellfun(@isempty,p)),2),'UniformOutput',false),...
    %     spCounts(1:end-1),baselineEnds,'UniformOutput',false);
    taskEnds = cellfun(@(a,pa) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,pw) ...
        findBins(ap(:,p)+pw,params.bins),{pa},taskWindow,'UniformOutput', false),ps,'UniformOutput',false),...
        a, 'UniformOutput', false),siteSegs,cellfun(@(t,cs) arrayfun(@(ti) find(strcmp(ti,cs)), t), ...
        taskAlignmentPoints, params.condSegMap.values,'UniformOutput', false),'UniformOutput',false);
    % spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(ha,pp)cellfun(@(h) (cellfun(@(tp,hp)...
    %     conv2(tp(:,max(1,hp(1)):max(1,hp(1))+200),ones(1,2),'valid').*...
    %     [repmat([1 NaN], size(tp,1), 100)],squeeze(num2cell(cell2mat(reshape(pp,size(pp,1),1,size(pp,2))),[1 2])),...
    %     num2cell(h,2),'UniformOutput',false))', ha, 'UniformOutput',false), a,tt,'UniformOutput',false),...
    %     p(~cellfun(@isempty,p)),(s(~cellfun(@isempty,p),:)),'UniformOutput',false),rawSpikes,taskEnds,'UniformOutput',false);
    spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(ha,pp)cellfun(@(h) cell2mat(cellfun(@(tp,hp)...
    sum(tp(:,max(1,hp(1)):min(size(tp,2),hp(end))),2)',squeeze(num2cell(cell2mat(reshape(pp,size(pp,1),1,size(pp,2))),[1 2])),...
    num2cell(h,2),'UniformOutput',false))', ha, 'UniformOutput',false), a,tt,'UniformOutput',false),...
    p(~cellfun(@isempty,p)),num2cell(s(~cellfun(@isempty,p),:),2),'UniformOutput',false),rawSpikes,psthPhaseEnds,'UniformOutput',false);
    %%
    spCounts = cellfun(@(cb) [vertcat(cb{:})] ,spCounts, 'UniformOutput',false);
    % spCounts = cellfun(@(c) cellfun(@(a) cellfun(@(s) cellfun(@(m)...
    %     m(:,1:2:min(cellfun(@(o) min(cellfun(@(n) size(n,2), o)),a))), s,'UniformOutput',false),...
    %     a, 'Uniformoutput', false),c),spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(c) cellfun(@(a) num2cell(cell2mat(reshape(cellfun(@(s) median(cell2mat(reshape(s,1,1,[])),3,'omitnan'),...
        num2cell(vertcat(a{:}),1), 'UniformOutput',false),1,1,[])),[2 3]),num2cell(c,2),'UniformOutput',false),spCounts,'UniformOutput',false);
    spCounts = cellfun(@(s) vertcat(s{:}), spCounts, 'UniformOutput',false);
    maxSitesPerCond = max(cellfun(@length, spCounts));
    spCounts(cellfun(@length,spCounts)~=maxSitesPerCond) = cellfun(@(s,r) vertcat(s,...
        cellfun(@(t) NaN(size(t)), r(end-(size(r,1)-size(s,1)-1):end), 'UniformOutput',false)),...
        spCounts(cellfun(@length,spCounts)~=maxSitesPerCond),spCounts(find(...
        cellfun(@length,spCounts)==maxSitesPerCond,1)),'UniformOutput',false);
    %%
    spCounts = cellfun(@(c1,c2,c3) cellfun(@(s1,s2,s3) [horzcat(squeeze(s1),NaN(size(s1,3),max(cellfun(@(z) size(z,2),{s1,s2,s3}))-size(s1,2))),...
        horzcat(squeeze(s2),NaN(size(s2,3),max(cellfun(@(z) size(z,1),{s1,s2,s3}))-size(s2,2))),horzcat(squeeze(s3),NaN(size(s3,3),...
        max(cellfun(@(z) size(z,2),{s1,s2,s3}))-size(s3,2)))],vertcat(c1,repmat(c1,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c1,1),1)),...
        vertcat(c2,repmat(c2,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c2,1),1)),vertcat(c3,repmat(c3,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c3,1),1)),'UniformOutput',false),...
        spCounts(1), spCounts(2), spCounts(3), 'UniformOutput',false);
    counts = vertcat(spCounts{:});
    allUnitsTrials = cellfun(@transpose,counts,'UniformOutput',false);
    clear rawCounts rawBCounts spCounts counts rawSpikes
end
allUnitsTrials = allUnitsTrials(tUnitsInd);
%% decoder setup
%somatotopyTrialLabs = cellfun(@(s,r) cellstr(repmat(string(r),size(s,1),1)), trialConds, unitSomatotopy,'UniformOutput',false);
savePath = loadPath + "Decoding\All\100_Bins\";
if(~exist(savePath,'dir'))
    mkdir(savePath);
end
ftLength = 1;
num_labels = 2;
num_cv_splits = 10;
testUnits = [1, 2, 5, 10, 15, 25, 35, 50, 75, 100];
binning_parameters = struct('sampling_interval', 1,'end_time', size(allUnitsTrials{1},2),'start_time', 1,'bin_width',ftLength);
binning_parameters.the_bin_start_times = binning_parameters.start_time:binning_parameters.sampling_interval:binning_parameters.end_time;
binning_parameters.the_bin_widths = repmat(binning_parameters,1,length(binning_parameters.the_bin_start_times));
somatotopicLabs = unique(unitSomatotopy);
somatotopicLabs(end+1) = {'Arm Hand'};
somatotopicLabs(end+1) = {'Arm Hand Face Trunk'};
if(strcmpi(countsOrFR,"counts"))
    fp = {};
    cle = poisson_naive_bayes_CL;
else
    fp = {};%zscore_normalize_FP;
    cle = max_correlation_coefficient_CL;
end
goodInds = cellfun(@(p) ~all(isnan(p(:,end-3:end)),2), allUnitsTrials, 'UniformOutput',false);
trialUnits = cellfun(@(n,b) n(b,:), allUnitsTrials, goodInds, 'UniformOutput',false);
trialUnitLabs = cellfun(@(c,b) c(b,:), trialConds, goodInds, 'UniformOutput', false);
if(strcmpi(countsOrFR,"FR"))
    if(ftLength==2)
        trialUnits = cellfun(@(t,s) cell2mat(arrayfun(@(tt) [tt*s(1),tt*s(end)], t(:,9:12), 'UniformOutput',false)), trialUnits, xyLoc,'UniformOutput',false);
        %trialUnits = cellfun(@(t)t(:,1:2:8), trialUnits, 'UniformOutput',false);
    else
        trialUnits = cellfun(@(t)t(:,9:12), trialUnits, 'UniformOutput',false);
    end
    goodUnitsF = cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
    trialUnits = cellfun(@(t,b) t(b,:), trialUnits, goodUnitsF, 'UniformOutput',false);
    trialUnitLabs = cellfun(@(t,b) t(b,:), trialUnitLabs, goodUnitsF, 'UniformOutput',false);
end
subpopulations = cellfun(@(s) cellfun(@(b) contains(unitSomatotopy,strsplit(s))...
    & sum(b,2)>0, fUnits,'UniformOutput',false), somatotopicLabs, 'UniformOutput',false);

goodUnits = find_sites_with_k_label_repetitions(trialUnitLabs,num_cv_splits*num_labels,{'E','L','P'});
trialC = trialUnits(goodUnits);
trialLabs = trialUnitLabs(goodUnits);
unitMonkey = monkeyUnitInd(goodUnits);
currSomatotopy = unitSomatotopy(goodUnits);
trialInds =  cellfun(@(t) true(size(t,1),1), trialC, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialC, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs, trialInds, 'UniformOutput',false);
if(ftLength==1)
    dsr = avg_DS(trainingSet,trainingLabs,num_cv_splits,num_labels);
    dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
else
    dsr = basic_DS(trainingSet,trainingLabs,num_cv_splits);
    dsr.binned_site_info.binning_parameters = binning_parameters;
end
dsr.num_times_to_repeat_each_label_per_cv_split = num_labels;
dsr.time_periods_to_get_data_from = arrayfun(@(a) unique([a,a+binning_parameters.bin_width-1]),...
    binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time,'UniformOutput',false);
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;%min(nUnits,length(siteIDs));
dsr.sample_sites_with_replacement = 0;%length(siteIDs)<nUnits;cv = standard_resample_CV(dsr,cle);
%%
CVAcc = NaN(num_cv_splits,length(dsr.time_periods_to_get_data_from),1);%length(dsr.time_periods_to_get_data_from));
somaUnits = repmat({repmat({CVAcc},length(fUnits),length(testUnits),length(somatotopicLabs))},numRuns,1);
continueIterations = true;
for iter = 1:numRuns
    %clT = repmat({cle},1,length(dsr.time_periods_to_get_data_from));
    iterAcc = CVAcc;
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    for s = 1:length(somatotopicLabs)
        for f = 1:length(fTypes)
            for n = 1:length(testUnits)
                popInds = find(ismember(goodUnits,find(subpopulations{s}{f})));
                unitSample = popInds(randperm(length(popInds),min(testUnits(n), length(popInds))));
                for iCV = 1:num_cv_splits
                    parfor iTrainingInterval = 1:size(all_XTr,2)
                        if(isempty(fp))
                            XTrF = all_XTr{iTrainingInterval}{iCV};
                        else
                            [~,XTrF] = fp.set_properties_with_training_data(all_XTr{iTrainingInterval}{iCV});
                        end
                        if(strcmp(countsOrFR,"counts"))
                            XTrF = fix(XTrF);
                        end
                        XTrF = XTrF(unitSample,:);
                        if(ftLength==1)
                            %cl = clT{iTrainingInterval};
                            %clT{iTrainingInterval} = cl.train(XTrF, all_YTr);
                            clT = cle.train(XTrF,all_YTr);
                        else
                            clT = fitcknn(XTrF',all_YTr,'NSMethod','kdtree');
                            clT.Distance = 'eucledian';
                            clT.NumNeighbors = 3;
                        end
                        XTrFt = fix(all_XTrt{iTrainingInterval}{iCV});
                        XTrFt = XTrFt(unitSample,:);
                        [predicted_labels decision_values] = clT.test(XTrFt);
                        acc(iTrainingInterval) = sum(predicted_labels==all_Ytrt)/length(all_Ytrt);
                        % for iTestInterval = 1:size(all_XTrt,2)
                        %     if(isempty(fp))
                        %         XTrFt = all_XTrt{iTestInterval}{iCV};
                        %     else
                        %         [~,XTrFt] = fp.set_properties_with_training_data(all_XTrt{iTestInterval}{iCV});
                        %     end
                        %     if(strcmp(countsOrFR,"counts"))
                        %         XTrFt = fix(XTrFt);
                        %     end
                        %     XTrFt = XTrFt(unitSample,:);
                        %     if(ftLength==1)
                        %         [predicted_labels decision_values] = clT{iTrainingInterval}.test(XTrFt);
                        %     else
                        %         predicted_labels = clT.predict(XTrFt');
                        %     end
                        %     iterAcc(iCV, iTrainingInterval,iTestInterval) = sum(predicted_labels==all_Ytrt)/length(all_Ytrt);
                        % end
                    end
                    iterAcc(iCV,:) = acc;
                end
                somaUnits{iter}{f,n,s} =iterAcc;               
            end
        end
    end
    if(mod(iter,25)==0)
        disp(num2str(100*(iter/numRuns))+"%");
    end
    if(iter>numRuns)
        continueIterations = false;
        %elseif(convergence_values < .001)
        %    continueIterations = false;
    end
    %[absoluteDiff, percChange] = cv.get_convergence_values(squeeze(mean(CVAcc(iter,:,:),2,'omitnan')));
    %convergence_values = [absoluteDiff, percChange];
    %convergence_values = convergence_values(1).* 100;
    %iter = iter+1;
end
%%
close all;
cls = distinguishable_colors(length(somatotopicLabs));
ls = {'-','-.'};
thresh = [];
for s = 1:length(somatotopicLabs)
    for f = 1:length(fTypes)
        for itrain =1:length(phaseNames)
            figure(f);
            subplot(2,2,itrain);
            hold on;
            title(phaseNames(itrain))%(f) + " Units");
            allAcc = cellfun(@(u) cell2mat(cellfun(@(v) v(:,itrain,itrain), u(f,:,s), 'UniformOutput',false)), somaUnits, 'UniformOutput', false);
            allAcc = cell2mat(reshape(allAcc,1,1,[]));
            errorbar(testUnits,mean(allAcc,[1,3],'omitnan'),std(allAcc,0,[1,3],'omitnan')./sqrt(numRuns*num_cv_splits),'color',cls(s,:),'Linestyle',ls(1+(s>4)));
            inter = polyxpoly(testUnits,mean(allAcc,[1,3],'omitnan'),testUnits,repmat(0.9,size(testUnits)));
            if(~isempty(inter))
                thresh(f,s,itrain) = fix(inter(1));
            else
                thresh(f,s,itrain) = NaN;
            end
            if(s==length(somatotopicLabs))
                if(itrain==1)
                    legend(somatotopicLabs,'Location','northwest','AutoUpdate','off')
                end
                ylim([0.3 1])
                plot([1 max(testUnits)],[0.90 0.90],'k--','LineWidth',1);
                plot([1 max(testUnits)],[0.33 0.33],'k:','LineWidth',1.5);
            end
        end
    end
end
%%
for t = 1:4
    saveFigures(figure(t),savePath+"Summary\",phaseNames(t),[]);
end
save(savePath+"allIters",'somaUnits','somatotopicLabs','fTypes','phaseNames','testUnits','thresh');

% for ut = 1:length(somatotopicLabs)
        %     dst = basic_DS(testingSet,testingLabs,1);
        %     dst.num_times_to_repeat_each_label_per_cv_split = num_labels;
        %     dst.binned_site_info.binning_parameters = binning_parameters;
        %     dst.num_resample_sites = dsr.num_resample_sites;
        %     dst.sample_sites_with_replacement = 1;
        %     if(ut==length(somatotopicLabs))
        %         siteIds = 1:length(currSomatotopy);
        %     else
        %         siteIds = find(strcmp(currSomatotopy,somatotopicLabs{ut}));
        %     end
        %     dst.sites_to_use = siteIds;
        %     iter = 1;
        %     CVAcc = NaN(numTestSamples-1,iCV,binning_parameters.end_time,binning_parameters.end_time);
        %     while(iter<numTestSamples)
        %         [~, ~, all_XTe, all_YTe] = dst.get_data;
        %         % [~, ~, all_XTe, all_YTe] = deal({});
        %         % for cc = 1:num_cv_splits
        %         %     currLabs = randperm(length(siteIds));
        %         %     labInds = cellfun(@(l) currLabs(strcmp(testingLabs(currLabs),l)), condLabels,'UniformOutput',false);
        %         %     for r = 1:binning_parameters.end_time
        %         %         rTrialSeeds = cellfun(@(lm) randi(length(lm),nLabels*num_cv_splits,1),labInds,'UniformOutput',false);
        %         %         all_XTe{r}{cc} = cell2mat(cellfun(@(li,s)reshape(testingSet(li(s),r),[],nLabels),...
        %         %             labInds,rTrialSeeds,'UniformOutput',false)')';
        %         %         all_YTe = cat(1,repmat(1,num_cv_splits,1),repmat(2,num_cv_splits,1),...
        %         %         repmat(3,num_cv_splits,1));
        %         %     end
        %         % end
        %         for iCV = 1:length(all_XTe{1})
        %             for iTrainingInterval = 1:binning_parameters.end_time
        %                 cl = clT{iTrainingInterval};
        %                 for iTestInterval = 1:binning_parameters.end_time
        %                     YTe = all_YTe;
        %                     if(isempty(fp))
        %                         XTe = all_XTe{iTestInterval}{iCV};
        %                     else
        %                         [~,XTe] = fp.set_properties_with_training_data(all_XTe{iTestInterval}{iCV});
        %                     end
        %                     if(iter>0)
        %                         [predicted_labels decision_values] = cl.test(XTe);
        %                         CVAcc(iter, iCV,iTrainingInterval, iTestInterval) = sum(predicted_labels==YTe)/length(YTe);
        %                     end
        %                 end
        %             end
        %         end
        %         iter = iter+1;
        %     end
        %
        %end
% binning_parameters.sampling_interval = params.binSize*1000;
% binning_parameters.bin_width = params.binSize*1000;
% binning_parameters.start_time = -2000;
% binning_parameters.end_time = 2000;
% binning_parameters.the_bin_start_times = findBins(binning_parameters.start_time/1000,params.bins);
% binning_parameters.the_bin_widths = 10;
%%
function avgCounts = bootstrapBaseline(b,ps)
BOOTLINEIN = 12;
for c = 1:BOOTLINEIN
    cb = cellfun(@(bb) randi([bb(1) bb(end)-20],1),num2cell(b(all(~isnan(b),2),:),2),'UniformOutput',false);
    avgCounts(:,:,c) = cell2mat(cellfun(@(tp,hp)sum(tp(:,hp(1):hp(1)+20),2), ...
        squeeze(num2cell(cell2mat(reshape(ps{1}(:,all(~isnan(b),2)),size(ps{1},1),1,length(cb))),[1 2])),...
        cb,'UniformOutput',false)');
end
avgCounts =fix(mean(avgCounts,3));
end
%%
function trialCondTable =  getTrialPhaseTable(phaseFR,phaseNames,condLabs,monkeyLabs)
unitTrialPhase = cellfun(@(c) cellfun(@(s) median(cell2mat(reshape(cellfun(@(a) cat(3,a{:}), s, 'UniformOutput', false),...
    [1,1,1,size(s,2)])),4,'omitnan'),c, 'UniformOutput',false),phaseFR(1:length(condLabs)),'UniformOutput',false);
unitTrialPhase = cellfun(@(cs) cellfun(@(c) cat(2,cat(1,c,...
    NaN([max(cellfun(@(s) size(s,1), cs))-size(c,1),size(c,[2,3])])),...
    NaN([size(c,1),max(cellfun(@(s) size(s,2),cs))-size(c,2),size(c,3)])),...
    cs, 'UniformOutput', false), num2cell([unitTrialPhase{:}],2), 'UniformOutput',false);
unitTrialPhase = cellfun(@(s,un) cellfun(@(c) cellfun(@(u,n) array2table(...
    [permute(u,[2 3 1]),repmat(n,size(u,2),1)],'VariableNames',[phaseNames, "Unit"]),...
    num2cell(c,[2 3]),num2cell(un:un+size(c,1)-1)','UniformOutput',false),s,'UniformOutput',false),unitTrialPhase,...
    num2cell(cumsum([1;cellfun(@(c) max(cellfun(@(r) size(r,1), c)), unitTrialPhase(1:end-1))])),'UniformOutput',false);
trialCondTable = cellfun(@(t) cellfun(@(c,cn) addvars(vertcat(c{:}),repmat(cn,height(vertcat(c{:})),1),...
    'NewVariableNames','Condition'),t,condLabs,'UniformOutput',false), unitTrialPhase, 'UniformOutput', false);
trialCondTable = vertcat(trialCondTable{:});
trialCondTable = cellfun(@(t,c) cellfun(@(u) addvars(u,repmat(c,height(u),1),'NewVariableNames', 'Condition'),...
    t, 'UniformOutput',false),vertcat(unitTrialPhase{:}),repmat(condLabs,length(trialCondTable),1),'UniformOutput',false);
trialCondTable = cellfun(@(a,c) cellfun(@(s) addvars(s,repmat(c,height(s),1),'NewVariableNames','Monkey'), a, ...
    'UniformOutput',false),trialCondTable,repmat(num2cell(monkeyLabs.Monkey),1,size(trialCondTable,2)),'UniformOutput',false);
trialCondTable = vertcat(trialCondTable{:});
trialCondTable = vertcat(trialCondTable{:});
end
%%
function [phaseSubVals, sigs] = ttestTrials(dist1,dist2,taskPhase,paired,pVal)
if(paired)
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)>=max(1,floor(length(dist1)/2));
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}./...
        max(1,d1{taskPhase})),dist1,dist2,'UniformOutput',false)',1);
else
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest2(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)'>=floor(length(dist1)/2);
    smallestTrialCount = min(size(dist1{1}{taskPhase},2),size(dist2{1}{taskPhase},2));
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}(:,1:smallestTrialCount)-...
        d1{taskPhase}(:,1:smallestTrialCount)),dist1,dist2,'UniformOutput',false)',1);
end
phaseSubVals = cellfun(@(m) median(cat(3,m{:}),3), phaseSubVals, 'UniformOutput', false);
sigs = double(sigs);
end
%% OLD %%
%% average PSTH Units
%% averageFR
% avgTable = readtable(loadPath+ "Task_Units_FR_AVGPSTH.xlsx");
% unit2CondTable = table();
% for a = 1:height(avgTable)
%     unit2CondTable(end+1:end+3,:) = repmat(avgTable(a,[1:8,12]),3,1);
% end
% unit2CondTable = addvars(unit2CondTable,repmat(['E';'L';'P'],height(avgTable),1),'NewVariableNames',"Condition");
% baselineAvgs = m.avgBaseline;
% baselineAvgs = cellfun(@(c) cellfun(@(s) median(cell2mat(s),2,'omitnan'),c,'UniformOutput', false),...
%     baselineAvgs(~cellfun(@isempty,baselineAvgs)),'UniformOutput', false);
% baselineAvgs = array2table(cell2mat(cellfun(@(c) vertcat(c{:}), baselineAvgs, 'UniformOutput',false)),'VariableNames',cellfun(@(s) ['Baseline_',s],condLabels,'UniformOutput', false));
% baselineAvgs = addvars(baselineAvgs,transpose(1:height(baselineAvgs)), 'NewVariableNames','Unit');
% avgTable = join(avgTable,baselineAvgs,"Keys",'Unit');
% stackVars = [arrayfun(@(p) find(contains(avgTable.Properties.VariableNames,p)), phaseNames,'Uniformoutput', false),...
%     find(~strcmp(avgTable.Properties.VariableNames,"TaskUnits") & contains(avgTable.Properties.VariableNames,"TaskUnits"))];
% avgFR = stack(avgTable,stackVars,'ConstantVariables',~ismember(1:width(avgTable),[stackVars{:}]),...
%     'NewDataVariableName',[phaseNames,"TaskUnit"],'IndexVariableName',"Condition");
% avgFR = arrayfun(@(t) table2array(avgFR(avgFR.Unit==t,cellstr(phaseNames))),unique(avgFR.Unit),'UniformOutput',false);
% avgFR = avgFR(tUnitsInd);
%%
% nonModUnits = unique(join(cellfun(@string,table2cell(avg(find(table2array(avg(:,'TaskUnit')~=1)),["Unit", "Condition"])))));
% trialFRS = trialFRS(~contains(join([string(trialFRS.Unit),string(trialConds)]), nonModUnits),:);
% trialFRS = stack(trialFRS,phaseNames,'ConstantVariables', ~contains(trialFRS.Properties.VariableNames,phaseNames), ...
%     "NewDataVariableName","FR", 'IndexVariableName',"Phase");
% FR =arrayfun(@(u) trialFRS(trialFRS.Unit==u,{'Unit','Condition','Phase','FR'}), unique(trialFRS{:,'Unit'}), 'UniformOutput',false);
% FR = cellfun(@(u) unstack(addvars(u,cumsum(strcmp(string(u{:,'Phase'}),phaseNames(1))),...
%     'NewVariableNames','Trial'),'FR','Phase'), FR, 'UniformOutput',false);
% condLabels = cellfun(@(u) cellstr(table2cell([u(:,'Condition')])), FR, 'UniformOutput',false);
% FR = cellfun(@(u) table2array(u(:,contains(u.Properties.VariableNames,phaseNames))), FR, 'UniformOutput',false)';
% missingTrials = cellfun(@(n) any(any(isnan(n))), FR);
% goodSamples = ~missingTrials;
% missingTrials = find(missingTrials);
% for s = 1:length(missingTrials)
%     sNanTrials = find(any(isnan(FR{missingTrials(s)}),2));
%     for c = 1:length(sNanTrials)
%         condTrialNaN = condLabels{missingTrials(s)}(sNanTrials(c));
%         rSample = goodSamples(randsample(sum(goodSamples),fix(sum(goodSamples)/10)));
%         allTrials = vertcat(FR{rSample});
%         FR{missingTrials(s)}(sNanTrials(c),:) = mean(allTrials(...
%             cellfun(@(ls) strcmp(ls,condTrialNaN),vertcat(condLabels{rSample})),:),1,'omitnan');
%     end
% end
%% load even odd avgs
% even = readtable(loadPath+"Task_Units_FR_EVENAVGPSTH.xlsx");
% odd = readtable(loadPath+"Task_Units_FR_ODDAVGPSTH.xlsx");
% varTableNames = even.Properties.VariableNames;
% even = stack(even,arrayfun(@(p) find(contains(varTableNames,p)), phaseNames,'Uniformoutput', false),...
%    'ConstantVariables',1,'NewDataVariableName',phaseNames, 'IndexVariableName',"Condition");
% even.Condition = arrayfun(@(p) varTableNames{p}(end), even.Condition);
% odd = stack(odd,arrayfun(@(p) find(contains(varTableNames,p)), phaseNames,'Uniformoutput', false),...
%    'ConstantVariables',1,'NewDataVariableName',phaseNames, 'IndexVariableName',"Condition");
% odd.Condition = arrayfun(@(p) varTableNames{p}(end), odd.Condition);
%% clean data
% tbs = {even,odd,avgPSTH};
% for a = 1:length(tbs)
%     cleanTable = tbs{a};
%     varNames = cleanTable.Properties.VariableNames;
%     for c = 1:length(unique(cleanTable.Condition))
%         currConds = cleanTable.Condition == condNames{c};
%         fillUnits = any(isnan(table2array(cleanTable(:,contains(varNames,phaseNames))))& currConds,2);
%         fillUnitsInd = find(fillUnits);
%         for s = 1:length(fillUnitsInd)
%             cleanTable(fillUnitsInd(s),contains(varNames,phaseNames)) = mean(cleanTable(~fillUnits(...'
%                 randsample(sum(~fillUnits),fix(sum(~fillUnits)/10))),contains(varNames,phaseNames)),1,'omitnan');
%         end
%         cleanTable{any(isnan(table2array(cleanTable(:,contains(varNames,["rSI","gSI"])))),2),contains(varNames,["rSI","gSI"])} = 0;
%     end
%     cleanTable = filloutliers(cleanTable,"clip","gesd","DataVariables",'FR');
%     if(any(~contains(unit2CondTable.Properties.VariableNames,cleanTable.Properties.VariableNames)))
%         cleanTable = join(unit2CondTable,cleanTable,"Keys",["Unit","Condition"]);
%     end
%     cleanTable = stack(cleanTable,["Go","Reach","Grasp","Withdraw"],'ConstantVariables',~contains(cleanTable.Properties.VariableNames,phaseNames),'NewDataVariableName','FR','IndexVariableName','Phase');
%     tbs{a} = cleanTable;
% end
% even = tbs{1};
% odd = tbs{2};
% avg = tbs{3};
% avg = unstack(avg,'FR',{'Phase'});
% avg = unstack(avg,contains(avg.Properties.VariableNames,[phaseNames,"TaskUnit"]),'Condition');
%% Trial PSTHS
% phaseFR = m.phaseFR;
% condXrep = cellfun(@(c) cellfun(@(s) median(cell2mat(reshape(cellfun(@(a) cat(3,a{:}), s, 'UniformOutput', false),...
%      [1,1,1,size(s,2)])),4,'omitnan'),c, 'UniformOutput',false),phaseFR(1:3),'UniformOutput',false);
% unitTrialPhase = cellfun(@(cs,ms,mt) cellfun(@(c) cat(2,cat(1,c,NaN([ms-size(c,1),size(c,[2,3])])),...
%     NaN([size(c,1),mt-size(c,2),size(c,3)])), cs, 'UniformOutput', false),...
%     num2cell([condXrep{:}],2), num2cell(unitCount), num2cell(max(cell2mat(cellfun(@(c) cellfun(@(s) size(s,2),...
%     c), condXrep,'UniformOutput',false)),[],2)),'UniformOutput',false);
% unitNums = num2cell([[1;1+cumsum(cellfun(@(s) size(s{1},1),unitTrialPhase(1:end-1,1)))],...
%     cumsum(cellfun(@(s) size(s{1},1), unitTrialPhase))],2);
% trialTable = cellfun(@(s,un) cellfun(@(c) cellfun(@(u,n) array2table(...
%     [permute(u,[2 3 1]),repmat(n,size(u,2),1)],'VariableNames',[phaseNames, "Unit"]),...
%     num2cell(c,[2 3]),num2cell(un(1):un(end))','UniformOutput',false),s,'UniformOutput',false),...
%     unitTrialPhase,unitNums,'UniformOutput',false);
% trialCondTable = cellfun(@(c) varfun(@(s) vertcat(s{:}),cell2table([c{:}],...
%     "VariableNames",["E","L","P"])), trialTable, 'UniformOutput', false);
% trialCondTable = varfun(@(t) vertcat(t{:,:}), vertcat(trialCondTable{:}));
% frs = stack(trialCondTable,1:width(trialCondTable));
% trialConds = table2array(frs(:,1));
% trialConds = arrayfun(@(s) s{1}(end), string(trialConds));
% trialFRS = array2table(table2array(frs(:,end)),'VariableNames',[phaseNames,"Unit"]);
% trialFRS = addvars(trialFRS,trialConds,'NewVariableNames',"Condition");
% % nonModUnits = unique(join(cellfun(@string,table2cell(avg(find(table2array(avg(:,'TaskUnit')~=1)),["Unit", "Condition"])))));
% % trialFRS = trialFRS(~contains(join([string(trialFRS.Unit),string(trialConds)]), nonModUnits),:);
% trialFRS = stack(trialFRS,phaseNames,'ConstantVariables', ~contains(trialFRS.Properties.VariableNames,phaseNames), ...
%     "NewDataVariableName","FR", 'IndexVariableName',"Phase");
% cellTrialFRS =arrayfun(@(u) trialFRS(trialFRS.Unit==u,{'Unit','Condition','Phase','FR'}), unique(trialFRS{:,'Unit'}), 'UniformOutput',false);
% cellTrialFRS = cellfun(@(u) unstack(addvars(u,cumsum(strcmp(string(u{:,'Phase'}),phaseNames(1))),...
%     'NewVariableNames','Trial'),'FR','Phase'), cellTrialFRS, 'UniformOutput',false);
% condLabels = cellfun(@(u) cellstr(table2cell([u(:,'Condition')])), cellTrialFRS, 'UniformOutput',false);
% cellTrialFRS = cellfun(@(u) table2array(u(:,contains(u.Properties.VariableNames,phaseNames))), cellTrialFRS, 'UniformOutput',false)';
% missingTrials = cellfun(@(n) any(any(isnan(n))), cellTrialFRS);
% goodSamples = ~missingTrials;
% missingTrials = find(missingTrials);
% for s = 1:length(missingTrials)
%     sNanTrials = find(any(isnan(cellTrialFRS{missingTrials(s)}),2));
%     for c = 1:length(sNanTrials)
%         condTrialNaN = condLabels{missingTrials(s)}(sNanTrials(c));
%         rSample = goodSamples(randsample(sum(goodSamples),fix(sum(goodSamples)/10)));
%         allTrials = vertcat(cellTrialFRS{rSample});
%         cellTrialFRS{missingTrials(s)}(sNanTrials(c),:) = mean(allTrials(...
%             cellfun(@(ls) strcmp(ls,condTrialNaN),vertcat(condLabels{rSample})),:),1,'omitnan');
%     end
% end
%% avgPSTHS clean
% avgPSTHS = cellfun(@(r) transpose(squeeze(r)),num2cell(mean(normFR,4,'omitnan'),[2 3]),'UniformOutput',false);
% for c = 1:length(unique(conditions))
%     fillUnits = cellfun(@(r) find(all(isnan(r),1)), avgPSTHS, 'UniformOutput', false);
%     fillUnitsInd = find(~cellfun(@isempty,fillUnits));
%     for s = 1:length(fillUnitsInd)
%         condInd = find(all(isnan(avgPSTHS{fillUnitsInd(s)}),1));
%         goodInds = cellfun(@isempty,fillUnits);
%         avgPSTHS{fillUnitsInd(s)}(:,condInd) = mean([avgPSTHS{goodInds(...'
%             randsample(sum(goodInds),fix(sum(goodInds/10))),condInd)}],2,'omitnan');
%     end
% end
% avgPSTHS = cellfun(@transpose,squeeze(num2cell(filloutliers(cell2mat(reshape(avgPSTHS,[1, 1, length(avgPSTHS)])),"clip","gesd",1),[1 2])),'UniformOutput',false);
%% dPCA
% firingRates = arrayfun(@(s) firingRates(strcmp(unitSomatotopy,s),:,:,:), somatotopicLabs, 'UniformOutput',false);
% numSomatotopy = cellfun(@(s) size(s,1), firingRates);
% firingRates = cellfun(@(s) cat(1,s,NaN(max(numSomatotopy)-size(s,1), size(s,2), size(s,3), size(s,4))),firingRates, 'UniformOutput',false);
% firingRates = permute(cat(5,firingRates{:}), [1 2 5 3 4]);
% trialNum = repmat(maxTrials,[1 1 length(unique(unitSomatotopy))]);
% for r = 1:length(somatotopicLabs)
%     trialNum(strcmp(unitSomatotopy,somatotopicLabs(r)),:,r) = 0;
% end
% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% margNames = {'Condition', 'Somatotopy', 'Independent', 'Condition/Somatotopy Interaction'};
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
% firingRatesAverage = mean(normFR,5,'omitnan');
% X = firingRatesAverage(:,:);
% X = bsxfun(@minus, X, mean(X,2));
% [W,~,~] = svd(X, 'econ');
% W = W(:,1:20);
% minimal plotting
% dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);
% computing explained variance
% explVar = dpca_explainedVariance(firingRatesAverage, W, W, 'combinedParams', combinedParams);
% a bit more informative plotting
% dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, 'explainedVar', explVar, ...
%     'time', time, 'timeEvents', timeEvents, 'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours);