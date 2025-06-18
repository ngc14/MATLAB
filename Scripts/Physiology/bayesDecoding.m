%% load data
loadPath = "S:\Lab\ngc14\Working\";%"E:\Data\"
phaseNames = ["Go", "Reach", "Grasp", "Withdraw"];
fTypes = ["Reach","Grasp","Both","Shallow","Deep"];
condLabels = {'E','L','P'};
taskPhase = false;
numRuns = 500;
mtform = load('S:\Lab\ngc14\Working\Both\S2G_Reference_tform.mat');
phaseWinSz = 0.2;
m = matfile("S:\Lab\ngc14\Working\Both\Baseline_FR\phaseAnalysis_Face.mat");

rawSpikes = m.rawSpikes;
avgPhase = m.avgPhase;
siteSegs = m.siteSegs;
siteChannels = m.siteChannels;
siteChannels = siteChannels{2};
phaseFR = m.phaseFR;
taskBaseline = m.taskBaseline;
sdm = m.siteDateMap;
taskFR = m.taskFR;
siteImaging = m.siteActiveInd;
chUnitMap = cell(height(sdm),1);
chUnitMap(strcmp([sdm.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([sdm.Monkey],"Skipper")) = {[32:-1:1]};
mappedChannels = cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap,siteChannels, 'Uniformoutput', false);
mappedChannels = [mappedChannels{:}];
params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5);
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), m.siteRep'));
monkeyUnitInd = cellstr(mapSites2Units(cellfun(@length, siteChannels), sdm.Monkey));
siteImaging = cellfun(@(si) mapSites2Units(cellfun(@length, siteChannels), si), siteImaging, 'UniformOutput',false);
xyLoc = num2cell(mapSites2Units(cellfun(@length,siteChannels),num2cell([sdm.x,sdm.y],2)),2);
xyLoc(strcmp(monkeyUnitInd,"Skipper")) = cellfun(@(s) round(transformPointsForward(mtform.tform,s)), xyLoc(strcmp(monkeyUnitInd,"Skipper")),'UniformOutput',false);
add_ndt_paths_and_init_rand_generator
%% trial by trial condition labels organized by unit
trialCondTable = getTrialPhaseTable(phaseFR(1:length(condLabels)),phaseNames,condLabels,sdm);
trialConds = trialCondTable(:,strcmp(trialCondTable.Properties.VariableNames,'Condition'));
trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t,:)),unique(trialCondTable.Unit),'UniformOutput',false);
%% only use units that are task modulated according to the task phase window
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(condLabels)),taskFR(1:length(condLabels)),'UniformOutput',false);
taskUnits = cell2mat(cellfun(@cell2mat, taskUnits,'UniformOutput',false));
tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(a) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
    a,'UniformOutput',false), s,'UniformOutput', false),phaseFR(1:length(condLabels)),'UniformOutput', false), phaseNames, 'UniformOutput',false);
[~,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials(r,g,1,true,0.05),cr,cg,'UniformOutput',false),...
    tPhase{strcmp(phaseNames,"Reach")},tPhase{strcmp(phaseNames,"Grasp")}, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
AUCVals = cellfun(@(c) cell2mat(c), avgPhase(1:length(condLabels)),'UniformOutput',false);
rAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Reach")), AUCVals, 'UniformOutput',false);
gAUC = cellfun(@(a) a(:,strcmp(phaseNames,"Grasp")), AUCVals, 'UniformOutput',false);
rUnits = cellfun(@(rg,r,g) rg==1 & r > g, rgInds, rAUC, gAUC, 'UniformOutput',false);
gUnits = cellfun(@(rg,r,g) rg==1 & r < g, rgInds, rAUC, gAUC, 'UniformOutput',false);
bothUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, num2cell(taskUnits,1), 'UniformOutput',false);
fUnits = {cell2mat(rUnits), cell2mat(gUnits), cell2mat(bothUnits),repmat(mappedChannels'<16,1,3), repmat(mappedChannels'>=16,1,3)};
%fUnits{end+1} = taskUnits;
tUnitsInd = sum(taskUnits,2)>0;

% trialConds = trialConds(tUnitsInd);
unitSomatotopy = unitSomatotopy(tUnitsInd);
monkeyUnitInd = string(monkeyUnitInd(tUnitsInd));
xyLoc  = xyLoc(tUnitsInd);
fUnits = cellfun(@(f) f(tUnitsInd,:), fUnits, 'UniformOutput',false);
%fUnits = cellfun(@(f) cellfun(@(l) f & l, lUnits, 'UniformOutput', false), fUnits, 'UniformOutput',false);
clear trialCondTable tPhase phaseFR phaseBaseline taskBaseline
%% get trial spike counts
if(taskPhase)
    taskWindow = {[phaseWinSz, 0]};
    windowPad = 200;
    taskAlignmentPoints = {["GoSignal" "StartLift"],["GoSignal","StartLift"],...
        ["GoSignal","StartHold"]};
    taskEnds = cellfun(@(a,pa) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,pw) ...
        findBins(ap(:,p)+pw,params.bins),{pa},taskWindow,'UniformOutput', false),ps,'UniformOutput',false),...
        a, 'UniformOutput', false),siteSegs,cellfun(@(t,cs) arrayfun(@(ti) find(strcmp(ti,cs)), t), ...
        taskAlignmentPoints, params.condSegMap.values,'UniformOutput', false),'UniformOutput',false);
    spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(ha,pp)cellfun(@(h) (cellfun(@(tp,hp)...
        conv2(tp(:,max(1,hp(1)):max(1,hp(1))+windowPad),ones(1,2),'valid').*...
        [repmat([1 NaN], size(tp,1), windowPad/2)],squeeze(num2cell(cell2mat(reshape(pp,size(pp,1),1,size(pp,2))),[1 2])),...
        num2cell(h,2),'UniformOutput',false))', ha, 'UniformOutput',false), a,tt,'UniformOutput',false),...
        p(~cellfun(@isempty,p)),(s(~cellfun(@isempty,p),:)),'UniformOutput',false)',rawSpikes,taskEnds,'UniformOutput',false);
    spCounts = cellfun(@(cb) cellfun(@(a) cell2mat(reshape(cat(1,a{:}),1,1,[],size(a,2))), ...
        cellfun(@(c) [c{:}],cb,'UniformOutput',false), 'UniformOutput',false),spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(c) cellfun(@(m) median(m(:,1:2:end,:,:),4,'omitnan'),c,'UniformOutput',false),spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(s)cellfun(@(n) num2cell(n,[2 3]), s, 'UniformOutput', false), spCounts, 'UniformOutput',false);
else
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
    spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(ha,pp)cellfun(@(h) cell2mat(cellfun(@(tp,hp)...
        sum(tp(:,max(1,hp(1)):min(size(tp,2),hp(end))),2)',squeeze(num2cell(cell2mat(reshape(pp,size(pp,1),1,size(pp,2))),[1 2])),...
        num2cell(h,2),'UniformOutput',false))', ha, 'UniformOutput',false), a,tt,'UniformOutput',false),...
        p(~cellfun(@isempty,p)),num2cell(s(~cellfun(@isempty,p),:),2),'UniformOutput',false),rawSpikes,psthPhaseEnds,'UniformOutput',false);
    spCounts = cellfun(@(cb) [vertcat(cb{:})] ,spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(c) cellfun(@(a) num2cell(cell2mat(reshape(cellfun(@(s) median(cell2mat(reshape(s,1,1,[])),3,'omitnan'),...
        num2cell(vertcat(a{:}),1), 'UniformOutput',false),1,1,[])),[2 3]),num2cell(c,2),'UniformOutput',false),spCounts,'UniformOutput',false);
end
spCounts = cellfun(@(s) squeeze(vertcat(s{:})), spCounts, 'UniformOutput',false);
maxSitesPerCond = max(cellfun(@length, spCounts));
spCounts(cellfun(@length,spCounts)~=maxSitesPerCond) = cellfun(@(s,r) vertcat(s,...
    repmat(cellfun(@(t) NaN(size(t)), r(end-(size(r,1)-size(s,1)-1):end), 'UniformOutput',false),1,size(s,2))),...
    spCounts(cellfun(@length,spCounts)~=maxSitesPerCond),spCounts(find(...
    cellfun(@length,spCounts)==maxSitesPerCond,1)),'UniformOutput',false);
spCounts = cellfun(@(c1,c2,c3) cellfun(@(s1,s2,s3) [horzcat(squeeze(s1)',NaN(size(s1,3),...
    max(cellfun(@(z) size(z,2),{s1,s2,s3}))-size(s1,2))),horzcat(squeeze(s2)',NaN(size(s2,3),...
    max(cellfun(@(z) size(z,1),{s1,s2,s3}))-size(s2,2))),horzcat(squeeze(s3)',NaN(size(s3,3),...
    max(cellfun(@(z) size(z,2),{s1,s2,s3}))-size(s3,2)))],vertcat(c1,repmat(c1,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c1,1),1)),...
    vertcat(c2,repmat(c2,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c2,1),1)),vertcat(c3,repmat(c3,max(cellfun(@(m) size(m,1),{c1,c2,c3}))-size(c3,1),1)),'UniformOutput',false),...
    spCounts(1), spCounts(2), spCounts(3), 'UniformOutput',false);
counts = vertcat(spCounts{:});
allUnitsTrials = cellfun(@transpose,counts,'UniformOutput',false);
clear rawCounts rawBCounts spCounts counts rawSpikes
%% decoder setup
savePath = loadPath + "Decoding\All\Bayes\Train_Test\";
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
somatotopicLabs = somatotopicLabs(end-1:end);
ds = cellfun(@(p) ~all(isnan(p(:,end-3:end)),2), allUnitsTrials, 'UniformOutput',false);
fp = {};
goodInds = cellfun(@(p) ~all(isnan(p(:,end-3:end)),2), allUnitsTrials, 'UniformOutput',false);
trialUnits = cellfun(@(n,b) n(b,:), allUnitsTrials, goodInds, 'UniformOutput',false);
trialUnitLabs = cellfun(@(c,b) c(b,:), trialConds, goodInds, 'UniformOutput', false);
subpopulations = cellfun(@(s) cellfun(@(b) contains(unitSomatotopy,strsplit(s))...
    & sum(b,2)>0, fUnits,'UniformOutput',false), somatotopicLabs, 'UniformOutput',false);

goodUnits = find_sites_with_k_label_repetitions(trialUnitLabs,num_cv_splits*num_labels,{'E','L','P'});
trialUnits = trialUnits(goodUnits);
trialLabs = trialUnitLabs(goodUnits);
unitMonkey = monkeyUnitInd(goodUnits);
currSomatotopy = unitSomatotopy(goodUnits);
trialInds =  cellfun(@(t) true(size(t,1),1), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs, trialInds, 'UniformOutput',false);

dsr = avg_DS(trainingSet,trainingLabs,num_cv_splits,num_labels);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.num_times_to_repeat_each_label_per_cv_split = num_labels;
dsr.time_periods_to_get_data_from = arrayfun(@(a) unique([a,a+binning_parameters.bin_width-1]),...
    binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time,'UniformOutput',false);
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
dsr.sample_sites_with_replacement = 0;
%%
timepoints = length(binning_parameters.the_bin_start_times);
iterAcc = NaN(num_cv_splits,timepoints,timepoints);
unitAcc = cell(num_cv_splits,timepoints,timepoints);
somaUnits = repmat({repmat({iterAcc},length(fTypes),1,length(somatotopicLabs))},numRuns,1);
uUnits = repmat({repmat({iterAcc},length(fTypes),1,length(somatotopicLabs))},numRuns,1);
cl = poisson_naive_bayes_CL;
parfor iter = 1:numRuns
    %%
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    tic;
    for s = 1:length(somatotopicLabs)
        for f = 1:length(fTypes)
            for n = 1:length(testUnits)
                popInds = find(ismember(goodUnits,find(subpopulations{s}{f})));
                %unitSample = [1:length(goodUnits)]';
                unitSample = popInds(randperm(length(popInds),min(testUnits(n), length(popInds))));
                iterAcc = NaN(num_cv_splits,timepoints,timepoints);
                unitAcc = cell(num_cv_splits,timepoints,timepoints);
                for iCV = 1:num_cv_splits
                    for iTrainingInterval = 1:timepoints
                        XTrF = all_XTr{iTrainingInterval};
                        for iTestingInterval = 1:timepoints
                            XTst = all_XTrt{iTestingInterval};
                            if(isempty(fp))
                                tr = XTrF{iCV};
                                XTst = XTst{iCV};
                            else
                                [~,tr] = fp.set_properties_with_training_data(tr{iCV});
                                [~,XTst] = fp.set_properties_with_training_data(XTst{iCV});
                            end
                            XTst = fix(XTst);
                            clT = cl.train(fix(tr), all_YTr);
                            [ia,ui] = fastTest(clT.lambdas(unitSample,:),XTst(unitSample,:));
                            iterAcc(iCV,iTrainingInterval,iTestingInterval) = sum(ia-all_Ytrt==0)/length(all_Ytrt);
                            ui = cellfun(@(u) sum((u-all_Ytrt==0))/length(all_Ytrt), ui, 'UniformOutput',false);
                            unitAcc{iCV,iTrainingInterval,iTestingInterval} = cell2mat(ui);
                        end
                    end
                end
                uUnits{iter} = unitAcc;
                somaUnits{iter}{f,n,s} = iterAcc;
            end
        end
    end
    toc
end
% unitAccPhase = mean(cell2mat(reshape(cellfun(@(n) squeeze(mean(n,2)), cellfun(@(p) reshape([p{:}],[1410,10,4]), ...
%     uUnits, 'UniformOutput',false), 'UniformOutput',false),1,1,[])),3);
% somas = strsplit(string(somatotopicLabs(end)));
% close all;
% for s = 1:length(somas)
%     subplot(3,length(somas),s);
%     hold on;
%     title(somas(s));
%     histogram(unitAccPhase(fUnit{s},2),0.3:.02:.7,"Normalization","percentage");
%     h = histogram(unitAccPhase(strcmp(currSomatotopy,somas(s)),2),0.3:.02:.7,"Normalization","cdf","Visible","off");
%     scatter(h.BinEdges(find(h.Values>.5,1)),0,100,'red','|',"LineWidth",2)
% end
% fUnit = cellfun(@(f) sum(f(goodUnits,:),2)>0, fUnits, 'UniformOutput', false);
% for s = 1:length(fUnit)
%     subplot(3,length(somas),s+3);
%     hold on;
%     title(fTypes(s));
%     histogram(unitAccPhase(fUnit{s},2),0.3:.02:.7,"Normalization","percentage");
%     h = histogram(unitAccPhase(fUnit{s},2),0.3:.02:.7,"Normalization","cdf","Visible","off");
%     scatter(h.BinEdges(find(h.Values>.5,1)),0,100,'red','|',"LineWidth",2)
% end
% close all
clear clA clT all_XTr all_XTrt taskFR AUCVals trialsUnits trainingSet
%%
ss = cellfun(@(c) cell2mat(cellfun(@(n) n{1}, c,'UniformOutput', false)), siteSegs, 'UniformOutput',false);
allAcc = cell2mat(reshape(cellfun(@(s) cell2mat(cellfun(@(t) ...
    reshape(mean(t,1,'omitnan'),[ones(1,length(size(s))),size(t,2)]), s, 'UniformOutput', false)), somaUnits, 'UniformOutput',false),...
    [ones(1,length(size(somaUnits{1}))+1),size(somaUnits,1)]));
avgSegs = cellfun(@(c) discretize([0, cumsum(abs(diff(mean(cell2mat(cellfun(@(n) n{1}, c,'UniformOutput', false)),1,'omitnan'))))],...
    0:params.binSize:length(params.bins)),siteSegs, 'UniformOutput',false);
avgAlignments = arrayfun(@(n) cellfun(@(c,a) (c(strcmp(a, n))-c(1))/(windowPad*params.binSize), avgSegs, ...
    params.condSegMap.values(cellstr(params.condNames))), ["StartReach","StartHold","StartWithdraw"],'UniformOutput',false);
cl = {'k','r'};
scl = distinguishable_colors(4,[cl,'w']);
for n = 1:length(testUnits)
    figure('Units','normalized','Position',[0 0 1 1]);
    for f = 1:length(fTypes)
        subplot(2,2,f);
        hold on;
        title(fTypes(f))
        for s = [1,3,4]
            e =errorbar(1:length(scl),100*squeeze(mean(allAcc(f,n,s,:,:),5)), squeeze(std(allAcc(f,n,s,:,:),0,5)).*100,...
                'Color',scl(s,:),'LineWidth',2-double(mod(s,2)==0));
            e.CapSize = 0;
        end
        ylim([30 100])
        % for a = 1:length(avgAlignments)-1
        %     if(a<2)
        %         avgP = mean(avgAlignments{a});
        %     else
        %         avgP = mean(avgAlignments{a}(1:end-1));
        %         avgP(end+1) = avgAlignments{a}(end);
        %     end
        %     for s = 1:length(avgP)
        %         plot(repmat(fix(avgP(s)),1,2), [30 100], ['-.',cl{s}],'LineWidth',1.25);
        %     end
        % end
        % xticklabels((windowPad.*(get(gca,'xtick')./params.binSize))./1000)
        xlim([0 60]);
    end
    legend(somatotopicLabs([1,3,4]))
    saveFigures(figure(n),savePath+"Summary\R+G\",num2str(testUnits(n)),[]);
end
%%
close all;
cls = distinguishable_colors(length(phaseNames));
ls = {'-','-.'};
thresh = [];
for s = 1:length(somatotopicLabs)
    for f = 1:length(fTypes)
        for itrain =1:length(phaseNames)
            for itest=1:length(phaseNames)
            figure(f);
            subplot(2,2,itrain);
            hold on;
            title(phaseNames(itrain))%(f) + " Units");
            allAcc = cellfun(@(u) cell2mat(cellfun(@(v) v(:,itrain,itest), u(f,:,s), 'UniformOutput',false)), somaUnits, 'UniformOutput', false);
            allAcc = cell2mat(reshape(allAcc,1,1,[]));
            errorbar(testUnits,mean(allAcc,[1,3],'omitnan'),std(allAcc,0,[1,3],'omitnan')./sqrt(numRuns*num_cv_splits),'color',cls(itest,:),'Linestyle',ls(1+(s>4)));
            inter = polyxpoly(testUnits,mean(allAcc,[1,3],'omitnan'),testUnits,repmat(0.9,size(testUnits)));
            if(~isempty(inter))
                thresh(f,s,itrain) = fix(inter(1));
            else
                thresh(f,s,itrain) = NaN;
            end
            if(s==length(somatotopicLabs) && itest==length(phaseNames))
                if(itrain==1)
                    legend(phaseNames,'Location','northwest','AutoUpdate','off')
                end
                ylim([0.3 1])
                plot([1 max(testUnits)],[0.90 0.90],'k--','LineWidth',1);
                plot([1 max(testUnits)],[0.33 0.33],'k:','LineWidth',1.5);
            end
            end
        end
    end
end
%%
for t = 1:4
    saveFigures(figure(t),savePath+"Summary\",phaseNames(t),[]);
end
save(savePath+"allIters",'somaUnits','somatotopicLabs','fTypes','phaseNames','testUnits');
%%
function [predVals,units] = fastTest(cl,tst)
curr_lambdas = repmat(cl, [1, 1, size(tst, 2)]);
XTe_repmat_for_all_classes = permute(repmat(tst, [1, 1, size(cl, 2)]), [1 3 2]);
unitLiklihood = -curr_lambdas + XTe_repmat_for_all_classes  .* log(curr_lambdas) - gammaln(XTe_repmat_for_all_classes  + 1);
log_likelihoods = sum(unitLiklihood, 1);
[vals inds] = randmax(permute(log_likelihoods, [2 3 1]));
[~,units] = cellfun(@(f) randmax(permute(f,[2 3 1])), num2cell(unitLiklihood,[2 3]), 'UniformOutput',false);
units = cellfun(@transpose, units, 'UniformOutput',false);
predVals = inds';
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