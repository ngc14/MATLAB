%% load data
phaseNames = ["Go", "Reach", "Hold", "Withdraw"];
fTypes = ["Reach","Grasp","Both","Shallow","Deep"];
conditions = ["Extra Small Sphere","Large Sphere", "Photocell"];
condLabels = {'E','L','P'};
taskPhase = false;
numRuns = 500;
windowPad = 200;
phaseWinSz = 0.2;
m = matfile_m("S:\Lab\ngc14\Working\M1_Physiology\Baseline_FR_OLD\phaseAnalysis_Face.mat");
add_ndt_paths_and_init_rand_generator

params = m.params;
rawSpikes = m.rawSpikes;
siteSegs = m.siteSegs;
siteChannels = m.siteChannels;
siteChannels = siteChannels{2};
phaseFR = m.phaseFR;
taskBaseline = m.taskBaseline;
taskFR = m.taskFR;
sdm = m.siteDateMap;
%%
siteRep = sdm.SiteRep;
simpRep = string(cellfun(@(r,t,f) r{find((f.*t)==(min((f./f).*t)),1)},...
    siteRep,sdm.Thresh, cellfun(@(r)  ~strcmp(r,""), siteRep,'UniformOutput',false),'UniformOutput', false))';
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
chUnitMap = cell(height(sdm),1);
chUnitMap(strcmp([sdm.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([sdm.Monkey],"Skipper")) = {[32:-1:1]};
mappedChannels =  cell2mat(cellfun(@(ch,l) ch(l(~isnan(l)))', chUnitMap,siteChannels, 'Uniformoutput', false));
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep'));
monkeyUnitInd = cellstr(mapSites2Units(cellfun(@length, siteChannels)', sdm.Monkey'));
replaceVals = cell2mat(cellfun(@(c) cellfun(@(s) any(size(s)==[0 0]),c),rawSpikes,'Uniformoutput',false));
for cc = 1:length(rawSpikes)
    rv = replaceVals(:,cc);
    if(any(rv))
        fillSize = cellfun(@(m) max(cell2mat(m')),num2cell(cellfun(@(s) size(s{1}),...
            [rawSpikes{~ismember(1:size(replaceVals,2),cc)}], 'UniformOutput', false),2),'UniformOutput',false);
    end
    rawSpikes{cc}(rv) = cellfun(@(f) repmat({repmat({NaN(1,length(params.bins))},[f(1),f(end)])},1,...
        length(cell2mat(params.PSTHAlignments.values(cellstr(conditions(cc)))))),fillSize(rv),'UniformOutput',false);
    siteSegs{cc}(rv) = cellfun(@(f) repmat({NaN(f(end),length(params.condSegMap(conditions(cc))))},1,length(cell2mat(params.PSTHAlignments.values(cellstr(conditions(cc)))))),...
        fillSize(rv),'UniformOutput',false);
end
xyLoc = num2cell([tPhys.X,tPhys.Y],2);
%params = PhysRecording(conditions,.01,.15,-6,5,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
% siteImaging = cellfun(@(si) mapSites2Units(cellfun(@length, siteChannels), si), siteImaging, 'UniformOutput',false);
% xyLoc = num2cell(mapSites2Units(cellfun(@length,siteChannels),num2cell([sdm.x,sdm.y],2)),2);
% xyLoc(strcmp(monkeyUnitInd,"Skipper")) = cellfun(@(s) round(transformPointsForward(mtform.tform,s)), xyLoc(strcmp(monkeyUnitInd,"Skipper")),'UniformOutput',false);
%% trial by trial condition labels organized by unit
trialCondTable = getTrialPhaseTable(phaseFR(1:length(condLabels)),phaseNames,condLabels,sdm);
trialConds = trialCondTable(:,strcmp(trialCondTable.Properties.VariableNames,'Condition'));
trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t & ~all(isnan(...
    trialCondTable{:,contains(trialCondTable.Properties.VariableNames,phaseNames)}),2),:)),unique(trialCondTable.Unit),'UniformOutput',false);
%% only use units that are task modulated according to the task phase window
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(condLabels)),taskFR(1:length(condLabels)),'UniformOutput',false);
taskUnits = cell2mat(cellfun(@cell2mat, taskUnits,'UniformOutput',false));
tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(a) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
    a,'UniformOutput',false), s,'UniformOutput', false),phaseFR(1:length(condLabels)),'UniformOutput', false), phaseNames, 'UniformOutput',false);
[rgVals,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials(r,g,1,true,0.05),cr,cg,'UniformOutput',false),...
    tPhase{strcmp(phaseNames,"Reach")},tPhase{strcmp(phaseNames,"Hold")}, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
AUCVals = cellfun(@(c) cell2mat(cellfun(@(d) mean(cell2mat(d),2,'omitnan'),c,'UniformOutput',false)), rgVals,'UniformOutput',false);
rUnits = cellfun(@(rg,r) rg==1 & r < 1, rgInds, AUCVals, 'UniformOutput',false);
gUnits = cellfun(@(rg,r) rg==1 & r > 1, rgInds, AUCVals, 'UniformOutput',false);
bothUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, num2cell(taskUnits,1), 'UniformOutput',false);
fUnits = {cell2mat(rUnits), cell2mat(gUnits), cell2mat(bothUnits),repmat(mappedChannels<=16,1,length(condLabels)), repmat(mappedChannels>16,1,length(condLabels))};
fUnits = cellfun(@(f) f(:,1), fUnits, 'Uniformoutput',false);
%fUnits{end+1} = taskUnits;
tUnitsInd = sum(taskUnits,2)>0;
goodInds = ~any(cell2mat(cellfun(@(n) all(isnan(n),2),AUCVals,'UniformOutput',false)),2) & tUnitsInd;
% trialConds = trialConds(tUnitsInd);
unitSomatotopy = unitSomatotopy(goodInds);
fUnits = cellfun(@(f) f(goodInds,:), fUnits, 'UniformOutput',false);
trialConds = trialConds(goodInds);
monkeyUnitInd = string(monkeyUnitInd(goodInds));
xyLoc  = xyLoc(goodInds);
%fUnits = cellfun(@(f) cellfun(@(l) f & l, lUnits, 'UniformOutput', false), fUnits, 'UniformOutput',false);
clear trialCondTable tPhase phaseFR phaseBaseline taskBaseline
%% get trial spike counts
if(taskPhase)
    taskWindow = {[phaseWinSz, 0]};
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
    phaseAlignmentPoints = {["GoSignal","StartReach","StartHold","StartWithdraw"],...
        ["GoSignal","StartReach","StartHold","StartWithdraw"],...
        ["GoSignal","StartReach","StartHold","StartWithdraw"]};
    phaseWindows = repmat({{[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
        [-phaseWinSz*(5/4), -phaseWinSz*(1/4)],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions));
    phaseWindows{end}{3} = [-.1 0];
    phaseInds = cellfun(@(p,cc) arrayfun(@(b)find(strcmp(string(cc),b)),...
        p),phaseAlignmentPoints,values(params.condSegMap,...
        params.condSegMap.keys()),'UniformOutput',false);
    psthPhaseEnds = cellfun(@(a,pa,pw) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,pw) ...
        findBins(ap(:,p)+pw,params.bins),num2cell(pa),pw,'UniformOutput', false),ps,'UniformOutput',false),...
        a, 'UniformOutput', false),siteSegs,phaseInds,phaseWindows,'UniformOutput', false);
    spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(h) cell2mat(cellfun(@(tp,hp)...
        cellfun(@(i) sum(hp(max(1,i(1)):max(1,i(end)))>0),tp), cellfun(@(s) num2cell(squeeze(s),1),...
        num2cell(cat(3,a{1}{:}),[2 3]),'UniformOutput',false)',h,'UniformOutput',false)'), ...
        num2cell(tt{1},2), 'UniformOutput',false),p,s,'UniformOutput',false),rawSpikes,psthPhaseEnds,'UniformOutput',false);
    spCounts = cellfun(@(cb) vertcat(cb{:}) ,spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(c) cellfun(@(a)num2cell(vertcat(a{:}),1),num2cell(c(goodInds),2),'UniformOutput',false),spCounts,'UniformOutput',false);
end
spCounts = cellfun(@(s) cellfun(@(t) squeeze(horzcat(t{:})),s,'UniformOutput',false), spCounts, 'UniformOutput',false);
maxSitesPerCond = max(cellfun(@length, spCounts));
spCounts = cellfun(@(c) cellfun(@(s) s(:,strcmp(phaseNames,"Reach")), c, 'UniformOutput',false),spCounts,'UniformOutput',false);
allUnitsTrials = horzcat(spCounts{:});
clear spCounts
%% decoder setup
savePath = "S:\Lab\ngc14\Working\PSTHS\Decoding\";
somatotopicLabs = unique(unitSomatotopy);
if(~exist(savePath,'dir'))
    mkdir(savePath);
end
fp = {};
ftLength = 1;
num_repeated_labels = 3;
num_cv_splits = 10;
testUnits = 35; [1, 2, 5, 10, 15, 25, 35, 50, 75, 100];
binning_parameters = struct('sampling_interval', 1,'end_time', size(allUnitsTrials{1},2),'start_time', 1,'bin_width',ftLength);
binning_parameters.the_bin_start_times = binning_parameters.start_time:binning_parameters.sampling_interval:binning_parameters.end_time;
binning_parameters.the_bin_widths = repmat(diff([0,binning_parameters.the_bin_start_times]),1,length(binning_parameters.the_bin_start_times));
trialUnits = cellfun(@(n,b) vertcat(n{:}), num2cell(allUnitsTrials,2), 'UniformOutput',false);
trialUnitLabs = trialConds;
subpopulations = cellfun(@(s) cellfun(@(b) contains(unitSomatotopy,strsplit(s))...
    & b(:,1)>0, fUnits,'UniformOutput',false), somatotopicLabs, 'UniformOutput',false);

goodUnits = find_sites_with_k_label_repetitions(trialUnitLabs,num_cv_splits*num_repeated_labels,{'E','L','P'});
trialUnits = trialUnits(goodUnits);
trialLabs = trialUnitLabs(goodUnits);
unitMonkey = monkeyUnitInd(goodUnits);
currSomatotopy = unitSomatotopy(goodUnits);
trialInds =  cellfun(@(t) ~isnan(t), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs, trialInds, 'UniformOutput',false);
fUnit = cellfun(@(f) sum(f(goodUnits,:),2)>0, fUnits, 'UniformOutput', false);

dsr = avg_DS(trainingSet,trainingLabs,num_cv_splits,num_repeated_labels);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.num_times_to_repeat_each_label_per_cv_split = num_repeated_labels;
dsr.time_periods_to_get_data_from = arrayfun(@(a) unique([a,a+binning_parameters.bin_width-1]),...
    binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time,'UniformOutput',false);
dsr.sample_sites_with_replacement = 0;
%%
timepoints = length(binning_parameters.the_bin_start_times);
iterAcc = NaN(num_cv_splits,timepoints,timepoints);
unitAcc = cell(num_cv_splits,timepoints,timepoints);
somaUnits = repmat({repmat({iterAcc},length(fTypes),length(testUnits),length(somatotopicLabs))},numRuns,1);
uUnits = repmat({repmat({iterAcc},length(fTypes),length(testUnits),length(somatotopicLabs))},numRuns,1);
allUnits = cell(numRuns,1);
dsr.sites_to_use = find(goodUnits);
dsr.num_resample_sites = -1;%,testUnits(n);
cl = poisson_naive_bayes_CL;
hbar = parforProgress(numRuns);
parfor iter = 1:numRuns
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    for n = 1:length(testUnits)
        for s = 1:length(somatotopicLabs)
            for f = 1:length(fTypes)
                popInds = find(ismember(find(goodUnits),find(subpopulations{s}{f})));
                unitSample = popInds(randperm(length(popInds),min(testUnits(n), length(popInds))));%unitSample = transpose(1:length(goodUnits));
                iterAcc = NaN(num_cv_splits,timepoints,timepoints);%unitAcc = cell(num_cv_splits,timepoints,timepoints);
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
                            XTst = fix(XTst(unitSample,:));
                            clT= cl.train(fix(tr(unitSample,:)), all_YTr);
                            [ia,~] = clT.test(XTst);
                            iterAcc(iCV,iTrainingInterval,iTestingInterval) = sum(ia-all_Ytrt==0)/length(all_Ytrt);
                            % [~,ui] = fastTest(clT.lambdas,XTst);
                            % unitAcc{iCV,iTrainingInterval,iTestingInterval} = cell2mat(cellfun(@(u) sum((u-all_Ytrt==0))/length(all_Ytrt), ui, 'UniformOutput',false));
                        end
                    end
                end
                uUnits{iter}{f,n,s} = unitSample;
                somaUnits{iter}{f,n,s} = iterAcc;
            end
        end
    end
    send(hbar, iter);
end
delete(gcp('nocreate'));
uUnits = cellfun(@(p) cellfun(@(r) resize(r,[testUnits,num_cv_splits],'FillValue',NaN),p,'UniformOutput',false),uUnits,'UniformOutput',false);
%%
unitAccPhase = mean(cell2mat(reshape(cellfun(@(n) squeeze(mean(n,[1 2],'omitnan')), cellfun(@(p)  cell2mat(permute(p,[2 4 3 1])), ...
    uUnits, 'UniformOutput',false), 'UniformOutput',false),1,1,[])),3);
somas = (string(somatotopicLabs));
close all;
for s = 1:length(somas)
    subplot(3,length(somas),s);
    hold on;
    title(somas(s));
    histogram(unitAccPhase(fUnit{s}),0.3:.02:.7,"Normalization","percentage");
    h = histogram(unitAccPhase(strcmp(currSomatotopy,somas(s))),0.3:.02:.7,"Normalization","cdf","Visible","off");
    scatter(h.BinEdges(find(h.Values>.5,1)),0,100,'red','|',"LineWidth",2)
end
for s = 1:length(fUnit)
    subplot(3,length(somas),s+length(somas));
    hold on;
    title(fTypes(s));
    histogram(unitAccPhase(fUnit{s}),0.3:.02:.7,"Normalization","percentage");
    h = histogram(unitAccPhase(fUnit{s}),0.3:.02:.7,"Normalization","cdf","Visible","off");
    scatter(h.BinEdges(find(h.Values>.5,1)),0,100,'red','|',"LineWidth",2)
end
close all
clear clA clT all_XTr all_XTrt taskFR AUCVals trialsUnits trainingSet
%%
allAcc = cell2mat(reshape(cellfun(@(s) cell2mat(cellfun(@(t) ...
    reshape(mean(t,1,'omitnan'),[ones(1,length(size(s))),size(t,2)]), s, 'UniformOutput', false)), somaUnits, 'UniformOutput',false),...
    [ones(1,length(size(somaUnits{1}))+1),size(somaUnits,1)]));
avgSegs = cellfun(@(c) discretize([0, cumsum(abs(diff(mean(cell2mat(cellfun(@(n) n{1}, c,'UniformOutput', false)),1,'omitnan'))),'omitnan')],...
    0:params.binSize:length(params.bins)),siteSegs, 'UniformOutput',false);
avgAlignments = arrayfun(@(n) cellfun(@(c) (c(strcmp(maxSegL, n))-c(1))/(windowPad*params.binSize), avgSegs), ...
    ["StartReach","StartHold","StartWithdraw"],'UniformOutput',false);
clr = {'k','r'};
scl = distinguishable_colors(4,[clr,'w']);
for n = 1:length(testUnits)
    figure('Units','normalized','Position',[0 0 1 1]);
    for f = 1:length(fTypes)
        subplot(2,3,f);
        hold on;
        title(fTypes(f))
        s = [1,3,4,2];
        e =boxchart(100*squeeze(allAcc(f,n,s,:))','Notch','on','MarkerStyle','none');
            %'Color',scl(s,:),'LineWidth',2-double(mod(s,2)==0));
         set(e,'MarkerStyle','none');

        ylim([30 100]); 
        %xlim([.5 4.5]);xticks(1:4);
        xticklabels(somas(s));
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
        %xlim([0 60]);
    end
    %legend(somatotopicLabs([1,3,4]))
    saveFigures(figure(n),savePath,"Units",[]);
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
            allAcc = cellfun(@(u) cell2mat(cellfun(@(v) v(:,itrain,itest), u{1}(f,:,s), 'UniformOutput',false)), somaUnits, 'UniformOutput', false);
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