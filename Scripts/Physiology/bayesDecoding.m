%% load data
conditions = ["Extra Small Sphere","Large Sphere", "Photocell"];
phaseNames = ["Go", "Reach", "Hold", "Withdraw"];
phaseWinSz = 0.2;
phaseWindows = repmat({{[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz*(5/4), -phaseWinSz*(1/4)],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions));
phaseWindows{end}{3} = [-.1 0];
taskAlign = containers.Map(conditions,repmat({{["GoSignal" "StartHold"]}},1,length(conditions)));
phaseAlign = containers.Map(conditions,cellfun(@num2cell,repmat({["GoSignal","StartReach","StartHold","StartWithdraw"]},1,length(conditions)),'UniformOutput',false));
fTypes = ["Reach","Grasp","Both","Shallow","Deep","Task"];
taskPhase = false;
numRuns = 500;
savePath = "S:\Lab\ngc14\Working\PSTHS\Decoding\";
params = PhysRecording(conditions,.01,.15,-6,5,containers.Map(conditions,...
    repmat({"StartReach"},1,length(conditions))));
%%
[sdm,siteSegs,siteTrialPSTHS,rawSpikes,siteChannels,~,simpRep,siteLocation,~,monkeys,~,...
    conditions,chUnitMap,siteTrialInfo] = getAllSessions(params,"Single","M1","");
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[phaseWinSz, 0]}},1,length(conditions)),siteSegs,siteTrialPSTHS,false,true);
[phaseBaseline,phaseFR] = calculatePhases(params,phaseAlign,phaseWindows,siteSegs,siteTrialPSTHS,false,true);
clear siteTrialPSTHS siteTrialInfo
%%
mappedChannels =  cell2mat(cellfun(@(ch,l) ch{end}(l(~isnan(l)))', chUnitMap,siteChannels, 'Uniformoutput', false)');
replaceVals = cell2mat(cellfun(@(c) cellfun(@(s) ~all(cellfun(@iscell,s)) | any(size(s)==[0 0]),c),rawSpikes,'Uniformoutput',false));
for cc = 1:length(rawSpikes)
    rv = replaceVals(:,cc);
    if(any(rv))
        fillSize = num2cell(cell2mat(cellfun(@(m) max(cell2mat(m)),[num2cell(cellfun(@(s) size(s), ...
            [rawSpikes{~ismember(1:size(replaceVals,2),cc)}],'Uniformoutput',false),2),num2cell(cellfun(@(s) size(s{1}),...
            [rawSpikes{~ismember(1:size(replaceVals,2),cc)}], 'UniformOutput', false),2),],'Uniformoutput',false)),2);
        rawSpikes{cc}(rv) = cellfun(@(f) repmat(repmat({repmat({NaN},[1,f(end)])},[1,f(1)]),1,...
            length(cell2mat(params.PSTHAlignments.values(cellstr(conditions(cc)))))),fillSize(rv),'UniformOutput',false);
        siteSegs{cc}(rv) = cellfun(@(f) repmat({NaN(f(end),length(params.condSegMap(string(conditions{cc}))))},1,length(cell2mat(...
            params.PSTHAlignments.values(cellstr(conditions(cc)))))),fillSize(rv),'UniformOutput',false);
    end
end
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep));
% xyLoc = num2cell(mapSites2Units(cellfun(@length,siteChannels),num2cell([sdm.x,sdm.y],2)),2);
%% trial by trial condition labels organized by unit
trialCondTable = getTrialPhaseTable(phaseFR(1:length(conditions)),phaseNames,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false),sdm);
trialConds = trialCondTable(:,strcmp(trialCondTable.Properties.VariableNames,'Condition'));
trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t & ~all(isnan(trialCondTable{:,...
    contains(trialCondTable.Properties.VariableNames,phaseNames)}),2),:)),unique(trialCondTable.Unit),'UniformOutput',false);
%% only use units that are task modulated according to the task phase window
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(conditions)),taskFR(1:length(conditions)),'UniformOutput',false);
taskUnits = cellfun(@cell2mat, taskUnits,'UniformOutput',false);
tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(a) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
    a,'UniformOutput',false), s,'UniformOutput', false),phaseFR(1:length(conditions)),'UniformOutput', false), phaseNames, 'UniformOutput',false);
[rgVals,rgInds] = cellfun(@(cr,cg) cellfun(@(r,g) ttestTrials(r,g,1,true,0.05),cr,cg,'UniformOutput',false),...
    tPhase{strcmp(phaseNames,"Reach")},tPhase{strcmp(phaseNames,"Hold")}, 'UniformOutput',false);
rgInds = cellfun(@cell2mat, rgInds, 'UniformOutput',false);
AUCVals = cellfun(@(c) cell2mat(cellfun(@(d) mean(cell2mat(d),2,'omitnan'),c,'UniformOutput',false)), rgVals,'UniformOutput',false);
rUnits = cellfun(@(rg,r) rg==1 & r < 1, rgInds, AUCVals, 'UniformOutput',false);
gUnits = cellfun(@(rg,r) rg==1 & r > 1, rgInds, AUCVals, 'UniformOutput',false);
bothUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, taskUnits, 'UniformOutput',false);
fUnits = {cell2mat(rUnits), cell2mat(gUnits), cell2mat(bothUnits),repmat(mappedChannels<=16,1,length(conditions)), ...
    repmat(mappedChannels>16,1,length(conditions)),cell2mat(taskUnits)};
goodInds = ~any(cell2mat(cellfun(@(n) all(isnan(n),2),AUCVals,'UniformOutput',false)),2) & any(cell2mat(taskUnits),2);
unitSomatotopy = unitSomatotopy(goodInds);
mappedChannels = mappedChannels(goodInds);
fUnits = cellfun(@(f) f(goodInds,:), fUnits, 'UniformOutput',false);
trialLabs = trialConds(goodInds);
clear trialCondTable tPhase phaseFR taskBaseline
%% get trial spike counts
if(taskPhase)
    windowPad = 200;
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
        ap(:,p)+pw,num2cell(pa),pw,'UniformOutput', false),ps,'UniformOutput',false),...
        a, 'UniformOutput', false),siteSegs,phaseInds,phaseWindows,'UniformOutput', false);
    spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(h) cell2mat(cellfun(@(tp,hp)...
        cellfun(@(i) sum(hp>=i(1) & hp<=i(end)),tp), cellfun(@(s) num2cell(squeeze(s),1),...
        num2cell(cat(3,a{1}{:}),[2 3]),'UniformOutput',false)',h,'UniformOutput',false)'), ...
        tt, 'UniformOutput',false),p,s,'UniformOutput',false),rawSpikes,psthPhaseEnds,'UniformOutput',false);
    spCounts = cellfun(@(cb) [cb{:}]' ,spCounts, 'UniformOutput',false);
    spCounts = cellfun(@(c) cellfun(@(a)num2cell(vertcat(a{:}),1),num2cell(c(goodInds),2),'UniformOutput',false),spCounts,'UniformOutput',false);
end
spCounts = cellfun(@(s) cellfun(@(t) squeeze(horzcat(t{:})),s,'UniformOutput',false), spCounts, 'UniformOutput',false);
allUnitsTrials = cellfun(@round,horzcat(spCounts{:}),'UniformOutput',false);
clear spCounts
%% decoder setup
somatotopicLabs = unique(unitSomatotopy);
if(~exist(savePath,'dir'))
    mkdir(savePath);
end
fp = {};
num_repeated_labels = 3;
num_cv_splits = 10;
testUnits = [1, 2, 5, 10, 15, 25, 35, 50, 75, 100];
binning_parameters = struct('sampling_interval', 1,'end_time', size(allUnitsTrials{1},2),'start_time', 1,'bin_width',1);
binning_parameters.the_bin_start_times = binning_parameters.start_time:binning_parameters.sampling_interval:binning_parameters.end_time;
binning_parameters.the_bin_widths = repmat(diff([0,binning_parameters.the_bin_start_times]),1,length(binning_parameters.the_bin_start_times));
timepoints = length(binning_parameters.the_bin_start_times);
subpopulations = cellfun(@(s) cellfun(@(b) contains(unitSomatotopy,strsplit(s))...
    & b(:,1)>0, fUnits,'UniformOutput',false), somatotopicLabs, 'UniformOutput',false);

goodUnits = find_sites_with_k_label_repetitions(trialConds(goodInds),num_cv_splits*num_repeated_labels,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false));
trialUnits = cellfun(@(n,b) vertcat(n{:}), num2cell(allUnitsTrials(goodUnits,:),2), 'UniformOutput',false);
trialLabs = trialLabs(goodUnits);
trialInds =  cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs, trialInds, 'UniformOutput',false);

dsr = avg_DS(trainingSet,trainingLabs,num_cv_splits,num_repeated_labels);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.num_times_to_repeat_each_label_per_cv_split = num_repeated_labels;
dsr.time_periods_to_get_data_from = arrayfun(@(a) unique([a,a+binning_parameters.bin_width-1]),...
    binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time,'UniformOutput',false);
dsr.nAvg = num_repeated_labels;
dsr.sample_sites_with_replacement = 0;
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
%%
somaUnits = repmat({repmat({NaN(num_cv_splits,timepoints,timepoints)},length(fTypes),length(testUnits),length(somatotopicLabs))},numRuns,1);
uUnits = repmat({cell(length(fTypes),length(testUnits),length(somatotopicLabs))},numRuns,1);
rng('shuffle');
cl = poisson_naive_bayes_CL;
hbar = parforProgress(numRuns);
parfor iter = 1:numRuns
    rs = RandStream.create('threefry4x64_20','NumStreams',numRuns,'StreamIndices',iter,'Seed','shuffle');
    RandStream.setGlobalStream(rs);
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    for s = 1:length(somatotopicLabs)
        for f = 1:length(fTypes)
            for n = 1:length(testUnits)
                popUnits = find(ismember(find(goodUnits),find(subpopulations{s}{f})));
                units =  popUnits(randperm(length(popUnits),min(length(popUnits),testUnits(n))));
                iterAcc = NaN(num_cv_splits,timepoints,timepoints);
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
                            clT= cl.train(fix(tr(units,:)), all_YTr);
                            [ia,~] = clT.test(fix(XTst(units,:)));
                            iterAcc(iCV,iTrainingInterval,iTestingInterval) = sum(ia-all_Ytrt==0)/length(all_Ytrt);
                        end
                    end
                end
                somaUnits{iter}{f,n,s} = iterAcc;
                uUnits{iter}{f,n,s} = units;
                %clT= cl.train(fix(tr),all_YTr);[~,ui]= fastTest(clT.lambdas,XTst);cell2mat(cellfun(@(u) sum((u-all_Ytrt==0))/length(all_Ytrt),ui,'UniformOutput',false));
            end
        end
    end
    send(hbar, iter);
end
delete(gcp('nocreate'));
unitAccPhase = cellfun(@(p) cell2mat(permute(cellfun(@(m) mean(m,1,'omitnan'),p,'Uniformoutput',false),...
    (length(size(p{1}))+ length(size(p))):-1:1)),somaUnits,'UniformOutput',false);
%uUnits = cellfun(@(p) cellfun(@(r) resize(r,[1,max(testUnits)],'FillValue',NaN),p,'UniformOutput',false),uUnits,'UniformOutput',false);
%%
trCl = [0 .7 0; 1 .5 0; .6 0 .2; 0 1 1;];
lnStyle = {'--',':','-.','-'};
figure();hold on;
for tr = 1:timepoints
    lnStyle = circshift(lnStyle,1);
    for ts = 1:timepoints
        distCols = 100.*cell2mat(cellfun(@(i) squeeze(mean(i(:,tr,ts,:,:,:),[4 6],'omitnan')),unitAccPhase,'UniformOutput',false)');
        errorbar(testUnits,mean(distCols,2,'omitnan')',std(distCols,0,2)','LineStyle',lnStyle{ts},'Color',trCl(tr,:));
    end
end
plot(get(gca,'XLim'),[85 85],'LineWidth',1,'Color','k'); ylim([25 100]);
saveFigures(gcf,savePath,"NUnits_"+num2str(numRuns),[]);
%%
plotPhase = ["Go","Reach","Hold"];
nUnits = 50;
typeGroup = cellfun(@(s) squeeze(s)',num2cell(100.*cell2mat(reshape(cellfun(@(i) mean(cell2mat(reshape(arrayfun(@(p) squeeze(i(:,p,...
p,:,testUnits==nUnits,:)),find(contains(string(phaseNames),plotPhase)),'UniformOutput',false),1,1,[])),3,'omitnan'),unitAccPhase,'UniformOutput',false),...
1,1,[])),[1 3]),'UniformOutput',false);
accTable = cellfun(@(t,n) array2table(t,'VariableNames',n+"_"+string(somatotopicLabs)'),typeGroup,fTypes,UniformOutput=false);
accTable = [accTable{:}];
writetable(accTable,savePath+"Decoding_"+num2str(nUnits)+"_"+num2str(numRuns),'FileType','spreadsheet','UseExcel',true);

figure(); hold on;
bx=boxchart(100.*cell2mat(cellfun(@(i) cell2mat(arrayfun(@(n) squeeze(mean(i(:,n,n,:,testUnits==nUnits,:),[4 6],'omitnan')),...
    1:timepoints,'UniformOutput',false)),unitAccPhase,'UniformOutput',false)),'Notch','on','MarkerStyle','none');
xticklabels(string(phaseNames));plot(get(gca,'XLim'),[85 85],'LineWidth',1,'Color','k');ylim([25 100]);
saveFigures(gcf,savePath,"Phases_"+num2str(numRuns)+"_"+num2str(nUnits),[]);


figure(); hold on;
bx=boxchart(100.*mean(cell2mat(cellfun(@(i) permute(mean(cell2mat(reshape(arrayfun(@(p) squeeze(i(:,p,p,:,testUnits==nUnits,:)),...
    find(contains(string(phaseNames),plotPhase)),'UniformOutput',false),1,1,[])),3,'omitnan'),[3 1 2]),unitAccPhase,'UniformOutput',false)),3,'omitnan'),'Notch','on','MarkerStyle','none');
xticklabels(somatotopicLabs);plot(get(gca,'XLim'),[85 85],'LineWidth',1,'Color','k');ylim([25 100]);
saveFigures(gcf,savePath,"Somatotopy_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

figure(); hold on;
bx=boxplotGroup(typeGroup([1 3 2]),'Notch','on','Symbol','','secondarylabels',somatotopicLabs,'primarylabels',["","",""],'interGroupSpace',2,'Colors',[1 0 0; 1 1 0; 0 0 1]);
plot(get(gca,'XLim'),[85 85],'LineWidth',1,'Color','k');ylim([25 100]);
allBxs = vertcat(bx.boxplotGroup(1:end).Children);
set(allBxs(contains(string({allBxs.Tag}),"Whisker")),'LineStyle','-')
saveFigures(gcf,savePath,"Types_"+num2str(numRuns)+"_"+num2str(nUnits),[]);


figure(); hold on;
bx=boxplotGroup(typeGroup(4:end-1),'Notch','on','Symbol','','secondarylabels',somatotopicLabs,'primarylabels',["",""],'interGroupSpace',2,'Colors',[.8 .8 .8; .1 .1 .1]);
plot(get(gca,'XLim'),[85 85],'LineWidth',1,'Color','k');ylim([25 100]);
allBxs = vertcat(bx.boxplotGroup(1:end).Children);
set(allBxs(contains(string({allBxs.Tag}),"Whisker")),'LineStyle','-');
saveFigures(gcf,savePath,"Laminar_"+num2str(numRuns)+"_"+num2str(nUnits),[]);
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
function trialCondTable =  getTrialPhaseTable(phaseFR,phaseNames,condLabs,siteTable)
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
    'UniformOutput',false),trialCondTable,repmat(num2cell(siteTable.Monkey),1,size(trialCondTable,2)),'UniformOutput',false);
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