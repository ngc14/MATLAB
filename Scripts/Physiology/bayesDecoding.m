%% load data
winSz = 0.2;
fTypes = ["Reach","Grasp","Both","Shallow","Deep","Task"];
conditions = ["Extra Small Sphere","Large Sphere", "Photocell"];
phaseNames = ["Go","Reach","Hold","Withdraw"];
params = PhysRecording(conditions,.01,.15,-6,5,containers.Map(conditions,...
    repmat({"StartReach"},1,length(conditions))));
taskAlign = containers.Map(conditions,repmat({{["GoSignal" "StartHold"]}},1,length(conditions)));
phaseAlign = containers.Map(conditions,cellfun(@(c) num2cell(string(c)),repmat({{"GoSignal","StartReach","StartHold","StartWithdraw"}},...
    1,length(conditions)),'UniformOutput',false));
phaseWin = repmat({{[0, winSz],[-winSz*(3/4),winSz*(1/4)],[-winSz*(5/4), -winSz*(1/4)],[-winSz*(1/4),winSz*(3/4)]}},1,length(conditions));
phaseWin{end}{3} = [-winSz/2 0];
savePath = "S:\Lab\ngc14\Working\Decoding\NTrials\";
%%
[siteDateMap,siteSegs,siteTrialPSTHS,rawSpikes,siteChannels,chMaps,~,~] = ...
    getAllSessions(params,"Single","M1","");
%%
simpRep =  cellfun(@(r,t) r(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh,'UniformOutput', true)';
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[winSz, 0]}},1,length(conditions)),siteSegs,siteTrialPSTHS,false,true);
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(conditions)),taskFR(1:length(conditions)),'UniformOutput',false);
taskUnits = cellfun(@cell2mat, taskUnits,'UniformOutput',false);
[~,phaseFR] = calculatePhases(params,phaseAlign,phaseWin,siteSegs,siteTrialPSTHS,false,true);
clear siteTrialPSTHS
%%
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
% xyLoc = num2cell(mapSites2Units(cellfun(@length,siteChannels),num2cell([sdm.x,sdm.y],2)),2);
%% only use units that are task modulated according to the task phase window
mappedChannels =  cell2mat(cellfun(@(ch,l) ch{end}(l(~isnan(l)))', chMaps,siteChannels, 'Uniformoutput', false)');
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep));
[rgVals,rgInds] = cellfun(@(crg) cellfun(@(a) cellfun(@(rg)ttestTrials({rg(1)},{rg(end)},1,true,0.05),...
    a,'UniformOutput',false),crg,'UniformOutput',false),phaseFR, 'UniformOutput',false);
rgVals = cellfun(@(v) vertcat(v{:}), rgVals, 'UniformOutput',false);
rgInds = cellfun(@(v) cell2mat(vertcat(v{:})), rgInds, 'UniformOutput',false);
rgVals = cellfun(@(c) cell2mat(cellfun(@(d) mean(cell2mat(d),2,'omitnan'),c,'UniformOutput',false)), rgVals,'UniformOutput',false);
rUnits = cellfun(@(rg,r) rg==1 & r < 1, rgInds, rgVals, 'UniformOutput',false);
gUnits = cellfun(@(rg,r) rg==1 & r > 1, rgInds, rgVals, 'UniformOutput',false);
bothUnits = cellfun(@(r,g,t) ~r & ~g & t==1, rUnits, gUnits, taskUnits, 'UniformOutput',false);
fUnits = {cell2mat(rUnits), cell2mat(gUnits), cell2mat(bothUnits),repmat(mappedChannels<=16,1,length(conditions)), ...
    repmat(mappedChannels>16,1,length(conditions)),cell2mat(taskUnits)};
goodInds = ~any(cell2mat(cellfun(@(n) all(isnan(n),2),rgVals,'UniformOutput',false)),2) & any(cell2mat(taskUnits),2);
unitSomatotopy = unitSomatotopy(goodInds);
mappedChannels = mappedChannels(goodInds);
fUnits = cellfun(@(f) f(goodInds,:), fUnits, 'UniformOutput',false);
somatotopicLabs = string(unique(unitSomatotopy)');
fTypes = [fTypes,somatotopicLabs];
%% trial by trial condition labels organized by unit
% trialCondTable = getTrialPhaseTable(phaseFR(1:length(conditions)),phaseNames,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false),sdm);
% trialConds = trialCondTable(:,strcmp(trialCondTable.Properties.VariableNames,'Condition'));
% trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t & ~all(isnan(trialCondTable{:,...
%     contains(trialCondTable.Properties.VariableNames,phaseNames)}),2),:)),unique(trialCondTable.Unit),'UniformOutput',false);
trialConds = cellfun(@(s) num2cell(cell2mat(cellfun(@(c,t) repmat(c(1),size(t{1},1),1),params.condAbbrev.values,s,...
    'UniformOutput',false)')),num2cell([siteSegs{:}],2),'UniformOutput',false);
trialConds = mapSites2Units(cellfun(@length,siteChannels),trialConds');
%% get trial spike counts
xl = linspace(-.75,1.25,1+range([-.75,1.25])./(winSz/2));
phaseInds = cellfun(@(p,cc) cellfun(@(w) arrayfun(@(b)find(strcmp(string(cc),b)),w),p,'UniformOutput',false),phaseAlign,...
    values(params.condSegMap,params.condSegMap.keys()),'UniformOutput',false);
phaseEnds = cellfun(@(a,pa,pw) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,n) ...
    ap(:,p)+n,pa,pw,'UniformOutput', false),ps,'UniformOutput',false),...
    a, 'UniformOutput', false),siteSegs,phaseInds,phaseWin,'UniformOutput', false);
phaseEnds = cellfun(@(a) cellfun(@(ps) cellfun(@(p) num2cell(cell2mat(arrayfun(@(t) t:(winSz/2):max(xl),p(:,3)+min(xl),'UniformOutput',false)),1),...
    ps,'UniformOutput',false),a,'UniformOutput',false),siteSegs,'Uniformoutput',false);
spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(h) cellfun(@(ap) cellfun(@(tp,hp)... sum(hp>=tp(1) & hp<=tp(end))
    sum(hp>=tp(1)-(winSz/2)/2 & hp<=tp(1)+(winSz/2)/2),num2cell(ap(~all(isnan(ap),2),:),2),h(1,~all(isnan(ap),2))'),[a{:}],'UniformOutput',false),tt,'UniformOutput',false), ...
    p,s,'UniformOutput',false),rawSpikes,phaseEnds,'UniformOutput',false);
spCounts = cellfun(@(cb) [cb{:}]' ,spCounts, 'UniformOutput',false);
spCounts = cellfun(@(c) cellfun(@(a) [a{:}], c,'UniformOutput',false),spCounts,'UniformOutput',false);
spCounts = cellfun(@(s) vertcat(s{:}), num2cell([spCounts{:}],2), 'UniformOutput',false);
allUnitsTrials = cellfun(@round,spCounts(goodInds),'UniformOutput',false);
clear spCounts
%% decoder setup
fp = {};
numRuns = 500;
num_repeated_labels = 5;
num_cv_splits =10;
num_fold = 30;
testUnits = [1, 2, 5, 10, 15, 25, 35, 50, 75, 100];
binning_parameters = struct('sampling_interval', 1,'end_time', size(allUnitsTrials{1},2),'start_time', 1,'bin_width',1);
binning_parameters.the_bin_start_times = binning_parameters.start_time:binning_parameters.sampling_interval:binning_parameters.end_time;
binning_parameters.the_bin_widths = repmat(diff([0,binning_parameters.the_bin_start_times]),1,length(binning_parameters.the_bin_start_times));
timepoints = length(binning_parameters.the_bin_start_times);
subpopulations =  cellfun(@(b) contains(unitSomatotopy,["Arm","Hand"]) & b(:,1)>0, fUnits,'UniformOutput',false); %cellfun(@(s) contains(unitSomatotopy,strsplit(s)) somatotopicLabs, 'UniformOutput',false);
subpopulations(end+1:end+length(somatotopicLabs)) = cellfun(@(s) strcmp(unitSomatotopy,s),somatotopicLabs,'UniformOutput',false);
trialLabs = trialConds(goodInds);

goodUnits = find_sites_with_k_label_repetitions(trialLabs,num_fold,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false));
trialUnits = cellfun(@(n,b) vertcat(n{:}), num2cell(allUnitsTrials(goodUnits,:),2), 'UniformOutput',false);
trialLabs = trialLabs(goodUnits);
trialInds =  cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs, trialInds, 'UniformOutput',false);

dsr = basic_DS(trainingSet,trainingLabs,num_cv_splits,num_fold);
dsr.binned_site_info.binning_parameters = binning_parameters;
dsr.num_times_to_repeat_each_label_per_cv_split = num_repeated_labels;
dsr.time_periods_to_get_data_from = arrayfun(@(a) unique([a,a+binning_parameters.bin_width-1]),...
    binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time,'UniformOutput',false);
%dsr.nAvg = 2;
dsr.sample_sites_with_replacement = 0;
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
%%
somaUnits = repmat({repmat({NaN(num_cv_splits,timepoints)},length(fTypes),length(testUnits))},numRuns,1);%,length(somatotopicLabs)
uUnits = repmat({cell(length(fTypes),length(testUnits))},numRuns,1); %,length(somatotopicLabs)
rng('shuffle');
cl = poisson_naive_bayes_CL;
hbar = parforProgress(numRuns);
parfor iter = 1:numRuns
    rs = RandStream.create('threefry4x64_20','NumStreams',numRuns,'StreamIndices',iter,'Seed','shuffle');
    RandStream.setGlobalStream(rs);
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    % for s = 1:length(somatotopicLabs)
       for f = 1:length(fTypes)
           for n = 1:length(testUnits)
               popUnits = find(ismember(find(goodUnits),find(subpopulations{f})));
               units =  popUnits(randperm(length(popUnits),min(length(popUnits),testUnits(n))));
                iterAcc = NaN(num_cv_splits,timepoints);
                for iCV = 1:num_cv_splits
                    for iTrainingInterval = 1:timepoints
                        XTrF = all_XTr{iTrainingInterval};
                        for iTestingInterval = iTrainingInterval%1:timepoints
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
                            iterAcc(iCV,iTrainingInterval) = sum(ia-all_Ytrt==0)/length(all_Ytrt);
                            % clT= cl.train(fix(tr),all_YTr);
                            % [~,ui]= fastTest(clT.lambdas,XTst);
                            % iterAcc{iTrainingInterval,iTestingInterval,iCV}=cellfun(@(u) sum((u-all_Ytrt==0))/length(all_Ytrt),ui);
                        end
                    end
                end
                somaUnits{iter}{f,n} = mean(iterAcc,1,'omitnan');
                uUnits{iter}{f,n} = units;
            end
        end
    send(hbar, iter);
end
delete(gcp('nocreate'));
unitAccPhase = cellfun(@(p) cell2mat(permute(cellfun(@(m) squeeze(mean(m,1,'omitnan')),p,'Uniformoutput',false),...
    (length(size(p{1}))+ length(size(p))):-1:1)),somaUnits,'UniformOutput',false);
%uUnits = cellfun(@(p) cellfun(@(r) resize(r,[1,max(testUnits)],'FillValue',NaN),p,'UniformOutput',false),uUnits,'UniformOutput',false);
%%
nUnits = 35; accLine = 90; trCl = [0 .7 0; 1 .5 0; 0 .2 .8; .4 .4 .4; .6 0 .2;];
accTUnit = cellfun(@(s) s(:,testUnits==nUnits),somaUnits,'UniformOutput',false);
em=cell2mat(cellfun(@(a) a([2,3,6]),cellfun(@(c) mean(cell2mat(cellfun(@(a) cell2mat(cellfun(@(n) mean(n,1,'omitnan'),a,'UniformOutput',false)),c,'UniformOutput',false)),1,'omitnan'),siteSegs,'UniformOutput',false)','UniformOutput',false));
em=[mean(em(:,1:2),1,'omitnan'),min(em(:,end)),max(em(:,end))]; 
[~,ei] = arrayfun(@(a) min(abs(xl-a)),em);
figure(); tiledlayout(1,3);
for d = 0:2
nexttile(); hold on; pll = cellfun(@(p,c) shadedErrorBar(xl,mean(p,1,'omitnan'),.5*std(p,0,1,'omitnan'),'lineProps',{'Color',trCl(c-(3*d),:),'LineWidth',2}),...
    arrayfun(@(a) cell2mat(cellfun(@(u) u{a},accTUnit,'UniformOutput',false)),(3*d)+(1:3+(d==2)),'UniformOutput',false),num2cell((3*d)+(1:3+(d==2))));
arrayfun(@(a,s) plot([a a], [0,1], s),em,["k--","k--","b--","m--"]);
legend([pll.mainLine],fTypes((3*d)+(1:3+(d==2)))); ylim([0 1]);xlim([min(xl),max(xl)]);
end
saveFigures(gcf,savePath,"Temporal_"+num2str(nUnits)+"_"+num2str(num_repeated_labels),[]);
%%
plotPhase = ["Reach","Hold"];
lnStyle = {'--',':','-.','-'};
figure(); hold on;
for tr = 1:timepoints
    lnStyle = circshift(lnStyle,1);
    for ts = 1:timepoints
        trtsPhase = 100.*cell2mat(cellfun(@(i) mean(squeeze(i(tr,ts,:,:,:,:)),2,'omitnan'),unitAccPhase,'UniformOutput',false)');
        errorbar(testUnits,mean(trtsPhase,2,'omitnan')',std(trtsPhase,0,2)','LineStyle',lnStyle{ts},'Color',trCl(tr,:));
        if(tr==ts)
            distCols(:,tr) = trtsPhase(testUnits==nUnits,:);
        end
    end
end
plot([0,max(testUnits)]+[0 1],[accLine accLine],'LineWidth',1,'Color','k'); ylim([0 100]);
%saveFigures(gcf,savePath,"NUnits_"+num2str(numRuns),[]);
%%
typeGroup = cell2mat(reshape(cellfun(@(s) cell2mat(arrayfun(@(e) s(:,e),[ei(1:2),fix(mean(ei(end-1:end))),ei(end)],'UniformOutput',false)),[accTUnit{:}],'UniformOutput',false)',numRuns,1,[]));
accTable = array2table(100.*squeeze(mean(typeGroup(:,2:3,:),2,'omitnan')),'VariableNames',fTypes); % +"_"+string(somatotopicLabs)'
accTable = [accTable,array2table(100.*typeGroup(:,:,strcmp(fTypes,"Task")),'VariableNames',phaseNames+"-Phase")];
%writetable(accTable,savePath+"Decoding_"+num2str(nUnits)+"_"+num2str(numRuns),'FileType','spreadsheet','UseExcel',true);

figure(); hold on;
bx=boxchart(accTable,phaseNames(1:end-1)+"-Phase",'Notch','on','MarkerStyle','none');
xticks(1:length(phaseNames)-1);xticklabels(string(phaseNames(1:end-1)));plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 100]);
saveFigures(gcf,savePath,"Phases_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

figure(); hold on;
bx=boxchart(accTable,["Arm" "Hand" "Trunk" "Face"],'Notch','on','MarkerStyle','none');
xticks(1:length(somatotopicLabs));xticklabels(["Arm" "Hand" "Trunk" "Face"]);
plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 100]);
saveFigures(gcf,savePath,"Somatotopy_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

% figure(); hold on;
% bx=boxchart(100.*mean(cell2mat(cellfun(@(i) permute(mean(cell2mat(reshape(arrayfun(@(p) squeeze(i(p,p,:,:,testUnits==nUnits,:)),...
%     find(contains(string(phaseNames),plotPhase)),'UniformOutput',false),1,1,[])),3,'omitnan'),[3 1 2]),unitAccPhase,'UniformOutput',false)),3,'omitnan'),'Notch','on','MarkerStyle','none');
% xticklabels(somatotopicLabs);plot(get(gca,'XLim'),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 100]);
% saveFigures(gcf,savePath,"Somatotopy_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

figure(); hold on;
bx=boxchart(accTable,["Reach","Both","Grasp"],'Notch','on','MarkerStyle','none');
xticks(1:3);xticklabels(fTypes([1 3 2]));
plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 100]);
% xticklabels(fTypes([1 3 2]));
% allBxs = vertcat(bx.boxplotGroup(1:end).Children);
% set(allBxs(contains(string({allBxs.Tag}),"Whisker")),'LineStyle','-')
saveFigures(gcf,savePath,"Types_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

figure(); hold on;%,'MarkerStyle','none','secondarylabels',fTypes(4:end-1),'primarylabels',["",""],'interGroupSpace',2,'Colors',[.8 .8 .8; .1 .1 .1]);
bx=boxchart(accTable,["Shallow","Deep"],'Notch','on','MarkerStyle','none');
xticks(1:2);xticklabels(["Shallow","Deep"]);
plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 100]);
% allBxs = vertcat(bx.boxplotGroup(1:end).Children);
% set(allBxs(contains(string({allBxs.Tag}),"Whisker")),'LineStyle','-');
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