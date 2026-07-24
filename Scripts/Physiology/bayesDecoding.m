%% load data
winSz = 0.2;
fTypes = ["Reach","Grasp","Both","Shallow","Deep","Task"];
conditions = ["Extra Small Sphere","Large Sphere", "Photocell"];
phaseNames = ["Go","Reach","Hold","Withdraw"];
winT = [-.5 .5];
params = PhysRecording(conditions,.01,.15,-6,5,containers.Map(conditions,repmat({"StartReach"},1,length(conditions))));
taskAlign = containers.Map(conditions,repmat({{["GoSignal" "StartHold"]}},1,length(conditions)));
phaseAlign = containers.Map(conditions,cellfun(@(c) num2cell(string(c)),repmat({{"GoSignal","StartReach","StartHold","StartWithdraw"}},...
    1,length(conditions)),'UniformOutput',false));
phaseWin = repmat({{[0, winSz],[-winSz*(3/4),winSz*(1/4)],[-winSz*(5/4), -winSz*(1/4)],[-winSz*(1/4),winSz*(3/4)]}},1,length(conditions));
phaseWin{end}{3} = [-winSz/2 0];
smoothKernel = (winSz-.05)/params.binSize;
savePath = "S:\Lab\ngc14\Working\Revisions\Decoding\PSTH\";
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
%%
goSegs = cellfun(@(c,p) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(p,"GoSignal"))-.5,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,params.condSegMap.values,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
normBaseline = cellfun(@(cc) vertcat(cc{:}),cellfun(@(c) cellfun(@(n) num2cell(n,2), c,'UniformOutput',false),normBaseline,'UniformOutput',false),'UniformOutput',false);
normBaseline = cellfun(@(d) horzcat(d{:}), num2cell(horzcat(normBaseline{:}),2),'UniformOutput',false);
clear goSegs
%%
siteTrialPSTHS = cellfun(@(cp) cellfun(@(s)cellfun(@(r) squeeze(num2cell(r,[1,2]))',s,'UniformOutput',false)',cellfun(@(p)...
    num2cell(permute(permute(p,[1 3 2]),[1 3 2]),[2,3]),vertcat(cp(:)),'UniformOutput',false),'UniformOutput',false),siteTrialPSTHS,'Uniformoutput', false);
siteTrialPSTHS = cellfun(@(n) [n{:}]' ,siteTrialPSTHS, 'UniformOutput',false);
siteTrialPSTHS = cellfun(@(c) cellfun(@(a) vertcat(a{:}),c,'UniformOutput',false),siteTrialPSTHS,'UniformOutput',false);
siteTrialPSTHS = cellfun(@(s) vertcat(s{:}), num2cell([siteTrialPSTHS{:}],2),'UniformOutput',false);
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
subpopulations =  cellfun(@(b) contains(unitSomatotopy,["Arm","Hand"]) & b(:,1)>0, fUnits,'UniformOutput',false); %cellfun(@(s) contains(unitSomatotopy,strsplit(s)) somatotopicLabs, 'UniformOutput',false);
subpopulations(end+1:end+length(somatotopicLabs)) = cellfun(@(s) strcmp(unitSomatotopy,s),somatotopicLabs,'UniformOutput',false);
%% trial by trial condition labels organized by unit
% trialCondTable = getTrialPhaseTable(phaseFR(1:length(conditions)),phaseNames,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false),sdm);
% trialConds = arrayfun(@(t) table2cell(trialConds(trialCondTable.Unit==t & ~all(isnan(trialCondTable{:,contains(trialCondTable.Properties.VariableNames,phaseNames)}),2),:)),unique(trialCondTable.Unit),'UniformOutput',false);
trialConds = cellfun(@(s) num2cell(cell2mat(cellfun(@(c,t) repmat(c(1),size(t{1},1),1),params.condAbbrev.values,s,'UniformOutput',false)')),num2cell([siteSegs{:}],2),'UniformOutput',false);
trialConds = mapSites2Units(cellfun(@length,siteChannels),trialConds');
trialLabs = trialConds(goodInds);
%% get trial spike counts
xl = [winT(1):(winSz/2):winT(end),winT(end)];
phaseInds = cellfun(@(p,cc) cellfun(@(w) arrayfun(@(b)find(strcmp(string(cc),b)),w),p,'UniformOutput',false),phaseAlign.values,...
    values(params.condSegMap,params.condSegMap.keys()),'UniformOutput',false);
phaseEnds = cellfun(@(a,pa,pw) cellfun(@(ps) cellfun(@(ap)cellfun(@(p,n) ...
    ap(:,p)+n,pa,pw,'UniformOutput', false),ps,'UniformOutput',false),a, 'UniformOutput', false),siteSegs,phaseInds,phaseWin,'UniformOutput', false);
phaseEnds = cellfun(@(a) cellfun(@(ps) cellfun(@(ap) cellfun(@(p) p(:,3)+ap,ps,'UniformOutput',false),...
    num2cell(xl),'UniformOutput',false),a,'UniformOutput',false),siteSegs,'Uniformoutput',false);
spCounts = cellfun(@(p,s) cellfun(@(tt,a)cellfun(@(h) cellfun(@(ap) cellfun(@(tp,hp) arrayfun(@(t) ... sum(hp>=tp(1) & hp<=tp(end)),...
    sum(hp>=t-(winSz/2) & hp<=t+(winSz/2)),tp),num2cell(ap(~all(isnan(ap),2),:),2),h(1,~all(isnan(ap),2))','UniformOutput',false),...
    [a{:}],'UniformOutput',false),tt,'UniformOutput',false),p,s,'UniformOutput',false),rawSpikes,phaseEnds,'UniformOutput',false);
spCounts = cellfun(@(cb) [cb{:}]' ,spCounts, 'UniformOutput',false);
spCounts = cellfun(@(c) cellfun(@(a) [a{:}], c,'UniformOutput',false),spCounts,'UniformOutput',false);
spCounts = cellfun(@(s) cell2mat(vertcat(s{:})), num2cell([spCounts{:}],2), 'UniformOutput',false);
allUnitsTrials = spCounts(goodInds);
clear spCounts
%% get trial PSTHs
xl = winT(1):params.binSize:winT(end);
allUnitsTrials = cellfun(@(a) sqrt(conv2(resize(a(:,findBins(winT(1),params.bins):findBins(winT(end),params.bins)),[size(a,1),smoothKernel+...
    range(findBins(winT,params.bins))]),gausswin(smoothKernel)'./sum(gausswin(smoothKernel)),'valid')),siteTrialPSTHS(goodInds), 'UniformOutput',false);
%cellfun(@(u) cell2mat(cellfun(@(t)accumarray(groupBins',t,[],@mean)',num2cell(u(:,findBins(xl(1),params.bins):findBins(xl(end),params.bins)),2),'UniformOutput',false)),normPSTH(goodInds),'UniformOutput',false);
% trialInds = cellfun(@(l) ~ismember(1:length(l),cell2mat(arrayfun(@(a) randsample(find(strcmp(l,a)),5),allLabs,'UniformOutput',false)))',trialLabs,'UniformOutput',false);
% testingSet = cellfun(@(t,i) t(~i,:), trialUnits,trialInds,'UniformOutput',false);
% testingLabs =  cellfun(@(t,i) t(~i,:), trialLabs, trialInds, 'UniformOutput',false);
%% decoder setup
testUnits = [1, 5, 10, 20, 35, 50, 100];
num_repeated_labels = 3;
num_cv_splits = 10;
nAvg = 2;
numRuns = 50;
fp = {};
cl = max_correlation_coefficient_CL;

goodUnits = find_sites_with_k_label_repetitions(trialLabs,num_repeated_labels*num_cv_splits,arrayfun(@(a) a{1}(1),conditions,'UniformOutput',false));
trialUnits = cellfun(@(n,b) vertcat(n{:}), num2cell(allUnitsTrials(goodUnits,:),2), 'UniformOutput',false);
trialInds =  cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialLabs(goodUnits), trialInds, 'UniformOutput',false);

dsr = avg_DS(trainingSet,trainingLabs,num_cv_splits,nAvg);
binning_parameters = struct('end_time', size(allUnitsTrials{1},2),'start_time', 1,'bin_width',1);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.time_periods_to_get_data_from = num2cell(binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time);
dsr.num_times_to_repeat_each_label_per_cv_split = num_repeated_labels;
dsr.sample_sites_with_replacement = 0;
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
dsr.randomly_shuffle_labels_before_running = 0;
dsr.nAvg = nAvg;
timepoints = length(dsr.time_periods_to_get_data_from);
%%
somaUnits = repmat({repmat({NaN(num_cv_splits,timepoints)},length(fTypes),length(testUnits))},numRuns,1);
uUnits = repmat({cell(length(fTypes),length(testUnits))},numRuns,1);
hbar = parforProgress(numRuns);
rng('shuffle');
parfor iter = 1:numRuns
    rs = RandStream.create('threefry4x64_20','NumStreams',numRuns,'StreamIndices',iter,'Seed','shuffle');
    RandStream.setGlobalStream(rs);
    for f = 1:length(fTypes)
        for n = 1:length(testUnits)
            popUnits = find(ismember(find(goodUnits),find(subpopulations{f})));
            units =  popUnits(randperm(length(popUnits),min(length(popUnits),testUnits(n))));
            [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
            iterAcc = NaN(num_cv_splits,length(timepoints));
            for iCV = 1:num_cv_splits
                for iTrainingInterval = 1:timepoints
                    XTrF = all_XTr{iTrainingInterval};%cell2mat(reshape(cellfun(@(i)i{iCV}(units,:),all_XTr,'UniformOutput',false),1,1,[]));%
                    for iTestingInterval = iTrainingInterval
                        XTsF = all_XTrt{iTestingInterval};%cell2mat(reshape(cellfun(@(i)i{iCV}(units,:),all_XTrt,'UniformOutput',false),1,1,[]))%
                        if(isempty(fp))
                            tr = XTrF{iCV}(units,:);
                            XTst = XTsF{iCV}(units,:);
                        else
                            [~,tr] =  fp.set_properties_with_training_data(XTrF{iCV}(units,:));
                            [~,XTst] =  fp.set_properties_with_training_data(XTsF{iCV}(units,:));
                            %[loadings,scores,~,~,~,pcaMean] = pca(reshape(permute(tr(units,:,:),[2 1 3]),size(tr,2),[]),'Algorithm','svd','NumComponents',10);testScores = (reshape(permute(XTst(units,:,:),[2 1 3]),size(XTst,2),[]) - pcaMean) * loadings(:, 1);
                        end
                        clT= cl.train(tr, all_YTr);
                        [ia,~] = clT.test(XTst);
                        iterAcc(iCV,iTestingInterval)= sum(ia-all_Ytrt==0)/length(all_Ytrt);
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
unitAccPhase = vertcat(unitAccPhase{:});
%uUnits = cellfun(@(p) cellfun(@(r) resize(r,[1,max(testUnits)],'FillValue',NaN),p,'UniformOutput',false),uUnits,'UniformOutput',false);
%%
accLine = 90; trCl = [0 .7 0; 1 .5 0; 0 .2 .8; .4 .4 .4; .6 0 .2;];
em=cell2mat(cellfun(@(a) a([2,3,6]),cellfun(@(c) mean(cell2mat(cellfun(@(a) cell2mat(cellfun(@(n) mean(n,1,'omitnan'),a,'UniformOutput',false)),c,'UniformOutput',false)),1,'omitnan'),siteSegs,'UniformOutput',false)','UniformOutput',false));
em=[mean(em(:,1:2),1,'omitnan'),min(em(:,end)),max(em(:,end))];
[~,ei] = arrayfun(@(a) min(abs(xl-a)),em);
phaseNames = arrayfun(@(s) string(num2str(round(s,2))),xl);
figure(); hold on;
trtsPhase = 100.*squeeze(mean(unitAccPhase(:,1:end,:,strcmp(fTypes,"Task")),2,'omitnan'));
errorbar(testUnits,squeeze(mean(trtsPhase,1,'omitnan')),squeeze(std(trtsPhase,0,1)),'LineStyle','-','LineWidth',2);
plot(get(gca,'XLim'),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]); legend("FR");
saveFigures(gcf,savePath,"_"+num2str(numRuns),[]);
for t = 1:length(testUnits)
    nUnits=testUnits(t);
    accTUnit = cellfun(@(s) permute(vertcat(s{:,testUnits==nUnits}),[3 2 1]),somaUnits,'UniformOutput',false);
    accTUnit =  cell2mat(accTUnit);
    accTable = array2table(100.*squeeze(mean(unitAccPhase(:,1:end,testUnits==nUnits,:),2,'omitnan')),'VariableNames',fTypes);
    accTable = [accTable,array2table(100.*unitAccPhase(:,1:end,testUnits==nUnits,strcmp(fTypes,"Task")),'VariableNames',"Time-"+phaseNames)];

    if(length(phaseNames)>1)
        figure(); nt=tiledlayout(1,3);
        for d = 0:2
            nexttile(nt); hold on; pll = cellfun(@(p,c) shadedErrorBar(xl,mean(p,1,'omitnan'),.5*std(p,0,1,'omitnan'),'lineProps',...
                {'Color',trCl(c-(3*d),:),'LineWidth',2},'patchSaturation',.1),arrayfun(@(a) accTUnit(:,:,a),(3*d)+(1:3+(d==2)),'UniformOutput',false),num2cell((3*d)+(1:3+(d==2))));
            arrayfun(@(a,s) plot([a a], [0,1], s),em,["k--","k--","b--","m--"]);
            legend([pll.mainLine],fTypes((3*d)+(1:3+(d==2)))); ylim([0 1]);xlim([min(xl),max(xl)]);
        end
        saveFigures(gcf,savePath,"Temporal_"+num2str(numRuns)+"_"+num2str(nUnits),[]);
    end
    writetable(accTable,savePath+"Decoding_"+num2str(nUnits)+"_"+num2str(numRuns),'FileType','spreadsheet','UseExcel',true);

    figure(); hold on;
    bx=boxchart(accTable,"Time-"+phaseNames(1:end),'Notch','on','MarkerStyle','none');
    xticks(1:5:length(phaseNames));xticklabels(string(phaseNames(1:5:end)));newXLim = [-.1 .1]+double(get(gca,'XLim'));
    plot(newXLim,[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);xlim(newXLim);
    saveFigures(gcf,savePath,"Phases_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,["Arm" "Hand" "Trunk" "Face"],'Notch','on','MarkerStyle','none');
    xticks(1:length(somatotopicLabs));xticklabels(["Arm" "Hand" "Trunk" "Face"]);
    plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
    saveFigures(gcf,savePath,"Somatotopy_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,fTypes([1 3 2]),'Notch','on','MarkerStyle','none');
    xticks(1:3);xticklabels(fTypes([1 3 2]));
    plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
    saveFigures(gcf,savePath,"Types_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,["Shallow","Deep"],'Notch','on','MarkerStyle','none');
    xticks(1:2);xticklabels(["Shallow","Deep"]);
    plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
    saveFigures(gcf,savePath,"Laminar_"+num2str(numRuns)+"_"+num2str(nUnits),[]);
end
% allBxs = vertcat(bx.boxplotGroup(1:end).Children);
% set(allBxs(contains(string({allBxs.Tag}),"Whisker")),'LineStyle','-');
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