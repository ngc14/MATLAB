%% load data
winSz = 0.2;
fTypes = ["Reach","Grasp","Both","Shallow","Deep","Task"];
phaseNames = ["GoSignal","StartReach","StartHold","StartWithdraw"];
winT = [-.5 .5];
num_repeated_labels = 3;
num_cv_splits = 10;
nAvg = 2;
accLine = 80;
trCl = [0 .7 0; 1 .5 0; 0 .2 .8; .4 .4 .4; .6 0 .2;];
savePath = "S:\Lab\ngc14\Working\Revisions\Decoding\MaxCorr\Fr_Sqrt\";

params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5,...
    containers.Map(["Extra Small Sphere","Large Sphere", "Photocell"],repmat({"StartReach"},1,3)));
xl = winT(1):params.binSize:winT(end);
timepoints = arrayfun(@(s) string(num2str(round(s,2))),xl);
smoothKernel = fix((winSz-.05)/params.binSize);

taskAlign = containers.Map(params.condNames,repmat({{["GoSignal" "StartHold"]}},1,length(params.condNames)));
phaseAlign = containers.Map(params.condNames,cellfun(@(c) num2cell(cell2mat(c)),...
    repmat({{phaseNames}},1,length(params.condNames)),'UniformOutput',false));
phaseWin = repmat({{[0, winSz],[-winSz*(3/4),winSz*(1/4)],[-winSz*(5/4), -winSz*(1/4)],[-winSz*(1/4),winSz*(3/4)]}},1,length(params.condNames));
phaseWin{strcmp(params.condNames,"Photocell")}{contains(phaseNames,"Hold")} = [-winSz/2 0];
%%
[siteDateMap,siteSegs,siteTrialPSTHS,rawSpikes,siteChannels,chMaps,~,~] = ...
    getAllSessions(params,"Single","M1","");
somatotopicLabs = unique(cell2mat([siteDateMap.SiteRep]'));
fTypes = [fTypes,somatotopicLabs];
simpRep =  cellfun(@(r,t) r(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh,'UniformOutput', true)';
mappedChannels =  cell2mat(cellfun(@(ch,l) ch{end}(l(~isnan(l)))', chMaps,siteChannels, 'Uniformoutput', false)');
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep));
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[winSz, 0]}},1,length(params.condNames)),siteSegs,siteTrialPSTHS,false,true);
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(params.condNames)),taskFR(1:length(params.condNames)),'UniformOutput',false);
taskUnits = cellfun(@cell2mat, taskUnits,'UniformOutput',false);
[~,phaseFR] = calculatePhases(params,phaseAlign,phaseWin,siteSegs,siteTrialPSTHS,false,true);
[rgVals,rgInds] = cellfun(@(crg) cellfun(@(a) cellfun(@(rg)ttestTrials({rg(contains(phaseNames,"Reach"))},...
    {rg(contains(phaseNames,"Hold"))},1,true,0.05),a,'UniformOutput',false),crg,'UniformOutput',false),phaseFR, 'UniformOutput',false);
rgInds = cellfun(@(v) cell2mat(vertcat(v{:})), rgInds, 'UniformOutput',false);
rgVals = cellfun(@(v) vertcat(v{:}), rgVals, 'UniformOutput',false);
rgVals = cellfun(@(c) cell2mat(cellfun(@(d) mean(cell2mat(d),2,'omitnan'),c,'UniformOutput',false)), rgVals,'UniformOutput',false);
clear phaseFR rawSpikes taskFR taskBaseline
%% optional baseline normalization
goSegs = cellfun(@(c,p) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(p,"GoSignal"))-4,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,params.condSegMap.values,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(3/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
normBaseline = cellfun(@(cc) vertcat(cc{:}),cellfun(@(c) cellfun(@(n) num2cell(n,2), c,'UniformOutput',false),normBaseline,'UniformOutput',false),'UniformOutput',false);
clear goSegs
%%
siteTrialPSTHS = cellfun(@(cp) cellfun(@(s)cellfun(@(r) squeeze(num2cell(r,[1,2]))',s,'UniformOutput',false)',cellfun(@(p)...
    num2cell(permute(permute(p,[1 3 2]),[1 3 2]),[2,3]),vertcat(cp(:)),'UniformOutput',false),'UniformOutput',false),siteTrialPSTHS,'Uniformoutput', false);
siteTrialPSTHS = cellfun(@(n) [n{:}]' ,siteTrialPSTHS, 'UniformOutput',false);
siteTrialPSTHS = cellfun(@(c) cellfun(@(a) vertcat(a{:}),c,'UniformOutput',false),siteTrialPSTHS,'UniformOutput',false);
%avgUnits = cellfun(@(cs,cb) cellfun(@(s,b) mean(s,1,'omitnan')./mean(max(1,b),2,'omitnan'),cs,cb,'UniformOutput',false),siteTrialPSTHS,normBaseline,'UniformOutput',false)
siteTrialPSTHS = cellfun(@(s) vertcat(s{:}), num2cell([siteTrialPSTHS{:}],2),'UniformOutput',false);
normBaseline = cellfun(@(d) horzcat(d{:}), num2cell(horzcat(normBaseline{:}),2),'UniformOutput',false);
%% only use units that are task modulated according to the task phase window
rUnits = cellfun(@(rg,r,t) rg==1 & r < 1 & t==1, rgInds, rgVals,taskUnits, 'UniformOutput',false);
gUnits = cellfun(@(rg,r,t) rg==1 & r >= 1 & t==1, rgInds, rgVals,taskUnits, 'UniformOutput',false);
bothUnits = cellfun(@(rg,t) ~rg & t==1, rgInds, taskUnits, 'UniformOutput',false);
goodUnits = any(cell2mat(taskUnits),2)';
unitsConds = [sum(cell2mat(rUnits(contains(params.condNames,"Sphere"))),2),sum(cell2mat(gUnits(contains(params.condNames,"Sphere"))),2),sum(cell2mat(bothUnits(contains(params.condNames,"Sphere"))),2)];
[unitTypeV,unitTypeT] = max(unitsConds,[],2);
unitTypeT(unitTypeV==0) = 0;
isMaxTie = unitsConds==unitTypeV;
unitCondAssign = find(sum(isMaxTie,2)>1 & unitTypeV~=0);
unitTypeT(unitCondAssign(sum(isMaxTie(unitCondAssign,:),2)==3)) = 3;
unitTypeT(unitCondAssign(sum(isMaxTie(unitCondAssign,:),2)==2 & ~isMaxTie(unitCondAssign,end))) = 3;
unitTypeT = rUnits{1}.*1 + gUnits{1}.*2 + bothUnits{1}.*3;
fUnits = {unitTypeT==1, unitTypeT==2, unitTypeT==3,mappedChannels<=16,mappedChannels>16,goodUnits};
clear isMaxTie unitsConds unitTypeT unitTypeV
%% get trial PSTHs organized by unit
trialUnits = cellfun(@(a) sqrt(conv2(resize(a(:,findBins(winT(1),params.bins):findBins(winT(end),params.bins)),[size(a,1),smoothKernel+...
    range(findBins(winT,params.bins))]),gausswin(smoothKernel)'./sum(gausswin(smoothKernel)),'valid')),siteTrialPSTHS(goodUnits), 'UniformOutput',false);
trialUnits = cellfun(@(n) vertcat(n{:}), num2cell(trialUnits,2), 'UniformOutput',false);
trialConds = mapSites2Units(cellfun(@length,siteChannels),cellfun(@(s) num2cell(cell2mat(cellfun(@(c,t) repmat(c(1),...
    size(t{1},1),1),params.condAbbrev.values,s,'UniformOutput',false)')),num2cell([siteSegs{:}],2),'UniformOutput',false)');

trialInds =  cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialConds(goodUnits), trialInds, 'UniformOutput',false);

subpopulations =  cellfun(@(b) contains(unitSomatotopy(goodUnits),["Arm","Hand"]) & b(goodUnits), fUnits(1:end-1),'UniformOutput',false);
subpopulations(end+1) = cellfun(@(b) b(goodUnits), fUnits(end),'UniformOutput',false);
subpopulations(end+1:end+length(somatotopicLabs)) = cellfun(@(s) strcmp(unitSomatotopy(goodUnits),s),somatotopicLabs,'UniformOutput',false);
%subpopulations = arrayfun(@(s) cellfun(@(f) f(goodUnits) & strcmp(unitSomatotopy(goodUnits),s),fUnits,'UniformOutput',false),["Arm","Hand"],'UniformOutput',false);
%subpopulations = [subpopulations{:}];
%fTypes = reshape(fTypes(1:5)'+["-Arm","-Hand"],[],1)';

goodUnitsTrials = find_sites_with_k_label_repetitions(trainingLabs,num_repeated_labels*num_cv_splits,arrayfun(@(a) a{1}(1),params.condNames,'UniformOutput',false))';
subpopulations = cellfun(@(s) ismember(goodUnitsTrials,find(s)), subpopulations,'UniformOutput',false);
dsr = avg_DS(trainingSet(goodUnitsTrials),trainingLabs(goodUnitsTrials),num_cv_splits,nAvg);
clear trialUnits trialInds trialConds trainingSet trainingLabs
%% decoder setup
testUnits = [1,10,20,30,40,50,60,70,80,90,100];
testUnits = 60;
numRuns = 1000;
fp = {};
cl = max_correlation_coefficient_CL;

binning_parameters = struct('end_time', length(timepoints),'start_time', 1,'bin_width',1);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.time_periods_to_get_data_from = num2cell(binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time);
dsr.num_times_to_repeat_each_label_per_cv_split = num_repeated_labels;
dsr.sample_sites_with_replacement = 0;
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
dsr.randomly_shuffle_labels_before_running = 0;
dsr.nAvg = nAvg;
%% get subpopulation decoding
somaUnits= repmat({repmat({NaN(num_cv_splits,length(timepoints))},length(fTypes),length(testUnits))},numRuns,1);
uUnits = repmat({cell(length(fTypes),length(testUnits))},numRuns,1);
hbar = parforProgress(numRuns);
rng('default');
parfor iter = 1:numRuns
    rs = RandStream.create('threefry4x64_20','NumStreams',numRuns,'StreamIndices',iter,'Seed','shuffle');
    RandStream.setGlobalStream(rs);
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    for n = 1:length(testUnits)
        for f = 1:length(fTypes)
            popUnits = find(subpopulations{f});
            units =  popUnits(randperm(length(popUnits),min(length(popUnits),testUnits(n))));
            iterAcc= NaN(num_cv_splits,length(timepoints));
            for iCV = 1:num_cv_splits
                for iTrainingInterval = 1:length(timepoints)
                    XTrF = all_XTr{iTrainingInterval};
                    %for iTestingInterval = iTrainingInterval
                    XTsF = all_XTrt{iTrainingInterval};
                    if(isempty(fp))
                        tr = XTrF{iCV}(units,:);
                        XTst = XTsF{iCV}(units,:);
                    else
                        [~,tr] =  fp.set_properties_with_training_data(XTrF{iCV}(units,:));
                        [~,XTst] =  fp.set_properties_with_training_data(XTsF{iCV}(units,:));
                    end
                    clT= cl.train(tr, all_YTr);
                    [ia,~] = clT.test(XTst);
                    iterAcc(iCV,iTrainingInterval)= sum(ia-all_Ytrt==0)/length(all_Ytrt);
                    %end
                end
            end
            somaUnits{iter}{f,n} = mean(iterAcc,1,'omitnan');
            uUnits{iter}{f,n} = units;
        end
    end
    send(hbar, iter);
end
unitAccPhase = cellfun(@(p) cell2mat(permute(cellfun(@(m) squeeze(mean(m,1,'omitnan')),p,'Uniformoutput',false),...
    (length(size(p{1}))+ length(size(p))):-1:1)),somaUnits,'UniformOutput',false);
unitAccPhase = vertcat(unitAccPhase{:});
uUnits = cellfun(@(p) cellfun(@(r) resize(r,[1,max(testUnits)],'FillValue',NaN),p,'UniformOutput',false),uUnits,'UniformOutput',false);
%% permutation testing
permutations = 10000;
allAcc = squeeze(mean(unitAccPhase(:,:,testUnits==nUnits,:),2,'omitnan'));
shuffleGroups = {fTypes(1:3), fTypes(4:5), fTypes(6:8), fTypes(9:10)};
shuffleGroups = cellfun(@(i) find(contains(fTypes,i)),shuffleGroups,'UniformOutput',false);

F_perm = cell(permutations,length(shuffleGroups));
diffRatio = NaN(length(shuffleGroups),max(cellfun(@(s) nchoosek(length(s),2),shuffleGroups)));
hbar = parforProgress(permutations);
parfor p = 1:(permutations)
    rt = RandStream.create('threefry4x64_20','NumStreams',permutations,'StreamIndices',p,'Seed','shuffle');
    RandStream.setGlobalStream(rt);
    shuffleUnits = cell(1,length(shuffleGroups));
    for s = 1:length(shuffleGroups)
        shuffledAcc = allAcc(:,shuffleGroups{s});
        shuffledAcc = shuffledAcc(:);
        shuffleUnits{s} = shuffledAcc(randperm(length(shuffledAcc)));
    end
    F_perm(p,:) = shuffleUnits;
    send(hbar,p);
end
for s = 1:length(shuffleGroups)
    groupIDs = repelem(fTypes(shuffleGroups{s}),size(allAcc,1));
    fo = anova(100.*allAcc(:,shuffleGroups{s}));
    fo =fo.stats.F(1);
    fs = NaN(length(F_perm(:,s)),1);
    hbar = parforProgress(length(F_perm));
    parfor p = 1:length(F_perm)
        an = anova(groupIDs,100.*F_perm{p,s});
        fs(p) = an.stats.F(1);
        send(hbar,p);
    end
    if((1+sum(fs>=fo))/(length(F_perm)+1)<0.001)
        pairs = nchoosek(1:length(shuffleGroups{s}),2);
        diffO = 100.*mean(allAcc(:,shuffleGroups{s}(pairs(:,1)))-allAcc(:,shuffleGroups{s}(pairs(:,2))),1,'omitnan');
        for p = 1:length(diffO)
            allDiffs = NaN(length(F_perm),1);
            currPair = ismember(groupIDs,fTypes(shuffleGroups{s}(pairs(p,:))));
            parfor f = 1:length(F_perm)
                allDiffs(f)= diff(groupsummary(100.*F_perm{f,s}(currPair),groupIDs(currPair)',"mean"));
            end
            diffRatio(s,p) = (1+sum(abs(allDiffs)>=abs(diffO(p))))/(length(F_perm)+1);
        end
    end
end
%% Figures
if(length(testUnits)>1)
    trtsPhase = 100.*squeeze(mean(unitAccPhase(:,1:end,:,:),2,'omitnan'));
    somaOrder = [find(contains(fTypes,"Arm")), find(contains(fTypes,"Hand"))];
    testTicks = linspace(1,550,length(testUnits));
    nAx = 2; nLs = length(fTypes)/nAx; nL = mod(length(fTypes),nAx);
    figure(); colororder(trCl); nt = tiledlayout(1,nAx);
    arrayfun(@(d) errorbar(nexttile(nt,d),squeeze(mean(trtsPhase(:,:,somaOrder((nLs*(d-1))+(1:(nLs+nL)))),1,'omitnan')),squeeze(std(...
        trtsPhase(:,:,somaOrder((nLs*(d-1))+(1:(nLs+nL)))),0,1)),'LineWidth',2),1:nAx);
    arrayfun(@(d) legend(nexttile(nt,d),fTypes(somaOrder((nLs*(d-1))+(1:(nLs+nL)))),'location','southeast'),1:nAx);
    arrayfun(@(d) set(nexttile(nt,d),'ylim',[30 100],'xtick',1:length(testTicks),'xticklabels',testUnits),1:nAx);
    saveFigures(gcf,savePath,"Errorbars-"+num2str(numRuns),[]);

    nAx = 3; nLs = floor(length(fTypes)/nAx); nL = mod(length(fTypes),nAx);
    figure(); colororder(trCl); nt=tiledlayout(1,nAx);
    arrayfun(@(d) boxchart(nexttile(nt,d),reshape(repmat(testTicks,size(trtsPhase,1),1,nLs+(nL*(d==nAx))),1,[]),reshape(...
        trtsPhase(:,:,ismember(fTypes,fTypes((nLs*(d-1))+(1:(nLs+(nL*(d==nAx))))))),1,[]),'GroupByColor',reshape(permute(...
        repmat(1:(nLs+(nL*(d==nAx))),size(trtsPhase,1),1,length(testTicks)),[1 3 2]),1,[]),'Notch','on','MarkerStyle','none','BoxWidth',30,'BoxFaceAlpha',1),1:nAx);
    arrayfun(@(d) legend(nexttile(nt,d),fTypes((nLs*(d-1))+(1:(nLs+(nL*(d==nAx))))),'location','southeast'),1:nAx);
    arrayfun(@(d) set(nexttile(nt,d),'ylim',[30 100],'xtick',testTicks,'xticklabels',testUnits),1:nAx);
    saveFigures(gcf,savePath,"Boxplots-"+num2str(numRuns),[]);

    figure(); hold on;
    [~,maxT] = max(100.*unitAccPhase(:,:,:,end-(length(somatotopicLabs)-1):end)>=accLine,[],2);
    maxT(maxT==1) = length(xl);
    figure(); errorbar(testUnits,squeeze(mean(xl(maxT),1,'omitnan')),squeeze(std(xl(maxT),0,1)));
    legend(somatotopicLabs);
    saveFigures(gcf,savePath,"Latency_Somatotopy-"+num2str(numRuns),[]);
end

for t = 1:length(testUnits)
    nUnits=testUnits(t);

    accTUnit = cell2mat(cellfun(@(s) permute(vertcat(s{:,testUnits==nUnits}),[3 2 1]),somaUnits,'UniformOutput',false));
    accTable = array2table(100.*squeeze(mean(unitAccPhase(:,:,testUnits==nUnits,:),2,'omitnan')),'VariableNames',fTypes);
    if(any(strcmp(fTypes,"Task")))
        accTable = [accTable,array2table(100.*unitAccPhase(:,:,testUnits==nUnits,strcmp(fTypes,"Task")),'VariableNames',"Time-"+timepoints)];
    end
    writetable(accTable,savePath+"Decoding_"+num2str(nUnits)+"_"+num2str(numRuns),'FileType','spreadsheet','UseExcel',true);

    if(length(timepoints)>1)
        em=cell2mat(cellfun(@(a) a([2,3,6]),cellfun(@(c) mean(cell2mat(cellfun(@(a) cell2mat(cellfun(@(n) mean(n,1,'omitnan'),a,'UniformOutput',false)),c,'UniformOutput',false)),1,'omitnan'),siteSegs,'UniformOutput',false)','UniformOutput',false));
        em=[mean(em(:,1:2),1,'omitnan'),min(em(:,end)),max(em(:,end))];
        figure(); nt=tiledlayout(1,nAx);
        accAx = arrayfun(@(d) arrayfun(@(a) accTUnit(:,:,a),(nLs*(d-1))+(1:(nLs+(nL*(d==nAx)))),'UniformOutput',false),1:nAx,'UniformOutput',false);
        d = cellfun(@(a,d) cellfun(@(p,c) shadedErrorBar(xl,mean(p,1,'omitnan'),std(p,0,1),'lineProps',...
            {'Color',trCl(c,:),'LineWidth',2},'patchSaturation',.1,'ax',nexttile(nt,d)),a,num2cell(1:length(a))),accAx,num2cell(1:nAx),'UniformOutput',false);
        cellfun(@(x,n) legend(nexttile(nt,n),[x.mainLine],fTypes((nLs*(n-1))+(1:(nLs+(nL*(n==nAx))))),'location','southeast','AutoUpdate','off'),d,num2cell(1:nAx));
        arrayfun(@(d) hold(nexttile(nt,d), 'on'), 1:nAx);
        arrayfun(@(d) set(nexttile(nt,d),'ylim',[0 1],'xlim',[min(xl),max(xl)]),1:nAx);
        arrayfun(@(d) arrayfun(@(a,s) plot(nexttile(nt,d),[a a], [0,1], s),em,["k--","k--","b--","m--"]),1:nAx);
        saveFigures(gcf,savePath,"Temporal_"+num2str(numRuns)+"_"+num2str(nUnits),[]);
    end
    figure(); hold on;
    boxchart(accTable,fTypes,"Notch",'on','MarkerStyle','none');
    xticklabels(fTypes);
    saveFigures(gcf,savePath,"Groups_"+num2str(numRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,"Time-"+timepoints(1:end),'Notch','on','MarkerStyle','none');
    xticks(1:5:length(timepoints));xticklabels(string(timepoints(1:5:end)));newXLim = [-.1 .1]+double(get(gca,'XLim'));
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
%allBxs= bx.boxplotGroup(1:end).Children; set(allBxs(contains(string({allBxs.Tag}),"Whisker")));
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