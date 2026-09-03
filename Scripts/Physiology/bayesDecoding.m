%% load data
phaseWindowSz = 0.2;
fTypes = ["Reach","Grasp","Both","Shallow","Deep","Task"];
phaseNames = ["GoSignal","StartReach","StartHold","StartWithdraw"];
eventAlign = {"StartReach"};
winAroundEvent = [-.5 .5];
nRepeatedCondTests = 3;
nCVFolds = 10;
nAvg = 2;
nRuns = 40;
permutations = 1000;
nInputUnits = 60; %[1,10,20,30,40,50,60,70,80,90,100];
pVal = .01;
fp = {};
classifier = max_correlation_coefficient_CL;
PERMUTATIONTEST=true;
savePath = "S:\Lab\ngc14\Working\Revisions\Decoding\New_Figure\";

params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5,...
    containers.Map(["Extra Small Sphere","Large Sphere", "Photocell"],repmat(eventAlign,1,3)));
xl = winAroundEvent(1):params.binSize:winAroundEvent(end);
timepoints = arrayfun(@(s) string(num2str(round(s,2))),xl);
smoothKernel = fix((phaseWindowSz-.05)/params.binSize);

taskAlign = containers.Map(params.condNames,repmat({{["GoSignal" "StartHold"]}},1,length(params.condNames)));
phaseAlign = containers.Map(params.condNames,cellfun(@(c) num2cell(cell2mat(c)),...
    repmat({{phaseNames}},1,length(params.condNames)),'UniformOutput',false));
phaseWin = repmat({{[0, phaseWindowSz],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)],[-phaseWindowSz*(5/4), -phaseWindowSz*(1/4)],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)]}},1,length(params.condNames));
phaseWin{strcmp(params.condNames,"Photocell")}{contains(phaseNames,"Hold")} = [-phaseWindowSz/2 0];
%%
[siteDateMap,siteSegs,siteTrialPSTHS,~,siteChannels,chMaps,~,~] = ...
    getAllSessions(params,"Single","M1","");
somatotopicLabs = unique(cell2mat([siteDateMap.SiteRep]'));
simpRep =  cellfun(@(r,t) r(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh,'UniformOutput', true)';
mappedChannels =  cell2mat(cellfun(@(ch,l) ch{end}(l(~isnan(l)))', chMaps,siteChannels, 'Uniformoutput', false)');
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep));
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[phaseWindowSz, 0]}},1,length(params.condNames)),siteSegs,siteTrialPSTHS,false,true);
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(params.condNames)),taskFR(1:length(params.condNames)),'UniformOutput',false);
taskUnits = cellfun(@cell2mat, taskUnits,'UniformOutput',false);
[~,phaseFR] = calculatePhases(params,phaseAlign,phaseWin,siteSegs,siteTrialPSTHS,false,true);
[rgVals,rgInds] = cellfun(@(crg) cellfun(@(a) cellfun(@(rg)ttestTrials({rg(contains(phaseNames,"Reach"))},...
    {rg(contains(phaseNames,"Hold"))},1,true,0.05),a,'UniformOutput',false),crg,'UniformOutput',false),phaseFR, 'UniformOutput',false);
rgInds = cellfun(@(v) cell2mat(vertcat(v{:})), rgInds, 'UniformOutput',false);
rgVals = cellfun(@(v) vertcat(v{:}), rgVals, 'UniformOutput',false);
rgVals = cellfun(@(c) cell2mat(cellfun(@(d) mean(cell2mat(d),2,'omitnan'),c,'UniformOutput',false)), rgVals,'UniformOutput',false);
%% optional baseline normalization
goSegs = cellfun(@(c,p) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(p,"GoSignal"))-4,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,params.condSegMap.values,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(3/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
normBaseline = cellfun(@(cc) vertcat(cc{:}),cellfun(@(c) cellfun(@(n) num2cell(n,2), c,'UniformOutput',false),normBaseline,'UniformOutput',false),'UniformOutput',false);
%normBaseline = cellfun(@(d) horzcat(d{:}), num2cell(horzcat(normBaseline{:}),2),'UniformOutput',false);
clear goSegs phaseFR taskFR taskBaseline normBaseline chMaps siteDateMap simpRep
%% combined conditions of PSTHS to get trial PSTHS organized by units
siteTrialPSTHS = cellfun(@(cp) cellfun(@(s)cellfun(@(r) squeeze(num2cell(r,[1,2]))',s,'UniformOutput',false)',cellfun(@(p)...
    num2cell(permute(permute(p,[1 3 2]),[1 3 2]),[2,3]),vertcat(cp(:)),'UniformOutput',false),'UniformOutput',false),siteTrialPSTHS,'Uniformoutput', false);
siteTrialPSTHS = cellfun(@(n) [n{:}]' ,siteTrialPSTHS, 'UniformOutput',false);
siteTrialPSTHS = cellfun(@(c) cellfun(@(a) vertcat(a{:}),c,'UniformOutput',false),siteTrialPSTHS,'UniformOutput',false);
%avgUnits = cellfun(@(cs,cb) cellfun(@(s,b) mean(s,1,'omitnan')./mean(max(1,b),2,'omitnan'),cs,cb,'UniformOutput',false),siteTrialPSTHS,normBaseline,'UniformOutput',false)
siteTrialPSTHS = cellfun(@(s) vertcat(s{:}), num2cell([siteTrialPSTHS{:}],2),'UniformOutput',false);
%% only use units that are task modulated according to the task phase window
rUnits = cellfun(@(rg,r,t) rg==1 & r < 1 & t==1, rgInds, rgVals,taskUnits, 'UniformOutput',false);
gUnits = cellfun(@(rg,r,t) rg==1 & r >= 1 & t==1, rgInds, rgVals,taskUnits, 'UniformOutput',false);
bothUnits = cellfun(@(rg,t) ~rg & t==1, rgInds, taskUnits, 'UniformOutput',false);
goodUnits = any(cell2mat(taskUnits),2)';
% unitsConds = [sum(cell2mat(rUnits(contains(params.condNames,"Sphere"))),2),...
%     sum(cell2mat(gUnits(contains(params.condNames,"Sphere"))),2),...
%     sum(cell2mat(bothUnits(contains(params.condNames,"Sphere"))),2)];
% [unitTypeV,unitTypeT] = max(unitsConds,[],2);
% unitTypeT(unitTypeV==0) = 0;
% isMaxTie = unitsConds==unitTypeV;
% unitCondAssign = find(sum(isMaxTie,2)>1 & unitTypeV~=0);
% unitTypeT(unitCondAssign(sum(isMaxTie(unitCondAssign,:),2)==3)) = 3;
% unitTypeT(unitCondAssign(sum(isMaxTie(unitCondAssign,:),2)==2 & ~isMaxTie(unitCondAssign,end))) = 3;
fUnits = { rUnits{1}, gUnits{1}, bothUnits{1},mappedChannels<=16,mappedChannels>16,goodUnits};
clear rUnits gUnits bothUnits taskUnits rgInds rgVals mappedChannels 
%% organize training PSTHS and condition labels and subpopulation indicies
trialUnits = cellfun(@(a) sqrt(conv2(resize(a(:,findBins(winAroundEvent(1),params.bins):findBins(winAroundEvent(end),params.bins)),[size(a,1),smoothKernel+...
    range(findBins(winAroundEvent,params.bins))]),gausswin(smoothKernel)'./sum(gausswin(smoothKernel)),'valid')),siteTrialPSTHS(goodUnits), 'UniformOutput',false);
trialUnits = cellfun(@(n) vertcat(n{:}), num2cell(trialUnits,2), 'UniformOutput',false);
trialConds = mapSites2Units(cellfun(@length,siteChannels),cellfun(@(s) num2cell(cell2mat(cellfun(@(c,t) repmat(c(1),...
    size(t{1},1),1),params.condAbbrev.values,s,'UniformOutput',false)')),num2cell([siteSegs{:}],2),'UniformOutput',false)');
em=cell2mat(cellfun(@(a) a(contains(params.condSegMap(params.condNames(1)),phaseNames)),cellfun(@(c) mean(cell2mat(cellfun(@(a) ...
    cell2mat(cellfun(@(n) mean(n,1,'omitnan'),a,'UniformOutput',false)),c,'UniformOutput',false)),1,'omitnan'),siteSegs,'UniformOutput',false)','UniformOutput',false));

trialInds =  cellfun(@(t) ~any(isnan(t),2), trialUnits, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialUnits, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) t(i,:), trialConds(goodUnits), trialInds, 'UniformOutput',false);

subpopulations =  cellfun(@(b) contains(unitSomatotopy(goodUnits),["Arm","Hand"]) & b(goodUnits), fUnits(1:end-1),'UniformOutput',false);
subpopulations(end+1) = cellfun(@(b) b(goodUnits), fUnits(end),'UniformOutput',false);
armHandIntersectionTypes = arrayfun(@(s) cellfun(@(f) f(goodUnits) & strcmp(unitSomatotopy(goodUnits),s),...
    fUnits(1:end-1),'UniformOutput',false),["Arm","Hand"],'UniformOutput',false);
subpopulations(end+1:end+length(somatotopicLabs)) = cellfun(@(s) strcmp(unitSomatotopy(goodUnits),s),somatotopicLabs,'UniformOutput',false);
subpopulations(end+1:end+length(fUnits(1:end-1))*2) = [armHandIntersectionTypes{:}];
subGroups = [fTypes,somatotopicLabs,reshape(fTypes(1:end-1)'+["-Arm","-Hand"],1,[])];

goodUnitsTrials = find_sites_with_k_label_repetitions(trainingLabs,nRepeatedCondTests*nCVFolds,arrayfun(@(a) a{1}(1),params.condNames,'UniformOutput',false))';
subpopulations = cellfun(@(s) ismember(goodUnitsTrials,find(s)), subpopulations,'UniformOutput',false);
dsr = avg_DS(trainingSet(goodUnitsTrials),trainingLabs(goodUnitsTrials),nCVFolds,nAvg);
clear trialUnits trialInds trialConds trainingSet trainingLabs siteSegs  goodUnitsTrials goodUnits siteTrialPSTHS unitSomatotopy siteChannels
%% decoder setup
binning_parameters = struct('end_time', length(timepoints),'start_time', 1,'bin_width',1);
dsr.the_basic_DS.binned_site_info.binning_parameters = binning_parameters;
dsr.time_periods_to_get_data_from = num2cell(binning_parameters.start_time:binning_parameters.bin_width:binning_parameters.end_time);
dsr.num_times_to_repeat_each_label_per_cv_split = nRepeatedCondTests;
dsr.sample_sites_with_replacement = 0;
dsr.sites_to_use = -1;
dsr.num_resample_sites = -1;
dsr.randomly_shuffle_labels_before_running = 0;
dsr.nAvg = nAvg;
%% bootstrap iterations of classifier training/testing
somaUnits= NaN(nRuns,length(timepoints),length(subGroups),length(nInputUnits));
delete(gcp('nocreate'));
parpool('Threads',min(maxNumCompThreads,6));
hbar = parforProgress(nRuns);
rng('shuffle', 'threefry4x64_20');
[~, all_YTr, ~, all_YTrt] = dsr.get_data;
if(PERMUTATIONTEST)
    shuffleGroups = cellfun(@(i) find(ismember(subGroups,i)),...
        {somatotopicLabs,subGroups(1:3)+"-Arm",subGroups(1:3)+"-Hand", subGroups(4:5)+"-Arm", subGroups(4:5)+"-Hand"},'UniformOutput',false);
    fullSubPop = cellfun(@(i) find(any([subpopulations{i}],2)), shuffleGroups,'UniformOutput',false);
    allGroups = cell2mat(cellfun(@(n,i) repmat(i,1,length(n)), shuffleGroups,num2cell(1:length(shuffleGroups)), 'UniformOutput',false));
    subGroupInd = cellfun(@(s) [1,cumsum(sum([subpopulations{s}],1))+1], shuffleGroups, 'UniformOutput',false);
    groupPerms = cell(permutations,1);
    for p = 1:permutations
        permPopulation = cellfun(@(s) s(randperm(length(s))), fullSubPop, 'UniformOutput',false);
        gp = cellfun(@(p,s) arrayfun(@(i) p(s(i):s(i+1)-1),1:length(s)-1,'UniformOutput',false), permPopulation, subGroupInd, 'UniformOutput',false);
        groupPerms{p} = [gp{:}];
    end
    [latPerm,somaPerm] = deal(NaN(nRuns,permutations,length(allGroups),length(nInputUnits))); 
end
for iter = 1:nRuns
    [all_XTr, ~, all_XTrt, ~] = dsr.get_data;
    for n = 1:length(nInputUnits)
        for f = 1:max(length(allGroups),length(subGroups))
            popUnits = find(subpopulations{f});
            unitInds = randperm(length(popUnits),length(popUnits));
            if(PERMUTATIONTEST)
                permAcc = NaN(permutations,length(timepoints));
            end
            parfor iTrainingInterval=1:length(timepoints)
                somaUnits(iter,iTrainingInterval,f,n) = classifierCV(all_XTr{iTrainingInterval},...
                    all_XTrt{iTrainingInterval},all_YTr,all_YTrt,popUnits(unitInds(1:min(length(popUnits),nInputUnits(n)))),nCVFolds,fp,classifier);
                if(PERMUTATIONTEST)
                    if(f<=length(allGroups))
                        for p = 1:permutations
                            permAcc(p,iTrainingInterval) = classifierCV(all_XTr{iTrainingInterval},...
                                all_XTrt{iTrainingInterval},all_YTr,all_YTrt,groupPerms{p}{f}(1:min(length(groupPerms{p}{f}),nInputUnits(n))),nCVFolds,fp,classifier);
                        end
                    end
                end
            end
            tic
            if(PERMUTATIONTEST)
                for p = 1:permutations
                    somaPerm(iter,p,f,n) = mean(permAcc(p,:),2,'omitnan'); 
                    latPerm(iter,p,f,n) = max(permAcc(p,:),[],2);
                end
            end
            toc;
        end
    end
    send(hbar, iter);
end
if(PERMUTATIONTEST)
    F_mean = squeeze(mean(somaPerm,1));
    SS_between = arrayfun(@(g) nRuns*sum(100.*(mean(somaPerm(:,allGroups==g),1)-mean(somaPerm(:,allGroups==g),'all')).^2),unique(allGroups));
    SS_within = arrayfun(@(g) sum(100.*(somaPerm(:,allGroups==g)-mean(somaPerm(:,allGroups==g),1)).^2,'all'),unique(allGroups));
    F_perm = arrayfun(@(g) (SS_between(g)/sum(allGroups==g)-1) /(SS_within(g)/(nRuns*sum(allGroups==g))-sum(allGroups==g)),unique(allGroups)); %max(0,var( - singleDrawVar)
    Lat_mean = squeeze(mean(latPerm,1));
    SS_between = arrayfun(@(g) nRuns*sum((mean(latPerm(:,allGroups==g),1)-mean(latPerm(:,allGroups==g),'all')).^2),unique(allGroups));
    SS_within = arrayfun(@(g) sum((latPerm(:,allGroups==g)-mean(latPerm(:,allGroups==g),1)).^2,'all'),unique(allGroups));
    Lat_perm = arrayfun(@(g) (SS_between(g)/sum(allGroups==g)-1) /(SS_within(g)/(nRuns*sum(allGroups==g))-sum(allGroups==g)),unique(allGroups));
end
allAcc = squeeze(mean(somaUnits(:,:,:,nInputUnits==60),2,'omitnan'));
%% permutation testing
rng('shuffle', 'threefry4x64_20');
hbar = parforProgress(permutations);
for p = 1:permutations
    [latPerm,somaPerm] = deal(NaN(nRuns,length(allGroups),1)); 
    permPopulation = cellfun(@(s) s(randperm(length(s))), fullSubPop, 'UniformOutput',false);
    groupPerms = cellfun(@(p,s) arrayfun(@(i) p(s(i):s(i+1)-1),1:length(s)-1,'UniformOutput',false), permPopulation, subGroupInd, 'UniformOutput',false);
    groupPerms = [groupPerms{:}];
    parfor iter = 1:nRuns
        [sp,lp] = deal(NaN(1,length(allGroups)));
        for f = 1:length(allGroups)
            units = groupPerms{f}(uUnits{iter}{groupInds(f)});
            iterAcc= zeros(1,length(timepoints));
            for iTrainingInterval = 1:length(timepoints)
                for iCV = 1:nCVFolds
                    if(isempty(fp))
                        tr = cacheTrain{iter}{iTrainingInterval}{iCV}(units,:);
                        XTst = CTs.cacheTest{iter}{iTrainingInterval}{iCV}(units,:);
                    else
                        [~,tr] =  fp.set_properties_with_training_data(cacheTrain{iter}{iTrainingInterval}{iCV}(units,:));
                        [~,XTst] =  fp.set_properties_with_training_data(cacheTest{iter}{iTrainingInterval}{iCV}(units,:));
                    end
                    clT= classifier.train(tr, all_YTr);
                    [ia,~] = clT.test(XTst);
                    iterAcc(iTrainingInterval)=  iterAcc(iTrainingInterval) +sum(ia-all_YTrt==0)/length(all_YTrt);
                end
            end
            iterAcc = iterAcc / nCVFolds;
            sp(f) = mean(iterAcc);
            [~,lp(f)] = max(iterAcc);
        end
        somaPerm(iter,:) = sp;
        latPerm(iter,:) = lp;
    end
    F_mean(p,:) = mean(somaPerm,1);
    SS_between = arrayfun(@(g) nRuns*sum(100.*(mean(somaPerm(:,allGroups==g),1)-mean(somaPerm(:,allGroups==g),'all')).^2),unique(allGroups));
    SS_within = arrayfun(@(g) sum(100.*(somaPerm(:,allGroups==g)-mean(somaPerm(:,allGroups==g),1)).^2,'all'),unique(allGroups));
    F_perm(p,:) = arrayfun(@(g) (SS_between(g)/sum(allGroups==g)-1) /(SS_within(g)/(nRuns*sum(allGroups==g))-sum(allGroups==g)),unique(allGroups)); %max(0,var( - singleDrawVar)
    Lat_mean(p,:) = mean(latPerm,1);
    SS_between = arrayfun(@(g) nRuns*sum((mean(latPerm(:,allGroups==g),1)-mean(latPerm(:,allGroups==g),'all')).^2),unique(allGroups));
    SS_within = arrayfun(@(g) sum((latPerm(:,allGroups==g)-mean(latPerm(:,allGroups==g),1)).^2,'all'),unique(allGroups));
    Lat_perm(p,:) = arrayfun(@(g) (SS_between(g)/sum(allGroups==g)-1) /(SS_within(g)/(nRuns*sum(allGroups==g))-sum(allGroups==g)),unique(allGroups));
    send(hbar, p);
end
save(savePath+"Permutation_Results.mat",'F_mean','F_perm','Lat_mean','Lat_perm', 'allAcc');
%%
for s = 1:length(shuffleGroups)
    diffRatio = table(Size=[1,nchoosek(length(shuffleGroups{s}),2)],VariableNames=join(subGroups(nchoosek(shuffleGroups{s},2)),"-",2)',...
        VariableTypes=repmat("doublenan",1,nchoosek(length(shuffleGroups{s}),2)));
    pairs = nchoosek(1:length(shuffleGroups{s}),2);
    fo = anova(100.*allAcc(:,shuffleGroups{s}));
    omnibusP = (1+sum(F_perm(:,s) >= fo.stats.F(1)))/(permutations+1);
    diffO = 100.*(allAcc(:,shuffleGroups{s}(pairs(:,1)))-allAcc(:,shuffleGroups{s}(pairs(:,2))));
    if(omnibusP<pVal)
        Find = find(allGroups==s);
        for p = 1:size(pairs,1)
            allDiffs = diff(F_mean(:,Find(pairs(p,:))),1,2);
            diffRatio{1,p}  = (1+sum(abs(allDiffs)>=abs(mean(diffO(:,p),1,'omitnan'))))/(permutations+1);
        end
    end
    diffRatio = addvars(diffRatio,omnibusP,Before=1);
    diffRatio{end+1,:} = [fo.stats.MeanSquares(1),mean(diffO,1,'omitnan')];
    diffRatio{end+1,:} = [fo.stats.F(1),std(diffO,0,1,'omitnan')];
    diffRatio{end+1,:} = [omnibusP<pVal,diffRatio{1,2:end}<(pVal/size(pairs,1))];
    diffRatio = addvars(diffRatio,["pVal";"MeanSquares/Mean";"F-stat/SD";"p<"+num2str(pVal)],Before=1,NewVariableNames="Value");
    diffRatio{:, 2:end} = compose("%.5f", diffRatio{:, 2:end});
    writetable(diffRatio, savePath+strjoin(subGroups(shuffleGroups{s}),"_")+"-stats.txt", 'FileType', 'text','Delimiter','\t');
end
%%
rng('shuffle', 'threefry4x64_20');
parfor p = 1:(permutations)
    shuffleUnits = cell(1,length(shuffleGroups));
    for s = 1:length(shuffleGroups)
        shuffledAcc = allAcc(:,shuffleGroups{s});
        for r = 1:size(shuffledAcc,1)
            shuffledAcc(r,:) = shuffledAcc(r,randperm(size(shuffledAcc,2)));
        end
        shuffleUnits{s} =  shuffledAcc;
    end
    F_perm(p,:) = shuffleUnits;
end
for s = 1:length(shuffleGroups)
    fo = anova(100.*allAcc(:,shuffleGroups{s}));
    omnibusP = (1+sum(fi>=fo.stats.F(1)))/(permutations+1);
    fi = cellfun(@(n) n.stats.F(1),cellfun(@(a) anova(100.*a),F_perm(:,s),'UniformOutput',false));
    diffRatio = table(Size=[1,nchoosek(length(shuffleGroups{s}),2)],VariableNames=join(subGroups(nchoosek(shuffleGroups{s},2)),"-",2)',...
        VariableTypes=repmat("double",1,nchoosek(length(shuffleGroups{s}),2)));
    pairs = nchoosek(1:length(shuffleGroups{s}),2);
    diffO = 100.*(allAcc(:,shuffleGroups{s}(pairs(:,1)))-allAcc(:,shuffleGroups{s}(pairs(:,2))));
    if(omnibusP<pVal)
        for p = 1:size(pairs,1)
            allDiffs = cellfun(@(f) f(:,pairs(p,:)), F_perm(:,s), 'UniformOutput',false);
            allDiffs = cellfun(@(f,i) diff(100.*mean([f(sub2ind(size(f),1:size(f,1),i'+1));f(sub2ind(size(f),1:size(f,1),~i'+1))],...
                2,'omitnan')),allDiffs,cellfun(@(r) randi(2,[permutations,1])-1, F_perm(:,s),'UniformOutput',false));
            diffRatio{1,p}  = (1+sum(abs(allDiffs)>=abs(mean(diffO(:,p),1,'omitnan'))))/(permutations+1);
        end
    end
    diffRatio = addvars(diffRatio,omnibusP,Before=1);
    diffRatio{end+1,:} = [fo.stats.MeanSquares(1),mean(diffO,1,'omitnan')];
    diffRatio{end+1,:} = [fo.stats.F(1),std(diffO,0,1,'omitnan')];
    diffRatio = addvars(diffRatio,["pVal";"MeanSquares/Mean";"F-stat/SD"],Before=1,NewVariableNames="Value");
    diffRatio{:, 2:end} = compose("%.5f", diffRatio{:, 2:end});
    writetable(diffRatio, savePath+strjoin(subGroups(shuffleGroups{s}),"_")+"-stats.txt", 'FileType', 'text','Delimiter','\t');
end
%% Figures
accLine = 80;
trCl = [0 .7 0; 1 .5 0; 0 .2 .8; .4 .4 .4; .6 0 .2;];
if(length(nInputUnits)>1)
    trtsPhase = 100.*squeeze(mean(unitAccPhase(:,1:end,:,:),2,'omitnan'));
    testTicks = linspace(1,550,length(nInputUnits));

    nAx = 3; nLs = floor(length(subGroups)/nAx); nL = mod(length(subGroups),nAx);
    figure(); colororder(trCl); nt=tiledlayout(1,nAx);
    %errorbar(nexttile(nt,d),squeeze(mean(trtsPhase(:,:,(nLs*(d-1))+(1:(nLs+nL))),1,'omitnan')),squeeze(std(trtsPhase(:,:,(nLs*(d-1))+(1:(nLs+nL))),0,1)))
    arrayfun(@(d) boxchart(nexttile(nt,d),reshape(repmat(testTicks,size(trtsPhase,1),1,nLs+(nL*(d==nAx))),1,[]),reshape(...
        trtsPhase(:,:,ismember(subGroups,subGroups((nLs*(d-1))+(1:(nLs+(nL*(d==nAx))))))),1,[]),'GroupByColor',reshape(permute(...
        repmat(1:(nLs+(nL*(d==nAx))),size(trtsPhase,1),1,length(testTicks)),[1 3 2]),1,[]),'Notch','on','MarkerStyle','none','BoxWidth',30,'BoxFaceAlpha',1),1:nAx);
    arrayfun(@(d) legend(nexttile(nt,d),subGroups((nLs*(d-1))+(1:(nLs+(nL*(d==nAx))))),'location','southeast'),1:nAx);
    arrayfun(@(d) set(nexttile(nt,d),'ylim',[30 100],'xtick',testTicks,'xticklabels',nInputUnits),1:nAx);
    saveFigures(gcf,savePath,"Boxplots-"+num2str(nRuns),[]);

    figure(); hold on;
    [~,maxT] = max(100.*unitAccPhase(:,:,:,arrayfun(@(s) find(strcmp(subGroups,s)), somatotopicLabs))>=accLine,[],2);
    maxT(maxT==1) = length(xl);
    figure(); boxchart(reshape(repmat(1:length(somatotopicLabs),nRuns,1),1,[]),reshape(squeeze((xl(maxT))),1,[]),'Notch','on');
    xticks(1:length(somatotopicLabs));xticklabels(somatotopicLabs);
    saveFigures(gcf,savePath,"Latency_Somatotopy-"+num2str(nRuns),[]);
end

for t = 1:length(nInputUnits)
    nUnits=nInputUnits(t);
    accTUnit = cell2mat(cellfun(@(s) permute(vertcat(s{:,nInputUnits==nUnits}),[3 2 1]),somaUnits,'UniformOutput',false));
    accTable = array2table(100.*squeeze(mean(unitAccPhase(:,:,nInputUnits==nUnits,:),2,'omitnan')),'VariableNames',subGroups);
    if(any(strcmp(subGroups,"Task")))
        accTable = [accTable,array2table(100.*unitAccPhase(:,:,nInputUnits==nUnits,strcmp(subGroups,"Task")),'VariableNames',"Time-"+timepoints)];
    end
    writetable(accTable,savePath+"Decoding_"+num2str(nUnits)+"_"+num2str(nRuns),'FileType','spreadsheet','UseExcel',true);

    if(length(timepoints)>1)
        em=[mean(em(:,1:2),1,'omitnan'),min(em(:,end)),max(em(:,end))];
        figure(); nAx=4; nLs = floor(length(subGroups)/nAx); nL = mod(length(subGroups),nAx); nt=tiledlayout(1,nAx);
        accAx = arrayfun(@(d) arrayfun(@(a) accTUnit(:,:,a),(nLs*(d-1))+(1:(nLs+(nL*(d==nAx)))),'UniformOutput',false),1:nAx,'UniformOutput',false);
        d = cellfun(@(a,d) cellfun(@(p,c) shadedErrorBar(xl,mean(p,1,'omitnan'),std(p,0,1),'lineProps',...
            {'Color',trCl(c,:),'LineWidth',2},'patchSaturation',.1,'ax',nexttile(nt,d)),a,num2cell(1:length(a))),accAx,num2cell(1:nAx),'UniformOutput',false);
        cellfun(@(x,n) legend(nexttile(nt,n),[x.mainLine],subGroups((nLs*(n-1))+(1:(nLs+(nL*(n==nAx))))),'location','southeast','AutoUpdate','off'),d,num2cell(1:nAx));
        arrayfun(@(d) hold(nexttile(nt,d), 'on'), 1:nAx);
        arrayfun(@(d) set(nexttile(nt,d),'ylim',[0 1],'xlim',[min(xl),max(xl)]),1:nAx);
        arrayfun(@(d) arrayfun(@(a,s) plot(nexttile(nt,d),[a a], [0,1], s),em,["k--","k--","b--","m--"]),1:nAx);
        saveFigures(gcf,savePath,"Temporal_"+num2str(nRuns)+"_"+num2str(nUnits),[]);
    end
    figure(); hold on;
    boxchart(accTable,subGroups,"Notch",'on','MarkerStyle','none');
    xticks(1:length(subGroups));xticklabels(subGroups);
    saveFigures(gcf,savePath,"Groups_"+num2str(nRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,"Time-"+timepoints(1:end),'Notch','on','MarkerStyle','none');
    xticks(1:5:length(timepoints));xticklabels(string(timepoints(1:5:end)));newXLim = [-.1 .1]+double(get(gca,'XLim'));
    plot(newXLim,[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);xlim(newXLim);
    saveFigures(gcf,savePath,"Phases_"+num2str(nRuns)+"_"+num2str(nUnits),[]);

    figure(); hold on;
    bx=boxchart(accTable,["Arm" "Hand" "Trunk" "Face"],'Notch','on','MarkerStyle','none');
    xticks(1:length(somatotopicLabs));xticklabels(["Arm" "Hand" "Trunk" "Face"]);
    plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
    saveFigures(gcf,savePath,"Somatotopy_"+num2str(nRuns)+"_"+num2str(nUnits),[]);

    for i = 1:3
        somaSet = "";
        if(i==2)
            somaSet="-Arm";
        elseif(i==3)
            somaSet="-Hand";
        end
        figure(); hold on;
        bx=boxchart(accTable,subGroups([1 3 2])+somaSet,'Notch','on','MarkerStyle','none');
        xticks(1:3);xticklabels(subGroups([1 3 2]));
        plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
        saveFigures(gcf,savePath,"Types_"+somaSet+num2str(nRuns)+"_"+num2str(nUnits),[]);

        figure(); hold on;
        bx=boxchart(accTable,["Shallow","Deep"]+somaSet,'Notch','on','MarkerStyle','none');
        xticks(1:2);xticklabels(["Shallow","Deep"]);
        plot([-1 1]+double(get(gca,'XLim')),[accLine accLine],'LineWidth',1,'Color','k');ylim([25 105]);
        saveFigures(gcf,savePath,"Laminar_"+somaSet+num2str(nRuns)+"_"+num2str(nUnits),[]);
    end
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

function iterAcc = classifierCV(XTrF,XTsF,yTr,yTs,units,nCVs,f,cl)
iterAcc = 0;
for iCV = 1:nCVs
    if(isempty(f))
        tr = XTrF{iCV}(units,:);
        XTst = XTsF{iCV}(units,:);
    else
        [~,tr] =  f.set_properties_with_training_data(XTrF{iCV}(units,:));
        [~,XTst] =  f.set_properties_with_training_data(XTsF{iCV}(units,:));
    end
    clT= cl.train(tr, yTr);
    [ia,~] = clT.test(XTst);
    iterAcc= iterAcc + (sum(ia-yTs==0)/length(yTs));
end
iterAcc = iterAcc/nCVs;
end