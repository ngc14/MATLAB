% load("D:\EMG.mat");
% load('S:\Lab\ngc14\Working\Both\EMG\EMG.mat','alignedSig', 'normTrials','sessionDates','rawActivity');
Fs = 2000;
num_cv_splits = 10;
num_labels = 2;
ftLength = 1;
testUnits = [1, 2, 5, 10, 15, 25, 35, 50, 75, 100];
alignments = [{'GoSignal'},{'GoSignal'},{'StartReach'},{'StartHold'},{'StartWithdraw'}];
alignWindows = {[-.85 .15],[-.15 .5],[-.5 .5],[-.20 .20],[-.15 .5]};
phaseWindows = {[-.6 -.1],[0.0, .2],[-.15 .05],[-.2 0],[-.15 .05]};
[alignWindows,phaseWindows] = cellfun(@(a,p) deal(a.*Fs,p.*Fs),alignWindows,phaseWindows,'UniformOutput',false);
numRuns = 500;
phaseNames = ["Baseline","Go", "Reach", "Grasp","Withdraw"];
condLabels = {'E', 'L', 'P'};
groupings = {{"Deltoid.mat","Biceps.mat","Triceps.mat"},...
    {"Wrist Extensor.mat","Wrist Flexor.mat","Digit Extensor.mat","Digit Flexor.mat"}};
groupNames = cellfun(@(g) g{1}(1:end-4), groupings, 'UniformOutput', false);%["Arm", "Hand"];
groupInds = cellfun(@(g) contains(string([groupings{:}]),string(g)), groupings, 'UniformOutput',false);
%%
% siteInds = cellfun(@(n) arrayfun(@(a) find(cellfun(@length,n)==a), unique(cellfun(@length,n)),'UniformOutput',false), normTrials, 'UniformOutput',false);
% for a = 1:size(siteInds,1)
%     currConds = siteInds(a,:);
%     currSig = alignedSig(a,~strcmp(phaseNames, "Baseline") & ~strcmp(phaseNames, "Withdraw"));
%     sessionEMGs{a} = cellfun(@(c) cellfun(@(t,m) cellfun(@(n) t(n,:),m,'UniformOutput',false)', c,currConds,'UniformOutput',false),currSig, 'UniformOutput',false);
% end
% sessionEMGs = cellfun(@(c) cellfun(@(m) cellfun(@(g) vertcat(m{g}), groupInds, 'UniformOutput',false), c, 'UniformOutput',false), sessionEMGs, 'UniformOutput',false);
% sessionEMGs = cellfun(@(c) vertcat(c{:}), sessionEMGs, 'UniformOutput',false);
% sessionEMGs = cat(3,sessionEMGs{:});
%%
alignedPhase = {};
for a = 1:length(phaseNames)
    p=fix(phaseWindows{a}-alignWindows{a});
    alignedPhase{a} = cellfun(@(a) cellfun(@(m) mean(m(:,p(1)+1:...
        min(arrayfun(@(a) a.*double((a>p(1))./(a>p(1))),[size(m,2),size(m,2)+p(end)-1]))),2,'omitnan'),a,...
        'UniformOutput',false),alignedSig(:,a), 'UniformOutput', false);
end
alignedPhase = cellfun(@(r) reshape(r,1,1,[]), [alignedPhase{:}], 'UniformOutput', false);
alignedPhase = reshape([alignedPhase{:}],size(alignedPhase,1),size(alignedPhase,2),[]);
rAc = arrayfun(@(g) num2cell(alignedPhase(:,:,g),2), 1:length(string([groupings{:}])), 'UniformOutput',false);
rAc = cellfun(@(g) cellfun(@cell2mat, g,'UniformOutput',false), rAc, 'UniformOutput',false);
allSessionSegs = horzcat(sessionSegs{:});
%%
maxSession = num2cell(horzcat(sessionDates{:}),2);
splitGroup = {};
trainingFeatures = {};
for f = 1:length(groupInds)
    currVal = num2cell(horzcat(rAc{groupInds{f}}),2);
    recordingSessions= cellfun(@(m) unique(string([m{groupInds{f}}])), maxSession, 'UniformOutput',false);
    currGroup= cellfun(@(r,m,s) cellfun(@(mg,ms)arrayfun(@(i) ms(strcmp(mg,i),:), s,'UniformOutput',false),...
        m(groupInds{f}),r,'UniformOutput',false), currVal, maxSession,recordingSessions,'UniformOutput',false);
    currGroup = cellfun(@(c)cellfun(@(m) cellfun(@(pa,mi) vertcat(pa,NaN(mi-size(pa,1),size(pa,2))),m,...
        num2cell(max(cell2mat(cellfun(@(a) cellfun(@(s) size(s,1),a),c,'UniformOutput',false)'),[],1)),'UniformOutput',false), c,'UniformOutput', false),...
        currGroup, 'UniformOutput', false);
    splitGroup{f} = cellfun(@(c) cellfun(@(a) mean(pagemtimes(reshape(cellfun(@(s) sum(all(~isnan(s),2)), a),1,1,[])....
        ./sum(cellfun(@(s) sum(all(~isnan(s),2)), a)),cell2mat(reshape(a,1,1,[]))),3,'omitnan'),num2cell(vertcat(c{:}),1),'UniformOutput',false),...
        currGroup, 'UniformOutput',false);
    trainingFeatures{f} = cellfun(@(s1,s2,s3) vertcat(s1,s2,s3), splitGroup{f}{1},splitGroup{f}{2},splitGroup{f}{3},'UniformOutput',false)';
end
%%
%    trainingLabs{f} = cellfun(@(l,p) cellfun(@(s) repmat(cellstr(l),size(s,1),1), p, 'UniformOutput', false), condLabels', splitGroup{f}, 'UniformOutput',false);
allLabs = {};
for g=1:length(groupings)
    currFeatures = splitGroup{g}';%squeeze(sessionEMGs(:,g,:));
    phaseLabs =  cellfun(@(pp,l) cellfun(@(lb) repmat(l,size(lb,1),1),pp, 'UniformOutput',false),...
        currFeatures,condLabels, 'UniformOutput',false);
    allLabs{g} = cellfun(@(s1,s2,s3) cellstr(cat(1,s1,s2,s3)), phaseLabs{1}, phaseLabs{2}, phaseLabs{3},'UniformOutput',false);
end
allLabs = horzcat(allLabs{:})';
allFeatures = vertcat(trainingFeatures{:});
subpopulations = cellfun(@(t,g) repmat(string(g),size(t,1),1), trainingFeatures, groupNames, 'UniformOutput',false);
subpopulations = vertcat(subpopulations{:});
goodUnits = find_sites_with_k_label_repetitions(allLabs,num_cv_splits*num_labels,{'E','L','P'});
trialC = allFeatures(goodUnits);
trialLabs = allLabs(goodUnits);
subpopulations = subpopulations(goodUnits);
subpopulations = cellfun(@(s) strcmp(subpopulations, s), unique(subpopulations), 'UniformOutput', false);
trialInds =  cellfun(@(t) true(size(t,1),1), trialC, 'UniformOutput',false);
trainingSet = cellfun(@(t,i) t(i,:), trialC, trialInds, 'UniformOutput',false);
trainingLabs =  cellfun(@(t,i) cellstr(t(i,:)), trialLabs, trialInds, 'UniformOutput',false);

binning_parameters = struct('sampling_interval', 1,'end_time', size(allFeatures{1},2)*ftLength,'start_time', 1,'bin_width',ftLength);
binning_parameters.the_bin_start_times = binning_parameters.start_time:binning_parameters.sampling_interval:binning_parameters.end_time;
binning_parameters.the_bin_widths = repmat(binning_parameters,1,length(binning_parameters.the_bin_start_times));
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
fp = {};
cle = max_correlation_coefficient_CL;
CVAcc = NaN(num_cv_splits,length(dsr.time_periods_to_get_data_from),length(dsr.time_periods_to_get_data_from));
somaUnits = repmat({repmat({CVAcc},length(testUnits),length(groupings))},numRuns,1);
continueIterations = true;
parfor iter = 1:numRuns
    clT = repmat({cle},1,length(dsr.time_periods_to_get_data_from));
    iterAcc = CVAcc;
    [all_XTr, all_YTr, all_XTrt, all_Ytrt] = dsr.get_data;
    for s = 1:length(groupings)
        for n = 1:length(testUnits)
            popInds = find(ismember(goodUnits,find(subpopulations{s})));
            unitSample = popInds(randperm(length(popInds),min(testUnits(n), length(popInds))));
            for iCV = 1:num_cv_splits
                for iTrainingInterval = 1:size(all_XTr,2)
                    if(isempty(fp))
                        XTrF = all_XTr{iTrainingInterval}{iCV};
                    else
                        [~,XTrF] = fp.set_properties_with_training_data(all_XTr{iTrainingInterval}{iCV});
                    end
                    XTrF = XTrF(unitSample,:);
                    if(ftLength==1)
                        cl = clT{iTrainingInterval};
                        clT{iTrainingInterval} = cl.train(XTrF, all_YTr);
                    else
                        clT = fitcknn(XTrF',all_YTr,'NSMethod','kdtree');
                        clT.Distance = 'eucledian';
                        clT.NumNeighbors = 3;
                    end
                    for iTestInterval = 1:size(all_XTrt,2)
                        if(isempty(fp))
                            XTrFt = all_XTrt{iTestInterval}{iCV};
                        else
                            [~,XTrFt] = fp.set_properties_with_training_data(all_XTrt{iTestInterval}{iCV});
                        end
                        XTrFt = XTrFt(unitSample,:);
                        if(ftLength==1)
                            [predicted_labels decision_values] = clT{iTrainingInterval}.test(XTrFt);
                        else
                            predicted_labels = clT.predict(XTrFt');
                        end
                        iterAcc(iCV, iTrainingInterval,iTestInterval) = sum(predicted_labels==all_Ytrt)/length(all_Ytrt);
                    end
                end
                somaUnits{iter}{n,s} =iterAcc;
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
%close all;
cls = distinguishable_colors(length(groupInds));
ls = {'-','-.'};
thresh = [];
for s = 1:length(groupInds)
    for itrain =1:length(phaseNames)
        subplot(2,3,itrain);
        hold on;
        title(phaseNames(itrain))
        allAcc = cellfun(@(u) cell2mat(cellfun(@(v) mean(v(:,itrain,itrain),2), u(:,s), 'UniformOutput',false)'), somaUnits, 'UniformOutput', false);
        allAcc = cell2mat(reshape(allAcc,1,1,[]));
        errorbar(testUnits,mean(allAcc,[1,3],'omitnan'),std(allAcc,0,[1,3],'omitnan')./sqrt(numRuns*num_cv_splits),'color',cls(s,:),'Linestyle',ls(1+(s>4)));
        inter = polyxpoly(testUnits,mean(allAcc,[1,3],'omitnan'),testUnits,repmat(0.9,size(testUnits)));
        if(~isempty(inter))
            thresh(s,itrain) = fix(inter(1));
        else
            thresh(s,itrain) = NaN;
        end
        if(s==length(groupInds))
            if(itrain==1)
                legend([groupNames],'Location','northwest','AutoUpdate','off')
            end
            ylim([0.3 1])
            plot([1 max(testUnits)],[0.90 0.90],'k--','LineWidth',1);
            plot([1 max(testUnits)],[0.33 0.33],'k:','LineWidth',1.5);
        end
    end
end