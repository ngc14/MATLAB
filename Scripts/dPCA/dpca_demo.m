% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of timepoints (note: all trials should have same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. For the neurons and 
% conditions with less,fill remaining entries in firingRates with 0 or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If it's filled up with zeros 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)
binWidth = 10; sTrials = 20; time=-.5:1/(1000/binWidth):1;
currD= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
    v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
currD =  cellfun(@(a)cellfun(@(u)max(0,u),a,'UniformOutput',false), currD,'UniformOutput',false);
%% Define parameter grouping
% firingRates array has [N S D T E] size; here we ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
% As explained in the paper, we group stimulus with stimulus/time interaction etc.:
mv = sum(cell2mat(cellfun(@(m)mean(cell2mat(m),2,'omitnan').*1000>1,num2cell([currD{:}],2),'UniformOutput',false)'),2)>sTrials/2 | ...
    mean(cell2mat(cellfun(@(n)mean(cat(3,n{:}),3,'omitnan'),currD,'UniformOutput',false)),2,'omitnan').*1000>1;
trialLength = floor(size(currD{1}{1}, 2) / binWidth);
trialPSTH = cell(1,length(currD));
for n = 1:length(currD)
    trialPSTH{n} = repmat({NaN(length(mv),trialLength)},max(1,sTrials),1);
    for s = 1:length(trialPSTH{n})
        for t = 1:trialLength
            iStart = binWidth * (t-1) + 1;
            iEnd   = binWidth *t;
            trialPSTH{n}{s}(:,t) = sum(currD{n}{s}(:,iStart:iEnd),2);%normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
        end
    end
end
firingRates = cellfun(@(c) cellfun(@(s) (conv2(resize(s(:,unique(round(ms_bins./binWidth))),[size(s,1),length(unique(round(ms_bins./binWidth)))+length(gausswin(ceil(smoothWin/binWidth)))-1],...
    'Pattern','edge','side','both'),transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),trialPSTH,'UniformOutput',false);
firingRates = cellfun(@(t) reshape(cellfun(@(r) reshape(r(mv,:),sum(mv),1,[]),t,'UniformOutput',false),1,[]),firingRates,'UniformOutput',false);
firingRates = cell2mat(permute(cat(4,firingRates{:}),[1 4 3 2]));
%firingRates = cell2mat(cellfun(@(r) circshift(r,randi([2*binWidth,size(r,3)-2*binWidth],1),3), num2cell(firingRates,3),'UniformOutput',false));
trialNum = ones(size(firingRates,1:2)).*size(firingRates,length(size(firingRates)));
combinedParams = {{1,[1 2]}, {2}};
margNames = {'Condition','Condition-Invariant'};
margColours = [23 100 171; 200 160 43; 150 150 150;]/256;
timeEventConds = cell2mat(cellfun(@(i) mean(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):timeBins(end)),1,'omitnan'),segInds,'UniformOutput',false));
timeEvents = time(round([mean(timeEventConds(:,1:2),1,'omitnan'),timeEventConds(3,3),mean([timeEventConds(1,3);timeEventConds(2,3)],1,'omitnan')]));
% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less and marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
%    combinedParams = {{1, [1 2]}, {2}}
% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,1:trialNum(n,s))), 1)), 'Something is wrong!')
        end
    end
end
% firingRatesAverage = cell2mat(cellfun(@(r) reshape(cell2mat(r),size(r{1},1),1,[]),cellfun(@(s) cellfun(@(t) cell2mat(cellfun(@(n)circshift(n,randi([2*binWidth,length(n)-2*binWidth],1)),n,num2cell(t(mv,unique(round(ms_bins./binWidth))),2),'UniformOutput',false)),s,'UniformOutput',false)',trialPSTH, 'UniformOutput',false),'UniformOutput',false));
firingRatesAverage = mean(firingRates(:,:,:,:), length(size(firingRates)),'omitnan');
%% Step 1: PCA of the dataset
X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));
[W,~,~] = svd(X, 'econ');
W = W(:,1:20);
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);
explVar = dpca_explainedVariance(firingRatesAverage, W, W,'combinedParams', combinedParams);
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, 'time', time,'timeEvents', timeEvents,...
    'marginalizationNames', margNames,'marginalizationColours', margColours);
% Step 2: PCA in each marginalization separately
dpca_perMarginalization(firingRatesAverage, @dpca_plot_default,'combinedParams', combinedParams);
% Step 3: dPCA without regularization and ignoring noise covariance
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which mar
[W,V,whichMarg] = dpca(firingRatesAverage, 5,'combinedParams', combinedParams);
explVar = dpca_explainedVariance(firingRatesAverage, W, V,'combinedParams', combinedParams);
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 2, ...
    'legendSubplot', 16);
%% Step 4: dPCA with regularization
% Once computed, you can simply load lambdas out of file:load('tmp_optimalLambdas.mat', 
% 'optimalLambda'). Note that this now includes noise covariance matrix 
% Cnoise which provides substantial regularization itself (even with lambda=0).
somaIndex = contains(string(somaTable(mv)),["Arm","Hand"]);
optimalLambda = dpca_optimizeLambda(firingRatesAverage(somaIndex,:,:), ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), ...
    'combinedParams', combinedParams, 'simultaneous', false,'numRep', 10);
Cnoise = dpca_getNoiseCovariance(firingRatesAverage(somaIndex,:,:), ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), 'simultaneous', false,'type','pooled');
[W,V,whichMarg] = dpca(firingRatesAverage(somaIndex,:,:,:),50,'combinedParams', combinedParams,'lambda', optimalLambda,'Cnoise', Cnoise);
explVar = dpca_explainedVariance(firingRatesAverage(somaIndex,:,:,:), W, V, 'combinedParams', combinedParams);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 2,           ...
    'legendSubplot', {16,params.condNames},'ylims',[]);
%
X = firingRatesAverage(somaIndex,:)';
Xcen = bsxfun(@minus, X, mean(X));
Z = Xcen * W;
%%
dataDim = size(firingRatesAverage);
first3Proj = reshape(Z(:,[find(whichMarg==2,3),find(whichMarg==1,3)])', [length([find(whichMarg==2,3),find(whichMarg==1,3)]) dataDim(2:end)]);
figure();
for c = 1:3
    subplot(1,3,c); hold on; title(conditions(c));
    plot3(squeeze(first3ProjArm(1,c,:)),squeeze(first3ProjArm(2,c,:)),squeeze(first3ProjArm(3,c,:)),'Color',[0 1 .2]);
    plot3(squeeze(first3ProjHand(1,c,:)),squeeze(first3ProjHand(2,c,:)),squeeze(first3ProjHand(3,c,:)),'Color',[.8 0 1])
    arrayfun(@(e) scatter3(first3ProjArm(1,c,e),first3ProjArm(2,c,e),first3ProjArm(3,c,e),'filled','MarkerFaceColor',[0 1 .2]), round(mean(segVals{c}(:,1:3),1,'omitnan')))
    arrayfun(@(e) scatter3(first3ProjHand(1,c,e),first3ProjHand(2,c,e),first3ProjHand(3,c,e),'filled','MarkerFaceColor',[.8 0 1]), round(mean(segVals{c}(:,1:3),1,'omitnan')))
    all3Dim = squeeze(mean(first3ProjArm,2,'omitnan')); scatter3(all3Dim(1,1),all3Dim(2,1),all3Dim(3,1),'black','*','sizeData',550)
    all3Dim = squeeze(mean(first3ProjHand,2,'omitnan')); scatter3(all3Dim(1,1),all3Dim(2,1),all3Dim(3,1),'black','*','sizeData',550)
    view(10,15);
end
saveFigures(gcf,"S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\","Traj-Condition-Invariant",[]);
arrayfun(@(a) view(a,90,0), gcf().Children);
saveFigures(gcf,"S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\","2DTraj-Condition-Invariant",[]);


figure();
for c = 1:3
    subplot(1,3,c); hold on; title(conditions(c));
    plot3(squeeze(first3ProjArm(4,c,:)),squeeze(first3ProjArm(5,c,:)),squeeze(first3ProjArm(6,c,:)),'Color',[0 1 .2]);
    plot3(squeeze(first3ProjHand(4,c,:)),squeeze(first3ProjHand(5,c,:)),squeeze(first3ProjHand(6,c,:)),'Color',[.8 0 1])
    arrayfun(@(e) scatter3(first3ProjArm(4,c,e),first3ProjArm(5,c,e),first3ProjArm(6,c,e),'filled','MarkerFaceColor',[0 1 .2]), round(mean(segVals{c}(:,1:3),1,'omitnan')))
    arrayfun(@(e) scatter3(first3ProjHand(4,c,e),first3ProjHand(5,c,e),first3ProjHand(6,c,e),'filled','MarkerFaceColor',[.8 0 1]), round(mean(segVals{c}(:,1:3),1,'omitnan')))
    all3Dim = squeeze(mean(first3ProjArm,2,'omitnan')); scatter3(all3Dim(4,1),all3Dim(5,1),all3Dim(6,1),'black','*','sizeData',550)
    all3Dim = squeeze(mean(first3ProjHand,2,'omitnan')); scatter3(all3Dim(4,1),all3Dim(5,1),all3Dim(6,1),'black','*','sizeData',550)
    view(10,15);
end
saveFigures(gcf,"S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\","Traj-Condition",[]);
arrayfun(@(a) view(a,90,0), gcf().Children);
saveFigures(gcf,"S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\","2DTraj-Condition",[]);

a = corr(Z(:,1:dims));
b = V(:,1:dims)'*V(:,1:dims);
[~, psp] = corr(V(:,1:dims), 'type', 'Kendall');
map = tril(a,-1)+triu(b);
r = [5 48 97]/256;w = [.95 .95 .95];b = [103 0 31]/256;c1 = zeros(128,3);c2 = zeros(128,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 128);c2(:,i) = linspace(w(i), b(i), 128);
end
figure();image(round(map*128)+128 + 2);colormap([c1;c2]); colorbar('Ticks',[1 128 256],'TickLabels',[-1 0 1]);
%% Optional: estimating "signal variance"
explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);
% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,'Cnoise',Cnoise);
%% Optional: decoding
decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
decodingClasses = {(1:3)',[1:3]'};
accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'filename', 'tmp_classification_accuracy.mat');
dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
    'filename', 'tmp_classification_accuracy.mat');
dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif);

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end