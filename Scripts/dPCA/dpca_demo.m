%% This section creates toy data.
%
% It should be replaced by actual experimental data. The data should be
% joined in three arrays of the following sizes (for the Romo-like task):
%
% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

%% Define parameter grouping
% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
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
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:
binWidth = 10; sTrials = 20;
currD= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
    v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
currD =  cellfun(@(a)cellfun(@(u)max(0,u),a,'UniformOutput',false), currD,'UniformOutput',false);
%%
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
%firingRates = cell2mat(cellfun(@(r) r(:,:,randperm(length(r))), num2cell(firingRates,3),'UniformOutput',false));
%firingRatesAverage = mean(firingRates, length(size(firingRates)),'omitnan');
trialNum = ones(size(firingRates,1:2)).*size(firingRates,length(size(firingRates)));
combinedParams = {{1, [1 2]}, {2}};
margNames = {'Condition','Condition-Independent'};
margColours = [23 100 171; 150 150 150;]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = time(round(mean(cell2mat(cellfun(@(i) findBins(params.bins(i),timeBins(1):1/(1000/binWidth):timeBins(end)),...
    segInds,'UniformOutput',false)),1,'omitnan')));
% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,s,1:trialNum(n,s))), 1)), 'Something is wrong!')
        end
    end
end
firingRatesAverage = cell2mat(cellfun(@(r) reshape(cell2mat(r),size(r{1},1),1,[]),cellfun(@(s) cellfun(@(t) cell2mat(cellfun(@(n)n,....(randperm(length(n))),...
    num2cell(t(mv,unique(round(ms_bins./binWidth))),2),'UniformOutput',false)),s,'UniformOutput',false)',smoothedPSTH, 'UniformOutput',false),'UniformOutput',false));
somaIndex = contains(string(somaTable(mv)),["Arm","Hand"]);
firingRatesAverage = firingRatesAverage(somaIndex,:,:);
%% Step 1: PCA of the dataset
X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 5, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

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

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), ...
    'combinedParams', combinedParams, ...
    'simultaneous', false, ...
    'numRep', 10, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), 'simultaneous', false,'type','pooled');

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 2,           ...
    'legendSubplot', 16);

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