% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of timepoints (note:trials should have same length in time)
% trialNum -- number of trials for each neuron in each S,D condition 
% firingRates -- all single-trial data together, massive array. Here
% If it's filled up with zeros 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)
binWidth = 10; sTrials = 20; time=-.5:1/(1000/binWidth):1;
combinedParams = {{1,[1 2]}, {2}};
margNames = {'Condition','Condition-Invariant'};
margColours = [23 100 171; 200 160 43; 150 150 150;]/256;
timeEventConds = cell2mat(cellfun(@(i) mean(findBins(params.bins(i),time(1):1/(1000/binWidth):time(end)),1,'omitnan'),segInds,'UniformOutput',false));
timeEvents = time(round([mean(timeEventConds(:,1:2),1,'omitnan'),timeEventConds(3,3),mean([timeEventConds(1,3);timeEventConds(2,3)],1,'omitnan')]));
currD= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
    v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
currD =  cellfun(@(a)cellfun(@(u)max(0,u),a,'UniformOutput',false), currD,'UniformOutput',false);
%% Define parameter grouping
% firingRates array has [N S D T E] size; ignore the 1st dimension (neurons)
% marginalizations: 1 - stimulus, 2 - decision, 3 - time
% 3 pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% 1 three-way interaction:
%    [1 2 3] - rest
% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
%
% For two parameters (stimulus and time), firingRates array of size [N S T E]
% marginalizations: 1 - stimulus, 2 - time, [1 2] - stimulus/time interaction
%    combinedParams = {{1, [1 2]}, {2}}
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
            trialPSTH{n}{s}(:,t) = sum(currD{n}{s}(:,iStart:iEnd),2);
        end
    end
end
firingRates = cellfun(@(c) cellfun(@(s) (conv2(resize(s(:,unique(round(...
    ms_bins./binWidth))),[size(s,1),length(unique(round(ms_bins./binWidth)))+length(gausswin(ceil(smoothWin/binWidth)))-1],...
    'Pattern','edge','side','both'),transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),trialPSTH,'UniformOutput',false);
firingRates = cellfun(@(t) reshape(cellfun(@(r) reshape(r(mv,:),sum(mv),1,[]),t,'UniformOutput',false),1,[]),firingRates,'UniformOutput',false);
firingRates = cell2mat(permute(cat(4,firingRates{:}),[1 4 3 2]));
%firingRates = cell2mat(cellfun(@(r) circshift(r,randi([2*binWidth,size(r,3)-2*binWidth],1),3), num2cell(firingRates,3),'UniformOutput',false));
trialNum = ones(size(firingRates,1:2)).*size(firingRates,length(size(firingRates)));
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
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);
explVar = dpca_explainedVariance(firingRatesAverage, W, W,'combinedParams', combinedParams);
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, 'time', time,'timeEvents', timeEvents,...
    'marginalizationNames', margNames,'marginalizationColours', margColours);
% Step 2: PCA in each marginalization separately
dpca_perMarginalization(firingRatesAverage, @dpca_plot_default,'combinedParams', combinedParams);
% Step 3: dPCA w/out regularization/noise covariance. W=decoder,V=encoder(ordered by explained variance)
[W,V,whichMarg] = dpca(firingRatesAverage, 5,'combinedParams', combinedParams);
explVar = dpca_explainedVariance(firingRatesAverage, W, V,'combinedParams', combinedParams);
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default,'explainedVar', explVar, ...
    'marginalizationNames', margNames,'marginalizationColours', margColours, 'whichMarg', whichMarg,...
    'time', time,'timeEvents', timeEvents,'timeMarginalization', 2,'legendSubplot', 16);
%% Step 4: dPCA with regularization
%load('optimalLambda'). Note that it includes noise covariance matrix Cnoise 
% which provides substantial regularization itself (even with lambda=0).
somaIndex = contains(string(somaTable(mv)),["Arm","Hand"]); numRep = 10; dims = 20;
optimalLambda = dpca_optimizeLambda(firingRatesAverage(somaIndex,:,:), ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), ...
    'combinedParams', combinedParams, 'simultaneous', false,'numRep', numRep);
Cnoise = dpca_getNoiseCovariance(firingRatesAverage(somaIndex,:,:), ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), 'simultaneous', false,'type','pooled');
[W,V,whichMarg] = dpca(firingRatesAverage(somaIndex,:,:,:),dims*(length(combinedParams)+1),...
    'combinedParams', combinedParams,'lambda', optimalLambda,'Cnoise', Cnoise);
explVar = dpca_explainedVariance(firingRatesAverage(somaIndex,:,:,:), W, V, 'combinedParams', combinedParams);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 2,...
    'legendSubplot', {16,params.condNames},'ylims',[]);
dataDim = size(firingRatesAverage);
X = firingRatesAverage(somaIndex,:)';
Xcen = bsxfun(@minus, X, mean(X));
Z = Xcen * W;
first3Proj = reshape(Z(:,[find(whichMarg==2,dims),find(whichMarg==1,dims)])',[2*dims dataDim(2:end)]);
%[W,~,~] = svd(Xcen, 'econ'); W = W(:,1:dims);
%%
saveFig = false; lineColor = [.8 0 1];% [0 1 .2];
savePath = "S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\";
figure(1);
for c = 1:length(conditions)
    subplot(1,3,c); hold on; title(conditions(c));
    plot3(squeeze(first3Proj(1,c,:)),squeeze(first3Proj(2,c,:)),squeeze(first3Proj(3,c,:)),'Color',lineColor);
    arrayfun(@(e) scatter3(first3Proj(1,c,e),first3Proj(2,c,e),first3Proj(3,c,e),'filled','MarkerFaceColor',lineColor), ...
        round(mean(segVals{c}(:,1:3),1,'omitnan')))
    all3Dim = scatter3(first3Proj(1,1),first3Proj(2,1),first3Proj(3,1),'black','*','sizeData',550);
    view(10,15);
end
% XY:view(0,90); XZ:view(0,0); YZ:view(90,0);
if(saveFig)
    saveFigures(figure(1),savePath,"Traj-CondInvariant",[]);arrayfun(@(a)view(a,0,90),gcf().Children); saveFigures(figure(1),savePath,"2DTraj-CondInvariant",[]);
end
figure(2);
for c = 1:length(conditions)
    subplot(1,3,c); hold on; title(conditions(c));
    plot3(squeeze(first3Proj(4,c,:)),squeeze(first3Proj(5,c,:)),squeeze(first3Proj(6,c,:)),'Color',lineColor);
    arrayfun(@(e) scatter3(first3Proj(4,c,e),first3Proj(5,c,e),first3Proj(6,c,e),'filled','MarkerFaceColor',lineColor),...
        round(mean(segVals{c}(:,1:3),1,'omitnan')))
    scatter3(first3Proj(4,1),first3Proj(5,1),first3Proj(6,1),'black','*','sizeData',550);
    view(10,15);
end
if(saveFig)
    saveFigures(figure(2),savePath,"Traj-Cond",[]);arrayfun(@(a)view(a,0,90),gcf().Children);saveFigures(figure(2),savePath,"2DTraj-Cond",[]);
end
%%
a = corr(Z);b = V'*V;[~, psp] = corr(V(:,1:20), 'type', 'Kendall');
map = tril(a,-1)+triu(b);
r = [5 48 97]/256;w = [.95 .95 .95];b = [103 0 31]/256;c1 = zeros(128,3);c2 = zeros(128,3);
for i=1:3;c1(:,i) = linspace(r(i), w(i), 128);c2(:,i) = linspace(w(i), b(i), 128);end
image(round(map*128)+128 + 2);colormap([c1;c2]); colorbar('Ticks',[1 128 256],'TickLabels',[-1 0 1]);
%% estimating "signal variance"
explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams,'Cnoise', Cnoise, 'numOfTrials', trialNum);
% Note how the pie chart changes relative to the previous figure.
% it is displaying percentages of (estimated) signal PSTH variances,not total.
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames,'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 3,           ...
    'legendSubplot', 16,'Cnoise',Cnoise);
%% decoding
decodingClasses = {(1:3)',[1:3]'};%{[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda,'combinedParams', combinedParams, 'decodingClasses', decodingClasses, ...
    'simultaneous', false, 'numRep', 5, 'filename', 'tmp_classification_accuracy.mat');
dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, 'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, 'simultaneous', false, ...
    'numRep', 5,'numShuffles', 20, 'filename', 'tmp_classification_accuracy.mat');
dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg',whichMarg,'time',time,'timeEvents',timeEvents,'timeMarginalization', 3,           ...
    'legendSubplot',16, 'componentsSignif', componentsSignif);

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end