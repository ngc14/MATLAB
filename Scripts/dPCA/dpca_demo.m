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
binWidth = 10; sTrials = 20; time=-.5:1/(1000/binWidth):1; dims = 10; 
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
[W,~,~] = svd(firingRatesAverage(:,:), 'econ');
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
somaIndex = cell2mat(arrayfun(@(a) find(somaTable(mv)==a,min(groupcounts(somaTable(mv)))),unique(somaTable(mv)),'UniformOutput',false));% channels>16; % %
optimalLambda = dpca_optimizeLambda(firingRatesAverage(somaIndex,:,:),firingRates(somaIndex,:,:,:),...
    trialNum(somaIndex,:),'combinedParams', combinedParams, 'simultaneous', false,'numRep', 10);
Cnoise = dpca_getNoiseCovariance(firingRatesAverage(somaIndex,:,:), ...
    firingRates(somaIndex,:,:,:), trialNum(somaIndex,:), 'simultaneous', false,'type','pooled');
[W,V,whichMarg] = dpca(firingRatesAverage(somaIndex,:,:,:),dims*(length(combinedParams)+1),...
    'combinedParams', combinedParams,'lambda', optimalLambda,'Cnoise', Cnoise);
explVar = dpca_explainedVariance(firingRatesAverage(somaIndex,:,:,:), W, V, 'combinedParams', combinedParams);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 2,...
    'legendSubplot', {16,params.condNames},'ylims',[]);
Xcen = bsxfun(@minus, firingRatesAverage(somaIndex,:)', mean(firingRatesAverage(somaIndex,:),2)');
Z = Xcen * W;
%%
figure();
tiledlayout(3,3);
gcounts = {};
for c = 1:length(combinedParams)+1
    mg = find(whichMarg==c,dims/2);
    for g = 1:3
        if(g==1)
            gcounts{g} = discretize(channels(somaIndex),[1,12,22,33]);
            gName = "Laminar";
        elseif(g==2)
            gcounts{g} = discretize(max(0,location(somaIndex,1)*ImagingParameters.px2mm),[0:4,max(location(somaIndex,1)*ImagingParameters.px2mm)]);
            gName = "CR";
        else
            gcounts{g} = discretize(min(7,max(1,location(somaIndex,2)*ImagingParameters.px2mm)),[1:6,max(location(somaIndex,2)*ImagingParameters.px2mm)]);
            gName = "ML";
        end
        nexttile(); hold on;
        if(c>length(combinedParams))
            title(gName + " Counts");
            bx = bar(groupcounts(gcounts{g}),'FaceColor','flat');
            co = num2cell(colororder,2);
            for l = 1:size(bx.CData,1)
                bx.CData(l,:) = co{l};
            end
            xticks(1:size(bx.CData,1)); xticklabels(1:size(bx.CData,1));
        else
            title(margNames{c} +"-"+ gName); ylim([0 1]);
            groupVals = cell2mat(arrayfun(@(a) sum(W(gcounts{g}==a,mg).^2)./sum(gcounts{g}==a),unique(gcounts{g}),'UniformOutput',false))';
            bx=bar(1:length(mg),cell2mat(cellfun(@(a,b) a./b,num2cell(groupVals,2),num2cell(sum(groupVals,2)),'UniformOutput',false)),'stacked');
            xticks(1:length(mg)); xticklabels(1:length(mg));
        end
    end
end
%saveFigures(gcf,"S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\Weights\","Hand-Laminar",[]);
%%
somaLabels = somaTable(mv);
[~,peakTimes] = max(firingRatesAverage(somaIndex,1,:),[],3);
[~,sI] = sort(location(somaIndex,1)*ImagingParameters.px2mm);
somaLabels = somaLabels(somaIndex);
mvv = [find(somaLabels(sI)=="Arm",min(groupcounts(somaLabels)));find(somaLabels(sI)=="Hand",min(groupcounts(somaLabels)))];
WS = W; WS(somaLabels(mvv)=="Arm",:) = W(somaLabels(mvv)=="Hand",:); WS(somaLabels(mvv)=="Hand",:) = W(somaLabels(mvv)=="Arm",:);
VS = V; VS(somaLabels(mvv)=="Arm",:) = V(somaLabels(mvv)=="Hand",:); VS(somaLabels(mvv)=="Hand",:) = V(somaLabels(mvv)=="Arm",:);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 2,...
    'legendSubplot', {16,params.condNames},'ylims',[]);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), WS, VS, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 2,...
    'legendSubplot', {16,params.condNames},'ylims',[]);
Xcen = bsxfun(@minus, firingRatesAverage(somaIndex,:)', mean(firingRatesAverage(somaIndex,:)'));
corr( Xcen * W, Xcen*WS);
%%
saveFig = false;
somaColors = containers.Map(["Arm","Hand"],{[0 1 .2],[.8 0 1]});
savePath = "S:\Lab\ngc14\Working\DataHigh\Centered\Demixed\";
plotSoma = "Arm"; plotLaminar = "";pDims = [1:10];
if(strcmp(plotSoma,"Hand"))
    ls = ':';Z =ZHand; whichMarg = whichMargHand;
    if(strcmp(plotLaminar,"Deep"))
        Z = ZHandD; whichMarg=whichMargHandD;
    elseif(strcmp(plotLaminar,"Superficial"))
        Z = ZHandS; whichMarg=whichMargHandS;
    end
else
    ls='-'; Z =ZArm; whichMarg = whichMargArm;
    if(strcmp(plotLaminar,"Deep"))
        Z = ZArmD; whichMarg=whichMargArmD;
    elseif(strcmp(plotLaminar,"Superficial"))
        Z = ZArmS; whichMarg=whichMargArmS;
    end
end
projT = Z(:,cell2mat(arrayfun(@(f) find(whichMarg==f,dims),1:length(combinedParams),'UniformOutput',false)));
%somaDist = reshape(transpose((projTArm - projTHand).^2),length(time),length(conditions),dims,[]);sqrt(sum(somaDist(:,c,:,nc),'all'))
projT = reshape(projT',dims*length(combinedParams),length(conditions),[]);
for nc = 1:length(combinedParams)
    figure(6+nc); st=(nc-1)*dims;
    for c = 1:length(conditions)
        lineColor = hsv2rgb(rgb2hsv(cell2mat(colors.values(cellstr("ArmHand_"+params.condAbbrev.values(cellstr(conditions(c)))))))-...
            [0 0 .4*strcmp(plotLaminar,"Deep")]);
        subplot(1,3,c); hold on; view(10,-15); xlim([-2 2.5]); ylim([-2 1.5]);
        [~,z,tr] = procrustes(squeeze(projTArm(st+[1:dims],c,:))',squeeze(projTHand(st+[1:dims],c,:))',"scaling",false,"reflection",'best');
        %zt = tr.b*squeeze(projTHand(st+[1:dims],c,:))'*tr.T + tr.c; zt=NaN(size(zt));z=zt;z(:,[pDims,find(~ismember(1:3,pDims))])=zt(:,1:3);
        title(params.condAbbrev(conditions(c))); %+ ", dist: "+ num2str(procrustes(squeeze(projTArm(st+pDims,c,:))',squeeze(projTHand(st+pDims,c,:))',"scaling",false,"reflection",'best'),'%.2f');
        %plot3(z(:,1),z(:,2),z(:,3),'Color',lineColor,'LineStyle',ls,'LineWidth',1.5); scatter3(z(1,1),z(1,2),z(1,3),'black','x','sizeData',150);
        plot3(squeeze(projT(st+1,c,:)),squeeze(projT(st+2,c,:)),squeeze(projT(st+3,c,:)),'Color',lineColor,'LineStyle',ls,'LineWidth',1.5);
        scatter3(projT(st+1,1),projT(st+2,1),projT(st+3,1),'black','*','sizeData',250);
        arrayfun(@(e) scatter3(projT(st+1,c,e),projT(st+2,c,e),projT(st+3,c,e),'filled','MarkerFaceColor',lineColor),round(mean(segVals{c}(:,1:3),1,'omitnan')))
    end
    % XY:view(0,90); XZ:view(0,0); YZ:view(90,0);
    if(saveFig)
        saveFigures(figure(nc),savePath,"Traj-"+margNames(nc),[]);
        %arrayfun(@(a)view(a,0,90),gcf().Children); saveFigures(figure(nc),savePath,"2DTraj-"+margNames(nc),[]);
    end
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