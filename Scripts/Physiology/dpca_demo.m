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
sTrials = 20;
binWidth = 10;
smoothWin = 150;
dims = 20;
time= -.5:binWidth/1000:2.5;
combinedParams = {{1,[1 3]},{2,[2,3]},{3},{[1,2],[1,2,3]}};
margNames = {'Condition','Somatotopy','Movement-Invariant','Cond/Soma Interaction'};
margColours = [23 100 171; 200 160 43; 150 150 150; 180 25 180]/256;
conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(conditions,.001,.001,-6,3,containers.Map(conditions,repmat({"StartReach"},1,length(conditions))));
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
tPhys = unitTable(conditions,params);
%%
tableInds = contains(string(tPhys.Somatotopy),["Arm","Hand"]);
somaTable = tPhys{tableInds,"Somatotopy"};
allLocations = tPhys{tableInds,["XT","YT"]};
allSegs= tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+params.condAbbrev.values)};%
tablePSTHD= tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+params.condAbbrev.values)};
clear tPhys
numUnits = all(cell2mat(cellfun(@(s) size(s,2), tablePSTHD,'UniformOutput',false))>=sTrials,2);%
somaTable = somaTable(numUnits);
allLocations = allLocations(numUnits,:);
allSegs = allSegs(numUnits,:);
%%
segInds = cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),["GoSignal","StartReach","StartHold","StartWithdraw"])),cellfun(@cell2mat,...
    num2cell(cellfun(@(aa) findBins(mean(aa,1,'omitnan'),params.bins),allSegs,'UniformOutput',false),1),'Uniformoutput',false),'UniformOutput',false);
currD= cellfun(@(v) permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(resize(max(0,d),[size(d,1),max(sTrials*2,size(d,2))],'FillValue',NaN),sTrials)',...%downsampleTrials(max(0,d),sTrials),...
    v,'Uniformoutput',false),1,1,[])),[3 1 2]),num2cell(tablePSTHD(numUnits,:),1), 'UniformOutput',false);
currD = squeeze(cellfun(@squeeze,num2cell(cellfun(@squeeze,num2cell(cat(4,currD{:}),[2,3]),'UniformOutput',false),4),'UniformOutput',false));
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
trialPSTH = cell(1,length(currD));
trialLength = floor(size(currD{1}{1}, 2) / binWidth);
for n = 1:length(currD)
    trialPSTH{n} = repmat({NaN(sTrials,trialLength)},length(params.condNames),1);
    for c = 1:length(trialPSTH{n})
        for t = 1:trialLength
            iStart = binWidth * (t-1) + 1;
            iEnd   = binWidth *t;
            trialPSTH{n}{c}(:,t) = sum(currD{n}{c}(:,iStart:iEnd),2);
        end
    end
end
clear currD;
%%
mv = sum(cell2mat(cellfun(@(m)mean(cell2mat(m'),2,'omitnan').*100>1,num2cell([trialPSTH{:}],1),'UniformOutput',false)),1)>sTrials/2 | ...
    mean(cell2mat(cellfun(@(n)mean(cat(2,n{:}),2,'omitnan'),trialPSTH,'UniformOutput',false)),1,'omitnan').*100>1;
firingRates = cellfun(@(c) cellfun(@(s) (conv2(resize(s(:,fix(findBins(time,params.bins)/binWidth)),...
    [size(s,1),length(time)+length(gausswin(ceil(smoothWin/binWidth)))-1],'Pattern','edge','side','both'),...
    transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),trialPSTH,'UniformOutput',false);
%+(rand(size(c)).*(std(c,0,2)*std(c,0,1)))
firingRates = cellfun(@(f) [f,cellfun(@(c) reshape(c(randperm(numel(c))),size(c)),f,'UniformOutput',false)],firingRates, 'UniformOutput',false);
firingRates(somaTable=="Hand") = cellfun(@(f) fliplr(f), firingRates(somaTable=="Hand"),'UniformOutput',false);

firingRates = reshape(firingRates(mv),[ones(1,sum(size(firingRates{1})~=1)),sum(mv)]);
firingRates = cat(length(size(firingRates))+sum(size(firingRates{1})~=1),firingRates{:});
firingRates = permute(cell2mat(permute(firingRates,[3 4 1 2 5])),[5 3 4 2 1]);
%firingRates = cell2mat(cellfun(@(r) circshift(r,randi([2*binWidth,size(r,3)-2*binWidth],1),3), num2cell(firingRates,3),'UniformOutput',false));
trialNum = ones(size(firingRates,[1:3])).*size(firingRates,length(size(firingRates)));
for n = 1:size(firingRates,1)
    for c = 1:size(firingRates,2)
        for d = 1:size(firingRates,3)
            assert(isempty(find(isnan(firingRates(n,c,d,:,1:trialNum(n,c))), 1)), 'Something is wrong!')
        end
    end
end
% firingRatesAverage = cell2mat(cellfun(@(r) reshape(cell2mat(r),size(r{1},1),1,[]),cellfun(@(s) cellfun(@(t) cell2mat(cellfun(@(n)circshift(n,randi([2*binWidth,length(n)-2*binWidth],1)),n,num2cell(t(mv,unique(round(ms_bins./binWidth))),2),'UniformOutput',false)),s,'UniformOutput',false)',trialPSTH, 'UniformOutput',false),'UniformOutput',false));
firingRatesAverage = mean(firingRates, length(size(firingRates)),'omitnan');
firingRatesAverage(isnan(firingRatesAverage)) = 0;
firingRates(isnan(firingRates)) = 0;
timeEventConds = cell2mat(cellfun(@(i) mean(findBins(params.bins(i),time),1,'omitnan'),...
    cellfun(@(s) fix(s(mv,~all(isnan(s),1))),segInds,'UniformOutput',false)','UniformOutput',false));
timeEvents = time(round([mean(timeEventConds(:,1:2),1,'omitnan'),median(timeEventConds(:,3:4),1,'omitnan')]));
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
% which provides substantial regularization itself (even with lambda=0). %
somaIndex = cell2mat(arrayfun(@(a) find(somaTable(mv)==a,min(groupcounts(somaTable(mv)))),unique(somaTable(mv)),'UniformOutput',false));
somaIndex = contains(string(somaTable(mv)),["Arm","Hand"]);
optimalLambda = dpca_optimizeLambda(firingRatesAverage(somaIndex,:,:,:),firingRates(somaIndex,:,:,:,:),...
    trialNum(somaIndex,:,:),'combinedParams', combinedParams, 'simultaneous', false,'numRep', 10);
Cnoise = dpca_getNoiseCovariance(firingRatesAverage(somaIndex,:,:,:), ...
    firingRates(somaIndex,:,:,:,:), trialNum(somaIndex,:,:), 'simultaneous', false,'type','averaged');
[W,V,whichMarg] = dpca(firingRatesAverage(somaIndex,:,:,:),dims*(length(combinedParams)+1),...
    'combinedParams', combinedParams,'lambda', optimalLambda,'Cnoise', Cnoise);
explVar = dpca_explainedVariance(firingRatesAverage(somaIndex,:,:,:), W, V, 'combinedParams', combinedParams);
dpca_plot(firingRatesAverage(somaIndex,:,:,:), W, V, @dpca_plot_default, ...
    'explainedVar', explVar,'marginalizationNames', margNames, 'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,'time', time,'timeEvents', timeEvents,'timeMarginalization', 3,...
    'legendSubplot', {16,params.condNames},'ylims',[]);
Z =  bsxfun(@minus, firingRatesAverage(somaIndex,:)', mean(firingRatesAverage(somaIndex,:),[2,3])')* W;
%% 2. Transform / Project TRIAL-AVERAGED PETHs (X) into dPCA space
targMarg =2; targNum =2;
colors = [1 0 0; 1 .7 0; 0 0 1];
phaseWindowSz = 0.2;
phaseAlign = num2cell(find(contains(maxSegL,["GoSignal","StartReach","StartHold","StartWithdraw"])));
phaseWin = {[0, phaseWindowSz],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)],[-phaseWindowSz*(5/4), -phaseWindowSz*(1/4)],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)]};

sizeX = size(firingRates(somaIndex,:,:,:)); %,
sizeX(1) = min(groupcounts(somaTable(mv)));
xfull_Arm = reshape(firingRates(find(somaTable(mv)=="Arm",sizeX(1)),:,:,:), sizeX(1), []);
Z_trials_Arm = reshape(WArm' * xfull_Arm, [],sizeX(2),length(time),sizeX(end)); % [Components, Stimulus, Time, Trials]
xfull_Hand = reshape(firingRates(find(somaTable(mv)=="Hand",sizeX(1)),:,:,:), sizeX(1), []);
Z_trials_Hand = reshape(WHand' * xfull_Hand, [],sizeX(2),length(time),sizeX(end)); % [Components, Stimulus, Time, Trials]

var_own = trace(cov([WArm(:,1:15)' * (xfull_Arm-mean(xfull_Arm,2,'omitnan'))]'));
var_cross= trace(cov([WHand(:,1:15)' * (xfull_Hand-mean(xfull_Hand,2,'omitnan'))]'));
alignmentIndex = var_cross / var_own;

figure; hold on;
armTarget=find(whichMargArm==targMarg,targNum);
handTarget=find(whichMargHand==targMarg,targNum);
armPhaseCond = [];
handPhaseCond = [];
for c = 1:sizeX(2)
    armSegs = mean(cell2mat(reshape(cellfun(@(s) downsampleTrials(s',sTrials)',allSegs(somaTable(mv)=="Arm",c),'UniformOutput',false),1,1,[])),3,'omitnan');
    handSegs = mean(cell2mat(reshape(cellfun(@(s) downsampleTrials(s',sTrials)',allSegs(somaTable(mv)=="Hand",c),'UniformOutput',false),1,1,[])),3,'omitnan');
    armPhaseCond(:,:,c) = cell2mat(cellfun(@(p,w) abs(squeeze(mean(Z_trials_Arm(armTarget(end),c,findBins(armSegs(:,p),time)+(w.*(1000/binWidth)),:),3))), phaseAlign, phaseWin,'UniformOutput',false));
    handPhaseCond(:,:,c) = cell2mat(cellfun(@(p,w) abs(squeeze(mean(Z_trials_Hand(handTarget(end),c,findBins(handSegs(:,p),time)+(w.*(1000/binWidth)),:),3))), phaseAlign, phaseWin,'UniformOutput',false));
    for t = 1:sTrials
        currArm = abs(squeeze(Z_trials_Arm(armTarget(end), c, :, t)));
        currHand = abs(squeeze(Z_trials_Hand(handTarget(end), c, :, t)));
        plot(time, currArm, 'Color', [colors(c, :), .5], 'LineWidth', 1,'LineStyle',':');
        plot(time, currHand, 'Color', [colors(c, :), .3], 'LineWidth', 1,'LineStyle','--');
    end
end
armTable = table(reshape(armPhaseCond,[],1),reshape(repmat("Arm",sTrials,length(phaseAlign),length(params.condNames)),[],1),...
    reshape(repmat(maxSegL(cell2mat(phaseAlign)),sTrials,length(params.condNames)),[],1),reshape(repelem(params.condAbbrev.values,sTrials,length(phaseAlign)),[],1),'VariableNames',["PCA","Somatotopy","Phase","Cond"]);
handTable = table(reshape(handPhaseCond,[],1),reshape(repmat("Hand",sTrials,length(phaseAlign),length(params.condNames)),[],1),...
    reshape(repmat(maxSegL(cell2mat(phaseAlign)),sTrials,length(params.condNames)),[],1),reshape(repelem(params.condAbbrev.values,sTrials,length(phaseAlign)),[],1),'VariableNames',["PCA","Somatotopy","Phase","Cond"]);
fullTable = [addvars(armTable,reshape(repelem([1:sTrials]',1,length(phaseAlign),length(params.condNames)),[],1),'NewVariableNames',"Trial");...
    addvars(handTable,reshape(repelem([1:sTrials]',1,length(phaseAlign),length(params.condNames)),[],1),'NewVariableNames',"Trial")];
fullTable = unstack(removevars(addvars(fullTable,string(fullTable.Somatotopy)+"-"+string(fullTable.Trial),string(fullTable.Phase)+"-"+string(fullTable.Cond),'NewVariableNames',{'Ind','P_C'}),{'Cond','Phase'}),"PCA",'P_C');
rmm = fitrm(fullTable,strjoin(fullTable.Properties.VariableNames(contains(fullTable.Properties.VariableNames,"_")),", ")+"~ 1+ Somatotopy",...
    'WithinDesign',table(repelem(maxSegL(cell2mat(phaseAlign))',length(params.condNames),1),repmat(params.condNames',length(phaseAlign),1),VariableNames={'Phase','Cond'}));
tbl = ranova(rmm, 'WithinModel', 'Phase*Cond');
tbl_posthoc = multcompare(rmm, 'Somatotopy', 'By', 'Phase','ComparisonType','Bonferroni')
%sizeX(1) = max(groupcounts(somaTable(mv)));
Z_avg_Arm = reshape(WArm' * reshape(mean(firingRates(find(somaTable(mv)=="Arm",min(groupcounts(somaTable(mv)))),:,:,:), 4, 'omitnan'), sizeX(1), []), [size(W,2), sizeX(2), length(time)]);
Z_avg_Hand = reshape(WHand' * reshape(mean(firingRates(find(somaTable(mv)=="Hand",min(groupcounts(somaTable(mv)))),:,:,:), 4, 'omitnan'), sizeX(1), []), [size(W,2), sizeX(2), length(time)]);
for c = 1:sizeX(2)
    h=plot(time, abs(squeeze(Z_avg_Arm(armTarget(end), c, :))), 'Color', max(colors(c, :)-[.01 .01 .01],0), 'LineWidth', 3);
    h=plot(time, abs(squeeze(Z_avg_Hand(handTarget(end), c, :))), 'Color', max(colors(c, :)-[.2 .2 .2],0), 'LineWidth', 3,'LineStyle','--');
end
yspan = get(gca,'YLim');
arrayfun(@(p) plot([p,p], [yspan(1),yspan(1) + range(yspan)], 'k--', 'LineWidth', 1), timeEvents);
xlabel('Time');
ylabel(["Z-score"]);
title(['Single-Trial Trajectories for ' margNames{targMarg}, ' Component ',  num2str(targNum)]);
legendEntries = cellfun(@(c) arrayfun(@(s) plot([NaN,NaN],[NaN,NaN],'LineWidth',3,'LineStyle',s,'Color',c),["-","--"]),num2cell(colors,2),'UniformOutput',false);
legend([legendEntries{:}], reshape(["Arm","Hand"]'+"="+params.condAbbrev.values,[],1), 'Location', 'best');
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
location = allLocations(mv,:);
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
somaLabels = somaTable(mv);
somaLabels = somaLabels(somaIndex);

if(0)
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
end

Z =  bsxfun(@minus, reshape(firingRatesAverage(:,:),length(somaLabels),[])', ...
    mean(firingRatesAverage(:,:)',1,'omitnan'))* W;

projT = Z(:,cell2mat(arrayfun(@(f) find(whichMarg==f,dims),unique(whichMarg),'UniformOutput',false)));
%somaDist = reshape(transpose((projTArm - projTHand).^2),length(time),length(conditions),dims,[]);sqrt(sum(somaDist(:,c,:,nc),'all'))
projT = reshape(projT',dims,length(unique(whichMarg)),length(conditions),[],length(time));
plotSeg = findBins(mean(cell2mat(reshape(cellfun(@(m) mean(m(:,contains(maxSegL,...
    ["GoSignal","StartReach","StartHold"])),1,'omitnan'), allSegs, 'UniformOutput',false),[],1)),1,'omitnan'),time);
for t = 2:length(plotSeg)
    figure();
    tiledlayout(3,length(combinedParams),'TileIndexing','columnmajor');
    for nc = 1:length(combinedParams)
        if(0)
            figure(6+nc); st=(nc-1)*dims;
            for c = 1:length(conditions)
                lineColor = hsv2rgb(rgb2hsv(colors(c,:)));
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
        for i = 1:3
            nexttile();
            hold on;
            if(i==1)
                view(0,90);
                title(margNames{nc});
            elseif(i==2)
                view(0,0);
            elseif(i==3)
                view(90,0);
            end
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            splitProj = squeeze(cellfun(@(s) squeeze(s),num2cell(squeeze(projT(1:3,nc,:,:,:)),[1 3 4]),'UniformOutput',false));
            for c = 1:length(splitProj)
                plot3(squeeze(splitProj{c}(1,1,1:plotSeg(t))),squeeze(splitProj{c}(2,1,1:plotSeg(t))),squeeze(splitProj{c}(3,1,1:plotSeg(t))),'Color',colors(c,:),'LineStyle','-','LineWidth',1.5);
                plot3(squeeze(splitProj{c}(1,2,1:plotSeg(t))),squeeze(splitProj{c}(2,2,1:plotSeg(t))),squeeze(splitProj{c}(3,2,1:plotSeg(t))),'Color',colors(c,:),'LineStyle','--','LineWidth',1.5);
            end
            for c = 1:length(splitProj)
                scatter3(squeeze(splitProj{c}(1,1,1)),squeeze(splitProj{c}(2,1,1)),squeeze(splitProj{c}(3,1,1)),'black','o','sizeData',150,'MarkerFaceColor','k');
                scatter3(squeeze(splitProj{c}(1,2,1)),squeeze(splitProj{c}(2,2,1)),squeeze(splitProj{c}(3,2,1)),'black','sq','sizeData',150,'MarkerFaceColor','k');
            end
        end
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