conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
pVal=0.05;
mainDir = "S:\Lab\";
monkey = "Gilligan";
sessionDate = "08_05_2019";
saveDir = strcat(mainDir,"ngc14\Working\Correlations\",sessionDate,"\");
saveFig = true;
evalWindow = [-1 2];
winSz = .2;
dirPath = mainDir+monkey+"\"+"All Data\"+monkey+"_"+sessionDate+"\Physiology\Results";
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
goInd = cellfun(@(as) find(contains(as,"Go")),allSegs,'UniformOutput',false);
reachInd = cellfun(@(as) find(contains(as,"Reach")),allSegs,'UniformOutput',false);
graspInd = cellfun(@(as) find(contains(as,"tHold"),1),allSegs,'UniformOutput',false);
binInds = findBins(evalWindow,params.bins);
close all;
%%
[condPSTHS,alignedSpikes,alignedTimes,allGoodTrials,corrTrialMatrix,corrMatrix] = deal(repmat({[]},1,length(conditions)));
[spikes,times,weights,allTrials,sessionConds,channels,events,labels,chMap] = getAllTrials(dirPath,"Single",true);
[~,~,unitNo] = unique(channels,'stable');
unitNames = arrayfun(@(s,u) strcat(num2str(s) ,"-",num2str(u)),chMap{end}(channels),1+[0;~diff(unitNo)]');
[unitNames,sIndex] = natsort(cellstr(unitNames));
spikes = spikes(sIndex);
NaNGraspInd = false(size(allTrials,1),1);
NaNGraspInd(ismember(allTrials(:,1),conditions)) = cellfun(@(s,t) sum(isnan(t))==1 & isempty(t(strcmp(params.condSegMap(s),"StartGrasp"))), ...
    allTrials(ismember(allTrials(:,1),conditions),1),times(ismember(allTrials(:,1),conditions))');
successfulTrials = NaNGraspInd | cellfun(@(b,t) ~isnan(str2double(b)) | isempty(b),allTrials(:,end-1));
allTrials = allTrials(successfulTrials,:);
times = times(successfulTrials);
for c = 1:length(conditions)
    currCond = conditions{c};
    taskCondInds = cell2mat(taskAlign(currCond));
    condInds = cellfun(@(a) contains(a,currCond),allTrials(:,1));
    condEvents = params.condSegMap(currCond);
    condEvents = condEvents(2:end-1);
    condAlign = cellfun(@(a) find(strcmp(condEvents,a)),...
        params.PSTHAlignments(currCond),'UniformOutput', false);
    siteTrialSegs{c} = num2cell(NaN(sum(condInds),length(maxSegL)),2);
    if(any(condInds) && ~isempty(spikes))
        for i = 1:sum(condInds)
            siteTrialSegs{c}{i}(condSegMappedInds{c}) = times{cumsum(condInds)==i & condInds==1}-...
                times{cumsum(condInds)==i & condInds==1}(condAlign{1});
        end
        % {alignedPSTHS}{units,trials}
        alignedSpikes(c) = cellfun(@(ap) cellfun(@(st) cellfun(@(s,t) ...
            s-t(ap),st(condInds),times(condInds),'UniformOutput',false),...
            spikes,'UniformOutput',false),condAlign,'UniformOutput',false);
        alignedTimes(c) = cellfun(@(ap) cellfun(@(t)[t(1:end-1)-t(ap), ...
            NaN(1,length(condEvents)-length(t)),t(end)-t(ap)], times(condInds),...
            'UniformOutput', false), condAlign, 'UniformOutput', false);
        trialHists = cellfun(@(ac) cellfun(@(a) histcounts(a,...
            [params.bins(1)-(params.sigmaSize/2):params.binSize:params.bins(end)+...
            (params.sigmaSize/2)]),ac,'UniformOutput',false),alignedSpikes{c},'UniformOutput', false);
        %{alignedPSTHS}{units,trials}{bins}
        smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,gausswin(params.sigma)/...
            sum(gausswin(params.sigma)),'valid')./params.binSize,ac,'UniformOutput', false),trialHists,'UniformOutput',false);
        taskCondInds = cellfun(@(a) findBins(a(ismember(allSegs{c},taskCondInds)),params.bins),alignedTimes{c},'UniformOutput',false);
        %{alignedPSTHS}{units,bins,trials}
        allGoodTrials{c} = cellfun(@(h) cellfun(@(u,s) mean(u(s(1):s(end))) > 1*params.binSize*(s(end)-s(1))...
            && mean(u(s(1):s(end))) < 200*params.binSize*(s(end)-s(1)),h,taskCondInds),smoothHists,'UniformOutput',false);
        unitTrialPSTH = cellfun(@(s,a) permute((a./a)'.*... %condWeights.*
            cat(1,s{:}),[3 2 1]),smoothHists,allGoodTrials{c},'UniformOutput', false);
        % PSTH(units X bins X trials) for each alignment
        condPSTHS{c} = cat(1,unitTrialPSTH{:});
    end
end
baselineAlign = cellfun(@(c,u) cell2mat(cellfun(@(t,i) max(1,mean(u(:,findBins(t(1)-3,params.bins):...
    findBins(t(1)-3,params.bins)+(1/params.binSize),i),2,'omitnan')),c,num2cell(1:size(u,3)),'UniformOutput',false)),alignedTimes,condPSTHS,'UniformOutput',false);
for b = 1:length(baselineAlign)
    condBaseline = baselineAlign{b};
    condBaseline(isnan(condBaseline)) = 1;
    baselineAlign{b} = condBaseline;
end
%%
normPSTH = cellfun(@(p,b) permute(permute(p,[1 3 2])./b,[1 3 2]),condPSTHS,baselineAlign,"UniformOutput",false);
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[0 0]}},length(conditions),1),...
   cellfun(@(a) {{cell2mat(a')}}, alignedTimes,'Uniformoutput',false),cellfun(@(c) {c},condPSTHS,'UniformOutput',false),false,true);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
taskUnits = any(cell2mat([tUnit{:}]),2);
channels = channels(taskUnits);
unitNo = unitNo(taskUnits);
unitNames = unitNames(taskUnits);
allGoodTrials = cellfun(@(a) a(taskUnits),allGoodTrials,'UniformOutput',false);
alignedSpikes = cellfun(@(a) a(taskUnits), alignedSpikes, 'UniformOutput',false);
condPSTHS = cellfun(@(a) a(taskUnits,:,:), condPSTHS, 'UniformOutput',false);
normPSTH = cellfun(@(c) c(taskUnits,:,:),normPSTH,'UniformOutput',false);
numTrials = min(cellfun(@(s) size(s,3), condPSTHS, 'UniformOutput',true));
figure();
condColors = {[1 0 0]; [.9 .7 0]; [0 .3 1]};
unitLabs = cell2mat(arrayfun(@(a) repmat(strcat("ChU_",replace(a,'-','_')),numTrials,1),unitNames,'UniformOutput',false)');
for c = 1:length(conditions)-1
    plotJointPSTHS(params.bins,{reshape(permute(normPSTH{c},[3 1 2]),[],length(params.bins),1)},...
        {repmat(cell2mat(vertcat(siteTrialSegs{1})),sum(taskUnits),1)},unitLabs,[any(diff(char(unitLabs)),2);1],...
        [],{evalWindow},[0 5],cell2struct(num2cell(repmat(condColors{c},sum(taskUnits),1),2),unique(unitLabs,'stable')));
end
if(saveFig)
    saveFigures(gcf,saveDir,"PSTHS",[]);
end
goSpk = cellfun(@(c,i,t,g) cellfun(@(a,gi) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(i)+winSz & s>at(i)),...
    a,t,num2cell(gi)),c,g,'UniformOutput',false),alignedSpikes(1:length(conditions)-1),...
    goInd(1:length(conditions)-1),alignedTimes(1:length(conditions)-1),allGoodTrials(1:length(conditions)-1),'UniformOutput',false);
reachSpk = cellfun(@(c,i,t,g) cellfun(@(a,gi) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(i) & s>at(i)-winSz),...
    a,t,num2cell(gi)),c,g,'UniformOutput',false),alignedSpikes(~cellfun(@isempty,reachInd)),...
    reachInd(~cellfun(@isempty,reachInd)),alignedTimes(~cellfun(@isempty,reachInd)),allGoodTrials(~cellfun(@isempty,reachInd)),'UniformOutput',false);
graspSpk = cellfun(@(c,i,t,g) cellfun(@(a,gi) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(i) & s>at(i)-winSz),...
    a,t,num2cell(gi)),c,g,'UniformOutput',false),alignedSpikes(~cellfun(@isempty,graspInd)),...
    graspInd(~cellfun(@isempty,graspInd)),alignedTimes(~cellfun(@isempty,graspInd)),allGoodTrials(~cellfun(@isempty,graspInd)),'UniformOutput',false);
goSpk = reshape(permute(cell2mat(permute(cat(3,goSpk{:}),[2 1 3])),[1 3 2]),sum(taskUnits),[]);
reachSpk = reshape(permute(cell2mat(permute(cat(3,reachSpk{:}),[2 1 3])),[1 3 2]),sum(taskUnits),[]);
graspSpk = reshape(permute(cell2mat(permute(cat(3,graspSpk{:}),[2 1 3])),[1 3 2]),sum(taskUnits),[]);
[goMatrix,reachMatrix,graspMatrix] = deal([]);
for u = 1:sum(taskUnits)
    figure();
    xVals = cellfun(@(m) min(m(u,:)):max(m(u,:)),{goSpk,reachSpk,graspSpk}, 'UniformOutput',false);
    for i = 1:sum(taskUnits)
        subplot(ceil(sum(taskUnits)/5),5,i);
        hold on;
        cm = cellfun(@(cp) corr(mean(cp(u,binInds(1):binInds(end),:),3,'omitnan')',...
            mean(cp(i,binInds(1):binInds(end),:),3,'omitnan')'),condPSTHS,'UniformOutput',false);
        ctm = {};
        for t = 1:numTrials
            ctm{t} = cellfun(@(cp) corr(cp(u,binInds(1):binInds(end),t)',...
                cp(i,binInds(1):binInds(end),t)'),condPSTHS,'UniformOutput',true);
        end
        ctm = cell2mat(ctm');
        for c = 1:length(conditions)
            corrMatrix{c}(u,i) = cm{c};
            corrTrialMatrix{c}(u,i,:) = ctm(:,c);
        end
        goMatrix(u,i) = corr(goSpk(u,:)',goSpk(i,:)',Rows='pairwise');
        reachMatrix(u,i) = corr(reachSpk(u,:)',reachSpk(i,:)',Rows='pairwise');
        graspMatrix(u,i) = corr(graspSpk(u,:)',graspSpk(i,:)',Rows='pairwise');
        s1=scatter(goSpk(u,:)',goSpk(i,:)','green','x','LineWidth',.5,'AlphaData',1,'SizeData',20);
        s2=scatter(reachSpk(u,:)',reachSpk(i,:)','magenta','+','LineWidth',.5,'AlphaData',1,'SizeData',20);
        s3=scatter(graspSpk(u,:)',graspSpk(i,:)','*','MarkerEdgeColor',[0 0 0],'LineWidth',.5,'AlphaData',1,'SizeData',20);
        plot(xVals{3},mean(graspSpk(i,isalmost(graspSpk(u,:),min(graspSpk(u,:)),1)),'omitnan')+(xVals{3}.*graspMatrix(u,i)),...
            'Color',[.1 .1 .1],'LineWidth',2.5);
        plot(xVals{2},mean(reachSpk(i,isalmost(reachSpk(u,:),min(reachSpk(u,:)),1)),'omitnan')+(xVals{2}.*reachMatrix(u,i)),...
            'Color',[.75 0 .75],'LineWidth',2.5);
        plot(xVals{1},mean(goSpk(i,isalmost(goSpk(u,:),min(goSpk(u,:)),1)),'omitnan')+(xVals{1}.*goMatrix(u,i)),...
            'Color',[0 .6 0],'LineWidth',2.5);
        arrayfun(@(s) set(s,'XJitter','rand'), [s1,s2,s3]);
        arrayfun(@(s) set(s,'YJitter','rand'), [s1,s2,s3]);
        arrayfun(@(s) set(s,'XJitterWidth',.85), [s1,s2,s3]);
        arrayfun(@(s) set(s,'YJitterWidth',.85), [s1,s2,s3]);
        title(strcat(unitNames(i),": ",num2str(min(sum(~isnan(goSpk(u,:))),sum(~isnan(goSpk(i,:))))),...
            " trials (",cell2mat(compose('%.2f; ',[goMatrix(u,i),reachMatrix(u,i),graspMatrix(u,i)])),")"));
        if(i==u)
            legend([s1,s2,s3],{'Go','Reach','Grasp'});
        end
    end
    if(saveFig)
        saveFigures(gcf,saveDir+"Unit_Correlations\",unitNames{u}+"_Correlations",[]);
    end
end
figure();
colormap('jet');
for c = 1:length(conditions)
    subplot(2,length(conditions)+2,c);hold on; axis image; axis ij;
    title(conditions(c));
    imagesc(corrMatrix{c});
    xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 1]);
end
avgMatrix = mean(cat(3,corrMatrix{1:length(conditions)-1}),3,'omitnan');
subplot(2,length(conditions)+2,length(conditions)+1); hold on; axis image; axis ij;
title("Average");
imagesc(avgMatrix);
xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 1]);
subplot(2,length(conditions)+2,length(conditions)+2); hold on; axis ij; axis tight;
title('Average marginals');
imagesc(mean(avgMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 1]);
% 2nd row
subplot(2,length(conditions)+2,(length(conditions)+2)+1);hold on; axis image; axis ij;
title('Go correlations');
imagesc(goMatrix);
xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
subplot(2,length(conditions)+2,(length(conditions)+2)+4); hold on; axis ij; axis tight;
title('Go marginals');
imagesc(mean(goMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim([0 .5]);
subplot(2,length(conditions)+2,(length(conditions)+2)+2);hold on; axis image; axis ij;
title('Reach correlations');
imagesc(reachMatrix);
xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
subplot(2,length(conditions)+2,(length(conditions)+2)+5); hold on; axis ij; axis tight;
title('Reach marginals');
imagesc(mean(reachMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim([0 .5]);
subplot(2,length(conditions)+2,(length(conditions)+2)+3); hold on; axis ij; axis image;
title('Grasp correlations');
imagesc(graspMatrix);
xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
subplot(2,length(conditions)+2,(length(conditions)+2)+6); hold on; axis ij; axis tight;
title('Grasp marginals');
imagesc(mean(graspMatrix.*(~logical(diag(true(sum(taskUnits),1)))./~logical(diag(true(sum(taskUnits),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 .5]);
fa = get(gcf,'Children');
cellfun(@(s) set(s,'FontSize',8,'TickLabelRotation',75), get(fa(strcmp(get(fa,'Type'),'axes')),'XAxis'));
if(saveFig)
    saveFigures(gcf,saveDir,"Correlations",[]);
end