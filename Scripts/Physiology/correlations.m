conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
plotUnits = false;
MIN_BLOCKS_FOR_UNIT = 13;
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
pVal=0.05;
mainDir = "S:\Lab\";
monkey = "Gilligan";
sessionDate = "08_05_2019";
evalWindow = [-1 2];
winSz = .1;
dirPath = mainDir+monkey+"\"+"All Data\"+monkey+"_"+sessionDate+"\Physiology\Results";
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
%%
[spikes,times,weights,currTrials,sessionConds,channels,events,labels,chMap] =...
    getAllTrials(dirPath,"Single",true);
[currTrialPSTHS,trialHists,alignedSpikes] = deal(repmat({[]},1,length(conditions)));
numUnits = size(spikes,2);
NaNGraspInd = false(size(currTrials,1),1);
NaNGraspInd(ismember(currTrials(:,1),conditions)) = cellfun(@(s,t) sum(isnan(t))==1 & isempty(t(strcmp(params.condSegMap(s),"StartGrasp"))), ...
    currTrials(ismember(currTrials(:,1),conditions),1),times(ismember(currTrials(:,1),conditions))');
successfulInds = cellfun(@(b,t) ~isnan(str2double(b)) | isempty(b),currTrials(:,end-1));
successfulTrials = successfulInds | NaNGraspInd;
currTrials = currTrials(successfulTrials,:);
times = times(successfulTrials);
[~,~,unitNo] = unique(channels,'stable');
unitNames = arrayfun(@(s,u) strcat(num2str(s) ,"_",num2str(u)),chMap{end}(channels),1+[0;~diff(unitNo)]');
[unitNames,sIndex] = natsort(cellstr(unitNames));
spikes = spikes(sIndex);
close all;
figure();
colormap('jet');
for c = 1:length(conditions)
    currCond = conditions{c};
    condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1));
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
        alignedTimes = cellfun(@(ap) cellfun(@(t)[t(1:end-1)-t(ap), ...
            NaN(1,length(condEvents)-length(t)),t(end)-t(ap)], times(condInds),...
            'UniformOutput', false), condAlign, 'UniformOutput', false);
        trialHists{c} = cellfun(@(ac) cellfun(@(a) histcounts(a,...
            [params.bins(1)-(params.sigmaSize/2):params.binSize:params.bins(end)+...
            (params.sigmaSize/2)]),ac,'UniformOutput',false),alignedSpikes{c},'UniformOutput', false);
        %{alignedPSTHS}{units,trials}{bins}
        smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,gausswin(params.sigma)/...
            sum(gausswin(params.sigma)),'valid')./params.binSize,ac,'UniformOutput', false),...
            trialHists{c},'UniformOutput',false);
        taskCondInds = taskAlign(conditions{c});
        taskCondInds = cellfun(@(n) cellfun(@(a) findBins(a(ismember(allSegs{c},taskCondInds{1})),...
            params.bins),n,'UniformOutput',false),alignedTimes,'UniformOutput',false);
        %{alignedPSTHS}{units,bins,trials}
        allGoodTrials = cellfun(@(h) cellfun(@(u,s) mean(u(s(1):s(end))) > 1*params.binSize*(s(end)-s(1))...
            && mean(u(s(1):s(end))) < 200*params.binSize*(s(end)-s(1)),h,taskCondInds{1}),smoothHists,'UniformOutput',false);
        unitTrialPSTH = cellfun(@(s,a) permute((a./a)'.*... %condWeights.*
            cat(1,s{:}),[3 2 1]),smoothHists,allGoodTrials,'UniformOutput', false);
        % PSTH(units X bins X trials) for each alignment
        currTrialPSTHS{c} = cat(1,unitTrialPSTH{:});
    else
        % pad stored info with empty arrays and NaN pad indicies for missing conditions
        currTrialPSTHS{c} = repmat(NaN(numUnits,length(bins),1),1,length(condAlign));
        alignedSpikes{c} = repmat(NaN(numUnits,1),1,length(condAlign));
    end
    reachInd = find(contains(allSegs{c},"Reach"));
    if(~isempty(reachInd))
        reachCorr = cellfun(@(a,g) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(reachInd) & s>at(reachInd)-winSz),...
            a,alignedTimes{1},num2cell(g)),alignedSpikes{c},allGoodTrials,'UniformOutput',false);
        graspInd = find(contains(allSegs{c},"Hold"),1);
        graspCorr = cellfun(@(a,g) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(graspInd) & s>at(graspInd)-winSz),...
            a,alignedTimes{1},num2cell(g)),alignedSpikes{c},allGoodTrials,'UniformOutput',false);
    else
        [reachCorr,graspCorr]= deal(cellfun(@(n) NaN(size(n)), allGoodTrials,'UniformOutput',false));
    end
    corrMatrix = []; corrTrialMatrix = [];
    binInds = findBins(evalWindow,params.bins);
    for u = 1:numUnits
        for i = 1:numUnits
            corrMatrix(u,i) = corr(mean(currTrialPSTHS{c}(u,binInds(1):binInds(end),:),3,'omitnan')',...
                mean(currTrialPSTHS{c}(i,binInds(1):binInds(end),:),3,'omitnan')');
            reachMatrix{c}(u,i) = corr(reachCorr{u}',reachCorr{i}',Rows='pairwise');
            graspMatrix{c}(u,i) = corr(graspCorr{u}',graspCorr{i}',Rows='pairwise');
            for t = 1:size(currTrialPSTHS{c},3)
                corrTrialMatrix(u,i,t) = corr(currTrialPSTHS{c}(u,binInds(1):binInds(end),t)',currTrialPSTHS{c}(i,binInds(1):binInds(end),t)');
            end
        end
    end
    condMatrix{c} = corrMatrix;
    subplot(3,length(conditions)+2,c);hold on; axis image; axis ij;
    title(conditions(c));
    imagesc(corrMatrix);
    xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 1]);
    subplot(3,length(conditions)+2,(length(conditions)+2)+c);hold on; axis image; axis ij;
    imagesc(reachMatrix{c});
    if(c==1)
        ylabel("Reach correlations");
    end
    xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
    subplot(3,length(conditions)+2,2*(length(conditions)+2)+c);hold on; axis image; axis ij;
    imagesc(graspMatrix{c});
    if(c==1)
        ylabel("Grasp correlations");
    end
    xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
    if(c==length(conditions))
        unitLabs = cell2mat(arrayfun(@(a) repmat(strcat("ChU_",a),size(currTrialPSTHS{c},3),1),unitNames,'UniformOutput',false)');
        avgMatrix = mean(cat(3,condMatrix{1:length(conditions)-1}),3,'omitnan');
        avgReach = mean(cat(3,reachMatrix{1:length(conditions)-1}),3,'omitnan');
        avgGrasp = mean(cat(3,graspMatrix{1:length(conditions)-1}),3,'omitnan');
        subplot(3,length(conditions)+2,c+1); hold on; axis image; axis ij;
        title("Average");
        imagesc(avgMatrix); 
        xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 1]);
        subplot(3,length(conditions)+2,(length(conditions)+2)+c+1);hold on; axis image; axis ij;
        imagesc(avgReach);
        xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
        subplot(3,length(conditions)+2,2*(length(conditions)+2)+c+1);hold on; axis image; axis ij;
        imagesc(avgGrasp);
        xticks(1:length(unitNames));xticklabels(unitNames);yticklabels([]);clim([0 .5]);
        subplot(3,length(conditions)+2,c+2); hold on; axis ij;
        imagesc(mean(avgMatrix.*(~logical(diag(true(numUnits,1)))./~logical(diag(true(numUnits,1)))),2,'omitnan'));
        xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 1]);
        subplot(3,length(conditions)+2,(length(conditions)+2)+c+2); hold on; axis ij;
        imagesc(mean(avgReach.*(~logical(diag(true(numUnits,1)))./~logical(diag(true(numUnits,1)))),2,'omitnan'));
        xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 .5]);
        subplot(3,length(conditions)+2,2*(length(conditions)+2)+c+2); hold on; axis ij;
        imagesc(mean(avgGrasp.*(~logical(diag(true(numUnits,1)))./~logical(diag(true(numUnits,1)))),2,'omitnan'));
        xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 .5]);
        plotJointPSTHS(params.bins,{reshape(permute(currTrialPSTHS{c},[3 1 2]),[],size(currTrialPSTHS{c},2),1)},...
            {repmat(mean(cell2mat(vertcat(siteTrialSegs{1:length(conditions)-1})),1,'omitnan'),numUnits*size(currTrialPSTHS{c},3),1)},...
            unitLabs,[any(diff(char(unitLabs)),2);1],[],{evalWindow},[0 15],...
            cell2struct(num2cell(distinguishable_colors(numUnits),2),unique(unitLabs,'stable')));
    end
end
