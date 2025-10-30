function [alignedSpikes, alignedTimes, allGoodTrials, condPSTHS] = getConditionTrialPSTHS(params,trialInfo,taskAlign,spikes,times)
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
conditions = params.condNames;
numTrials = min(cellfun(@(cn) sum(cellfun(@(ct) contains(ct,cn),trialInfo(:,1))),conditions));
for c = 1:length(conditions)
    currCond = conditions{c};
    taskCondInds = cell2mat(taskAlign(currCond));
    condInds = contains(trialInfo(:,1),currCond) & cumsum(contains(trialInfo(:,1),currCond))<=numTrials;
    condEvents = params.condSegMap(currCond);
    condEvents = condEvents(2:end-1);
    condAlign = cellfun(@(a) find(strcmp(condEvents,a)),...
        params.PSTHAlignments(currCond),'UniformOutput', false);
    if(any(condInds) && ~isempty(spikes))
        alignedSpikes(c) = cellfun(@(ap) cellfun(@(st) cellfun(@(s,t) ...
            s-t(ap),st(condInds),times(condInds),'UniformOutput',false),...
            spikes,'UniformOutput',false),condAlign,'UniformOutput',false);
        alignedTimes(c) = cellfun(@(ap) cellfun(@(t)[t(1:end-1)-t(ap), ...
            NaN(1,length(condEvents)-length(t)),t(end)-t(ap)], times(condInds),...
            'UniformOutput', false), condAlign, 'UniformOutput', false);
        trialHists = cellfun(@(ac) cellfun(@(a) histcounts(a,...
            [params.bins(1)-(params.sigmaSize/2):params.binSize:params.bins(end)+...
            (params.sigmaSize/2)]),ac,'UniformOutput',false),alignedSpikes{c},'UniformOutput', false);
        smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,gausswin(params.sigma)/...
            sum(gausswin(params.sigma)),'valid')./params.binSize,ac,'UniformOutput', false),trialHists,'UniformOutput',false);
        taskCondInds = cellfun(@(a) findBins(a(ismember(allSegs{c},taskCondInds)),params.bins),alignedTimes{c},'UniformOutput',false);
        allGoodTrials{c} = cellfun(@(h) cellfun(@(u,s) mean(u(s(1):s(end))) > 1*params.binSize*(s(end)-s(1))...
            && mean(u(s(1):s(end))) < 200*params.binSize*(s(end)-s(1)),h,taskCondInds),smoothHists,'UniformOutput',false);
        unitTrialPSTH = cellfun(@(s,a) permute((a./a)'.*... %condWeights.*
            cat(1,s{:}),[3 2 1]),smoothHists,allGoodTrials{c},'UniformOutput', false);
        condPSTHS{c} = cat(1,unitTrialPSTH{:});
    end
end
end