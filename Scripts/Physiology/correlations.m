conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[0, 0]}},1,length(conditions));
plotUnits = false;
MIN_BLOCKS_FOR_UNIT = 13;
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
pVal=0.05;
saveDir = "S:\Lab\ngc14\Working\PMd\Task_Units\";
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
%%
[spikes,times,weights,currTrials,sessionConds,channels,events,labels,chMap] =...
    getAllTrials("S:\Lab\Gilligan\All Data\Gilligan_08_30_2019\Physiology\Results","Single",true);
[currSeg,currTrialPSTHS,currActive,trialHists,alignedSpikes] = deal(repmat({[]},1,length(conditions)));
numUnits = size(spikes,2);
NaNGraspInd = false(size(currTrials,1),1);
NaNGraspInd(ismember(currTrials(:,1),conditions)) = cellfun(@(s,t) sum(isnan(t))==1 & isempty(t(strcmp(params.condSegMap(s),"StartGrasp"))), ...
    currTrials(ismember(currTrials(:,1),conditions),1),times(ismember(currTrials(:,1),conditions))');
successfulInds = cellfun(@(b,t) ~isnan(str2double(b)) | isempty(b),currTrials(:,end-1));
successfulTrials = successfulInds | NaNGraspInd;
currTrials = currTrials(successfulTrials,:);
times = times(successfulTrials);
for c = 1:length(conditions)
    currCond = conditions{c};
    condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1));
    condEvents = params.condSegMap(currCond);
    condAlign = cellfun(@(a) find(strcmp(condEvents,a)),...
        params.PSTHAlignments(currCond),'UniformOutput', false);
    % generate and smooth PSTHS  from current session
    if(any(condInds) && ~isempty(spikes))
        % {alignedPSTHS}{units,trials}
        alignedSpikes(c) = cellfun(@(ap) cellfun(@(st) cellfun(@(s,t) ...
            s-t(ap),st(condInds),times(condInds),'UniformOutput',false),...
            spikes,'UniformOutput',false),condAlign,'UniformOutput',false);
        alignedTimes = cellfun(@(ap) cell2mat(cellfun(@(t)[t(1:end-1)-t(ap), ...
            NaN(1,length(condEvents)-length(t)),t(end)-t(ap)], times(condInds),...
            'UniformOutput', false)'), condAlign, 'UniformOutput', false);
        trialHists{c} = cellfun(@(ac) cellfun(@(a) histcounts(a,...
            [params.bins(1)-(params.sigmaSize/2):params.binSize:params.bins(end)+...
            (params.sigmaSize/2)]),ac,'UniformOutput',false),alignedSpikes{c},'UniformOutput', false);
        %{alignedPSTHS}{units,trials}{bins}
        smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,gausswin(params.sigma)/...
            sum(gausswin(params.sigma)),'valid')./params.binSize,ac,'UniformOutput', false),...
            trialHists{c},'UniformOutput',false);
        %{alignedPSTHS}{units,bins,trials}
        unitTrialPSTH = cellfun(@(s) ... %condWeights.*
            cat(3,s{:}),smoothHists,'UniformOutput', false);
        % aligned trial times for the session for each (1) PSTH
        currSeg{c} = alignedTimes;
        % PSTH(units X bins X trials) for each alignment
        currTrialPSTHS{c} = cat(1,unitTrialPSTH{:});
    else
        % pad stored info with empty arrays and NaN pad indicies for missing conditions
        currSeg{c} = repmat({NaN(size(params.condSegMap(currCond)))},1,length(condAlign));
        currTrialPSTHS{c} = repmat(NaN(numUnits,length(bins),1),1,length(condAlign));
        alignedSpikes{c} = repmat(NaN(numUnits,1),1,length(condAlign));
    end
    corrMatrix = []; corrTrialMatrix = [];
    for u = 1:numUnits
        for i = 1:numUnits
            corrMatrix(u,i) = corr(mean(currTrialPSTHS{c}(u,:,:),3,'omitnan')',...
                mean(currTrialPSTHS{c}(i,:,:),3,'omitnan')');
            for t = 1:size(currTrialPSTHS{c},3)
                corrTrialMatrix(u,i,t) = corr(currTrialPSTHS{c}(u,:,t)',currTrialPSTHS{c}(i,:,t)');
            end
        end
    end
    figure();
    subplot(1,2,1);
    imagesc(corrMatrix);
    colormap('jet');
    colorbar;
    xticklabels(chMap{end}(channels));yticklabels(chMap{end}(channels));
    clim([0 1]);
    subplot(1,2,2);
    imagesc(mean(corrTrialMatrix,3,'omitnan'));
    colormap('jet');
    colorbar;
    xticklabels(chMap{end}(channels));yticklabels(chMap{end}(channels));
    clim([0 1]);
end