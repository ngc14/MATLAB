function [spikes,times,unitTrials,trials,conds,channel,eventNames,labels] = getSessionInfo2(folderName, singleOrAll)
sessionDir = dir(folderName+"\*.mat");
[~,sortInd] = natsort({sessionDir.name});
sessionDir = sessionDir(sortInd);
firstChannel = load([sessionDir(1).folder, '\', sessionDir(1).name]);
MIN_BLOCKS_FOR_UNIT = 15;
sessionCondSegs = firstChannel.sortedSpikeData.ConditionSegments;
conds = firstChannel.sortedSpikeData.Conditions;
events = containers.Map(conds, sessionCondSegs);
eventNames = events;
trials = firstChannel.sortedSpikeData.ArduinoData;
[allTimes, spikes] = deal({});
channel = [];
labels = string([]);
unitTrials = {};
for f = 1:length(sessionDir)
    l = load([sessionDir(f).folder, '\', sessionDir(f).name]);
    labeledData = isfield(l, 'label') || contains(fieldnames(l, '-full'),'label');
    if(labeledData)
        if(~isfield(l.sortedSpikeData, 'label'))
            labs = l.label;
        else
            labs = l.sortedSpikeData.label;
        end
        cSpikes = l.sortedSpikeData.SpikeTimes;
        segTimes = l.sortedSpikeData.SegTimes;
        unitTrials{f} = {};
        %exclude trials with neuron firing less than 1Hz or greater than 250 Hz
        if(any(~cellfun(@(a) length(a)==1,cSpikes(:))))
            if(size(cSpikes,1)<length(labs))
                labs = labs((length(labs)-size(cSpikes,1)+1):end);
            end
            allGoodTrials = ~(cellfun(@(a,b) length(a) <=(b(end)-b(1)) ...
                | length(a)>200*(b(end)-b(1)) | any(isnan(a)), cSpikes, segTimes));
% %             cSpikes(~allGoodTrials) = {NaN};
            blockInds = cumsum(mod(1:length(trials),length(conds))==1);
            unitTrials{f} = num2cell(allGoodTrials.*blockInds,2);
            goodUnitsOnChannel = (cellfun(@(tc) length(unique(tc))-1, unitTrials{f})>MIN_BLOCKS_FOR_UNIT)';
            if(strcmpi(singleOrAll, 'Single'))
                goodUnits = goodUnitsOnChannel & arrayfun(@(a) strcmp(a,'s'), labs);
            else
                goodUnits = goodUnitsOnChannel;
            end
            trialSegs = cellfun(@(e) values(events,{e}),repmat(trials(:,1)',sum(goodUnits),1), 'UniformOutput', false);
            labs = labs(goodUnits);

            segTimes = segTimes(goodUnits,:);
            trialSegs = cellfun(@(t) t{:}, trialSegs, 'UniformOutput', false);
            cSpikes = cSpikes(goodUnits,:);

            missGraspInds = cellfun(@(a,b) length(a)+1==length(b), segTimes,trialSegs);
            if(any(missGraspInds(:)))
                segTimes(missGraspInds) = cellfun(@(a,b) [b(1:find(strcmp(a,'StartGrasp'))-1),...
                    NaN, b(find(strcmp(a, 'StartGrasp')):end)],trialSegs(missGraspInds),...
                    segTimes(missGraspInds),'UniformOutput', false);
            end
            segTimes(cellfun(@length,trialSegs)~=cellfun(@length,segTimes)) = ...
                cellfun(@(n) NaN(1,length(n)), segTimes(cellfun(@length,trialSegs)...
                ~=cellfun(@length,segTimes)),'UniformOutput',false);

            allTimes(end+1:end+size(segTimes,1),:) = segTimes;
            channel(end+1:end+sum(goodUnits)) = f;
            labels(end+1:end+length(labs)) = labs;
            spikes(end+1:end+size(cSpikes,1),:) = cSpikes;
        end
    end
end
%%
times = {};
segInfoInds =  ~cell2mat(cellfun(@(a) all(arrayfun(@(b) isnan(b),a)) | length(a)==1, allTimes, 'UniformOutput', false));
[~, segInfoFirst] = sort(segInfoInds,1,'descend');
for segs=1:size(segInfoFirst,2)
    times{segs} = allTimes{segInfoFirst(1,segs),segs};
end
if(~isempty(times))
    condSegTrials = cellfun(@(e) values(events,{e}),trials(:,1)');
    badSegments = cellfun(@length, times)~=cellfun(@length,condSegTrials);
    if(any(badSegments))
        %disp(['NaNd segment times for trial(s): ', num2str(find(badSegments)), ' due to bad segmentation.']);
        times(badSegments) = cellfun(@(a) NaN(1,length(a)), condSegTrials(badSegments), 'UniformOutput', false);
        %         times(badSegments) = [];
        %         spikes(:,badSegments) = [];
        %         trials(badSegments,:) = [];
        %         % correct allTimes for weights correction
        %         allTimes(:,badSegments) = [];
    end
end
%% getting corresponding PSTH bins for grasp, reach, and task windows
weights = zeros(size(allTimes,1),length(conds));

for uN = 1:size(spikes,1)
    for c = 1:length(conds)
        condInds = strcmp(trials(:,1)', conds{c});
        weights(uN,c) = sum(cellfun(@(a) any(~isnan(a)), spikes(uN,condInds)));
    end
end
% ensure total number of trials per condition is calculated
% from correct trials
for c = 1:length(conds)
    weights(:,c) = weights(:,c)./max(max(weights(:,c),[],2));
end
