date = '05_17_2019';
binSize=.01; % bin size in seconds
secondsBeforeEvent = .2;
secondsAfterEvent = .5;
alignSegsConds = {{'StartReach', 'StartGrasp','StartLift'},{'StartReach', 'StartGrasp','StartLift'},...
    {'StartReach','StartGrasp', 'StartHold'},{'GoSignal','StartReplaceHold','StartReward'}};
colors = {'r', 'g', 'b', 'k', 'm', 'y', 'c', 'k'};
colors = {'r', 'g', 'b', 'k', 'r', 'g', 'b', 'k'};
maxAlign = max(cellfun(@length, alignSegsConds));
sigma = 5; % smoothing window in ms

Kernel = (-3*sigma:3*sigma);
BinSize = length(Kernel);
Kernel = (-BinSize/2:BinSize/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

sessionDir = dir(['S:\Lab\Gilligan\All Data\Gilligan_',date,'\Physiology\Results']);
if(~exist(['S:\Lab\Gilligan\All Data\Gilligan_',date,'\Physiology\Results\PSTH'], 'dir'))
    mkdir(['S:\Lab\Gilligan\All Data\Gilligan_',date,'\Physiology\Results\PSTH']);
end
sessionDir = sessionDir(~[sessionDir.isdir]);
[~,sortInd] = natsort({sessionDir.name});
sessionDir = sessionDir(sortInd);
arduinoData = [];
spikeTimes = [];
segTimes = [];
for f = 1:length(sessionDir)
    load([sessionDir(f).folder, '\', sessionDir(f).name]);
    conds = sortedSpikeData.Conditions;
    locs = sortedSpikeData.Locations;
    spikeTimes = [spikeTimes; sortedSpikeData.SpikeTimes];
    segTimes = [segTimes; sortedSpikeData.SegTimes];
end
locs = unique(locs);
conds = unique(conds);
arduinoData = sortedSpikeData.ArduinoData;
figure('units','normalized','outerposition',[0 0 1 1])
for l = 1:length(locs)
    currLoc = locs{l};
    locTrials = find(cellfun(@(a) strcmp(a,currLoc),arduinoData(:,8)));
    for c = 1:length(conds)
        currCond = conds{c};
        alignSegs =  alignSegsConds{c};
        condTrials = find(cellfun(@(a) strcmp(a,currCond),arduinoData(:,1)));
        currTrialInds = intersect(condTrials,  locTrials);
        currSpikes = spikeTimes(:,currTrialInds);
        currSegs = segTimes(:,currTrialInds);
        for a = 1:length(alignSegs)
            ax = subplot(length(locs), length(alignSegs), length(alignSegs)*(l-1)+a);
            hold on;
            title(alignSegsConds{1}{a});
            alignInd = find(strcmp(sortedSpikeData.ConditionSegments{c}, alignSegs(a)));
            alignTimes = cellfun(@(a) getAlignedTimes(a,alignInd), currSegs, 'UniformOutput', false);
            alignedSpikes = cellfun(@minus, currSpikes, alignTimes, 'UniformOutput', false);
            trialHists = cellfun(@(a) histcounts(a,...
                -secondsBeforeEvent-length(Kernel)/200:binSize:...
                secondsAfterEvent+length(Kernel)/200)./binSize, alignedSpikes,'UniformOutput', false);
            trialHistsSmooth = cellfun(@(a) conv(a,Kernel'), trialHists, 'UniformOutput', false);
            trialHistsSmooth = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
            totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
            p = plot(ax,-secondsBeforeEvent:binSize:secondsAfterEvent,nanmedian(totalHists,1),...
                'Color',colors{length(conds)*(l-1)+c}, 'LineWidth', 2);
            if(~strcmp(conds{c},'Rest'))
                otherPhases = find(~strcmp(alignSegs(a), alignSegs));
                for ao=1:length(otherPhases)
                    alignIndOther = find(strcmp(sortedSpikeData.ConditionSegments{c}, alignSegs(otherPhases(ao))));
                    nextPhase = cellfun(@(a) getAlignedTimes(a,alignIndOther), currSegs, 'UniformOutput', false);
                    alignVals = cellfun(@minus, nextPhase,alignTimes, 'UniformOutput', false);
                    segPlot = nanmedian([alignVals{1,:}]);
                    plot([segPlot,segPlot], [min(p.YData) max(p.YData)], 'Color',colors{length(conds)*(l-1)+c});
                end
            end
            %errorbar
            %                     mainLineColor=get(p, 'color');
            %                     edgeColor=mainLineColor+(1-mainLineColor)*0.55;
            %                     uE=mean(totalHists,1)+std(totalHists,1);
            %                     lE=mean(totalHists,1)-std(totalHists,1);
            %                     %Make the patch
            %                     yP=[lE,fliplr(uE)];
            %                     xP=[-secondsBeforeEvent:binSize:secondsAfterEvent,fliplr(-secondsBeforeEvent:binSize:secondsAfterEvent)];
            %                     xP(isnan(yP))=[];
            %                     yP(isnan(yP))=[];
            %                     d = patch(xP,yP,1);
            %                     set(d,'facecolor',edgeColor,'edgecolor','none','facealpha',.5)
        end
    end
end
saveas(gcf,['S:\Lab\Gilligan\All Data\Gilligan_',date,'\Physiology\Results\PSTH\All'], 'fig');
saveas(gcf,['S:\Lab\Gilligan\All Data\Gilligan_',date,'\Physiology\Results\PSTH\All'], 'png');
close gcf;

function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end