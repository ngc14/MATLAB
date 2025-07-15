function [siteDateMap,siteSegs,siteTrialPSTHS,rawSpikes,siteChannels,siteActiveInd, simpRep,...
    siteLocation, siteMasks, monkeys,vMask,conditions,channelMap] = getAllSessions(params,singleOrAllUnits,domain)
%  assign parameters
rawSpikes = [];
drivePath = "S:\Lab\";
monkeys = ["Gilligan", "Skipper"];
excludeRep = "Face";
% PSTH parameters: bin sizes, smoothing kernel, seconds prior to zero 
% alignment,seconds after zero alignment, alignment point(s) for each 
% condition (default value)
bins = params.bins;
conditions = cellstr(params.condNames);

siteDateMap = table();
allActivityMaps = containers.Map();
vMask = containers.Map();
for m = 1:length(monkeys)
    [monkeyTable, monkeyMask, monkeyActivityMaps] = getMonkeyInfo(drivePath,...
        monkeys(m),domain,true);
    siteDateMap = [siteDateMap; monkeyTable];
    allActivityMaps = [allActivityMaps; monkeyActivityMaps];
    vMask(monkeys(m)) = monkeyMask;
end
siteDateMap = siteDateMap(~cellfun(@isempty, siteDateMap.Date),:);
% load info from all sites
numSites = height(siteDateMap);
numSites=20;
[siteLocation, siteRep, siteThresh,siteSegs,siteChannels,...
    siteTrialPSTHS,siteActiveInd,rawSpikes,channelMap] = deal(cell(1,numSites));
hbar=parfor_progressbar(numSites,strcat("Iterating ", num2str(numSites), " instances..."));
parfor  i = 1:numSites
    currSession = siteDateMap(i,:);
    if(strcmp(currSession.Monkey,"Gilligan"))
        dateFormat = 'MM_dd_yyyy';
    else
        dateFormat = 'yyyy_MM_dd';
    end
    du = [currSession.Date{:}];
    du.Format = dateFormat;
    currSession.Date = du;
    physDir = strcat(drivePath,currSession.Monkey,"\All Data\", currSession.Monkey,...
        "_",string(currSession.Date),"\Physiology\Results_New\");
    if(~exist(physDir,'dir'))
        if(~ismember(currSession.Date,cellfun(@(d) datetime(d,'InputFormat','MM_dd_yyyy'),...
                {'05_02_2019'})))
            Spike_SortRawData(char(currSession.Date),char(currSession.Monkey));
            labelSingleUnits(char(currSession.Monkey),char(currSession.Date));
        end
    else
        firstChannelDir = dir(strcat(physDir,"*.mat"));
        firstChannelDir = load([firstChannelDir(1).folder,'\',firstChannelDir(1).name]);
        if(~isfield(firstChannelDir, 'label') && ~contains(fieldnames(firstChannelDir, '-full'),'label'))
            labelSingleUnits(char(currSession.Monkey),char(currSession.Date));
        end
    end
    [spikes,times,weights,currTrials,sessionConds,channels,~,~,chMap] =...
        getSessionInfo2(physDir,singleOrAllUnits);
    if(~isempty(spikes))
        [currSeg,currUnitChs,currTrialPSTHS,currActive,trialHists,alignedSpikes] = deal(repmat({[]},...
            1,length(conditions)));
        numUnits = size(spikes,1);
        for c = 1:length(conditions)
            currCond = conditions{c};
            condParamInd = cellfun(@(a) contains(a,conditions{c}),sessionConds);
            condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1))';
            condWeights = weights(:, condParamInd);
            condEvents = params.condSegMap(currCond);
            condAlign = cellfun(@(a) find(strcmp(condEvents,a)),...
                params.PSTHAlignments(currCond),'UniformOutput', false);
            monkeyImFile = strcat(currSession.Monkey,...
                    string(values(params.condAbbrev,{currCond})));
            if(allActivityMaps.isKey(monkeyImFile))
                condImgMap = values(allActivityMaps,{monkeyImFile});
                currActive{c} = condImgMap{1}(currSession.y, currSession.x);
            else
                currActive{c} = [NaN, NaN];
            end
            % generate and smooth PSTHS  from current session
            if(any(condInds) && ~isempty(spikes))
                % {alignedPSTHS}{units,trials}
                alignedSpikes{c} = cellfun(@(ap) cellfun(@(s,t) s-t(ap), ...
                    spikes(:,condInds),repmat(times(condInds),size(spikes,1),1),...
                    'UniformOutput',false),condAlign,'UniformOutput',false);
                alignedTimes = cellfun(@(ap) cell2mat(cellfun(@(t) t-t(ap),...
                    times(condInds),'UniformOutput', false)'), condAlign, 'UniformOutput', false);
                trialHists{c} = cellfun(@(ac) cellfun(@(a) histcounts(a,...
                    [bins,bins(end)+params.binSize])./(params.binSize),...
                    ac,'UniformOutput', false),alignedSpikes{c},'UniformOutput', false);
                %{alignedPSTHS}{units,trials}{bins}
                smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,...
                    gausswin(params.sigma)/sum(gausswin(params.sigma)),'same'),...
                    ac,'UniformOutput', false),trialHists{c},'UniformOutput',false);
                %cellfun(@(ac)  cellfun(@(a) conv(a,poisspdf(1:params.sigma,params.sigma)/...
                % sum(poisspdf(1:params.sigma,params.sigma)),'same'),ac,'UniformOutput', false),...
                % trialHists{c},'UniformOutput',false)
                %{alignedPSTHS}{units,bins,trials}
                unitTrialPSTH = cellfun(@(s) ...%condWeights.*...
                    reshape(cell2mat(s),[],length(bins),sum(condInds)),...
                    smoothHists,'UniformOutput', false);
                currUnitChs{c}= channels;
                % aligned trial times for the session for each (1) PSTH
                currSeg{c} = alignedTimes;
                % PSTH(units X bins X trials) for each alignment
                currTrialPSTHS{c} = unitTrialPSTH;
            else
                % pad stored info with empty arrays and NaN pad indicies for missing conditions
                currUnitChs{c} = NaN(1,numUnits);
                currSeg{c} = repmat({NaN(size(params.condSegMap(currCond)))},...
                    1,length(condAlign));
                currTrialPSTHS{c} = repmat({NaN(numUnits,...
                    length(bins),1)},1, ...
                    length(condAlign));
                alignedSpikes{c} = repmat({NaN(numUnits,1)},1,length(condAlign));
            end
           
        end
        % get current session joint label
        siteRep{i} = currSession.SiteRep{:};
        siteThresh{i} = currSession.Thresh{:};
        siteLocation{i} = [currSession.x, currSession.y];
        siteSegs{i} = currSeg;
        siteChannels{i} = currUnitChs;
        siteTrialPSTHS{i} = currTrialPSTHS;
        siteActiveInd{i} = currActive;
        rawSpikes{i} = alignedSpikes;
        channelMap{i} = chMap;
    end
     hbar.iterate(1);
end
hbar.close;hbar.delete;
%% remove sessions that had no trial information
emptyInds = cellfun(@isempty, siteLocation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludeRepInd = find(cellfun(@(f) length(f)==1 & sum(strcmp(f,...
    excludeRep))==length(f),siteRep(~emptyInds)));
emptyInds(arrayfun(@(f) find(cumsum(~emptyInds)==f,1),excludeRepInd)) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siteDateMap = siteDateMap(~emptyInds,:);
channelMap = channelMap(~emptyInds);
siteLocation = vertcat(siteLocation(~emptyInds));
siteRep = vertcat(siteRep(~emptyInds));
siteThresh = vertcat(siteThresh(~emptyInds));
siteSegs = num2cell(vertcat(siteSegs{~emptyInds}),1);
rawSpikes = num2cell(vertcat(rawSpikes{~emptyInds}),1);
siteChannels = num2cell(vertcat(siteChannels{~emptyInds}),1);
siteTrialPSTHS = num2cell(vertcat(siteTrialPSTHS{~emptyInds}),1);
siteActiveInd = num2cell(vertcat(siteActiveInd{~emptyInds}),1);
%%%%%%%%%% Replacing face sites with next best representation %%%%%%%%%%%%%
nextBest = cellfun(@(r)  ~strcmp(r,excludeRep), siteRep,'UniformOutput',false);
simpRep = string(cellfun(@(r,t,f) r{find((f.*t)==(min((f./f).*t)),1)},...
    siteRep,siteThresh, nextBest,'UniformOutput', false));
emptyInds = cellfun(@isempty, simpRep);
siteDateMap = siteDateMap(~emptyInds,:);
channelMap = channelMap(~emptyInds);
siteLocation = vertcat(siteLocation(~emptyInds))';
siteSegs = cellfun(@(ss) vertcat(ss(~emptyInds')),siteSegs,'UniformOutput',false);
siteChannels = cellfun(@(sc) sc(~emptyInds'),siteChannels,'UniformOutput',false);
siteTrialPSTHS = cellfun(@(stp) stp(~emptyInds'),siteTrialPSTHS,'UniformOutput',false);
siteActiveInd = cellfun(@(sa) sa(~emptyInds'),siteActiveInd,'UniformOutput',false);
rawSpikes = cellfun(@(rs) rs(~emptyInds'), rawSpikes, 'UniformOutput',false);
simpRep = simpRep(~emptyInds);
channelMap = channelMap(~emptyInds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% voronoi tiles for each monkey
siteMasks = repmat({},1,height(siteDateMap));
%%
for m = 1:length(monkeys)
    if(strcmp(monkeys(m),"Skipper"))
        mm = MotorMapping(50);
    else
        mm = MotorMapping(35);
    end
    mRefMask = vMask(monkeys(m));
    mRefMask = mRefMask{1};
    mMap = find(strcmp(siteDateMap.Monkey,monkeys(m)));
    if(~isempty(mMap))
        [verticies, vCells] = voronoin(fliplr([cell2mat(siteLocation(mMap)); ...
            [0 size(mRefMask,2); size(mRefMask,1) 0; 0 0;size(mRefMask,1) size(mRefMask,2)]]));
        for i = 1:length(mMap)
            currSite = siteLocation{mMap(i)};
            tempCircle = zeros(size(mRefMask)+2*mm.tileBuffer);
            tempCircle((currSite(2)-mm.siteRadius+mm.tileBuffer):...
                (currSite(2)+mm.siteRadius+mm.tileBuffer),...
                (currSite(1)-mm.siteRadius+mm.tileBuffer):...
                (currSite(1)+mm.siteRadius+mm.tileBuffer)) = ...
                mm.poolCircle;
            tempCircle = tempCircle(mm.tileBuffer:end-(...
                mm.tileBuffer+1),mm.tileBuffer:end-(...
                mm.tileBuffer+1));
            siteMasks{mMap(i)} = tempCircle & poly2mask(verticies(vCells{i},2),...
                verticies(vCells{i},1),size(tempCircle,1),size(tempCircle,2));
        end
    end
end
clear tempCircle verticies vCells poolCircle
% group multiple representation sites/units accordingly
remappedReps = cell(size(siteRep));
for sr = 1:length(siteRep)
    if(length(siteRep{sr})==1)
        remappedReps{sr} = siteRep{sr};
    else
        allSiteReps = unique(siteRep{sr}(siteThresh{sr} <= min(siteThresh{sr})*1.1));
        if(all(contains(allSiteReps, MotorMapping.forelimbRep)))
            remappedReps{sr} = "Forelimb";
        elseif(any(contains(allSiteReps, MotorMapping.forelimbRep)))
            remappedReps{sr} = "Mixed";
        else
            remappedReps{sr} = "Axial";
        end
    end
end
