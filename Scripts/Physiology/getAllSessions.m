function [siteDateMap,siteSegs,siteTrialPSTHS,rawSpikes,siteChannels,siteActiveInd, simpRep,...
    siteLocation, siteMasks, monkeys,vMask,conditions,channelMap,trialInfo] = ...
    getAllSessions(params,singleOrAllUnits,domain)
%  assign parameters
rawSpikes = [];
drivePath = "S:\Lab\";
monkeys = ["Gilligan", "Skipper"];
excludeRep = "Face";
condWindows = [0 0 0 .5];
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
if(strcmp(domain,"PMd"))
    siteDateMap = siteDateMap([2,4:17,19,20,21,23,24,25,27,28,31,32,38,43,46,47,48,49,51,53,56],:);
else
    %siteDateMap = siteDateMap([51,82,5,10,11,32,43,46,62],:);%,,33,37,39,45,49,53,60,61,67,69,71,72,73,77,84,85,87,93,97,19,47,59],:)
end
% load info from all sites
numSites = height(siteDateMap);
[siteLocation, siteRep, siteThresh,siteSegs,siteChannels,...
    siteTrialPSTHS,siteActiveInd,rawSpikes,channelMap] = deal(cell(1,numSites));
delete(gcp('nocreate'));
parpool('local'); 
parRun = size(gcp('nocreate'),1);
if(~parRun)
    hbar = waitbar(0, 'Processing...', 'Name',['Iterating ',num2str(numSites),' instances....']);
else
    hbar = parforProgress(numSites);
end
parfor  i = 1:numSites
    currSession = siteDateMap(i,:);
    if(strcmpi(currSession.Monkey,"Gilligan"))
        dateFormat = 'MM_dd_uuuu';
    else
        dateFormat = 'uuuu_MM_dd';
    end
    du = [currSession.Date{:}];
    du.Format = dateFormat;
    currSession.Date = char(du);
    physDir = strcat(drivePath,currSession.Monkey,"\All Data\", currSession.Monkey,...
        "_",string(currSession.Date),"\Physiology\");
    %delete(fullfile(fullfile(physDir,'*.cache')));
    physDir = strcat(physDir,"Results");
    if(isempty(dir(physDir+"\*.mat")))
        if(~ismember(currSession.Date,{'05_02_2019','11_11_2019','2021_09_20','2022_06_22','2022_06_28','2022_07_11','2022_07_12'}))
            disp(['Sorting and labeling session...']);
            Spike_SortRawData(currSession.Date,char(currSession.Monkey));
            labelSingleUnits(currSession.Date,char(currSession.Monkey));
        else
            disp(['Bad session: ', currSession.Date]);
        end
    end
    if(~isempty(dir(physDir+"\*.mat")))
        dirChannels = dir(physDir+"\*.mat");
        firstChannel = load([strcat(physDir,'\',dirChannels(1).name)]);
        if(~isfield(firstChannel,'label') && ~contains(fieldnames(firstChannel,'-full'),'label'))
            disp(['Labeling session...']);
            labelSingleUnits(currSession.Date,char(currSession.Monkey));
        end
    end
    [spikes,times,weights,currTrials,sessionConds,channels,~,~,chMap] =...
        getAllTrials(physDir,singleOrAllUnits,true);
    if(~isempty(spikes))
        [currSeg,currTrialPSTHS,currActive,trialHists,alignedSpikes] = deal(repmat({[]},1,length(conditions)));
        numUnits = size(spikes,1);
        NaNGraspInd = repmat({false},size(currTrials,1),1);
        NaNGraspInd(ismember(currTrials(:,1),conditions)) = cellfun(@(s,t) sum(isnan(t))==1 & isnan(t(strcmp(params.condSegMap(s),"StartGrasp"))), ...
            currTrials(ismember(currTrials(:,1),conditions),1),times(ismember(currTrials(:,1),conditions))','UniformOutput',false);
        successfulInds = cellfun(@(b,t) ~isnan(str2double(b)) | isempty(b),currTrials(:,end-1));
        successfulTrials = successfulInds | cell2mat(NaNGraspInd);
        currTrials = currTrials(successfulTrials,:);
        times = times(successfulTrials);
        for c = 1:length(conditions)
            currCond = conditions{c};
            condParamInd = cellfun(@(a) contains(a,conditions{c}),sessionConds);
            condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1));
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
                alignedSpikes(c) = cellfun(@(ap) cellfun(@(st) cellfun(@(s,t) ...
                    s-(t(ap)+condWindows(c)),st(condInds),times(condInds),'UniformOutput',false),...
                    spikes,'UniformOutput',false),condAlign,'UniformOutput',false);
                alignedTimes = cellfun(@(ap) cell2mat(cellfun(@(t)[t(1:end-1)-(t(ap)+condWindows(c)), ...
                    NaN(1,length(condEvents)-length(t)),t(end)-(t(ap)+condWindows(c))], times(condInds),...
                    'UniformOutput', false)'), condAlign, 'UniformOutput', false);
                trialHists{c} = cellfun(@(ac) cellfun(@(a) histcounts(vertcat(a(:)),...
                    [bins(1)-(params.sigmaSize/2):params.binSize:bins(end)+...
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
        end
        % get current session joint label
        siteRep{i} = currSession.SiteRep{:};
        siteThresh{i} = currSession.Thresh{:};
        siteLocation{i} = [currSession.x, currSession.y];
        siteChannels{i} = channels;
        siteActiveInd{i} = currActive;
        channelMap{i} = chMap;
        trialInfo{i} = currTrials(cellfun(@(a) ismember(a,conditions),currTrials(:,1)),:);
        siteSegs{i} = currSeg;
        siteTrialPSTHS{i} = currTrialPSTHS;
        rawSpikes{i} = alignedSpikes;
    end
    if(parRun)
        send(hbar, i);
    else
        waitbar(i/numSites,hbar,['Processed ' num2str(i),' of ', num2str(numSites), ' instances.']);
    end
end
if(~parRun)
    close(hbar);
end
%% remove sessions that had no trial information
emptyInds = cellfun(@isempty, siteLocation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludeRepInd = find(cellfun(@(f) length(f)==1 & sum(strcmp(f, excludeRep))==...
    length(f),siteRep(~emptyInds)));
emptyInds(arrayfun(@(f) find(cumsum(~emptyInds)==f,1),excludeRepInd)) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siteDateMap = siteDateMap(~emptyInds,:);
channelMap = channelMap(~emptyInds);
siteLocation = vertcat(siteLocation(~emptyInds));
siteRep = vertcat(siteRep(~emptyInds));
siteThresh = vertcat(siteThresh(~emptyInds));
siteSegs = num2cell(vertcat(siteSegs{~emptyInds}),1);
rawSpikes = num2cell(vertcat(rawSpikes{~emptyInds}),1);
trialInfo = vertcat(trialInfo(~emptyInds));
siteChannels = siteChannels(~emptyInds);
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
siteTrialPSTHS = cellfun(@(stp) stp(~emptyInds'),siteTrialPSTHS,'UniformOutput',false);
siteActiveInd = cellfun(@(sa) sa(~emptyInds'),siteActiveInd,'UniformOutput',false);
rawSpikes = cellfun(@(rs) rs(~emptyInds'), rawSpikes, 'UniformOutput',false);
simpRep = simpRep(~emptyInds);
channelMap = channelMap(~emptyInds);
siteChannels = siteChannels(~emptyInds');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% voronoi tiles for each monkey
siteMasks = repmat({},1,height(siteDateMap));
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
            tempCircle((currSite(2)-mm.siteRadius+mm.tileBuffer):(currSite(2)+...
                mm.siteRadius+mm.tileBuffer),(currSite(1)-mm.siteRadius+mm.tileBuffer):...
                (currSite(1)+mm.siteRadius+mm.tileBuffer))= mm.poolCircle;
            tempCircle = tempCircle(mm.tileBuffer:end-(mm.tileBuffer+1),...
                mm.tileBuffer:end-(mm.tileBuffer+1));
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
end