function getAllSessionsNew(monkeyIn, singleOrAllIn, sessionOrUnitIn,alignmentIn)
close all;
clear all;
clc;
%%
allEvents = containers.Map(["Extra Small Sphere", "Large Sphere", "Photocell", "Rest"], ...
    {["GoSignal","StartReach","StartGrasp","StartLift","StartHold",...
    "StartWithdraw","StartReplaceHold","StartReplaceSuccess","StartReward"],...
    ["GoSignal","StartReach","StartGrasp","StartLift","StartHold",...
    "StartWithdraw","StartReplaceHold","StartReplaceSuccess","StartReward"],...
    ["GoSignal","StartReach","StartGrasp","StartHold",...
    "StartWithdraw","StartReplaceHold","StartReplaceSuccess","StartReward"],...
    ["GoSignal","StartReplaceHold", "StartReplaceSuccess","StartReward"]});
%%%
if(~exist('monkeyIn', 'var'))
    monkey = ["Gilligan", "Skipper"];
else
    monkey = string(monkeyIn);
end
if(~exist('singleOrAllIn', 'var'))
    singleOrAll = "Single";
else
    singleOrAll = singleOrAllIn;
end
if(~exist('sessionOrUnitIn', 'var'))
    sessionOrUnitPSTHS = "Unit";
else
    sessionOrUnitPSTHS = sessionOrUnitIn;
end
if(~exist('alignmentIn', 'var'))
    PSTHAlignments = containers.Map(["Extra Small Sphere", "Large Sphere", "Photocell", "Rest"],...
        {["GoSignal","StartReach", "StartLift", "StartWithdraw"], ...
        ["GoSignal","StartReach", "StartLift", "StartWithdraw"],...
        ["GoSignal","StartReach", "StartHold", "StartWithdraw"],...
        ["GoSignal", "StartReplaceHold", "StartReplaceSuccess", "StartReward"]});
else
    PSTHAlignments = alignmentIn;
end

%  assign parameters
pVal = 0.05;
conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
sessionDirPrefix = "S:\Lab\";

% PSTH parameters and alignment
binSize=.01; % bin size in seconds
sigmaSize = .15; % smoothing window in seconds
sigma = round(sigmaSize/binSize);
secondsBeforePSTHAlignmentPoint = -6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 5; % time (in seconds) after alignement point to create PSTH
inclusiveBins = secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values
bins = inclusiveBins(1:end-1);

phaseWindowSize = 0.2;
baselineAlignments = ["GoSignal", "GoSignal"];
baselineWindow = [-4,-1];
phaseAlignments = {["GoSignal", "StartHold"],"GoSignal","StartReach","StartLift"};
phaseWindows = {[0 0],[0, phaseWindowSize],[-phaseWindowSize/2, phaseWindowSize/2], [-phaseWindowSize, 0]};

% get activity map for current condition.
spaceInds = cellfun(@(c) [0, regexp(c, '\s')],conditions, 'UniformOutput', false);
condAbbrev = cellfun(@(c,i) c(i+1),conditions, spaceInds, 'UniformOutput', false);
condAbbrev{3} = "PC";

% representation information from voronoi
squarePoolSize = 56; %28 = 0.5mm, 1.0mm^2; 42 = 0.75 mm, 2.2mm^2; 56 = 1mm, 3.4mm^2;
poolCircle = fspecial('disk',squarePoolSize);
erToName = containers.Map(["Digit", "Wrist", "Triceps", "Biceps", "Elbow", "Forearm", "Shoulder",...
    "Face","Ear", "Trunk", "Neck", "No response"],["Hand", "Hand", "Arm", "Arm","Arm", "Arm", "Arm", ...
    "Face", "Face", "Trunk", "Trunk", ""]);
erAbbrevs = erToName.keys;
repColors = struct("Arm", [.75 .75 .75], "Hand", [.35 .35 .35], "Forelimb",...
    [.5 .5 .5], "Mixed", [121 76 92]./255,"Trunk", [.5 .25 .1],...
    "Face",[155 9 144]./255,"Axial", [102 0 36]./255);
forelimbNames = ["Arm","Hand"];
repNames = fieldnames(repColors);
%%
[vMask, allActivityMaps] = deal(containers.Map());
siteDateMap = table();
for m = 1:length(monkey)
    %%% all mapped sites and x,y locations
    [~,~,rawRef]=xlsread(strcat(sessionDirPrefix,monkey{m},"\Mapping\",monkey{m},"_MM_Sites.xlsx"));
    heading = rawRef(1,:);
    rawRef(1,:) = [];
    [~,recordingInd] = find(strcmp(heading, "Recording"));
    [~,siteInd] = find(strcmp(heading, "Site #"));
    [~,threshInd] = find(strcmp(heading, "Threshold(s)"));
    [~,erInd] = find(strncmp(heading, "Evoked",length("Evoked")));
    [~,xInd] = find(strcmp(heading, "x"));
    [~,yInd] = find(strcmp(heading, "y"));
    recordingSiteInds = find(cellfun(@(a,b) strcmp(a,"Yes") & isnumeric(b),...
        rawRef(:,recordingInd),rawRef(:,siteInd)));
    % store information for each site
    currDateMap = array2table(cellfun(@str2double, string(rawRef(recordingSiteInds,...
        [siteInd,xInd,yInd]))),'VariableNames', ["Site", "x", "y"]);
    % mapped sites to determine identity of recorded sites
    mappedSiteInds = find(~cellfun(@(e) (length(e)==1 && isnan(e)) || ...
        strcmp(e, "No response"), rawRef(:,erInd)));
    allLocs = [cellfun(@double, rawRef(mappedSiteInds,xInd), 'UniformOutput', true),...
        cellfun(@double, rawRef(mappedSiteInds,yInd), 'UniformOutput', true)];
    [allReps, allThreshs] = deal(cell(length(recordingSiteInds),1));
    for r = 1:length(recordingSiteInds)
        currLoc = [double(rawRef{recordingSiteInds(r),xInd}), double(rawRef{recordingSiteInds(r),yInd})];
        % find mapped site closest to recorded site
        [~,minMappedInd] = min(sqrt(sum((allLocs-currLoc).^2,2)));
        % map evoked response names and current thresholds of closest mapped
        % site to current recorded site
        currERs = regexp(rawRef{mappedSiteInds(minMappedInd),erInd},',', 'split');
        thresh = rawRef{mappedSiteInds(minMappedInd),threshInd};
        mappedERs = cellfun(@(e) find(cellfun(@(ab) contains(strip(e),ab), ...
            erAbbrevs),1),currERs);
        allReps{r} = values(erToName,erAbbrevs(mappedERs));
        if(length(thresh)<length(mappedERs))
            disp(rawRef{mappedSiteInds(minMappedInd),1});
        end
        if(isnumeric(thresh))
            allThreshs{r} = thresh;
        else
            allThreshs{r} = cellfun(@(t) str2double(strip(t)),...
                regexp(thresh, ',', 'split'));
        end
    end
    currDateMap = addvars(currDateMap,allReps, allThreshs,repmat(...
        string(monkey{m}),length(allThreshs),1),'NewVariableNames',...
        ["SiteRep", "Thresh","Monkey"]);
    clear rawNum rawRef
    
    %%% confirm M1 sites based on assigned border
    MM = double(imread(strcat(sessionDirPrefix,monkey{m},"\Mapping\Motor Maps V3\MM_Simp_RGB-01.png")));
    M1border = bwareafilt(MM(:,:,1)>125 & MM(:,:,2)>125 & MM(:,:,3)<125, 1);
    % find border
    [row,col] = find(M1border);
    pointsX = unique(row);
    pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
    pointsX(end+1:end+2) = [size(MM,2), 1];
    pointsY(end+1:end+2) = [1 1];
    % create M1/PM mask
    M1mask = double(poly2mask(pointsY,pointsX, size(MM,1),size(MM,2)));
    % validate M1 sites based on location
    M1sites = logical(arrayfun(@(a,b) M1mask(b,a), currDateMap.x,currDateMap.y));
    currDateMap(~M1sites,:) = [];
    clear row col pointsX pointsY M1border M1mask
    %%% date and session information (double check M1 mask)
    [~,~,dateToPNValue] = xlsread(strcat(sessionDirPrefix, monkey{m},"\",monkey{m},"Session.xlsx"));
    [~, dateInd] = find(strcmp(dateToPNValue, "Date"));
    [~,domainInd] = find(strcmp(dateToPNValue, "Domain"));
    [~,PNInd] = find(strcmp(dateToPNValue, "Site #"));
    % clean up, remove: duplicate test sites (close/far), non-motor recordings,
    % improper logging (no site #)
    validInds = find(cellfun(@(a,b,c) ~any(isnan(a)) & any(strcmp(b, ...
        ["M1", "PMd", "PMv"])) & ~ismissing(string(c)),dateToPNValue(:,dateInd),...
        dateToPNValue(:,domainInd), dateToPNValue(:,PNInd), 'UniformOutput', true));
    % match session numbers from _MM_Sites.xlsx sheet to session.xlsx sheet
    [~,viPN,viSD] = intersect(cellfun(@str2double, string(dateToPNValue(...
        validInds,PNInd))),currDateMap.Site);
    validInds = validInds(viPN);
    % assign session number's date, add domain, duplicate session number from
    % session.xlsx sheet
    currDateMap = addvars(currDateMap(viSD,:), dateToPNValue(validInds,dateInd), ...
        dateToPNValue(validInds,domainInd), dateToPNValue(validInds,PNInd),...
        'NewVariableNames', ["Date", "Domain", "Site2"]);
    clear dateToPNValue
    %%% activity maps from imaging
    allActivityFiles = dir(strcat(sessionDirPrefix, monkey{m}, "\Mapping\tTests\*HSV.png"));
    allMaps = arrayfun(@(a) imresize(imread(strcat(a.folder, "\", a.name)),...
        [size(MM,1),size(MM,2)]), allActivityFiles, 'UniformOutput', false);
    % activity is red
    for ma = 1:length(allMaps)
        allActivityMaps(strcat(monkey{m},"_",allActivityFiles(ma).name(1:end-8)))=...
            allMaps{ma}(:,:,1)>100 & allMaps{ma}(:,:,2)<100 & allMaps{ma}(:,:,3)<100;
    end
    
    %%% vessel mask
    tempMask = imread(strcat(sessionDirPrefix,monkey{m}, "\Mapping\clean_mask_filled.bmp"));
    vMask(monkey{m}) = logical(tempMask(:,:,1)>200);
    
    siteDateMap = [siteDateMap;addvars(currDateMap)];
end
clear MM allMaps tempMask
%% main script: load info from all sites
numSites = height(siteDateMap);
[siteLocation, siteRep, siteThresh] = deal(cell(1,numSites));
[sessionSegs,sessionUnitChs,sessionTrialPSTHS,sessionActive,sessionPhasesAUC,baselineSampAUC] = deal({});
for  i = 1:numSites
    currSession = siteDateMap(i,:);
    if(strcmp(currSession.Monkey,"Gilligan"))
        dateFormat = 'mm_dd_yyyy';
    else
        dateFormat = 'yyyy_mm_dd';
    end
    if(~isempty(dir(strcat(sessionDirPrefix,currSession.Monkey,"\All Data\",...
            currSession.Monkey,"_",datestr(currSession.Date, dateFormat),"\Physiology\Results"))))
        [spikes,times,weights,currTrials,sessionConds,channels,events,~] =...
            getSessionInfo2(strcat(sessionDirPrefix,currSession.Monkey,...
            "\All Data\", currSession.Monkey,"_",datestr(currSession.Date,...
            dateFormat)),singleOrAll);
    else
        spikes = [];
    end
    if(~isempty(spikes))
        if(strcmp(sessionOrUnitPSTHS, "Unit"))
            numUnits = size(spikes,1);
        else
            numUnits = 1;
        end
        if(strcmp(conditions, "All"))
            conditions = sessionConds;
            %conditions(contains(conditions,"Rest")) = [];
        end
        if(length(events.keys)==length(conditions))
            eventsAll=events;
        end
        if(length(sessionSegs)<length(conditions))
            missingInitalization = length(conditions)-length(sessionSegs);
            sessionSegs(end+1:end+missingInitalization) = cell(1,missingInitalization);
            sessionTrialPSTHS(end+1:end+missingInitalization) = cell(1,missingInitalization);
            sessionUnitChs(end+1:end+missingInitalization) = cell(1,missingInitalization);
            sessionActive(end+1:end+missingInitalization) = cell(1,missingInitalization);
            sessionPhasesAUC(end+1:end+missingInitalization) = cell(1,missingInitalization);
            baselineSampAUC(end+1:end+missingInitalization) = cell(1,missingInitalization);
        end
        for c = 1:length(conditions)
            condParamInd = cellfun(@(a) contains(a,conditions{c}),sessionConds);
            condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1))';
            condWeights = weights(:, condParamInd);
            if(strcmp(conditions{c}, "Rest"))
                condPhaseAlignments = {{"GoSignal", "StartReplaceHold"},{"GoSignal"}};
                condPhaseWindows = {[0 0], [-1 1]};
            elseif(strcmp(conditions{c}, "Photocell"))
                condPhaseAlignments = cellfun(@(s) strrep(s,"StartLift",...
                    "StartHold"), phaseAlignments,'UniformOutput',false);
                condPhaseWindows = phaseWindows;
            else
                condPhaseAlignments = phaseAlignments;
                condPhaseWindows = phaseWindows;
            end
            % generate and smooth PSTHS  from current session
            if(any(condInds) & ~isempty(spikes))
                condEvents = events(conditions{c});
                condAlign = cellfun(@(a) find(strcmp(condEvents,a)),...
                    PSTHAlignments(conditions{c}),'UniformOutput', true);
                
                % {alignedPSTHS}{units,trials}
                alignedSpikes = arrayfun(@(ap) cellfun(@(s,t) s-t(ap), ...
                    spikes(:,condInds),repmat(times(condInds),size(spikes,1),1),...
                    'UniformOutput',false),condAlign,'UniformOutput',false);
                alignedTimes = arrayfun(@(ap) cell2mat(cellfun(@(t) t-t(ap), times(condInds), ...
                    'UniformOutput', false)'), condAlign, 'UniformOutput', false);
                
                trialHists = cellfun(@(ac) cellfun(@(a) ...
                    histcounts(a,inclusiveBins)./binSize, ac,'UniformOutput', false),...
                    alignedSpikes,'UniformOutput', false);
                %{alignedPSTHS}{units,trials}{bins}
                smoothHists = cellfun(@(ac)  cellfun(@(a) conv(a,gausswin(sigma)/...
                    sum(gausswin(sigma)),'same'),ac,'UniformOutput', false),...
                    trialHists, 'UniformOutput', false);
                %{alignedPSTHS}{units,bins,trials}
                unitPSTHS = cellfun(@(s) condWeights.*reshape(cell2mat(s),[size(s,1),...
                    length(bins),sum(condInds)]),smoothHists,'UniformOutput', false);
                
                baselineAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
                    baselineAlignments,'UniformOutput', true);
                psthPhaseAlignment = cellfun(@(pa) cellfun(@(a) find(strcmp(condEvents,a)),...
                    pa,'UniformOutput',true),condPhaseAlignments, 'UniformOutput', false);
                % find start and end index in (1) each PSTH for (2) each phase
                psthPhaseEnds = cellfun(@(at) cellfun(@(pa,pw) num2cell(...
                    findBins(at(:,pa)+pw,bins),2),psthPhaseAlignment,condPhaseWindows,'UniformOutput',false),...
                    alignedTimes,'UniformOutput', false);
                % create mask that NaNs out all values outside of interrogated phase
                % for each (1) PSTH for each (2) phase
                psthPhaseMask = cellfun(@(pr) cellfun(@(ps) NaNMaskBins(ps,length(bins)),...
                    pr,'UniformOutput', false), psthPhaseEnds, 'UniformOutput', false);
                % divide AUC by time if time is not fixed for the current
                % phase
                psthPhaseLength = cellfun(@(pp) cellfun(@(p) repmat(numel(unique(...
                    squeeze(sum(~isnan(p),2))))~=1,size(p,1),1).*...
                    squeeze(sum(~isnan(p),2))+1, pp,'UniformOutput',false),...
                    psthPhaseMask,'UniformOutput', false);
                % create equivalent sized start and end points in the baseline window
                % for each (1) PSTH for each (2) phase
                psthBaselineMask = cellfun(@(a,pp) cellfun(@(p) max(ones(size(p,1),2), ...
                    findBins(a(:,baselineAlignment)+baselineWindow,bins)...
                    -[zeros(size(p,1),1),nansum(p,2)]),pp,'UniformOutput', false),...
                    alignedTimes,psthPhaseMask, 'UniformOutput', false);
                % bootstrapped sample of AUC in baseline window for each
                % (1) PSTH for each (2) tested phase
                baselineSampAUC{c}(end+1,:) = cellfun(@(pb,hh,pp,pt) cellfun(@(p,w,t)...
                    AUCBaselineBootstrap(num2cell(p,2),pt,w)./t',pb,hh,pp, 'UniformOutput', false),...
                    psthBaselineMask,psthPhaseMask,psthPhaseLength,unitPSTHS,'UniformOutput',false);
                % AUC for each (1) PSTH for each tested phase
                sessionPhasesAUC{c}(end+1,:) = cellfun(@(udh,pp,ud)...
                    cellfun(@(ps,tl) permute(trapz(ud.*~isnan(permute(...
                    repmat(ps,[1,1,size(ud,1)]),[3 2 1])),2),[1 3 2])./tl'...
                    ,udh,pp,'UniformOutput',false),psthPhaseMask,psthPhaseLength,...
                    unitPSTHS,'UniformOutput', false);
                
                if(strcmp(sessionOrUnitPSTHS, 'Unit'))
                    sessionUnitChs{c}(end+1) = {channels};
                else
                    sessionUnitChs{c}(end+1) = {[]};
                end
                % aligned trial times for the session for each (1) PSTH
                sessionSegs{c}(end+1)= {alignedTimes};
                % PSTH(units X bins X trials) for each alignment
                sessionTrialPSTHS{c}(end+1) = {unitPSTHS};
            else
                % pad stored info with empty arrays and NaN pad indicies for missing conditions
                sessionUnitChs{c}(end+1) = {NaN(1,numUnits)};
                sessionSegs{c}(end+1) = {repmat({NaN(size(allEvents(conditions{c})))}...
                    ,1,length(PSTHAlignments))};
                sessionTrialPSTHS{c}{end+1} = repmat({NaN(numUnits,length(bins),1)},...
                    1, length(PSTHAlignments));
                baselineSampAUC{c}(end+1,:) = repmat({repmat({NaN(numUnits,1)},...
                    1,length(condPhaseAlignments))},1,length(PSTHAlignments));
                sessionPhasesAUC{c}(end+1,:) =  repmat({repmat({NaN(numUnits,1)},...
                    1,length(condPhaseAlignments))},1,length(PSTHAlignments));
            end
            condImgMap = allActivityMaps(strcat(currSession.Monkey,"_",condAbbrev{c}));
            sessionActive{c}(end+1) = condImgMap(currSession.y, currSession.x);
        end
        currSite = [currSession.x, currSession.y];
        % get current session joint label
        siteRep{i} = cellfun(@cellstr, currSession.SiteRep,'UniformOutput', false);
        siteThresh{i}= currSession.Thresh;
        siteLocation{i}= num2cell(currSite,2);
        % voronoi tiling
        tempCircle = zeros(size(vMask(currSession.Monkey))+400);
        tempCircle((currSite(2)-squarePoolSize+200):(currSite(2)+squarePoolSize+200),...
            (currSite(1)-squarePoolSize+200):(currSite(1)+squarePoolSize+200))...
            = poolCircle;
        tempCircle = tempCircle(200:end-201,200:end-201);
        siteMasks{i} = tempCircle>0;
    end
    disp(strcat(num2str(currSession.Site),", ",datestr(currSession.Date,dateFormat),...
        ".  Session ", num2str(i) ," of ", num2str(numSites)));
end
clear trialHists smoothHists unitPSTHS alignedSegs alignedSpikes currTrials alignedTimes times tempCircle
%% remove sessions that had no trial information
emptyInds = cellfun(@isempty, siteLocation);
siteLocation = siteLocation(~emptyInds);
siteRep = siteRep(~emptyInds);
siteThresh = siteThresh(~emptyInds);
siteMasks = siteMasks(~emptyInds);
siteLocation = horzcat(siteLocation{:});
siteRep = horzcat(siteRep{:});
siteThresh = horzcat(siteThresh{:});
allUnits = cellfun(@(sc) cellfun(@(si) size(si,2),sc, 'UniformOutput', true), ...
    sessionUnitChs, 'UniformOutput', false);

siteDateMap = siteDateMap(~emptyInds,:);
vXY = containers.Map();
for m = 1:length(monkey)
    mtInds = strcmp(siteDateMap.Monkey,monkey{m});
    [verticies, vCells] = voronoin([cell2mat(siteLocation(mtInds)'); ...
        [0 size(vMask(monkey{m}),2); size(vMask(monkey{m}),1) 0; 0 0; ...
        size(vMask(monkey{m}),1) size(vMask(monkey{m}),2)]]);
    vXY(monkey{m}) = cellfun(@(v) [verticies(v,2),verticies(v,1)], vCells, 'UniformOutput', false);
end
%% group multiple representation sites/units accordingly
remappedReps = cell(size(siteRep));
for s = 1:length(siteRep)
    if(length(siteRep{s})==1)
        remappedReps(s) = siteRep{s};
    else
        allSiteReps = unique(siteRep{s}(siteThresh{s} <= min(siteThresh{s})*1.1));
        if(all(contains(allSiteReps, forelimbNames)))
            remappedReps{s} = "Forelimb";
        elseif(any(contains(allSiteReps, forelimbNames)))
            remappedReps{s} = "Mixed";
        else
            remappedReps{s} = "Axial";
        end
    end
end
simpRep = string(cellfun(@(r,t) r{t==min(t)}, siteRep,siteThresh, 'UniformOutput', false));
sessionToUnitInds = cellfun(@(ac) cell2mat(arrayfun(@(m,n) ones(1,m)*n,ac,...
    1:length(ac),'UniformOutput',false)), allUnits,'UniformOutput', false);
%%
alignLimits = {[-.6,.15],[-.6, .15], [-.15, .6],[-.15,.6]};
phaseNames = ["Task", "Go", "Reach", "Grasp"];
saveDirPath = "S:\Lab\ngc14\Figures\Physiology\Results\NEW\";
taskUnitBase = {{}};
taskUnitPhase = {{}};
AUCRange = [0 250];
countRange = [0 20];
FRLim = [0 1];
%%
close all
restTaskUnits = false(1,max(cellfun(@nansum,allUnits,'UniformOutput',true)));
maxUnitFR = cellfun(@(c) cell2mat(cellfun(@(p) cell2mat(cellfun(@(u) ...
    max(nanmax(u,[],3),[],2), p,'UniformOutput', false)), c,...
    'UniformOutput', false)'), sessionTrialPSTHS, 'UniformOutput',false);
maxUnitFR = nanmax(cat(3,maxUnitFR{:}),[],3);
unitRep = cellfun(@(a) mapSessionInds2Units(a,cellstr(simpRep)), allUnits, 'UniformOutput', false);

for c = 1:length(conditions)
    condEvents = allEvents(conditions{c});
    PSTHAlignInd = cellfun(@(a) find(strcmp(condEvents,a)),PSTHAlignments(conditions{c}),'UniformOutput', true);
    taskTrialPSTH = sessionTrialPSTHS{c};
    condBaseline = baselineSampAUC{c};
    condPhases = sessionPhasesAUC{c};
    condSegs = sessionSegs{c};
    activityInds = sessionActive{c};
    condRep = unitRep{c};
    %
    [~,allTaskUnits] = ttestTrials(condBaseline,condPhases,[],...
        find(strcmp(phaseNames,"Task")),true,pVal);
    condTaskSessions = mapUnitInds2Session(sessionToUnitInds{c},allTaskUnits);
    taskTrialPSTH = cellfun(@(ti,st) cellfun(@(s) (ti./ti)'.*s, st, 'UniformOutput',false),...
        condTaskSessions,taskTrialPSTH,'UniformOutput',false);
    
    RGunits = false(1,length(allTaskUnits));
    if(~strcmp(conditions{c}, "Rest"))
        restTaskUnits(find([condTaskSessions{:}])) = 1;
        condPhaseNames = phaseNames;
    else
        condPhaseNames = phaseNames(1);
    end
    for m=1:length(monkey)+1
        if(m==length(monkey)+1)
            mInds = true(1,height(siteDateMap));
            saveName = strcat("BOTH_",sessionOrUnitPSTHS,"_",singleOrAll);
        else
            mInds = strcmp(siteDateMap.Monkey,monkey{m})';
            saveName = strcat(monkey{m},"_",sessionOrUnitPSTHS,"_",singleOrAll);
        end
        mtInds = mapSessionInds2Units(allUnits{c},mInds) & allTaskUnits;
        unitPSTHS = cellfun(@(p) nanmean(p,3),vertcat(taskTrialPSTH{:}), 'UniformOutput',false);
        unitPSTHS = cellfun(@cell2mat,num2cell(unitPSTHS,1),'UniformOutput',false);
        unitNormPSTHS = cellfun(@(u,n) u./n, unitPSTHS,...
            num2cell(maxUnitFR,1),'UniformOutput', false);
        trialSegs =  num2cell(vertcat(condSegs{mInds}),1);
        
        for pw = 1:length(condPhaseNames)
            [phaseBaseSubAUC, phaseMod] = ttestTrials(condBaseline,...
                condPhases,[],pw,true,pVal);
            phaseMod = (mtInds & phaseMod);
            phaseBaseSubAUC = nanmean((phaseMod./phaseMod)'.*phaseBaseSubAUC,2)';
            
            % unit counts
            currRep=string(condRep);
            if (strcmp(phaseNames(pw), "Task"))
                singleInds = mapSessionInds2Units(allUnits{c},mInds);
                currRep(~singleInds)="";
            else
                currRep(~mtInds) = "";
            end
            if(any(contains(["Reach", "Grasp"], phaseNames{pw})))
                RGunits(phaseMod) = 1;
                baseSubTrialPhase = cellfun(@(r,b) r{pw}-b{pw},condPhases,condBaseline,...
                    'UniformOutput', false);
                if(strcmp(phaseNames(pw),"Reach"))
                    rAUC = num2cell(baseSubTrialPhase,1);
                else
                    gAUC = num2cell(baseSubTrialPhase,1);
                end
            end
            unitCountFig = modulatedUnitsPerRep(currRep,phaseMod,...
                strcat("Task Unit Count: ", phaseNames(pw), " phase"),...
                strcat(phaseNames(pw)," Units"), repColors);
            saveFigures(unitCountFig, strcat(saveDirPath,"\Count\", condAbbrev{c},"\"),...
                strcat(saveName,"_", phaseNames(pw)),[]);
            if(m<length(monkey)+1)
                sitePhaseMod = mapUnitInds2Session(sessionToUnitInds{c},...
                    phaseMod);
                countMapFig = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
                    siteMasks(mInds),cellfun(@nansum,sitePhaseMod),false,5,countRange);
                saveFigures(countMapFig, strcat(saveDirPath,"\Maps\Count\", condAbbrev{c},"\"),...
                    strcat(saveName,"_", phaseNames(pw)),[]);
            end
            % unit AUC
            FRFig = modulatedUnitsPerRep(currRep,phaseBaseSubAUC,...
                strcat("AUC: ", phaseNames(pw), " phase"), ...
                strcat(phaseNames(pw)," Units"), repColors);
            saveFigures(FRFig, strcat(saveDirPath,"\AUC\", condAbbrev{c},"\"),...
                strcat(saveName,"_", phaseNames(pw)),[]);
            if(m<length(monkey)+1)
                siteAUC = mapUnitInds2Session(sessionToUnitInds{c},...
                    phaseBaseSubAUC);
                if(pw==1)
                    FRMapFig = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
                        siteMasks(mInds),cellfun(@nanmean,siteAUC),false,5,[1 20]);
                else
                    FRMapFig = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
                        siteMasks(mInds),cellfun(@nanmean,siteAUC),false,5,AUCRange);
                end
                saveFigures(FRMapFig, strcat(saveDirPath,"\Maps\AUC\", condAbbrev{c},"\"),....
                    strcat(saveName,"_", phaseNames(pw)),[]);
            end
            %% unit PSTHS
            sessionPhases = mapUnitInds2Session(sessionToUnitInds{c},phaseMod);

            phasePSTHS = cellfun(@(v) (phaseMod./phaseMod)'.*vertcat(v{:}),...
                num2cell(unitPSTHS,1), 'UniformOutput',false);
            phaseJointFig = plotJointPSTHS(bins, phasePSTHS,trialSegs,...
                currRep,simpRep(mInds),activityInds,PSTHAlignInd,alignLimits);
            saveFigures(phaseJointFig, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\"), strcat(saveName,"_", phaseNames(pw)),[]);
            
            normPhasePSTH = cellfun(@(p) p(phaseMod',:),unitNormPSTHS,'UniformOutput',false);
            psthCat = horzcat(normPhasePSTH{2});
            [~,sortVals] = max(psthCat(~any(isnan(psthCat),2),:),[],2);
            
            phaseHMFig = unitJointPSTH(bins,normPhasePSTH,currRep(phaseMod),...
                bins(sortVals),trialSegs,sessionToUnitInds{c}(phaseMod),...
                activityInds,PSTHAlignInd,alignLimits,FRLim);
            labelUnitPSTHS(phaseHMFig);
            cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\Units\Peak\",phaseNames(pw),"\"), strcat(saveName,"_", rh,"_",...
                "_Heatmaps"),[]),phaseHMFig,[repNames(contains(string(repNames),unique(simpRep)))',"All"]);
            
            unitLocation = cell2mat(mapSessionInds2Units(allUnits{c}, siteLocation)');
            unitLocation =(unitLocation-min(unitLocation)).*(18/1000);
            mlSort = unitLocation(phaseMod,2);
            rcSort = unitLocation(phaseMod,1);
            
            mlHMFig = unitJointPSTH(bins,normPhasePSTH,currRep(phaseMod),...
                mlSort,trialSegs,sessionToUnitInds{c}(phaseMod),...
                activityInds,PSTHAlignInd,alignLimits,FRLim);
            labelUnitPSTHS(mlHMFig);

            cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\Units\MedLat\",phaseNames(pw),"\"), strcat(saveName,"_", rh,"_",...
                "_Heatmaps"),[]),mlHMFig,[repNames(contains(string(repNames),unique(simpRep)))',"All"]);
            %%
            rcHMFig = unitJointPSTH(bins,normPhasePSTH,currRep(phaseMod),...
                rcSort,trialSegs,sessionToUnitInds{c}(phaseMod),...
                activityInds,PSTHAlignInd,alignLimits,FRLim);
            labelUnitPSTHS(rcHMFig);
             cellfun(@(f) camroll(gca(f),90),rcHMFig,'UniformOutput',false);
             cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\Units\RostCaud\",phaseNames(pw),"\"), strcat(saveName,"_", rh,"_",...
                "_Heatmaps"),[]),rcHMFig,[repNames(contains(string(repNames),unique(simpRep)))',"All"]);
            %close all;
        end
        if(c<length(conditions))
            currRep = string(condRep);
            currRep(~RGunits) = "";
            rgUnitsVals = cellfun(@(rp,gp) cell2mat(cellfun(@(r,g) [nanmean(r,2),...
                nanmean(g,2),ttest(r',g','Alpha', pVal/2,'Dim',1)'],rp,gp,...
                'UniformOutput',false)),rAUC,gAUC, 'UniformOutput',false);
            rgUnitsVals = cat(3,rgUnitsVals{:});
            
            sigRGunits = nansum(rgUnitsVals(:,3,:),3)>size(rgUnitsVals,3)/2;
            rgVals = nanmean(rgUnitsVals(:,[1,2],:),3);
            [~,mRGInd] = max(abs(rgVals),[],2);
            RUCounts = modulatedUnitsPerRep(currRep,[sigRGunits & mRGInd==1]',...
                strcat("R>G Units"), "R>G", repColors);
            GUCounts = modulatedUnitsPerRep(currRep,[sigRGunits & mRGInd==2]',...
                strcat("R<G Units"), "R<G", repColors);
            saveFigures(RUCounts, strcat(saveDirPath,"\Count\", condAbbrev{c},"\"),...
                strcat(saveName,"_Rx" ),[]);
            saveFigures(GUCounts, strcat(saveDirPath,"\Count\", condAbbrev{c},"\"),...
                strcat(saveName,"_Gx" ),[]);
            
            rPSTHS = cellfun(@(v) v(mRGInd==1,:),unitPSTHS, 'UniformOutput',false);
            rJointFig = plotJointPSTHS(bins,rPSTHS,trialSegs,currRep(mRGInd==1),...
                simpRep(mInds),activityInds,PSTHAlignInd,alignLimits);
            saveFigures(rJointFig, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\"), strcat(saveName,"_Rx"),[]);
            gPSTHS = cellfun(@(v) v(mRGInd==2,:),unitPSTHS, 'UniformOutput',false);
            gJointFig = plotJointPSTHS(bins,gPSTHS,trialSegs,currRep(mRGInd==2),...
                simpRep(mInds),activityInds,PSTHAlignInd,alignLimits);
            saveFigures(gJointFig, strcat(saveDirPath, "\PSTH\",...
                condAbbrev{c},"\"), strcat(saveName,"_Gx"),[]);
            
            if(m<length(monkey)+1)
                siteRG = mapUnitInds2Session(sessionToUnitInds{c},mRGInd);
                countMapR = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
                    siteMasks(mInds),cellfun(@(r) nansum(r==1),siteRG(mInds)),...
                    false,5,countRange./2);
                saveFigures(countMapR, strcat(saveDirPath,"\Maps\Count\", condAbbrev{c},"\"),...
                    strcat(saveName,"_Rx"),[]);
                countMapG = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
                    siteMasks(mInds),cellfun(@(g) nansum(g==2),siteRG(mInds)),...
                    false,5,countRange./2);
                saveFigures(countMapG, strcat(saveDirPath,"\Maps\Count\", condAbbrev{c},"\"),...
                    strcat(saveName,"_Gx"),[]);
            end
        end
    end
end

%%units modulating differentially for reach-to-grasp conditions
ssInds = find(strcmp(conditions,"Extra Small Sphere"));
lsInds = find(strcmp(conditions,"Large Sphere"));
testPhase = find(strcmp(phaseNames, "Task"));
for m=1:length(monkey)+1
    if(m==length(monkey)+1)
        mInds = true(1,height(siteDateMap));
        saveName = strcat("BOTH_",sessionOrUnitPSTHS,"_",singleOrAll);
    else
        mInds = strcmp(siteDateMap.Monkey,monkey{m})';
        saveName = strcat(monkey{m},"_",sessionOrUnitPSTHS,"_",singleOrAll);
    end
    currRep = cellfun(@(a) mapSessionInds2Units(a,cellstr(simpRep)), max(allUnits{ssInds},...
        allUnits{lsInds},[],2), 'UniformOutput', false);
    
    ssVals = cellfun(@(s,b) s{testPhase}-b{testPhase}, taskUnitPhase{ssInds},...
        taskUnitPhase{lsInds},'UniformOutput', false);
    sphereUnitsVals = cellfun(@(rp,gp) cell2mat(cellfun(@(r,g) [nanmean(r,2),...
            nanmean(g,2),ttest(r',g','Alpha', pVal/2,'Dim',1)'],rp,gp,...
            'UniformOutput',false)),rAUC,gAUC, 'UniformOutput',false);
      
            taskUnitBase{c}{pw} = condBaseline;
            taskUnitPhase{c}{pw} = condPhases;
            
    currRep(~RGunits) = "";
    rgUnitsVals = cellfun(@(rp,gp) cell2mat(cellfun(@(r,g) [nanmean(r,2),...
        nanmean(g,2),ttest(r',g','Alpha', pVal/2,'Dim',1)'],rp,gp,...
        'UniformOutput',false)),rAUC,gAUC, 'UniformOutput',false);
    rgUnitsVals = cat(3,rgUnitsVals{:})
    [~,allTaskUnits] = ttestTrials(condBaseline,condPhases,[],...
        find(strcmp(phaseNames,"Task")),true,pVal);
        condTaskSessions = mapUnitInds2Session(sessionToUnitInds{c},allTaskUnits);
    taskTrialPSTH = cellfun(@(ti,st) cellfun(@(s) (ti./ti)'.*s, st, 'UniformOutput',false),...
        condTaskSessions,taskTrialPSTH,'UniformOutput',false);


    mtInds = mapSessionInds2Units(allUnits{ss},mInds) & mapSessionInds2Units(...
        allUnits{ls},mInds);
    unitNormPSTHS = cellfun(@(p) nanmean(p,3),...
        vertcat(taskTrialPSTH{:}), 'UniformOutput',false);
    unitNormPSTHS = cellfun(@(u,n) u./n, unitNormPSTHS,...
        num2cell(maxUnitFR,1),'UniformOutput', false);
    trialSegs =  num2cell(vertcat(condSegs{mInds}),1);
    
    bothSphereSessions = allUnits{ssInds}(mInds)>0 & allUnits{lsInds}(mInds)>0;
    bothSphereUnits = mapSessionInds2Units(allUnits{ssInds}(bothSphereSessions),bothSphereSessions);
    currSessionInds = mapUnitInds2Session(sessionToUnitInds{ssInds}(bothSphereUnits),bothSphereUnits)';
   [~, spInds] = ttestTrials(taskUnitAUCS{ssInds,1},taskUnitAUCS{lsInds,1},...
        [],find(strcmp(phaseNames,"Task")),false,pVal);
    unitSPBoth = modulatedUnitsPerRep(unitRep(bothSphereUnits),spInds',"Sphere modulated units","Spheres",repColors);
    saveFigures(unitSPBoth,strcat(saveDirPath,"Count\Cond\"),strcat(saveName,"_Sphere_modulated_units"),[]);
   
    if(m<length(monkey)+1)
        siteGAUC = mapUnitInds2Session(sessionToUnitInds{c}(mtInds),...
            rgVals(gUnits));
        FRMapGAUCFig = mapUnitVals(vXY(monkey{m}),vMask(monkey{m}),...
            siteMasks(mInds),cellfun(@nanmean,siteGAUC),false,5,AUCRange);
        saveFigures(FRMapGAUCFig, strcat(saveDirPath,"\Maps\AUC\", condAbbrev{c},"\"),....
            strcat(saveName,"_GUnits"),[]);
    end
end
end
%%

function labelUnitPSTHS (cellFigs)
ax1 = cellfun(@gca, cellFigs,'UniformOutput',false);
ax2 = cellfun(@(x1) axes('Position', x1.Position), ax1,'UniformOutput',false);
cellfun(@(a1,a2) linkprop([a1,a2], {'XLim','YLim','Position','view'}), ax1,ax2,'UniformOutput',false);
cellfun(@(y1,y2) ylim(y2,[y1.YLim]), ax1,ax2,'UniformOutput',false);
cellfun(@(y1,y2)  yticks(y2,y1.YTick), ax1,ax2,'UniformOutput',false);
cellfun(@(a2) set(a2,'Visible', 'off'), ax2,'UniformOutput', false);
cellfun(@(a1,a2) text(a1,repmat(max(xlim(a1))*1.01,size(a2.YTick)), ...
    a2.YTick,arrayfun(@(s) sprintf('%d ',fix(s)), a2.YTick,'UniformOutput',false),...
    'HorizontalAlignment','left','VerticalAlignment','middle'),...
    ax1,ax2,'UniformOutput', false);
end
function [phaseSubVals, sigs] = ttestTrials(dist1,dist2,uInds,taskPhase,paired,pVal)
if(isempty(uInds))
    uInds = cellfun(@(u) true(1,size(u{1},1)),dist1,'UniformOutput',false);
end
if(paired)
    sigs = nansum(cell2mat(cellfun(@(d1,d2,b) ttest(d1{taskPhase}(b,:)',...
        d2{taskPhase}(b,:)','Alpha', pVal,'Dim',1),dist1,dist2,uInds,...
        'UniformOutput', false)'),1)>size(dist1,2)/2;
else
    sigs = nansum(cell2mat(cellfun(@(d1,d2,b) ttest2(d1{taskPhase}(b,:)',...
        d2{taskPhase}(b,:)','Alpha', pVal,'Dim',1),dist1,dist2,uInds,...
        'UniformOutput', false)'),1)>size(dist1,2)/2;
end
phaseSubVals = cell2mat(cellfun(@(d1,d2,b) ...
    nanmean(d2{taskPhase}(b,:)-d1{taskPhase}(b,:),2),...
    dist1,dist2,uInds,'UniformOutput',false));
end
function sessionInds = mapUnitInds2Session(session2Units, currInds)
sessionInds = arrayfun(@(a) currInds(session2Units==a),1:max(session2Units),'UniformOutput', false);
end
function unitInds = mapSessionInds2Units(unit2Session,currInds)
unitInds = arrayfun(@(sr,si) repmat(sr,si,1),currInds, unit2Session,'UniformOutput', false);
unitInds = vertcat(unitInds{:})';
end
function bins = findBins(timesT,allBins)
if(~any(~isnan(timesT)))
    bins =NaN(size(timesT));
else
    bins = discretize(timesT,allBins);
end
end
function nanMask = NaNMaskBins(ap,binLength)
if(all(cellfun(@(p) ~any(~isnan(p)), ap)))
    nanMask = NaN(size(ap,1),binLength);
else
    nanMask = cell2mat(cellfun(@(p) [NaN(1,max(1,p(1))-1),...
        ones(1,max(0,1+min(binLength,p(end))-max(1,p(1)))),...
        NaN(1,binLength-min(binLength,p(end)))],ap,'UniformOutput',false));
end
end

function AUCPSTH = AUCBaselineBootstrap(ap,psths,winds)
NUMSAMPLES = 5;
if(all(cellfun(@(p) ~any(~isnan(p)), ap)))
    AUCPSTH = NaN;
else
    sampStart = cellfun(@(ss) randi([ss], 1, NUMSAMPLES),ap, 'UniformOutput', false);
    phaseL = nansum(winds,2);
    medAUC = [];
    for n = 1:NUMSAMPLES
        sinds = repmat(cellfun(@(s) s(n), sampStart, 'UniformOutput',true),1,2);
        sinds(:,2) = sinds(:,2) + phaseL;
        nanMask = cell2mat(cellfun(@(p) [NaN(1,max(1,p(1))-1),...
            ones(1,max(0,1+min(length(winds),p(end))-max(1,p(1)))),...
            NaN(1,length(winds)-min(length(winds),p(end)))],num2cell(sinds,2),'UniformOutput',false));
        medAUC(:,:,n) = permute(trapz(psths.*permute(...
            repmat(nanMask==1,1,1,size(psths,1)),[3 2 1]), 2),[1 3 2]);
    end
    AUCPSTH = nanmedian(medAUC,3);
end
end