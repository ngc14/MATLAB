function getAllSessions(monkeyIn, singleOrAllIn, sessionOrUnitIn, conditionIn,alignmentIn)
close all;
%%
allEvents = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'}, ...
    {{'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartLift','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReach','StartGrasp','StartHold',...
    'StartWithdraw','StartReplaceHold','StartReplaceSuccess','StartReward'},...
    {'GoSignal','StartReplaceHold', 'StartReplaceSuccess','StartReward'}});
%%%
if(~exist('monkeyIn', 'var'))
    monkey = 'Skipper';
else
    monkey = monkeyIn;
    if(strcmp(monkey, 'Both'))
       disp('st'); 
    end
end
if(~exist('singleOrAllIn', 'var'))
    singleOrAll = 'Single';
else
    singleOrAll = singleOrAllIn;
end
if(~exist('sessionOrUnitIn', 'var'))
    sessionOrUnitPSTHS = 'Unit';
else
    sessionOrUnitPSTHS = sessionOrUnitIn;
end
if(~exist('conditionIn', 'var'))
    conditions = {'All'};
else
    conditions = conditionIn;
end
if(~exist('alignmentIn', 'var'))
    PSTHAlignments = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'},...
        {{'GoSignal','StartReach', 'StartLift', 'StartWithdraw'}, ...
        {'GoSignal','StartReach', 'StartLift', 'StartWithdraw'},...
        {'GoSignal','StartReach', 'StartHold', 'StartWithdraw'},...
        {'GoSignal', 'StartReplaceHold', 'StartReplaceSuccess', 'StartReward'}});
else
    PSTHAlignments = alignmentIn;
end
%%% assign parameters
pVal = 0.05;
movementConds = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};
% get activity map for current condition.
spaceInds = cellfun(@(c) [0, regexp(c, '\s')], [movementConds,'Rest'], 'UniformOutput', false);
condAbbrev = cellfun(@(c,i) c(i+1),[movementConds,'Rest'], spaceInds, 'UniformOutput', false);
condAbbrev{3} = 'PC';

sessionDirPrefix = ['S:\Lab\', monkey,'\'];
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'mm_dd_yyyy';
else
    dateFormat = 'yyyy_mm_dd';
end
% PSTH parameters and alignment
binSize=.01; % bin size in seconds
sigmaSize = .15; % smoothing window in seconds
sigma = round(sigmaSize/binSize);
secondsBeforePSTHAlignmentPoint = -6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 5; % time (in seconds) after alignement point to create PSTH
inclusiveBins = secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values
bins = inclusiveBins(1:end-1);

% phase window parameters
baselineAlignmentRange = {'GoSignal', 'GoSignal'};
taskAlignmentRange = {'GoSignal', 'StartHold'};
phaseAlignments = {'GoSignal','StartReach','StartLift'};
phaseWindowSize = .20; % size of window for phases in seconds
baselineWindow = [-3,-1];
taskWindow = [0, 0];
% representation information from voronoi
squarePoolSize = 56; %28 = 0.5mm, 1.0mm^2; 42 = 0.75 mm, 2.2mm^2; 56 = 1mm, 3.4mm^2;

erToName = containers.Map({'Digit', 'Wrist', 'Triceps', 'Biceps', 'Elbow', 'Forearm', 'Shoulder',...
    'Face','Ear', 'Trunk', 'Neck', 'No response'},{'Hand', 'Hand', 'Arm', 'Arm','Arm', 'Arm', 'Arm', ...
    'Face', 'Face', 'Trunk', 'Trunk', NaN});
erAbbrevs = erToName.keys;
repColors = struct('Arm', [.75 .75 .75], 'Hand', [.35 .35 .35], 'Forelimb',...
    [.5 .5 .5], 'Mixed', [121 76 92]./255,'Trunk', [.5 .25 .1],...
    'Face',[155 9 144]./255,'Axial', [102 0 36]./255);
forelimbNames = {'Arm','Hand'};
%%% all mapped sites and x,y locations
[~,~,rawRef]=xlsread(['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx']);
heading = rawRef(1,:);
rawRef(1,:) = [];
[~,recordingInd] = find(strcmp(heading, 'Recording'));
[~,siteInd] = find(strcmp(heading, 'Site #'));
[~,threshInd] = find(strcmp(heading, 'Threshold(s)'));
[~,erInd] = find(strncmp(heading, 'Evoked',length('Evoked')));
[~,xInd] = find(strcmp(heading, 'x'));
[~,yInd] = find(strcmp(heading, 'y'));
recordingSiteInds = find(cellfun(@(a,b) strcmp(a,'Yes') & isnumeric(b),...
    rawRef(:,recordingInd),rawRef(:,siteInd)));
% store information for each site
siteDateMap = array2table(cellfun(@str2double, string(rawRef(recordingSiteInds,...
    [siteInd,xInd,yInd]))),'VariableNames', {'Site', 'x', 'y'});
% mapped sites to determine identity of recorded sites
mappedSiteInds = find(~cellfun(@(e) (length(e)==1 && isnan(e)) || strcmp(e, 'No response'), rawRef(:,erInd)));
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
    mappedERs = cellfun(@(e) find(cellfun(@(ab) ~isempty(strfind(strip(e),ab)), erAbbrevs),1),currERs);
    allReps{r} = values(erToName,erAbbrevs(mappedERs));
    if(length(thresh)<length(mappedERs))
        disp(rawRef{mappedSiteInds(minMappedInd),1});
    end
    if(isnumeric(thresh))
        allThreshs{r} = thresh;
    else
        allThreshs{r} = cellfun(@(t) str2double(strip(t)), regexp(thresh, ',', 'split'));
    end
end
siteDateMap = addvars(siteDateMap,allReps, allThreshs,'NewVariableNames', {'SiteRep', 'Thresh'});
clear rawNum rawRef
%%% confirm M1 sites based on assigned border
MM = double(imread([sessionDirPrefix,'Mapping\Motor Maps V3\MM_Simp_RGB-01.png']));
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
M1sites = logical(arrayfun(@(a,b) M1mask(b,a), siteDateMap.x,siteDateMap.y));
siteDateMap(~M1sites,:) = [];
clear row col pointsX pointsY MM M1border M1mask
%%% date and session information (double check M1 mask)
[~,~,dateToPNValue] = xlsread([sessionDirPrefix, monkey, 'Session.xlsx']);
[~, dateInd] = find(strcmp(dateToPNValue, 'Date'));
[~,domainInd] = find(strcmp(dateToPNValue, 'Domain'));
[~,PNInd] = find(strcmp(dateToPNValue, 'Site #'));
% clean up, remove: duplicate test sites (close/far), non-motor recordings,
% improper logging (no site #)
validInds = find(cellfun(@(a,b,c) ~any(isnan(a)) & any(strcmp(b, {'M1', 'PMd', 'PMv'}))...
    & ~ismissing(string(c)),dateToPNValue(:,dateInd),dateToPNValue(:,domainInd),...
    dateToPNValue(:,PNInd), 'UniformOutput', true));
% match session numbers from _MM_Sites.xlsx sheet to session.xlsx sheet
[~,viPN,viSD] = intersect(cellfun(@str2double, string(dateToPNValue(...
    validInds,PNInd))),siteDateMap.Site);
validInds = validInds(viPN);
% assign session number's date, add domain, duplicate session number from
% session.xlsx sheet
siteDateMap = addvars(siteDateMap(viSD,:), dateToPNValue(validInds,dateInd), ...
    dateToPNValue(validInds,domainInd), dateToPNValue(validInds,PNInd),...
    'NewVariableNames', {'Date', 'Domain', 'Site2'});
clear dateToPNValue
%%% use colored voronoi for joint labels
MMColorReferenceName = [sessionDirPrefix,'Mapping\Motor Maps V3\', monkey,...
    ' MM All ', num2str(squarePoolSize),'px_RNR.mat'];
poolCircle = fspecial('disk',squarePoolSize);
MMColorReference = load(MMColorReferenceName);
if(strcmp(monkey,'Skipper'))
    MMColorReference = MMColorReference.MM_image_flat;    
else
    MMColorReference = MMColorReference.MM_joint_image_flat2;
end
%%% activity maps from imaging
allActivityFiles = dir([sessionDirPrefix, 'Mapping\tTests\*HSV.png']);
allActivityMaps = arrayfun(@(a) imresize(imread([a.folder, '\', a.name]),...
    [size(MMColorReference,1),size(MMColorReference,2)]), allActivityFiles, 'UniformOutput', false);
% activity is red
allActivityMaps = containers.Map(cellfun(@(a) a(1:end-8), {allActivityFiles.name}, 'UniformOutput', false),...
    cellfun(@(a) a(:,:,1)>100 & a(:,:,2)<100 & a(:,:,3)<100, allActivityMaps, 'UniformOutput', false));
%%% vessel mask
vMask = imread(['S:\Lab\',monkey, '\Mapping\clean_mask_filled.bmp']);
vMask = logical(vMask(:,:,1)>200);
imMaskDims = size(vMask);
%% main script: load info from all sites
numSites = height(siteDateMap);
siteMasks = zeros(imMaskDims(1),imMaskDims(2),numSites);
if(~(exist('siteLocation','var') | exist('siteRep', 'var') | exist('siteThresh', 'var')))
    [siteLocation, siteRep, siteThresh] = deal(cell(1,200));
end
if(~(exist('masterPSTH', 'var') | exist('masterSegs', 'var')))
    [masterPSTH,masterSegs,masterChannels,masterFano,masterXCorr,masterActivity,...
        baselineFR,phaseFR,unitPhaseInds,taskFR,reachFR,taskUnitInds,masterPhaseBins,allPSTHS] = deal({});
end
totalUnits = 0;
for i = 1:numSites
    currSession = siteDateMap(i,:);
    if(~isempty(dir([sessionDirPrefix,'All Data\', monkey,'_',...
            datestr(currSession.Date, dateFormat),'\Physiology\Results'])))
        [spikes,times,weights,currTrials,sessionConds,channels,events,~] =...
            getSessionInfo2([sessionDirPrefix,'All Data\', monkey,'_',...
            datestr(currSession.Date, dateFormat)],singleOrAll);
    else
        spikes = [];
    end
    if(~isempty(spikes))
        if(strcmp(sessionOrUnitPSTHS, 'Unit'))
            numUnits = size(spikes,1);
        else
            numUnits = 1;
        end
        totalUnits = totalUnits + numUnits;
        if(strcmp(conditions, 'All'))
            conditions = sessionConds;
            %conditions(contains(conditions,'Rest')) = [];
        end
        if(length(masterPSTH)<length(conditions))
            missingInitalization = length(conditions)-length(masterPSTH);
            masterPSTH(end+1:end+missingInitalization) = ...asrepmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            masterSegs(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            allPSTHS(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            masterFano(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            baselineFR(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            taskFR(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
                        reachFR(end+1:end+missingInitalization) = repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);

            masterPhaseBins(end+1:end+missingInitalization) =repmat({cell(1,length(PSTHAlignments))},1,missingInitalization);
            phaseFR(end+1:end+missingInitalization) = repmat({cell(1,length(phaseAlignments))},1,missingInitalization);
            unitPhaseInds(end+1:end+missingInitalization) = repmat({cell(1,length(phaseAlignments))},1,missingInitalization);
            
            taskUnitInds(end+1:end+missingInitalization) = cell(1,missingInitalization);
            masterActivity(end+1:end+missingInitalization) = cell(1,missingInitalization);
            masterChannels(end+1:end+missingInitalization) = cell(1,missingInitalization);
            masterXCorr(end+1:end+missingInitalization) = cell(1,missingInitalization);
        end
        for c = 1:length(conditions)
            activation = allActivityMaps(condAbbrev{c});
            condParamInd = cellfun(@(a) contains(a,conditions{c}),sessionConds);
            condInds = cellfun(@(a) contains(a,conditions{c}),currTrials(:,1))';
            % generate and smooth PSTHS  from current session
            if(any(condInds) & ~isempty(spikes))
                condEvents = events(conditions{c});
                condAlign = cellfun(@(a) find(strcmp(condEvents,a)),PSTHAlignments(conditions{c}),'UniformOutput', true);
                if(strcmp(conditions{c},'Rest'))
                    condPhaseAlignments = find(strcmp(condEvents, phaseAlignments{1}));
                    condWindows = {[-phaseWindowSize 0]};
                    taskAlignmentRange = {'GoSignal', 'StartReplaceSuccess'};
                                        reachAlignmentRange = taskAlignmentRange;

                else
                    if(strcmp(conditions{c},'Photocell'))
                        condPhaseAlignments = cellfun(@(e) find(strcmp(...
                            condEvents,replace(e,'Lift','Hold'))),...
                            phaseAlignments, 'UniformOutput', true);
                    else
                        condPhaseAlignments = cellfun(@(e) find(strcmp(...
                            condEvents,e)),phaseAlignments,'UniformOutput',true) ;
                    end
                    condWindows = phaseWindows;
                    taskAlignmentRange = {'GoSignal', 'StartHold'};
                    reachAlignmentRange = {'GoSignal', 'StartHold'};
                end
                taskAlignmentRange = {'GoSignal', 'StartReplaceHold'};
                condWeights = weights(:, condParamInd);
                %
                % {alignedPSTHS}{units,trials}
                alignedSpikes = arrayfun(@(ap) cellfun(@(s,t) s-t(ap), ...
                    spikes(:,condInds),repmat(times(condInds),size(spikes,1),1),...
                    'UniformOutput',false),condAlign,'UniformOutput',false);
                alignedTimes = arrayfun(@(ap) cellfun(@(t) t-t(ap), times(condInds), ...
                    'UniformOutput', false), condAlign, 'UniformOutput', false);
                
                baselineAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
                    baselineAlignmentRange,'UniformOutput', true);
                taskAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
                    taskAlignmentRange, 'UniformOutput', true);
                baselineRange = cellfun(@(at) cellfun(@(a) findBins(bins,...
                    a(baselineAlignment)+baselineWindow),at,'UniformOutput', false),...
                    alignedTimes, 'UniformOutput', false);
                taskRange = cellfun(@(at) cellfun(@(a) findBins(bins,...
                    a(taskAlignment)+taskWindow),at,'UniformOutput', false),...
                    alignedTimes,'UniformOutput', false);
                reachAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
                    reachAlignmentRange, 'UniformOutput', true);
                reachRange = cellfun(@(rt) cellfun(@(a) findBins(bins,...
                    a(reachAlignment)), rt, 'Uniformoutput', false),...
                    alignedTimes, 'UniformOutput', false);
                
                
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
                % determine unit modulation type(s)
                if(strcmp(sessionOrUnitPSTHS, 'Unit'))
                    %{alignedPSTHS}{trial #}{unit #, bins}
                    PSTHCalc = cellfun(@(p) squeeze(num2cell(p,...
                        [1 2]))', unitPSTHS,'UniformOutput', false);
                else
                    PSTHCalc = cellfun(@(u) squeeze(nanmean(u,1)),...
                        unitPSTHS,'UniformOutput',false);
                end
                
                %{alignedPSTHS}{trials}{units,bins}
                NUMSAMPLES = 10;
                baselineDist = cellfun(@(pc,bsr) cellfun(@(p,sr) cell2mat(arrayfun(@(s) ...
                    trapz(p(:,s:(s+(phaseWindowSize/binSize))),2),randi(max([...
                    find(sr,1,'first'), find(sr,1,'last')-(phaseWindowSize/binSize);1,1],...
                    [],1),[1,NUMSAMPLES]),'UniformOutput', false)),pc,bsr,...
                    'UniformOutput',false), PSTHCalc, baselineRange, 'UniformOutput', false);
                baselineSamp = cellfun(@(bd) num2cell(cell2mat(cellfun(@(b)...
                    nanmean(b,2), bd, 'Uniformoutput', false)),2),baselineDist ,'UniformOutput', false);
                
                %                 baselineSamp = cellfun(@(pc,bsr) num2cell(cell2mat(cellfun(@(p,sr)...
                %                     nanmedian(cell2mat(),pc,bsr,'UniformOutput',false)),2),PSTHCalc,baselineRange, 'UniformOutput', false);
                baselineTaskLength = cellfun(@(pc,bsr,tr) num2cell(cell2mat(cellfun(@(p,sr,t) ...
                    trapz(p(:,(max(1,find(sr,1,'last')...
                    -sum(t)):1:find(sr,1,'last'))),2),pc, bsr,tr,'UniformOutput',false)),2),...
                    PSTHCalc, baselineRange,taskRange,'UniformOutput',false);
                baselineReachLength = cellfun(@(pc,bsr,tr) num2cell(cell2mat(cellfun(@(p,sr,t) ...
                    trapz(p(:,(max(1,find(sr,1,'last')...
                    -sum(t)):1:find(sr,1,'last'))),2),pc, bsr,tr,'UniformOutput',false)),2),...
                    PSTHCalc, baselineRange,reachRange,'UniformOutput',false);
                %
                taskT = cellfun(@(up,ui,bi) num2cell(cell2mat(cellfun(@(p,i)...
                    trapz(p(:,i),2),up,ui,'UniformOutput',false))-cell2mat(bi),2),...
                    PSTHCalc,taskRange,baselineTaskLength, 'UniformOutput', false);
                reachT = cellfun(@(up,ui,bi) num2cell(cell2mat(cellfun(@(p,i)...
                    trapz(p(:,i),2),up,ui,'UniformOutput',false))-cell2mat(bi),2),...
                    PSTHCalc,reachRange,baselineReachLength, 'UniformOutput', false);
                if(numUnits==1)
                    taskT = cellfun(@(t) t', taskT,'UniformOutput',false);
                    reachT = cellfun(@(t) t', reachT, 'UniformOutput', false);
                end
                taskFR{c} = cellfun(@(tf,af) vertcat(tf, af), taskFR{c}, ...
                    [taskT,repmat({repmat({NaN},1,numUnits)},1,length(PSTHAlignments)-length(condAlign))], 'UniformOutput', false);
                reachFR{c} = cellfun(@(rf, af) vertcat(rf,af), reachFR{c},...
                    [reachT,repmat({repmat({NaN},1,numUnits)},1,length(PSTHAlignments)-length(condAlign))], 'UniformOutput', false);
                baselineFR{c} =  cellfun(@(bf, af) vertcat(bf,af), baselineFR{c}, ...
                    [baselineSamp, repmat({repmat({NaN},1,numUnits)},1,length(PSTHAlignments)-length(condAlign))],'UniformOutput', false);
                
                %%
                pTypes = repmat({[]},1,length(phaseAlignments));
                pFR = repmat({{}},1,length(phaseAlignments));
                pBins = repmat({{}},1,length(PSTHAlignments));
                %%
                for pw = 1:length(condPhaseAlignments)
                    pTFR = [];
                    % determine task modulated units
                    % unit assignment based on phase windows
                    phaseRange = cellfun(@(at) cellfun(@(a) findBins(bins,...
                        a(condPhaseAlignments(pw))+condWindows{pw}),at,'UniformOutput',...
                        false), alignedTimes, 'UniformOutput', false);
                    pTFR = cellfun(@(pc,pi) num2cell(cell2mat(cellfun(@(p,i) trapz(...
                        p(:,i),2),pc,pi,'UniformOutput',false)),2),...
                        PSTHCalc,phaseRange,'UniformOutput', false);
                    pFR{pw} = cellfun(@(pb,bb) num2cell(cell2mat(cellfun(@(p,b) p-b,...
                        pb,bb, 'UniformOutput', false)),2),pTFR,baselineSamp,'UniformOutput', false);
                    pBins =cellfun(@(p) num2cell(cell2mat(p'),2), phaseRange,'UniformOutput', false);
                    for u = 1:numUnits
                        if(pw==1)
                            taskBaseline = cell2mat(cellfun(@(b) b{u}, baselineTaskLength,'UniformOutput',false)')';
                            taskEncoding = cell2mat(cellfun(@(p) p{u}, taskT,'UniformOutput',false)')';
                            [~,pTTask] = ttest2(taskBaseline(~isnan(taskBaseline)), ...
                                taskEncoding(~isnan(taskEncoding)));
                            if(any(pTTask<pVal/length(alignedTimes)))
                                taskUnitInds{c}(end+1) = true;
                            else
                                taskUnitInds{c}(end+1) = false;
                            end
                        end
                        [~,closestAlignment] = min(abs(condAlign-condPhaseAlignments(pw)));
                        [~,pT] = ttest2(cell2mat(cellfun(@(b) b{u}, baselineSamp,'UniformOutput',false)')',...
                            cell2mat(cellfun(@(p) p{u}, pTFR,'UniformOutput',false)')');
                        pTypes{pw}(end+1) = pT(closestAlignment)<pVal/length(alignedTimes);
                        
                    end
                end
                %%
                masterPhaseBins{c} = cellfun(@(pb,p) vertcat(pb,p),...
                    masterPhaseBins{c}, pBins, 'UniformOutput', false);
                phaseFR{c}(end+1,:) =  pFR;
                %                 unitPhaseInds{c} = cellfun(@(mT,t) vertcat(mT,t),...
                %                     unitPhaseInds{c},num2cell(cell2mat(pTypes),1),'UniformOutput', false);
                unitPhaseInds{c}(end+1,:) = pTypes;
                % collect data from all sessions
                %                 siteSampleActivity = activation(max(1,currSession.y-squarePoolSize):...
                %                     min(currSession.y+squarePoolSize,size(activation,2)),...
                %                     max(1,currSession.x-squarePoolSize):...
                %                     min(currSession.x+squarePoolSize,size(activation,2)));
                if(strcmp(sessionOrUnitPSTHS, 'Unit'))
                    masterPSTH{c} = cellfun(@(mp, np) [mp,np], ...
                        masterPSTH{c}, [cellfun(@(u) num2cell(nanmean(u,3),[2,3])',...
                        unitPSTHS,'UniformOutput', false),repmat({NaN},1,...
                        length(PSTHAlignments)-length(condAlign))],'UniformOutput', false);
                    masterChannels{c}(end+1:end+numUnits) = num2cell(channels);
                    %                     masterFano{c}= cellfun(@(mf,nf) [mf,nf], masterFano{c}, cellfun(@(p) squeeze(...
                    %                         var(cell2mat(p),0,2))./(0.01 + squeeze(nanmean(cell2mat(p),2)))',...
                    %                         PSTHCalc, 'UniformOutput', false),'UniformOutput', false);
                else
                    masterPSTH{c}(end+1) = {cellfun(@(u) nanmean(nanmean(u,3),1), unitPSTHS, 'UniformOutput', false)};
                    masterChannels{c}(end+1) = {[]};
                    %                     masterFano{c}(i) =  {cellfun(@(p) var(cell2mat(p),0,2)./...
                    %                         (0.01 + nanmean(cell2mat(p),2)), PSTHCalc, 'UniformOutput', false)'};
                end
                masterSegs{c}= cellfun(@(at,nt) [at,nt], masterSegs{c},[...
                    cellfun(@(at) repmat({cell2mat(at')},1,numUnits),...
                    alignedTimes,'UniformOutput', false), repmat({repmat({NaN},1,numUnits)},...
                    1, length(PSTHAlignments)-length(condAlign))],'UniformOutput', false);
                allPSTHS{c}{end+1} = unitPSTHS;
            else
                % pad stored info with empty arrays and NaN pad indicies for missing conditions
                masterPSTH{c} = cellfun(@(mp,ep) [mp,ep], masterPSTH{c},...
                    repmat({repmat({NaN(size(bins))},1,numUnits)},1,...
                    length(PSTHAlignments)),'UniformOutput', false);
                masterChannels{c}(end+1:end+numUnits) = {missing};
                %                 masterFano{c}(end+1:end+numUnits) = {missing};
                masterSegs{c} = cellfun(@(mp,ep) [mp,ep], masterSegs{c},...
                    [repmat({repmat({NaN(sum(cellfun(@(ct) strcmp(ct,conditions{c}),....
                    currTrials(:,1))), length(allEvents(conditions{c})))},...
                    1,numUnits)},1,length(PSTHAlignments))],'UniformOutput', false);
                allPSTHS{c}{end+1} = repmat({NaN(numUnits,length(bins),length(currTrials))}, 1, length(PSTHAlignments));
                
                
                %masterXCorr{c}(end+1:end+numUnits) = {missing};
                taskUnitInds{c}(end+1:end+numUnits) = false;
                unitPhaseInds{c}(end+1,:) = repmat({NaN(1,numUnits)}, 1, length(phaseAlignments));
                
                reachFR{c} = cellfun(@(tf,af) vertcat(tf, af), taskFR{c}, ...
                    repmat({repmat({NaN(1,sum(cellfun(@(ct) strcmp(ct,conditions{c}),....
                    currTrials(:,1))))},numUnits,1)},1,length(PSTHAlignments)), 'UniformOutput', false);
              
                phaseFR{c}(end+1,:) = repmat({repmat({repmat({NaN(1,sum(cellfun(@(ct) strcmp(ct,conditions{c}),....
                    currTrials(:,1))),1)},numUnits,1)},1,length(PSTHAlignments))},1,length(phaseAlignments));
                masterPhaseBins{c} = cellfun(@(pb,p) vertcat(pb,p),...
                    masterPhaseBins{c}, repmat({NaN(2,length(phaseAlignments))}, 1, length(PSTHAlignments)),...
                    'UniformOutput', false);
            end
            masterActivity{c}(end+1:end+numUnits) = repmat(activation(...
                currSession.y,currSession.x),1,numUnits);
            %activityInd{c}(end+1:end+numUnits) = repmat(sum(...
            %siteSampleActivity(:))/length(siteSample(:)) > .25,1,numUnits);
            
            if(strcmp(datestr(currSession.Date, dateFormat),'05_17_2019'))
                disp('single');
            end
        end
        currSite = [currSession.x, currSession.y];
        % get current session joint label
        siteRep{i}(end+1:end+numUnits) = repmat(cellfun(@cellstr, currSession.SiteRep,'UniformOutput', false), 1, numUnits);
        siteThresh{i}(end+1:end+numUnits)= repmat(currSession.Thresh,1,numUnits);
        siteLocation{i}(end+1:end+numUnits) = num2cell(repmat(currSite,numUnits,1),2);
        % voronoi tiling
        tempCircle = zeros(imMaskDims(1)+400,imMaskDims(2)+400);
        tempCircle((currSite(2)-squarePoolSize+200):(currSite(2)+squarePoolSize+200),...
            (currSite(1)-squarePoolSize+200):(currSite(1)+squarePoolSize+200))...
            = poolCircle;
        tempCircle = tempCircle(200:end-201,200:end-201);
        siteMasks(:,:,i) = tempCircle>0;
    end
    disp([num2str(currSession.Site),', ',datestr(currSession.Date,dateFormat),...
        '.  Session ', num2str(i) ,' of ' num2str(numSites)]);
end
%% remove sessions that had no trial information
emptyInds = cellfun(@isempty, siteLocation);
siteLocation = siteLocation(~emptyInds);
siteRep = siteRep(~emptyInds);
siteThresh = siteThresh(~emptyInds);
siteLocation = horzcat(siteLocation{:});
siteRep = horzcat(siteRep{:});
siteThresh = horzcat(siteThresh{:});
% NAN pad conditions with fewer trials
%
emptyCondInds = cellfun(@(c) cellfun(@(in) isa(in,'missing'),c, 'UniformOutput', true),masterPSTH, 'UniformOutput', false);
masterPSTH = padMissingCells(masterPSTH,emptyCondInds);
% masterFano = padMissingCells(masterFano,emptyCondInds);
masterChannels = padMissingCells(masterChannels,emptyCondInds);
% baselineFR = padMissingCells(baselineFR(~strcmp(conditions,'Rest')),emptyCondInds(~strcmp(conditions,'Rest')));
% assign all units to corresponding session
smInd = false(1,size(siteMasks,3));
for sm = 1:size(siteMasks,3)
    smInd(sm) = any(any(siteMasks(:,:,sm)));
end
siteMasks = siteMasks(:,:,smInd);
sessionInds = logical([1; any(diff( cell2mat(siteLocation')),2)]);
unitSessionLabel = cumsum(sessionInds)';
[verticies, vCells] = voronoin([cell2mat(siteLocation(sessionInds)'); ...
    [0 size(vMask,2); size(vMask,1) 0; 0 0; size(vMask,1) size(vMask,2)]]);
vXY = cellfun(@(v) [verticies(v,2),verticies(v,1)], vCells, 'UniformOutput', false);
% group multiple representation sites/units accordingly
remappedReps = cell(size(siteRep));
for s = 1:length(siteRep)
    if(length(siteRep{s})==1)
        remappedReps(s) = siteRep{s};
    else
        allSiteReps = unique(siteRep{s}(siteThresh{s} <= min(siteThresh{s})*1.1));
        if(all(contains(allSiteReps, forelimbNames)))
            remappedReps{s} = 'Forelimb';
        elseif(any(contains(allSiteReps, forelimbNames)))
            remappedReps{s} = 'Mixed';
        else
            remappedReps{s} = 'Axial';
        end
    end
end
simpRep = cellfun(@(r,t) r{t==min(t)}, siteRep,siteThresh, 'UniformOutput', false);
%
saveDirPath = 'S:\Lab\ngc14\Figures\Physiology\Results\New\';
saveName = [monkey, '_',sessionOrUnitPSTHS,'_',singleOrAll];
% taskUnitInds = cellfun(@(m) true(1,length(m)), masterPSTH,'UniformOutput',false);

%% plot joint PSTHS
close all;



jointName = unique(simpRep);

%    figure('Units','normalized', 'Position',[0 0 1 1])
%
close all;
bothSphereCondInds = (~(cellfun(@(t) all(isnan(t)), taskFR{1}{2}) | ...
    cellfun(@(t) all(isnan(t)), taskFR{2}{2})));
bothSphereCondInds = bothSphereCondInds & cellfun(@(s,c) ~(all(s==0) & all(c==0)), taskFR{1}{2}, taskFR{2}{2});
spInds = false(1,length(bothSphereCondInds));
spInds(bothSphereCondInds) = cellfun(@(ain,bin) ttest2(ain(~isnan(ain)),bin(~isnan(bin)), 'alpha', pVal/2),...
    taskFR{1}{2}(bothSphereCondInds), taskFR{2}{2}(bothSphereCondInds));

phaseFR = cellfun(@(a) a(2:end,:), phaseFR, 'UniformOutput', false)
unitPhaseInds = cellfun(@(a) a(2:end,:), unitPhaseInds, 'UniformOutput', false)

currRep = simpRep;
currRep(~(taskUnitInds{1} | taskUnitInds{2})) = {''};
spInds = double(spInds);
spInds(~(taskUnitInds{1} | taskUnitInds{2})) = NaN;
unitSP = modulatedUnitsPerRep(simpRep,int32(spInds),"Sphere modulated units","Spheres",repColors);
figMapSpheres = mapUnitVals(vXY,vMask,siteMasks,sessionInds,int32(spInds),true,5,[1,10]);
saveFigures(unitSP,[saveDirPath,'Count\Cond\'],[saveName,'_Sphere_modulated_units'],[]);
saveFigures(figMapSpheres,[saveDirPath,'Maps\Count\Cond\'],[saveName,'_Sphere_modulated_units'],[]);
%%
for c = 1:length(conditions)
    condEvents = allEvents(conditions{c});
    alignments = cellfun(@(a) find(strcmp(condEvents,a)),PSTHAlignments(conditions{c}),'UniformOutput', true);
    currRep = simpRep;
    
    maskMatrix = ones(length(siteLocation),length(bins));
    if(c==4)
                            reachAlignmentRange = {'GoSignal', 'ReplaceSuccess'};

        tks = (taskUnitInds{1} ~= 0| taskUnitInds{2}~=0 | taskUnitInds{3} ~=0);
        currRep(~(tks)) = {''};
        maskMatrix(~tks,:) = NaN(sum(~tks),length(bins));
        NaNPSTHS{c} = cellfun(@(n) num2cell(vertcat(n{:}).*maskMatrix,2)', masterPSTH{c}, 'UniformOutput', false);
        plotJointPSTHS(bins,{masterPSTH{c}{1}},{masterSegs{c}{1}},currRep,masterActivity{c},...
            alignments,{[-3 3]});
        phaseAlignments = phaseAlignments(1);
        reachAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
            {'GoSignal', 'ReplaceSuccess'}, 'UniformOutput', true);
    else
                            reachAlignmentRange = {'GoSignal', 'StartHold'};
        tks = taskUnitInds{c};
        currRep(~tks) = {''};
        maskMatrix(~tks,:) = NaN(sum(~tks),length(bins));
                NaNPSTHS{c} = cellfun(@(n) num2cell(vertcat(n{:}).*maskMatrix,2)', masterPSTH{c}, 'UniformOutput', false);
                plotJointPSTHS(bins,NaNPSTHS{c},masterSegs{c},currRep,...
                    masterActivity{c}.*tks,alignments,alignLimits);
        reachAlignment = cellfun(@(a) find(strcmp(condEvents,a)),...
            reachAlignmentRange, 'UniformOutput', true);
    end
    saveFigures(gcf,[saveDirPath,'PSTH\',condAbbrev{c},'\'],saveName,[]);
    
%     
%     if(c<4)
%         rSpeed =cellfun(@(p) cellfun(@(a) 1./(abs(diff([a(:,reachAlignment)],1,2))), ...
%             p, 'UniformOutput', false), masterSegs{c}, 'UniformOutput', false);
%         rFR = cellfun(@(t) t, reachFR{c},'UniformOutput', false);
%         allTrialReps = vertcat(cellfun(@(r,s) repmat(string(r),length(s),1), simpRep,rSpeed{2},'UniformOutput',false));
%         allTrialReps = vertcat(allTrialReps{:});
%         [unitB,~,unitR,~,unitStats] = cellfun(@(us,uf) regress(uf', [ones(length(us),1), vertcat(us)]), rSpeed{2},reachFR{c}{2}', 'UniformOutput', false);
%         unitCorrs = cellfun(@(us,uf) corrcoef(us, uf'), rSpeed{2}, rFR{2}', 'UniformOutput', false);
%         unitCorrs = cellfun(@(r) r(2), unitCorrs, 'UniformOutput', false);
%         unitStats =cellfun(@(s) s(1), unitStats, 'UniformOutput', true);
%         
%         scatter(vertcat(rSpeed{2}{:}), [rFR{2}{:}]',36,cell2mat(cellfun(@(s) repColors.(s), allTrialReps, 'Uniformoutput', false)),'filled');
%         hold on;
%         ylabel('FR');
%         xlabel('Speed');
%         allCorr = corrcoef(vertcat(rSpeed{2}{:}), [rFR{2}{:}]');
%         title(['All units: r^2 = ', num2str(allCorr(2),3), '; mean r^2 = ',...
%             num2str(nanmean(unitStats),3)]);
% 
%         pAll = polyfit(vertcat(rSpeed{2}{:}),[(rFR{2}{:})]',1);
%         fAll = polyval(pAll,vertcat(rSpeed{2}{:}));
%         plot(vertcat(rSpeed{2}{:}),fAll,'k-')
%         saveFigures(gcf,[saveDirPath,'Reach_Speed\',condAbbrev{c},'\'],'All_Units',[]);
%         close all;
%         for j = 1:length(jointName)
%             figure(); hold on;
%             currJ = cellfun(@(s) strcmp(s, jointName{j}), simpRep, 'UniformOutput', true);
%             scatter(vertcat(rSpeed{2}{currJ}), [rFR{2}{currJ}]',36,repColors.(jointName{j}),'filled');
%             [b,~,r,~,stats] = regress([rFR{2}{currJ}]', ...
%                 [ones(sum(cellfun(@length,rSpeed{2}(currJ))),1), vertcat(rSpeed{2}{currJ})])
%             ylabel('FR');
%             xlabel('Speed');
%             title([jointName{j}, ' units: r^2 = ', num2str(stats(1),3),...
%                 '; mean r^2 = ',num2str(nanmean(unitStats(currJ)),3), '; corrcoef: = ', ...
%                 num2str(nanmean(unitStats(currJ)),2)]);
%             p = polyfit(vertcat(rSpeed{2}{currJ}),[(rFR{2}{currJ})]',1);
%             f = polyval(p,vertcat(rSpeed{2}{currJ}));
%             plot(vertcat(rSpeed{2}{currJ}),f,'k-')
%             saveFigures(gcf,[saveDirPath,'Reach_Speed\',condAbbrev{c},'\'],[jointName{j},'_Units'],[]);
%             close all
%         end
%         unitStats(~tks) = NaN;
%         rFig = modulatedUnitsPerRep(currRep,unitStats,...
%             "r-squared", "Task Units", repColors);
%         rMapFig = mapUnitVals(vXY,vMask,siteMasks,sessionInds,...
%             unitStats,false,5,[]);
%         saveFigures(rFig, [saveDirPath,'\Reach_Speed\', condAbbrev{c},'\'],[saveName],[])
%         saveFigures(rMapFig,[saveDirPath,'Maps\Reach_Speed\',condAbbrev{c},'\'],...
%             [saveName],[]);
%     end
% 
% 
%     figure('Units','normalized', 'Position',[0 0 1 1]);
%     hold on;
%     raSpeed = (vertcat(rSpeed{2}{:}));
%     raFR = horzcat(rFR{2}{:})';
%     for j = 1:length(jointName)
%             jointInds = cell2mat(cellfun(@(a,s) repmat(strcmp(a,jointName{j}), size(s,1),1), simpRep,masterSegs{1}{2},'UniformOutput', false)');
%             subplot(2,2,j);hold on;
%             scatter(raSpeed(jointInds), raFR(jointInds),36,'MarkerEdgeColor',...
%                 repColors.(jointName{j}));
%             rVal = corr(raSpeed(jointInds), raFR(jointInds));
% 
%     p = polyfit(raSpeed,raFR,1);
%     f = polyval(p,raSpeed);
%      plot(raSpeed,f,'k-')
% 
%             title([jointName{j}, ': r = ', num2str(rVal,3)]);
%     end
% 
% 
% 





if(strcmp(monkey,'Gilligan'))
    colorRange = [1 20];
else
    colorRange = [1 30];
end
for pw = 1:length(phaseAlignments)
    if pw==1
        phaseName = 'Go';
    elseif pw==2
        phaseName = 'Reach';
    else
        phaseName = 'Grasp';
    end
    % AUC by representation
    sessionVals = cell2mat(cellfun(@(pp) nanmean(cell2mat(pp),2),...
        cellfun(@(cl) cl{2}, phaseFR{c}(:,pw), 'UniformOutput',false),'UniformOutput', false));
    sessionVals(~tks) = NaN;
    FRFig = modulatedUnitsPerRep(currRep,sessionVals',...
        strcat("AUC: ", phaseName, " phase"), "Task Units", repColors);
    FRMapFig = mapUnitVals(vXY,vMask,siteMasks,sessionInds,...
        sessionVals,false,5,[1 300]);
    saveFigures(FRFig, [saveDirPath,'\AUC\', condAbbrev{c},'\'],[saveName,'_', phaseName],[]);
    
    saveFigures(FRMapFig, [saveDirPath,'\Maps\AUC\', condAbbrev{c},'\'],[saveName,'_', phaseName],[])
    
    % Count by representaiton
    phI = cell2mat(unitPhaseInds{c}(:,pw)');
    phI(tks==0) = NaN;
    unitFig = modulatedUnitsPerRep(currRep,int32(phI),...
        { string([phaseName, ' units'])}, string(condAbbrev{c}), repColors);
    mapFig = mapUnitVals(vXY,vMask,siteMasks,sessionInds,...
        int32(phI),true,5,colorRange);
    
    %Phase modulated PSTHS
    phaseJointFig =  plotJointPSTHS(bins, NaNPSTHS{c},masterSegs{c},...
        currRep,masterActivity{c},alignments,alignLimits);
    
    saveFigures(FRFig, [saveDirPath,'\AUC\', condAbbrev{c},'\'],[saveName,'_', phaseName],[]);
    saveFigures(phaseJointFig, [saveDirPath, '\PSTH\',...
        condAbbrev{c},'\'], [saveName,'_', phaseName],[]);
    saveFigures(unitFig,[saveDirPath,'Count\',condAbbrev{c},'\'],...
        [saveName,'_', phaseName, '_units'],[]);
    saveFigures(FRMapFig, [saveDirPath,'\Maps\AUC\', condAbbrev{c},'\'],[saveName,'_', phaseName],[])
    saveFigures(mapFig,[saveDirPath,'Maps\Count\',condAbbrev{c},'\'],...
        [saveName,'_', phaseName, '_units'],[]);
end
end
figMapG = mapUnitVals(vXY,vMask,vMask,sessionInds,int32(unitPhaseInds{2}.*tks),true,5,[1 10]);
saveFigures(figMapG,[saveDirPath,'Maps\Count\',condAbbrev{c},'\'],[saveName,'_Grasp_units'],[]);
figMapRGx = mapUnitVals(vXY,vMask,vMask,sessionInds,int32(unitPhaseInds{3}.*tks),true,5,[-5 5]);
saveFigures(figMapRGx,[saveDirPath,'Maps\Count\',condAbbrev{c},'\'],[saveName,'_Go_units'],[]);
saveFigures(gcf,[saveDirPath,condAbbrev{c},'\'],saveName,[]);
close all;
% phase analysis
% phaseAnalysis(repSave,sessionInds,phaseFR,taskFR,baselineFR,masterPSTH,...
%     masterPhaseBins,taskUnitInds,unitPhaseInds,condAbbrev,phaseNames,...
%     siteMasks,vXY,vMask,repColors,pVal,saveDirPath,saveName);
% close all;
end
%%

function bins = findBins(allBins, timesT)
if(~any(~isnan(timesT)))
    bins =false(size(allBins));
else
    bins = discretize(allBins,timesT)==1;
end
end

function paddedCell = padMissingCells(c,emptyIndex)
paddedCell = c;
for p = 1:length(c)
    sizeOfNonEmpty = cell2mat(cellfun(@size, paddedCell{p}(~emptyIndex{p}),'UniformOutput', false)');
    ty = cellfun(@class, paddedCell{p}(~emptyIndex{p}), 'UniformOutput', false);
    sameSize = ~any(diff(sizeOfNonEmpty));
    sizeOfNonEmpty = mode(sizeOfNonEmpty);
    sizeOfNonEmpty(~sameSize) = 1;
    if(strcmp(ty{1}, 'logical'))
        paddedCell{p}(emptyIndex{p}) = repmat(false(sizeOfNonEmpty),1,sum(emptyIndex{p}));
    else
        paddedCell{p}(emptyIndex{p}) = repmat({NaN(sizeOfNonEmpty)}, 1,sum(emptyIndex{p}));
    end
end
end