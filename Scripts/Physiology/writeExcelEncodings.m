clear all;
close all;

monkey = 'Skipper';
singleOrAll = 'Single';
PSTH_type = 'Unit'; % Unit or Trial
condsToAnalyze = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};

saveVals = 1;
saveFigs = 0;
pVal = 0.05; % p-value to test significance for units FR modulation in windows
binSize=.01; % bin size in seconds
sigma = .05; % smoothing window in seconds
sigma = round(sigma/binSize);
secondsBeforePSTHAlignmentPoint = 6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 7; % time (in seconds) after alignement point to create PSTH
bins = -secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values
phaseWindowSize = .15; % size of window for phases in seconds
baselineWindowSize = 2.4; % size of window for baseline calculation in seconds
baselineBufferTime = 1; % time after trial initiate to start calculating baseline activity in seconds
durationOfTrialInitiate = 5; % duration of time from trial initiate to go cue in seconds(check arduino files)
alignmentValue = 0; % for determining baseline (and rest)
alignmentPoint = 'StartReach';
savePath = ['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',singleOrAll,'\'];
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'mm_dd_yyyy';
else
    dateFormat = 'yyyy_mm_dd';
end
xCorrAligntments = {'GoSignal', 'StartReplaceHold'};
unitType = {'Go','Task', 'Reach', 'Grasp'};
unitAlignments = {{'GoSignal','GoSignal'}, {'StartReach', 'StartHold'}, ...
    {'StartReach', 'StartReach'}, {'StartLift', 'StartLift'}};
alignmentTimeRanges = {[0, phaseWindowSize], [-phaseWindowSize/2 0], ...
    [-phaseWindowSize/2, phaseWindowSize/2], [-phaseWindowSize, 0]};
phaseWindow = {'Go', 'Task', 'Reach', 'Grasp', 'Rest'};
encodingNames = cellfun(@(a) cellfun(@(b) [a, ' Units ', b, ' Phase FR Change'],...
    phaseWindow,'UniformOutput', false), [singleOrAll, unitType],'UniformOutput', false);
headingNames = ['Total Units', 'Single Units',cellfun(@(a) [a, ' Unit Counts'],unitType,'UniformOutput',false),...
    [encodingNames{:}]];
headingNames = num2cell(headingNames);
headingRow = 1;
excelFileNameRead = ['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx'];
[~,~,rawRef]=xlsread(excelFileNameRead);
numCols = size(rawRef,2);
[~,recordingInd] = find(strcmp(rawRef, 'Recording'));
[~,siteInd] = find(strcmp(rawRef, 'Site #'));
siteNumInds = find(cellfun(@(a,b) strcmp(a,'Yes'), rawRef(:,recordingInd)));
mappedSites = cell2mat(rawRef(headingRow+1:end,siteInd));
recordingSites = cell2mat(rawRef(siteNumInds,siteInd));
clear rawRef
[~,~,excelText] = xlsread(['S:\Lab\', monkey,'\', monkey, 'Session.xlsx']);
[~, dateInd] = find(strcmp(excelText, 'Date'));
[~,PNInd] = find(strcmp(excelText, 'Site #'));
[~,domainInd] = find(strcmp(excelText, 'Domain'));
cellInds = find(all(cell2mat(cellfun(@(a) cellfun(@(c) ~any(isnan(c)),a),...
    {excelText(:,dateInd),excelText(:,PNInd),excelText(:,domainInd)}, 'UniformOutput', false)),2));
excelText = excelText(cellInds,:);

validInds = find(cellfun(@(a) strcmp(a,'M1'), string(excelText(:,domainInd))) & ...
    cellfun(@isnumeric, excelText(:,PNInd)));
[allSites,indA,~] = intersect([excelText{validInds,PNInd}],recordingSites);
dates = excelText(validInds(indA),dateInd);

clear excelText excelSiteNums
%%
dispstat('','init');
dq = parallel.pool.DataQueue;
afterEach(dq, @updatedispstat);
updatedispstat(0,length(dates));
allVars={};
parfor s = 1:length(dates)
    %% load in PSTH variables
    currSessionFolder = ['S:\Lab\', monkey,'\All Data\', monkey,'_',datestr(dates{s}, dateFormat)];
    if(~isempty(dir([currSessionFolder,'\Physiology\Results\*.mat'])))
        [alignedSpikes,alignTimes,weights,sessionTrials,sessionConds,channels, events,labels] =...
            getSessionInfo(['S:\Lab\', monkey,'\All Data\', monkey,'_',datestr(dates{s}, dateFormat)], singleOrAll, alignmentPoint);
        if(~isempty(alignedSpikes))
            trialHists = cellfun(@(a) histcounts(a,bins)./binSize, alignedSpikes,'UniformOutput', false);
            trialHists = cellfun(@(a) conv(a,gausswin(sigma)/sum(gausswin(sigma)),'full'), trialHists, 'UniformOutput', false);
            trialHists = cellfun(@(a) a(1:end-(sigma-2)), trialHists, 'UniformOutput', false);
            allPSTHS=reshape(cell2mat(trialHists'), [size(trialHists,2), length(bins),size(trialHists,1)]);
            % movement conditions aligned to reach onset
            
            totalUnits = length(labels);
            singleUnits = sum(cellfun(@(a) strcmp(a,'s'), labels));
            numUnits = size(allPSTHS,3);
            
            unitMods =  repmat({false(1,numUnits)},length(condsToAnalyze), length(unitType));
            FRChanges = repmat({[]},length(condsToAnalyze), length(unitType));
            %% get start and end bins for baseline, rest, and phases to analyze
            baselineStartBin = cellfun(@(a) findBins(bins,alignmentValue-...
                (durationOfTrialInitiate - baselineBufferTime)),alignTimes);
            baselineEndBin = baselineStartBin + round(baselineWindowSize/binSize);
            
            restInds = cellfun(@(a) contains(a,'Rest'),sessionTrials(:,1));
            restStartBin = cellfun(@(a) findBins(bins,alignmentValue), alignTimes(restInds));
            restEndBin = restStartBin + round(baselineWindowSize/binSize);
            restCondInd = find(cellfun(@(a) strcmp(a,'Rest'),sessionConds));
            baselineFR = {};
            phaseFR = {};
            condFR = {};
            condModulated = {};
            % each condition in session to analyze
            for c = 1:length(condsToAnalyze)
                %% calculate number of unitType defined units in session
                conditionInd = find(cellfun(@(a) strcmp(a,condsToAnalyze{c}),sessionConds));
                condInds = cellfun(@(a) contains(a,condsToAnalyze{c}),sessionTrials(:,1));
                badSegInds = find(cellfun(@length, alignTimes(condInds))~=length(events(condsToAnalyze{c})));
                segInds = find(condInds);
                segInds(badSegInds) = [];
                avgPhaseBins = {};
                % plot and save PSTTHs
                if(saveFigs)
                    plotPSTH(['S:\Lab\', monkey,'All Data\', monkey,'_',datestr(dates{d}, dateFormat)],condsToAnalyze{c},...
                        allPSTHS(condInds,:,:), alignTimes(segInds), channels, bins);
                end
                % assign unit types for each unit in penetration
                for p = 1:length(unitAlignments)
                    startWindowInd = strcmp(events(condsToAnalyze{c}), unitAlignments{p}(1));
                    endWindowInd = strcmp(events(condsToAnalyze{c}),unitAlignments{p}(2));
                    phaseStartBins = cellfun(@(a) findBins(bins, a(startWindowInd)...
                        +alignmentTimeRanges{p}(1)), alignTimes(segInds));
                    phaseEndBins = cellfun(@(a) findBins(bins, a(endWindowInd)...
                        +alignmentTimeRanges{p}(2)), alignTimes(segInds));
                    [~, phaseFR{c,p}, baselineFR{c}] = calcFRChanges(allPSTHS(intersect(find(condInds),segInds),:,:),...
                        [baselineStartBin(intersect(find(condInds),segInds)); ...
                        baselineEndBin(intersect(find(condInds),segInds))]',...
                        [phaseStartBins; phaseEndBins]');
                    avgPhaseBins{c,p} = [round(nanmean(phaseStartBins)), round(nanmean(phaseEndBins))];
                    for u = 1:numUnits
                        unitEncoding = phaseFR{c,p}(:,u);
                        unitBaseline = baselineFR{c}(:,u);
                        [~,pT] = ttest2(unitBaseline(~isnan(unitBaseline)), unitEncoding(~isnan(unitEncoding)));
                        if(pT<pVal)
                            unitMods{c,p}(u) = 1;
                        end
                    end
                end
                avgBaselineWindow = [round(nanmean(baselineStartBin(segInds))), ...
                    round(nanmean(baselineEndBin(segInds)))];
                % cross correlation for units for all unit type definitions
                %         corrStartBin = cellfun(@(a) findBins(bins, a(strcmp(events{c},...
                %             xCorrAligntments{1}))), alignTimes(segInds));
                %         corrEndBin = cellfun(@(a) findBins(bins, a(strcmp(events{c},...
                %             xCorrAligntments{2}))), alignTimes(segInds));
                %         xCorrUnitTypes{c} = cellfun(@(a) corrcoef(permute(weights(a, conditionInd).*...
                %             permute(nanmean(allPSTHS(condInds,corrStartBin:corrEndBin,a),1),[3 2 1]),[2 1])),...
                %             unitMods(c,:),'UniformOutput', false);
                %         xCorrUnitTypes{c}([Inf,cellfun(@sum,unitMods(c,:))]<1) = deal({NaN});
                %% sesion caluclations
                unitMods(:,2:length(unitType)+1) = unitMods;
                unitMods(:,1) = deal({ones(1,numUnits)});
                % calculate average PSTH (depending on unit or trial PSTH)
                if(strcmp(PSTH_type, 'Unit'))
                    %  FRChange = cellfun(@(u) [cellfun(@(a) nanmean(nanmean(...
                    %  weights(u,condWeightInd)'.*a(:,u),1)), FRVals(c,:)), ...
                    %  nanmean(nanmean(weights(u,restWeightInd)'.*restFR(:,u),1))],...
                    %  unitMods(c,:), 'UniformOutput', false);
                    avgPSTHS = cellfun(@(u) squeeze(nanmean(weights(u, conditionInd).*...
                        permute(nanmean(allPSTHS(condInds,:,u),1),[3,1,2]),1))', unitMods(c,:), 'UniformOutput', false);
                    restPSTHS = cellfun(@(u) squeeze(nanmean(weights(u,restCondInd).*...
                        permute(nanmean(allPSTHS(restInds,:,u),1),[3,1,2]),1))', unitMods(c,:), 'UniformOutput', false);
                else
                    % FRChange = cellfun(@(u) [cellfun(@(a) nanmean(nanmean(...
                    % a(:,u),2)),FRVals(c,:)), nanmean(nanmean(restFR(:,u),2))],...
                    % unitMods(c,:),'UniformOutput', false);
                    avgPSTHS = cellfun(@(u) nanmean(squeeze(nanmean(allPSTHS(condInds,:,u),3)),1),...
                        unitMods(c,:), 'UniformOutput', false);
                    restPSTHS = cellfun(@(u) nanmean(squeeze(nanmean(allPSTHS(restInds,:,u),3)),1),...
                        unitMods(c,:), 'UniformOutput', false);
                end
                %         FRChange = cellfun(@(b) cellfun(@(a) nanmean(nanmean(b(a))), unitMods(c,:)),...
                %             [FRChanges,FRChangeRest], 'UniformOutput', false);
                [restFRChange,~,~,] = cellfun(@(p) calcFRChanges(p, ...
                    [nanmean(baselineStartBin(restInds)), nanmean(baselineEndBin(restInds))],...
                    [nanmean(restStartBin), nanmean(restEndBin)]),restPSTHS, 'UniformOutput', false);
                
                FRChange = cellfun(@(b,r) [cellfun(@(a) calcFRChanges(b,avgBaselineWindow,a),...
                    avgPhaseBins(c,:),'UniformOutput', false), r],avgPSTHS,restFRChange,'UniformOutput', false);
                condFR{c} = reshape(cell2mat([FRChange{:}]'),1,[]);
                condModulated{c} = cellfun(@sum,unitMods(c,2:end));
            end
            modulatedUnits = cell2mat(condModulated');
            sessionFRChanges = cell2mat(condFR');
            %     xCorrs = cell2mat(cellfun(@(a) cellfun(@(b) nanmedian(b(:)),a),xCorrUnitTypes,'UniformOutput', false)');
            allVars{s} = [repmat(totalUnits,length(condsToAnalyze),1), repmat(singleUnits,length(condsToAnalyze),1),...
                modulatedUnits,sessionFRChanges];%xCorrs];
        end
    end
    updatedispstat(s,length(dates));
    
end
unRunDates = cellfun(@isempty,allVars);
allVars = allVars(~unRunDates);
dates = dates(~unRunDates);
allSites = allSites(~unRunDates);
%% excel output
if(saveVals)
    dispstat('', 'init');
    startCol = char(numCols + 64 + 1);
    endCol = [char(floor((numCols + size(allVars{1},2))/26)+64),...
        char(mod((numCols + size(allVars{1},2)),26)+64)];
    for s = 1:length(dates)
        PN = allSites(s);
        writeRow = find(mappedSites==PN)+headingRow;
        if(isempty(writeRow))
            continue
        end
        for c = 1:length(condsToAnalyze)
            excelFileNameWrite{c} = [savePath, 'FRs\M1_' ,PSTH_type,...
                '_',condsToAnalyze{c} ,'_Summary.xlsx'];
            if(s==1)
                try
                    copyfile(excelFileNameRead, excelFileNameWrite{c});
                    writetable(table(headingNames{:}), excelFileNameWrite{c}, 'Range',...
                        [startCol,num2str(headingRow)],'WriteVariableNames', false, 'Sheet', 1);
                catch ME
                    disp('ERROR INITALIZING');
                end
            end
            tableVars = num2cell(allVars{s}(c,:));
            writeT = table(tableVars{:});
            try
                writetable(writeT,excelFileNameWrite{c},'Range',...
                    [startCol,num2str(writeRow),':',...
                    endCol,num2str(writeRow)],...
                    'WriteVariableNames', false, 'Sheet', 1);
            catch ME
                disp('ERROR WRITING TO FILE');
                continue;
            end
        end
        dispstat(['Writing session info... ', num2str(100*round(s/length(dates))), '%']);
    end
end

function bin = findBins(allBins, timePoint)
binSize = mode(diff(allBins));
bin = find(isalmost(allBins,timePoint,binSize/1.99),1);
if(isempty(bin))
    bin = NaN;
end
end

function updatedispstat(a,b)
persistent count total;
if(a==0)
    count = 1;
    total = b;
else
    count = count + 1;
end
dispstat(['Processing sessions... ',num2str(100*round(count/total,2)),'%']);
end