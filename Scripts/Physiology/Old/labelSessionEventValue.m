clear all;
close all;
monkey = 'Gilligan';
dateStart = datetime(2019,4,1, 'Format', 'MM_dd_y');
dateEnd = datetime(2019,12,6, 'Format', 'MM_dd_y');
dateArray = [dateStart:dateEnd]';
%dateArray = [datetime(2020,7,15, 'Format', 'MM_dd_y'),datetime(2020,7,10, 'Format', 'MM_dd_y'),datetime(2020,5,29, 'Format', 'MM_dd_y')];
singleUnitsOnly = 1;
unitPSTHVals = 1;

savePSTHS = 0;
saveVals = 1;

removeBadSegTrials = 0;

graspConds = {'Extra Small Sphere', 'Large Sphere'};
reachConds = {'Extra Small Sphere', 'Large Sphere','Photocell'};

binSize=.01; % bin size in seconds
sigma = .25; % smoothing window in seconds
sigma = round(sigma/binSize);
secondsBeforePSTHAlignmentPoint = 6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 4; % time (in seconds) after alignement point to create PSTH
PSTHDisplayLimits = [-1, 2.5]; % seconds before and afer reach onset to display the saved PSTHs
phaseWindowSize = .15; % size of window for phases in seconds
baselineWindowSize = 2.4; % size of window for baseline calculation in seconds
baselineBufferTime = 1; % time after trial initiate to start calculating baseline activity in seconds
durationOfTrialInitiate = 5; % duration of time from trial initiate to go cue in seconds(check arduino files)

pVal = 0.05; % p-value to test significance for units FR modulation in windows
correctedP = pVal/2;

graspConds = {'Extra Small Sphere'};
reachConds = {'Extra Small Sphere'};
correctedP = pVal;

excelFileNameRead = ['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx'];
[~,~,rawRef]=xlsread(excelFileNameRead);
startCell = char(size(rawRef,2) + 64 + 1);
headingRow = 1;
headingNames = {'Total single & multi units','Total single units',...
    'Task units', 'Go units', 'Reach units', 'Grasp units',	...
    'Task FR','Go FR','Reach FR','Grasp FR',...
    'Task W FR','Task Go FR','Task Reach FR','Task Grasp FR',};
rawRef = cell2mat(rawRef(headingRow+1:end,strcmp(rawRef(headingRow,:),'Site #')));

if(saveVals)
    excelFileNameWrite = ['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites - ESS'];
    if(unitPSTHVals)
        excelFileNameWrite = [excelFileNameWrite, ' - Unit.xlsx'];
    else
        excelFileNameWrite = [excelFileNameWrite, ' - Trial.xlsx'];
    end
    copyfile(excelFileNameRead, excelFileNameWrite);
    writetable(table(headingNames), excelFileNameWrite, 'Range',...
        [startCell,num2str(headingRow)],'WriteVariableNames', false, 'Sheet', 1);
end

%%
bins = -secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values

% Kernel = (-sigma/2:sigma/2);
% Factor = (1.00)/(std(Kernel)*sqrt(2*pi));
% Kernel = Factor*(exp(-(0.5*((Kernel.^2)/std(Kernel).^2))));


% armHandMask = imread(['S:\Lab\',monkey,'\Mapping\ArmHandMask.png']);
% armMask = armHandMask(:,:,1)>100 & armHandMask(:,:,1)<255;
% handMask = armHandMask(:,:,3)>100 & armHandMask(:,:,3)<255;
% ccArm = bwconncomp(armMask);
% ccHand = bwconncomp(handMask);
% sizeCCArm = cellfun(@length,ccArm.PixelIdxList);
% sizeCCHand = cellfun(@length,ccHand.PixelIdxList);
% [~,maxArm] = max(sizeCCArm);
% [~,maxHand] = max(sizeCCHand);
% armMask = zeros(size(armMask));
% handMask = zeros(size(handMask));
% armMask(ccArm.PixelIdxList{maxArm}) = 1;
% handMask(ccHand.PixelIdxList{maxHand}) = 1;

% MM = double(imread(['S:\Lab\',monkey,'\Mapping\Simplified_MM-01.png']));
% %border is yellow
% M1border = MM(:,:,1)>200 & MM(:,:,2) > 200 & MM(:,:,3) < 100;
% [row,col] = find(M1border);
% pointsX = unique(row);
% pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
% pointsX(end+1:end+2) = [size(MM,2), 1];
% pointsY(end+1:end+2) = [1 1];
% M1mask = double(poly2mask(pointsY,pointsX, size(MM,1),size(MM,2)));
% M1mask(M1mask==0) = NaN;
% M1mask = repmat(M1mask,[1,1,3]);
% MM = MM.*M1mask;
%
% %arm is light grey
% armMask = MM(:,:,1) > 200 & MM(:,:,2) > 200 & MM(:,:,3) > 200 & MM(:,:,1) < 250 & MM(:,:,2) < 250 & MM(:,:,3) < 250;
% armMask = imfilter(conv2(armMask, ones(5)/5^2, 'same')>.5,ones(8));
% %hand is dark grey
% handMask = MM(:,:,1) > 100 & MM (:,:,2) > 100 & MM(:,:,3) > 100 & MM(:,:,1) < 200 & MM(:,:,2) < 200 & MM(:,:,3) < 200;
% handMask = imfilter(conv2(handMask, ones(5)/5^2, 'same')>.5,ones(8));
%
% activationMask = imresize(imread(['S:\Lab\',monkey,'\Mapping\tTests\ESS_HSV.png']),[size(MM,1),size(MM,2)]);
% activationMask = activationMask(:,:,1)>200 & activationMask(:,:,2) < 100 & activationMask(:,:,3) < 100;
% ccAct = bwconncomp(activationMask);
% sizeCCAct = cellfun(@length,ccAct.PixelIdxList);
% smallDomains = sizeCCAct<200;
% pixels = cell2mat(ccAct.PixelIdxList(smallDomains)');
% activationMask(pixels) = 0;
% domains = regionprops(activationMask,'PixelIdxList');
% domainIm = zeros(size(MM,1),size(MM,2));
% for i = 1:length(domains)
%     domainIm(domains(i).PixelIdxList) = i;
% end
[~,~,dateToPNValue] = xlsread(['S:\Lab\', monkey,'\', monkey, 'Session.xlsx']);
[~, dateInd] = find(strcmp(dateToPNValue, 'Date'));
[~,domainInd] = find(strcmp(dateToPNValue, 'Domain'));
[~,PNInd] = find(strcmp(dateToPNValue, 'Site #'));
sessionDates = dateToPNValue(:,dateInd);
validInds = find(cellfun(@(a) any(isstrprop(a,'digit')), sessionDates));
inRangeDates = cellfun(@(a) datetime(a,'Format', 'MM_dd_y')>= dateStart & a <=dateEnd, sessionDates(validInds));
sessionDates = sessionDates(validInds(inRangeDates));
sessionDomains = dateToPNValue(validInds(inRangeDates),domainInd);
sessionSites = dateToPNValue(validInds(inRangeDates),PNInd);
clear dateToPNValue;

sessionDates = cellfun(@(a) datetime(a, 'Format', 'MM_dd_y'), sessionDates, 'UniformOutput', false);
%%
cd(['S:\Lab\',monkey,'\All Data\']);
for d = 1:length(dateArray)
    currName = [monkey,'_', datestr(dateArray(d), 'mm_dd_yyyy')];
    
    if(exist([currName, '\Physiology\Results\',currName,'_1.mat'], 'file'))
        disp(dateArray(d))
        sessionInd = find(cellfun(@(a) isequal(dateArray(d), a), sessionDates));
        PN = sessionSites{sessionInd};
        writeRow = find(rawRef==PN)+headingRow;
        if(strcmp(sessionDomains(sessionInd), 'M1') && ~isempty(writeRow))
            sessionDir = dir(['S:\Lab\Gilligan\All Data\', monkey,'_',datestr(dateArray(d), 'mm_dd_yyyy'),'\Physiology\Results\*.mat']);
            sessionDir = sessionDir(~[sessionDir.isdir]);
            [~,sortInd] = natsort({sessionDir.name});
            sessionDir = sessionDir(sortInd);
        else
            sessionDir = [];
        end
        %%
        trialInfo = [];
        segTimes = {};
        totalUnits = 0;
        singleUnits = 0;
        unitPSTHS = {};
        allUnitInd= 1;
        badSegs = [];
        for f = 1:length(sessionDir)
            load([sessionDir(f).folder, '\', sessionDir(f).name]);
            labeledData = isfield(sortedSpikeData, 'label') || any(contains(who('-file',[sessionDir(f).folder, '\', sessionDir(f).name]), 'label'));
            %% unit criteria checks
            if(labeledData)
                if(~isfield(sortedSpikeData, 'label'))
                    labels = label;
                else
                    labels = sortedSpikeData.label;
                end
                %exclude trials with neuron firing less than 1Hz or greater
                %than 250 Hz
                allGoodTrials = ~(cellfun(@(a,b) length(a) <(b(end)-b(1)) | length(a)>250*(b(end)-b(1)) | any(isnan(a)), sortedSpikeData.SpikeTimes, sortedSpikeData.SegTimes));
                missGraspInds = cellfun(@length, allGoodTrials)+1==cellfun(@(a)...
                    length(sortedSpikeData.ConditionSegments{strcmp(a,...
                    sortedSpikeData.Conditions)}), sortedSpikeData.ArduinoData(:,1)');
                allGoodTrials(missGraspInds) = cellfun(@(a,b) [b(1:2), NaN, b(...
                    find(cellfun(@(g) strcmp(g, 'StartGrasp'), ...
                    sortedSpikeData.ConditionSegments{strcmp(a,sessionConds)})):end)],...
                    sortedSpikeData.ArduinoData(missGraspInds,1)', allGoodTrials(missGraspInds), 'UniformOutput', false);
                goodUnitsOnChannel = sum(allGoodTrials,2)>length(allGoodTrials)/(length(sortedSpikeData.Conditions));
                if(singleUnitsOnly)
                    goodUnits = find(arrayfun(@(a) strcmp(a,'s'), labels) & goodUnitsOnChannel');
                else
                    goodUnits = find(goodUnitsOnChannel');
                end
                totalUnits = totalUnits + sum(goodUnitsOnChannel);
                singleUnits = singleUnits + length(goodUnits);
                for u = 1:length(goodUnits)                                       
                    trialInfo = sortedSpikeData.ArduinoData;
                    %remove trials with neuron firing less than 1Hz or greater
                    %than 250 Hz AND units with bad trial segmentation
                    goodTrialsInSession = allGoodTrials(goodUnits(u),:);
                    nPhases = cellfun(@length, sortedSpikeData.SegTimes(goodUnits(u),:));
                     missGraspInds = nPhases+1==cellfun(@(a)...
                        length(sortedSpikeData.ConditionSegments{strcmp(a,...
                        sortedSpikeData.Conditions)}), sortedSpikeData.ArduinoData(:,1)');
                    nPhases(missGraspInds) = nPhases(missGraspInds)+1;
                    badTrialSegs = ~ismember(nPhases,[1,cellfun(@length,sortedSpikeData.ConditionSegments)]);
                    badSegs(allUnitInd,:) = badTrialSegs;
                    if(removeBadSegTrials)
                        % bad trial segmentation (grasp phase missing)
                        % potentially throwing out data?
                        goodTrialsInSession = goodTrialsInSession & ~badTrialSegs;                        
                    end
                    %% unit PSTH generation
                    alignTimesUnit = cellfun(@(a) getAlignedTimes(a,2), sortedSpikeData.SegTimes(goodUnits(u),:), 'UniformOutput', false);
                    unitRestInds = cellfun(@(a) contains(a, 'Rest'), trialInfo(:,1));
                    % rest conditions aligned to cue (there is no reach
                    % onset for the rest condition)
                    alignTimesUnit(unitRestInds) =  cellfun(@(a) getAlignedTimes(a,1),...
                        sortedSpikeData.SegTimes(goodUnits(u),unitRestInds), 'UniformOutput', false);
                    alignedSpikes = cellfun(@minus, sortedSpikeData.SpikeTimes(goodUnits(u),:), alignTimesUnit, 'UniformOutput', false);
                    trialHists = cellfun(@(a) histcounts(a,bins)./binSize, alignedSpikes,'UniformOutput', false);
                    %                   trialHistsSmooth = cellfun(@(a) conv(a,Kernel, 'same'), trialHists, 'UniformOutput', false);
                    %                   trialHistsSmooth = cellfun(@(a) smoothdata(a,'gaussian', 25), trialHists, 'UniformOutput', false);
                    totalHists = cellfun(@(a) conv(a,gausswin(sigma)/sum(gausswin(sigma)),'full'), trialHists, 'UniformOutput', false);
                    totalHists = cellfun(@(a) a(1:end-(sigma-2)), totalHists, 'UniformOutput', false);
                    %                   totalHists = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
                    unitPSTHS{allUnitInd}(goodTrialsInSession) = totalHists(goodTrialsInSession);
                    unitPSTHS{allUnitInd}(~goodTrialsInSession) = {NaN(1,length(bins))};
                    segTimes{allUnitInd}(goodTrialsInSession) = sortedSpikeData.SegTimes(goodUnits(u),goodTrialsInSession);
                    currTrialSegs = sortedSpikeData.SegTimes(goodUnits(u),:);
                    if(allUnitInd==1)
                        allTrialSegs = currTrialSegs;
                    end
                    replaceInds = cellfun(@(a,b) any(isnan(a)) & ~any(isnan(b)), allTrialSegs,currTrialSegs);
                    allTrialSegs(replaceInds) = currTrialSegs(replaceInds);
                    unitLocation(allUnitInd,1) = f;
                    unitLocation(allUnitInd,2) = goodUnits(u);
                    allUnitInd = allUnitInd + 1;
                end
            end
        end
        
        %%PSTH generation by trials
        if(sum(any(~cellfun(@isempty,unitPSTHS))))                       
            % get rid of all bad segment trials (not every unit accounted
            % for each bad trial)
            if(removeBadSegTrials)
                allBadSegs = any(badSegs,1);
                badUnit = find(sum(badSegs,2)~=sum(allBadSegs));
                for b =1:length(badUnit)
                    wrongTrials = find(allBadSegs ~= badSegs(badUnit(b),:));
                    wrongTrials = wrongTrials - sum(badSegs(badUnit(b),:));
                    unitPSTHS{badUnit(b)}(wrongTrials) = {NaN(1,length(bins))};
                end
            end
            % get trial segment times for all trials (to save)
            alignedTrialSegs = allTrialSegs;
            restInds = cellfun(@(a) contains(a, 'Rest'), trialInfo(:,1));
            alignedTrialSegs(~restInds) =  cellfun(@(a) a-a(2), alignedTrialSegs(~restInds), 'UniformOutput', false);
            alignedTrialSegs(restInds) = cellfun(@(a) a-a(1), alignedTrialSegs(restInds), 'UniformOutput' ,false);
            
            %% getting corresponding PSTH bins for grasp, reach, and task windows
            allPSTHS=cellfun(@(a) cell2mat(a(:)), unitPSTHS, 'UniformOutput', false);
            allPSTHS=permute(reshape(cell2mat(allPSTHS')', [size(allPSTHS{1},2), size(allPSTHS{1},1),length(allPSTHS)]),[2,1,3]);
            % movement conditions aligned to reach onset
            alignTimes = cellfun(@(a) getAlignedTimes(a,2), allTrialSegs, 'UniformOutput', false);
            % rest condition aligned to cue onset
            alignTimes(restInds) = cellfun(@(a) getAlignedTimes(a,1), allTrialSegs(restInds), 'UniformOutput', false);
            
            % PSTH bin for the go cue onset
            goTimes = cellfun(@(a) getAlignedTimes(a,1), allTrialSegs, 'UniformOutput', false);
            alignedGoTimes = cellfun(@minus, goTimes, alignTimes, 'UniformOutput', false);
            goTimesBin = cellfun(@(a) find(isalmost(bins,a,binSize/1.99),1), alignedGoTimes, 'UniformOutput', false);
            % PSTH bin for the start of baseline calcuation
            baselineStartTimes = cellfun(@(a) a-(durationOfTrialInitiate - baselineBufferTime), alignedGoTimes,'UniformOutput', false);
            baselineStartBin = cellfun(@(a) find(isalmost(bins,a,binSize/1.99),1), baselineStartTimes, 'UniformOutput', false);
            % PSTH bin for the start of the hold
            holdTimes = cellfun(@(a) getAlignedTimes(a,length(a)-4), allTrialSegs(~restInds),'UniformOutput', false);
            alignedHoldTimes = cellfun(@minus, holdTimes, alignTimes(~restInds), 'UniformOutput', false);
            holdBin = cellfun(@(a) find(isalmost(bins,a,binSize/1.99),1), alignedHoldTimes, 'UniformOutput', false);
            
            % PSTH bin for the start and end of the reach
            reachTimes = cellfun(@(a) getAlignedTimes(a,2), allTrialSegs(~restInds), 'UniformOutput', false);
            alignedReachTimes = cellfun(@minus, reachTimes, alignTimes(~restInds), 'UniformOutput', false);
            % reach phase: window size centered on reach onset
            reachStartBin = cellfun(@(a) find(isalmost(bins,a-(phaseWindowSize/2),binSize/1.99),1), alignedReachTimes, 'UniformOutput', false);
            reachEndBin = cellfun(@(a) find(isalmost(bins,a+(phaseWindowSize/2),binSize/1.99),1), alignedReachTimes, 'UniformOutput', false);
            
            % PSTH bin for the start and end of the grasp
            liftTimes = cellfun(@(a) getAlignedTimes(a,4), allTrialSegs(~restInds), 'UniformOutput', false);
            alignedLiftTimes = cellfun(@minus, liftTimes, alignTimes(~restInds), 'UniformOutput', false);
            % grasp phase: window size backwards from start of lift
            graspStartBin = cellfun(@(a) find(isalmost(bins, a-phaseWindowSize,binSize/1.99),1), alignedLiftTimes, 'UniformOutput', false);
            graspEndBin = cellfun(@(a) find(isalmost(bins,a,binSize/1.99),1), alignedLiftTimes, 'UniformOutput', false);
            %% get unitPSTHs or trial PSTHs to save calculated values
            taskUnits = zeros(1,size(allPSTHS,3));
            reachUnits = zeros(1,size(allPSTHS,3));
            graspUnits = zeros(1,size(allPSTHS,3));
            goUnits = zeros(1,size(allPSTHS,3));
            taskFR = [];
            maxTaskFR = [];
            goFR = [];
            reachFR = [];
            graspFR = [];
            baseAllFR = [];
            maxBaseAllFR = [];
            baseFR = [];
            wFR = [];
            wBaseFR = [];
            maxBaseFR= [];
            maxWFR = [];
            wMaxBaseFR = [];
            
            for unit = 1:size(allPSTHS,3)              
                currPSTH = num2cell(allPSTHS(:,:,unit),2);
                % task phase: beginning of reach phase to start of lift
                taskFR(:,unit) = cell2mat(cellfun(@(a,b,c) nanmean(a(b:c)), currPSTH(~restInds)', reachStartBin, graspEndBin, 'UniformOutput', false));
                maxTaskFR(:,unit) = cell2mat(cellfun(@(a,b,c) max(a(b:c)), currPSTH(~restInds)', reachStartBin, graspEndBin, 'UniformOutput', false));
                reachFR(:,unit) = cell2mat(cellfun(@(a,b,c) nanmean(a(b:c)), currPSTH(~restInds)', reachStartBin, reachEndBin,'UniformOutput', false));
                graspFR(:,unit) = cell2mat(cellfun(@(a,b,c) nanmean(a(b:c)), currPSTH(~restInds)', graspStartBin, graspEndBin,'UniformOutput', false));
                % go phase: window size forward from go cue
                goFR(:,unit) = cell2mat(cellfun(@(a,b) nanmean(a(b:b+(phaseWindowSize/binSize))), currPSTH', goTimesBin,'UniformOutput', false));
                % baseline: baseline window size forward from baseline
                % start time (1 second after trial start)
                baseAllFR(:,unit) = cell2mat(cellfun(@(a,b) nanmean(a(b:b+(baselineWindowSize/binSize))), currPSTH', baselineStartBin,'UniformOutput', false));
                maxBaseAllFR(:,unit)= cell2mat(cellfun(@(a,b) max(a(b:b+(baselineWindowSize/binSize))), currPSTH', baselineStartBin,'UniformOutput', false));
                % withhold: baseline window size forwards from go cue for
                % rest trials
                wFR(:,unit) = cell2mat(cellfun(@(a,b) mean(a(b:b+(baselineWindowSize/binSize))), currPSTH(restInds)', goTimesBin(restInds),'UniformOutput', false));
                maxWFR(:,unit) = cell2mat(cellfun(@(a,b) max(a(b:b+(baselineWindowSize/binSize))), currPSTH(restInds)', goTimesBin(restInds),'UniformOutput', false));
                % store rest trials separately from movement trials
                wBaseFR(:,unit) = baseAllFR(restInds,unit);
                wMaxBaseFR(:,unit) = maxBaseAllFR(restInds,unit);
                baseFR(:,unit) = baseAllFR(~restInds,unit);
                maxBaseFR(:,unit) = maxBaseAllFR(~restInds,unit);
                
                %% unit counting and FR calculation
                alreadyFlagged = 0;
                for c = 1:length(reachConds)
                    condInds = cellfun(@(a) contains(a,reachConds{c}),trialInfo(~restInds,1));
                    [~, pT] = ttest2(baseFR(condInds,unit), goFR(condInds,unit));
                    if(pT<correctedP && ~alreadyFlagged)
                        alreadyFlagged = 1;
                        goUnits(unit) = 1;
                        %FRTaskUnits(end+1) = nanmean(taskFR(condInds))./nanmean(baseFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),1) = (baseFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),2) = (taskFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),1) = taskFR(condInds);
                        %FRTaskUnits(end+1-sum(condInds):end,2) = baseFR(condInds);
                    end
                end
                
                alreadyFlagged = 0;
                for c = 1:length(reachConds)
                    condInds = cellfun(@(a) contains(a,reachConds{c}),trialInfo(~restInds,1));
                    [~, pT] = ttest2(baseFR(condInds,unit), taskFR(condInds,unit));
                    if(pT<correctedP && ~alreadyFlagged)
                        alreadyFlagged = 1;
                        taskUnits(unit) = 1;
                        %FRTaskUnits(end+1) = nanmean(taskFR(condInds))./nanmean(baseFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),1) = (baseFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),2) = (taskFR(condInds));
                        %FRTaskUnits(end+1:end+sum(condInds),1) = taskFR(condInds);
                        %FRTaskUnits(end+1-sum(condInds):end,2) = baseFR(condInds);
                    end
                end
                
                alreadyFlagged = 0;
                for c = 1:length(graspConds)
                    condInds = cellfun(@(a) contains(a,graspConds{c}),trialInfo(~restInds,1));
                    %increase = mean(baseGraspSpikes(condInds)) < mean(graspSpikeCount(condInds));
                    [~, pG] = ttest2(baseFR(condInds), graspFR(condInds));
                    if(pG<correctedP && ~alreadyFlagged)
                        alreadyFlagged=1;
                        graspUnits(unit) = 1;
                        %FRGraspUnits(end+1) = nanmean(graspFR(condInds))./nanmean(baseFR(condInds));
                        %[~,pGE] = ttest2(reachSpikeCount(~reachOnlyInds),graspSpikeCount);
                        %if(pGE<correctedP)
                        %    if((increase && mean(graspSpikeCount) > mean(baseGraspSpikes)) || ...
                        %            (~increase && mean(graspSpikeCount) < mean(baseGraspSpikes)))
                        %        graspEx = graspEx + 1;
                        %    end
                        %end
                    end
                end
                
                alreadyFlagged = 0;
                for c = 1:length(reachConds)
                    condInds =cellfun(@(a) contains(a,reachConds{c}),trialInfo(~restInds,1));
                    % increase = mean(baseReachSpikes(condInds)) < mean(reachSpikeCount(condInds));
                    [~, pR] = ttest2(baseFR(condInds), reachFR(condInds));
                    if(pR<correctedP && ~alreadyFlagged)
                        alreadyFlagged=1;
                        reachUnits(unit) = 1;
                        %FRReachUnits(end+1) = nanmean(reachFR(condInds))./nanmean(baseFR(condInds));
                        %[~,pRE] = ttest2(reachSpikeCount(~reachOnlyInds), graspSpikeCount);
                        %if(pRE<correctedP)
                        %    if((increase && mean(reachSpikeCount) > mean(baseReachSpikes)) || ...
                        %             (~increase && mean(reachSpikeCount) < mean(baseReachSpikes)))
                        %        reachEx = reachEx + 1;
                        %    end
                        %end
                    end
                end
            end
        end
        
        %% saving figures and values
        if(sum(any(~cellfun(@isempty,unitPSTHS))))
            %% weight unit PSTHS by how much each trial contributes to summary PSTH (each unit weighted as a ratio of the number of trials present to the total number of trials for each condition)
            weights = [];
            for uN = 1:size(allPSTHS,3)
                for c = 1:length(sortedSpikeData.Conditions)
                    condInds = strcmp(trialInfo(:,1)', sortedSpikeData.Conditions{c});
                    weights(uN,c) = sum(cellfun(@(a) any(~isnan(a)), unitPSTHS{uN}(condInds)));
                end
            end
            % ensure total number of trials per condition is calculated
            % from correct trials 
            for c = 1:length(sortedSpikeData.Conditions)
                weights(:,c) = weights(:,c)./max(max(weights(:,c),[],2));
            end
            
            taskModulatedInds = taskUnits~=0;
            restInds = cellfun(@(a) contains(a, 'Rest'), trialInfo(:,1));
            condInds = cellfun(@(a) strcmp(a, 'Extra Small Sphere'), trialInfo(:,1))';
            condsNoRest = cellfun(@(a) strcmp(a, 'Extra Small Sphere'),trialInfo(~restInds,1))';
            condWeightInd = find(cellfun(@(a) strcmp(a, 'Extra Small Sphere'), sortedSpikeData.Conditions));
            restWeightInd = find(cellfun(@(a) strcmp(a, 'Rest'), sortedSpikeData.Conditions));
            if(removeBadSegTrials)
                condInds = condInds & ~allBadSegs';
                restInds = restInds & ~allBadSegs';
            end
            
            if(savePSTHS)
                figure('units', 'normalize', 'position', [0 0 1 1], 'Visible', 'off');
                plot(bins,nanmean(nanmean(bsxfun(@times,allPSTHS(condInds,:,:),...
                    reshape(weights(:,condWeightInd),1,1,[])),3),1),'r--', 'LineWidth', 2);
                hold on;
                plot(bins, nanmean(nanmean(allPSTHS(condInds,:,:),3),1), 'r-', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(allPSTHS(condInds,:,:),1),3),'r:', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(bsxfun(@times,allPSTHS(condInds,:,taskModulatedInds),...
                    reshape(weights(taskModulatedInds,condWeightInd),1,1,[])),3),1),'m--', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(allPSTHS(condInds,:,taskModulatedInds),1),3),'m:', 'LineWidth', 2);
                
                plot(bins,nanmean(nanmean(bsxfun(@times,allPSTHS(restInds,:,:),...
                    reshape(weights(:,restWeightInd),1,1,[])),3),1),'b--', 'LineWidth', 2);
                hold on;
                plot(bins, nanmean(nanmean(allPSTHS(restInds,:,:),3),1), 'g-', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(allPSTHS(restInds,:,:),1),3),'b:', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(bsxfun(@times,allPSTHS(restInds,:,taskModulatedInds),...
                    reshape(weights(taskModulatedInds,restWeightInd),1,1,[])),3),1),'c--', 'LineWidth', 2);
                plot(bins, nanmean(nanmean(allPSTHS(restInds,:,taskModulatedInds),1),3),'c:', 'LineWidth', 2);
                xlim(PSTHDisplayLimits);
                xlabel('Time in relation to reach onset (s)');
                ylabel('Firing Rate');
                
                %if(sum(goodUnitInds)>1)
                %    plot(cutBins(2:end),nanmean(weighting.*squeeze(nanmean(trialPSTHsA(:,:,goodUnitInds),1))',1),'y--', 'LineWidth', 2);
                %    plot(cutBins(2:end),nanmean(squeeze(nanmean(trialPSTHsA(:,:,goodUnitInds),1))',1),'y', 'LineWidth', 2);
                %else
                %    plot(cutBins(2:end),nanmean(weighting.*(nanmean(trialPSTHsA(:,:,goodUnitInds),1)),1),'y--', 'LineWidth', 2);
                %    plot(cutBins(2:end),nanmean(nanmean(trialPSTHsA(:,:,goodUnitInds),1),1),'y', 'LineWidth', 2);
                %end
                %if(sum(goodUnitInds)>1)
                %    plot(cutBins(2:end),nanmean(weighting.*squeeze(nanmean(trialPSTHsWA(:,:,goodUnitInds),1))',1),'c', 'LineWidth', 2);
                %    plot(cutBins(2:end),nanmean(squeeze(nanmean(trialPSTHsWA(:,:,goodUnitInds),1))',1),'c--', 'LineWidth', 2);
                %else
                %    plot(cutBins(2:end),nanmean(weighting.*squeeze((nanmean(trialPSTHsWA(:,:,goodUnitInds),1))),1),'c', 'LineWidth', 2);
                %    plot(cutBins(2:end),nanmean(nanmean(trialPSTHsWA(:,:,goodUnitInds),1),1),'c--', 'LineWidth', 2);
                %end
                legend('ESS - WUnits','ESS - UWUnits', 'ESS - Trials',...
                    'ESS - TM WUnits', 'ESS - TM Trials', ...
                    'Rest - WUnits','Rest - UWUnits','Rest - Trials',...
                    'Rest - TM WUnits','Rest - TM Trials',...
                    'Location', 'northwest','AutoUpdate','off');               
                segBins =  nanmean(cell2mat(alignedTrialSegs(condInds & ~any(badSegs,1))'),1);
                segBins = segBins(segBins> -secondsBeforePSTHAlignmentPoint & segBins<secondsAfterPSTHAlignmentPoint);
                currYLim = get(gca,'ylim');
                currYLim = [0, currYLim(end)];
                arrayfun(@(a) line([a a], currYLim, 'Color', 'black', 'LineStyle', '--'),segBins);
                ylim(currYLim);
                saveas(gcf,['S:\Lab\', monkey,'\Mapping\Encoding Maps\PSTHs\',num2str(PN),'.png' ]);
                saveas(gcf,['S:\Lab\', monkey,'\Mapping\Encoding Maps\PSTHs\',num2str(PN),'.fig' ]);
                pause(0.1);
                close all;
            end
            
            if(saveVals)
                if(unitPSTHVals)
                    FRTaskUnits = weights(:,condWeightInd).*(nanmean(taskFR(condsNoRest,:),1)./nanmean(baseFR(condsNoRest,:),1))';
                    FRTaskUnitsW = weights(:,restWeightInd).*(nanmean(wFR,1)./nanmean(wBaseFR,1))';
                    FRReachUnits = weights(:,condWeightInd).*(nanmean(reachFR(condsNoRest,:),1)./nanmean(baseFR(condsNoRest,:),1))';
                    FRGraspUnits = weights(:,condWeightInd).*(nanmean(graspFR(condsNoRest,:),1)./nanmean(baseFR(condsNoRest,:),1))';
                    FRGoUnits = weights(:,condWeightInd).*(nanmean(goFR(condsNoRest,:),1)./nanmean(baseFR(condsNoRest,:),1))';
                    
                    FRTReachUnits = weights(taskModulatedInds,condWeightInd).*...
                        (nanmean(reachFR(condsNoRest,taskModulatedInds),1)./nanmean(baseFR(condsNoRest,taskModulatedInds),1))';
                    FRTGraspUnits = weights(taskModulatedInds,condWeightInd).*...
                        (nanmean(graspFR(condsNoRest,taskModulatedInds),1)./nanmean(baseFR(condsNoRest,taskModulatedInds),1))';
                    FRTGoUnits = weights(taskModulatedInds,condWeightInd).*...
                        (nanmean(goFR(condsNoRest,taskModulatedInds),1)./nanmean(baseFR(condsNoRest,taskModulatedInds),1))';                    
                else                   
                    FRTaskUnits = nanmean(taskFR(condsNoRest,:),2)./nanmean(baseFR(condsNoRest,:),2);
                    FRTaskUnitsW = nanmean(wFR,2)./nanmean(wBaseFR,2);
                    FRReachUnits = nanmean(reachFR(condsNoRest,:),2)./nanmean(baseFR(condsNoRest,:),2);
                    FRGraspUnits = nanmean(graspFR(condsNoRest,:),2)./nanmean(baseFR(condsNoRest,:),2);
                    FRGoUnits = nanmean(goFR(condsNoRest,:),2)./nanmean(baseFR(condsNoRest,:),2);
                    
                    FRTReachUnits = nanmean(reachFR(condsNoRest,taskModulatedInds),2)./nanmean(baseFR(condsNoRest,taskModulatedInds),2);
                    FRTGraspUnits = nanmean(graspFR(condsNoRest,taskModulatedInds),2)./nanmean(baseFR(condsNoRest,taskModulatedInds),2);
                    FRTGoUnits = nanmean(goFR(condsNoRest,taskModulatedInds),2)./nanmean(baseFR(condsNoRest,taskModulatedInds),2);
                end
                save(['S:\Lab\', monkey,'\Mapping\Encoding Maps\PSTHs\',num2str(PN),'.mat'],...
                    'allPSTHS', 'alignedTrialSegs','weights','trialInfo','unitLocation',...
                    'taskUnits', 'goUnits', 'reachUnits', 'graspUnits', 'bins');
                status = 0;
                while(status==0)
                    writeT = table(totalUnits, singleUnits, sum(taskUnits),...
                        sum(goUnits), sum(reachUnits),sum(graspUnits),...
                        nanmean(FRTaskUnits), nanmean(FRGoUnits),...
                        nanmean(FRReachUnits), nanmean(FRGraspUnits),...
                        nanmean(FRTaskUnitsW), nanmean(FRTGoUnits),...
                        nanmean(FRTReachUnits), nanmean(FRTGraspUnits));
                    try
                        writetable(writeT,excelFileNameWrite,'Range',...
                            [startCell,num2str(writeRow),':',...
                            char(startCell+size(writeT,2)-1),num2str(writeRow)],...
                            'WriteVariableNames', false, 'Sheet', 1);
                        status = 1;
                         pause(.1);
                    catch ME
                        status = 0;
                    end
                    %
                end
            end
            %             reachVals = [reachVals; FRReachUnits'];
            %             graspVals = [graspVals; FRGraspUnits'];
            %             domainInds(end+1:end+length(FRReachUnits)) = domainIm(yLoc,xLoc);
            %             if(armMask(yLoc,xLoc)==1)
            %                 armHandVal(end+1:end+length(FRReachUnits)) = 1;
            %             elseif(handMask(yLoc,xLoc)==1)
            %                 armHandVal(end+1:end+length(FRReachUnits)) = 2;
            %             else
            %                 armHandVal(end+1:end+length(FRReachUnits)) = 0;
            %             end
        end
    end
end


%%
function maxVal = getMaxDeviation(hists,startInd,endInd)
currHist = hists(startInd:endInd);
[~,maxInd] = max(abs(currHist));
maxVal = currHist(maxInd);
end
function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end
function plotHists(ax, valsArm, valsHand,minVal, maxVal,minCutoff, maxCutoff,titleInfo,legendFlag)
stepVal = abs((minVal-maxVal)/9);
[bArm,~] = hist(valsArm,minVal:stepVal:maxVal);
[bHand,e] = hist(valsHand,minVal:stepVal:maxVal);
b = bar(ax,e, [bArm./sum(bArm);bHand./sum(bHand)]');
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 0 1];
YZ = ax.YLim(2);
line(ax,[minCutoff, minCuto+ff], [0, YZ],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
line(ax,[maxCutoff, maxCutoff], [0, YZ],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
ax.YLim(2) = YZ;
if(exist('legendFlag','var'))
    legend(ax,sprintf('Arm (n=%d)',length(valsArm)),sprintf('Hand (n=%d)',length(valsHand)),'Location', 'best');
    title(ax,titleInfo);
else
    title(ax,[titleInfo, sprintf(', Arm (n=%d), Hand (n=%d)',length(valsArm),length(valsHand))]);
end
end