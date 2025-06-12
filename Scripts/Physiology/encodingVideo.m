close all;
clear all;

PSTH_type = 'Trial';
monkey = 'Gilligan';
cond = 'ESS';
singleOrAll = 'All';
frameResolution = .01; %frame rate to capture snap shots of motor map encoding in seconds
baselineBufferTime = 1; % time after trial initiate to start calculating baseline activity in seconds
durationOfTrialInitiate = 5; % duration of time from trial initiate to go cue in seconds(check arduino files)
baselineWindowSize = 2.4; % size of window for baseline calculation in seconds
startVideo = .7;
endVideo = 7;
playBackFR = 6;

sizeMultiplier = 4;
FRChangeClip = 5;
topRightPoint = [690, 120];
timerPosition = [405,650];
timeFontSize = 40;
legendFontSize = 40;
bufferLegendSpace = 10;
numFRLegend = 5;

refMask = imread('S:\Lab\ngc14\Figures\Physiology\M1_Mask-01.png');
MMColorReferenceName = ['S:\Lab\',monkey,'\Mapping\Motor Maps V3\Gilligan MM All 42px_RNR.mat'];
PSTHPath = ['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHs\',singleOrAll,'\'];
sessionInfoPath = ['S:\Lab\',monkey,'\Mapping\Encoding Maps\PSTHS\', singleOrAll,'\FRs\','M1_',cond,'_',PSTH_type,'_Summary.xlsx'];
jointName = {'Hand', 'Arm', 'Face', 'Trunk'};
remap_Joint = [1;2;2;5;5;6;7];
jointVals = [2,5,6,7];
labelNames = {'Cue', 'Reach', 'Grasp', 'Lift', 'Hold', 'Replace'};
excelVariablesToRead = {'site','x','y','total'};
legendVals = linspace(1,FRChangeClip,numFRLegend);

MMColorReference = load(MMColorReferenceName);
cmap_Joint = MMColorReference.cmap_Joint;
MMColorReference = MMColorReference.MM_image_flat;
for r = 1:length(remap_Joint)
    cmap_Joint(r,:) = cmap_Joint(remap_Joint(r),:);
end
cmap_Joint(2:3,:) = [.55 0 0; .55 0 0];
cmap_Joint(4:5,:) = [0 .54 .54; 0 .54 .54];
cmap_Joint(6,:) = [.55 0 .5];
cmap_Joint(end,:) = [.55 .50 .37];
remap_Joint = unique(remap_Joint);
spsd = spreadsheetDatastore(sessionInfoPath,'Sheet',1);
spsd.ReadSize = 'sheet';
sessionInfo = read(spsd);
sessionVarNames = sessionInfo.Properties.VariableNames;
for s = 1:length(excelVariablesToRead)
    [matched, splitT] = regexp(sessionVarNames, excelVariablesToRead{s}, 'ignorecase', 'match','split','forcecelloutput');
    contendors = ~cellfun(@isempty,matched);
    if(any(contendors))
        contendors = find(contendors);
        if(length(contendors)==1)
            varInd(s) = contendors;
        else
            [~, varInd(s)] = min(cellfun(@(a) sum(cellfun(@length,a)), splitT(contendors)));
            varInd(s) = contendors(varInd(s));
        end
    else
        disp('No Matching Variables')
        return;
    end
end
containsUnits = ~isnan(sessionInfo.(sessionVarNames{varInd(4)}));
siteNum = sessionInfo.(sessionVarNames{varInd(1)})(containsUnits);
siteLoc = [sessionInfo.(sessionVarNames{varInd(2)})(containsUnits), ...
    sessionInfo.(sessionVarNames{varInd(3)})(containsUnits)];
%% PSTHs
allFiles = dir(PSTHPath);
allFiles = allFiles(cellfun(@(a) ~isnan(str2double(a)),{allFiles.name}));
[~,sortInd] = natsort({allFiles.name});
allFiles = allFiles(sortInd);
numSites = length(allFiles);
%%
allPenetrationSegs = [];
for i = 1:numSites
    folderInd = find(strcmp({allFiles.name}, num2str(siteNum(i))));
    loaded = load([allFiles(folderInd).folder, '\', allFiles(folderInd).name,...
        '\',allFiles(folderInd).name,'.mat']);
    if(i==1)
        bins = loaded.bins;
        binStrLengths = round(nanmean(arrayfun(@(a)numel(num2str(a)),bins)),0);
        binSize= mode(diff(bins));
        numBinsPerFrame = round(frameResolution/binSize,0);
        startBin = find(isalmost(bins,-startVideo,binSize/1.99),1);
        endBin = find(isalmost(bins, endVideo, binSize/1.99),1)-1;
        totalFrames = round(length(startBin:endBin)/numBinsPerFrame,0);
        allFrames = repmat(refMask, [1 1 1 totalFrames]);
    end
    ESSInds = cellfun(@(a) strcmp(a, 'Extra Small Sphere'), loaded.sessionTrials(:,1));
    condWeightInd = cellfun(@(a) contains(a,'Extra Small Sphere'),loaded.sessionConds);
    segs = loaded.alignedTrialSegs(ESSInds);
    trialSegs = cellfun(@length, segs);
    correctSegLength = trialSegs==mode(trialSegs);
    allPenetrationSegs(i,:) = nanmean(reshape(cell2mat(segs(correctSegLength)'),...
        [sum(correctSegLength),mode(trialSegs)]));
    trialStartBin = cellfun(@(a) find(isalmost(bins,a(1)-...
        (durationOfTrialInitiate - baselineBufferTime),binSize/1.99),1),...
        segs, 'UniformOutput', false);
    avgBaselineStartBin = find(isalmost(bins,nanmean(...
        cellfun(@(a) a(1)-(durationOfTrialInitiate - baselineBufferTime), segs)),binSize/1.99),1);
    sitePSTHS = loaded.allPSTHS;
    if(strcmp(PSTH_type, 'Unit'))
        weights =  loaded.weights;
    else
        weights = []; 
    end
    
    parfor t = 1:totalFrames
%         baseFR = [];
%         currFRS = [];
%         for unit = 1:size(sitePSTHS,3)
%             currPSTH = num2cell(sitePSTHS(ESSInds,:,unit),2);
%             baseFR(:,unit) = cell2mat(cellfun(@(a,b) nanmean(a(b:b+round(baselineWindowSize/binSize,0))), ...
%                 currPSTH', trialStartBin,'UniformOutput', false));
%             currBin = ((t-1)*numBinsPerFrame)+startBin;
%             currFRS(:,unit) = cell2mat(cellfun(@(a) nanmean(a(currBin:...
%                 currBin+numBinsPerFrame)),currPSTH', 'UniformOutput', false));
%         end
%         
%         if(strcmp(PSTH_type, 'Units'))
%             currFR = nanmean(nanmean(currFRS.*weights(:,ESSInds)',2)./nanmean(baseFR,2));
%         else
%             currFR = nanmean(nanmean(currFRS,1)./nanmean(baseFR,1));
%         end
        currFR = NaN;
        if(strcmp(PSTH_type, 'Unit'))
            % calculate average PSTH (depending on unit or trial PSTH)
            avgPSTHS = squeeze(nanmean(weights(:, condWeightInd).*permute(nanmean(sitePSTHS(ESSInds,:,:),1),[3,1,2]),1));
        else
            avgPSTHS = nanmean(squeeze(nanmean(sitePSTHS(ESSInds,:,:),3)),1);
        end
        currBin = ((t-1)*numBinsPerFrame)+startBin;
        currFR = nanmean(avgPSTHS(currBin:currBin+numBinsPerFrame))/nanmean(avgPSTHS...
            (avgBaselineStartBin:avgBaselineStartBin+round(baselineWindowSize/binSize,0)));
        jointInds = squeeze(MMColorReference([max(1,siteLoc(i,2)-15):...
            min(siteLoc(i,2)+15,size(MMColorReference,1))], [max(1,siteLoc(i,1)-15):...
            min(siteLoc(i,1)+15,size(MMColorReference,2))],:));
        jointInds = jointInds(~isnan(jointInds));
        colorInd = mode(jointInds);
        allFrames(:,:,:,t) = insertShape(allFrames(:,:,:,t), 'FilledCircle', ...
            [siteLoc(i,1), siteLoc(i,2), sizeMultiplier*min(FRChangeClip,currFR)],...
            'Color', [255.*cmap_Joint(colorInd,:)]);
        
        currBinString = num2str(bins(currBin), '%.1f');
        %% write legends onto image
        if (i==1)
            % time
            allFrames(:,:,:,t) = insertText(allFrames(:,:,:,t),...
                timerPosition,[currBinString, 's'],...
                'FontSize', timeFontSize, 'AnchorPoint','Center',...
                'BoxOpacity', 0, 'TextColor', 'black');
            
            % account for largest circle radius and buffer to offset
            FRCircleX = topRightPoint(1) -(FRChangeClip*sizeMultiplier);
            % align  vertically along existing legend
            FRCircleY = topRightPoint(2);
            % calculate step size to span existing legend
            FRStepSize = legendFontSize + bufferLegendSpace;
            % initalize with a half step
            FRCircleY = round(FRCircleY + FRStepSize/2,0);
            allFrames(:,:,:,t) = insertText(allFrames(:,:,:,t), ...
                [FRCircleX-bufferLegendSpace-sizeMultiplier*FRChangeClip,...
                FRCircleY-bufferLegendSpace], [char(hex2dec('0394')),' FR'],...
                'AnchorPoint','CenterBottom', 'BoxOpacity', 0, 'TextColor', ...
                'black', 'FontSize',legendFontSize);
            for FRL = 1:numFRLegend
                fStep = (FRL-1)*FRStepSize;
                legendName = num2str(legendVals(FRL),'%.0f');
                if(FRL==numFRLegend)
                    legendName = [legendName, '+'];
                end
                % FRChange circle
                allFrames(:,:,:,t) = insertShape(allFrames(:,:,:,t), 'FilledCircle', ...
                    [FRCircleX, FRCircleY+fStep,sizeMultiplier*legendVals(FRL)],...
                    'Color', 'black');
                
                % FRChange name(value)
                allFrames(:,:,:,t) = insertText(allFrames(:,:,:,t), ...
                    [FRCircleX-bufferLegendSpace-legendFontSize+...
                    (abs(2-numel(legendName))*(legendFontSize/2)),FRCircleY+fStep],...
                    legendName,'FontSize', legendFontSize,'AnchorPoint',...
                    'Center','BoxOpacity', 0,'TextColor', 'black');
            end
            
            jstepSize = round((FRCircleY+((numFRLegend-1)*FRStepSize)-topRightPoint(2)-bufferLegendSpace)/length(jointName));
            jX = topRightPoint(1)-(FRChangeClip*sizeMultiplier)-(legendFontSize*...
                (.5+max(cellfun(@numel,jointName)/2)))-bufferLegendSpace;
            jY = topRightPoint(end)-bufferLegendSpace + jstepSize/2;
            for j = 1:length(jointName)
                % joint name
                allFrames(:,:,:,t) = insertText(allFrames(:,:,:,t), ...
                    [jX,jY+((j-1)*jstepSize)],jointName{j},'FontSize', legendFontSize,...
                    'AnchorPoint', 'Center','BoxOpacity', 0,'TextColor',...
                    [255.*cmap_Joint(jointVals(j),:)]);
            end
        end
    end
end
labelFrames = arrayfun(@(a) floor((find(isalmost(bins, a, binSize/1.99),1)-startBin)/numBinsPerFrame),...
    nanmean(allPenetrationSegs), 'UniformOutput', false);
videoFileName = ['S:\Lab\', monkey, '\Mapping\Encoding Maps\Videos\', singleOrAll, '_',PSTH_type, '_PSTHs.avi'];

vid = VideoWriter(videoFileName);
vid.FrameRate = playBackFR;
open(vid);
for l = 1:length(labelNames)
    range = labelFrames{l}:labelFrames{l+1}-1;
    for r = 1:length(range)
        allFrames(:,:,:,range(r)) = insertText(allFrames(:,:,:,range(r)), ...
            [timerPosition(1), timerPosition(2)+timeFontSize],labelNames{l},...
            'FontSize', timeFontSize, 'AnchorPoint', 'CenterTop', 'TextColor',...
            255.*[.05  .05 .05], 'BoxOpacity', 0);
    end
end
writeVideo(vid,allFrames);
close(vid);