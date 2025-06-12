close all;
clear all;

useCorners = 1;
unitVals = 1;
monkey = 'Gilligan';
frameResolution = .1; %frame rate to capture snap shots of motor map encoding in seconds
baselineBufferTime = 1; % time after trial initiate to start calculating baseline activity in seconds
durationOfTrialInitiate = 5; % duration of time from trial initiate to go cue in seconds(check arduino files)
binSize=.01; % bin size in seconds
baselineWindowSize = 2.4; % size of window for baseline calculation in seconds
secondsBeforePSTH = 6;
secondsAfterPSTH = 3;
startVideo = 4;
endVideo = 2.7;
playBackFR = 3;
videoFileName = ['S:\Lab\', monkey, '\Mapping\Encoding Maps\Videos\SingleUnits_UnitPSTHs.avi'];
labelNames = {'Cue', 'Start Reach', 'Start Grasp', 'Start Lift', 'Start Hold', 'Start Withdrawl', 'Start Replace Hold'};


bins = -secondsBeforePSTH:binSize:secondsAfterPSTH;
PSTHPath = ['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\'];
allFiles = dir([PSTHPath, '*.mat']);
refMask = imread('S:\Lab\ngc14\Figures\Physiology\M1_Mask-01.png');
refMaskInd = repmat(sum(refMask,3)>1,[1, 1, 3]);


encodingColoring = [0 0 0; 0 0 .6; linspace(1,1,4)', linspace(.8,0,4)', linspace(.8,0,4)'];
encodingColoring(end+1,:) = [.5 .5 .5];

[~,sortInd] = natsort({allFiles.name});
allFiles = allFiles(sortInd);
numSites = length(allFiles);
tiledMap = load(['S:\Lab\', monkey, '\Mapping\Motor Maps V3\MM_Tiles.mat']);
v = tiledMap.v;
c = tiledMap.c;
siteMask = tiledMap.siteMask;
siteNum = tiledMap.siteNum;
siteLoc = tiledMap.siteLoc;

numBinsPerFrame = frameResolution/binSize;
startBin = find(isalmost(bins,-startVideo,binSize/1.99),1);
endBin = find(isalmost(bins, endVideo, binSize/1.99),1);
totalFrames = round(length(startBin:endBin)/numBinsPerFrame);

allMasks = {};
for i =1:numSites
    currMask = false(size(refMask,1), size(refMask,2));
    [x,y] = ind2sub([size(refMask,1), size(refMask,2)], find(siteMask(:,:,i)));
    [in,on] = inpolygon(x,y,v(c{i},2),v(c{i},1));
    currMask(sub2ind([size(refMask,1), size(refMask,2)], x(in|on),y(in|on))) = true;
    allMasks{i} = currMask;
end

allFrames = NaN(size(refMask,1), size(refMask,2), totalFrames);
allFrames = squeeze(num2cell(allFrames, [1 2]));

allPenetrationSegs = [];
for i = 1:numSites
    fileInd = find(strcmp({allFiles.name}, [num2str(siteNum(i)), '.mat']));
    loaded = load([allFiles(fileInd).folder, '\', allFiles(fileInd).name]);
    ESSInds = cellfun(@(a) strcmp(a, 'Extra Small Sphere'), loaded.trialInfo(:,1));
    segs = loaded.alignedTrialSegs(ESSInds);
    trialSegs = cellfun(@length, segs);
    correctSegLength = trialSegs==mode(trialSegs);
    allPenetrationSegs(end+1,:) = nanmean(reshape(cell2mat(segs(correctSegLength)'),[sum(correctSegLength),mode(trialSegs)]));
    trialStartBin = cellfun(@(a) find(isalmost(bins,a(1)-...
        (durationOfTrialInitiate - baselineBufferTime),binSize/1.99),1), segs, 'UniformOutput', false);
    sitePSTHS = loaded.allPSTHS;
    
    if(unitVals)
        weights =  loaded.weights;
    end
    
    siteValsOverTime = [];
    for t = 1:totalFrames
        baseFR = [];
        currFRS = [];
        for unit = 1:size(sitePSTHS,3)
            currPSTH = [];
            if(unitVals)
                currPSTH = num2cell(sitePSTHS(ESSInds,:,unit).*weights(unit)',2);
            else
                currPSTH = num2cell(sitePSTHS(ESSInds,:,unit),2);
                %nanmean(nanmean(sitePSTHS(ESSInds,:,:),1),3)
            end
            
            baseFR(:,unit) = cell2mat(cellfun(@(a,b) nanmean(a(b:b+(baselineWindowSize/binSize))), ...
                currPSTH', trialStartBin,'UniformOutput', false));
            currBin = ((t-1)*numBinsPerFrame)+startBin;
            currFRS(:,unit) = cell2mat(cellfun(@(a) nanmean(a(currBin:...
                currBin+numBinsPerFrame)),currPSTH', 'UniformOutput', false));
        end
        
        currFR = NaN;
        if(unitVals)
            currFR = nanmean(nanmean(currFRS,2)./nanmean(baseFR,2));
        else
            currFR = nanmean(nanmean(currFRS,1)./nanmean(baseFR,1));
        end
        
        if(currFR ==0)
            allFrames{t}(allMasks{i})= 1;
        elseif(currFR>0 && currFR<=.67)
            allFrames{t}(allMasks{i}) = 2;
        elseif(currFR>.67 && currFR<=1.3)
            allFrames{t}(allMasks{i}) = 8;
        elseif(currFR>1.3 && currFR<=2)
            allFrames{t}(allMasks{i}) = 3;
        elseif(currFR>2 && currFR<=2.6)
            allFrames{t}(allMasks{i}) = 4;
        elseif(currFR>2.6 && currFR<=3.3)
            allFrames{t}(allMasks{i}) = 5;
        elseif(currFR>3.3)
            allFrames{t}(allMasks{i}) = 6;
        end
    end
end

vid = VideoWriter(videoFileName);
vid.FrameRate = playBackFR;
open(vid);
rgbIm = cellfun(@(a) ind2rgb(a, encodingColoring), allFrames, 'UniformOutput', false);
allIms = [];
for d = 1:length(rgbIm)
    currIm = rgbIm{d};
    currIm(refMaskInd) = im2double(refMask(refMaskInd));
    allIms(:,:,:,d) = currIm;
end
allIms=permute(reshape(cell2mat(rgbIm'), [size(refMask,1), size(refMask,2),totalFrames ,3]),[1,2,4,3]);
allIms(repmat(refMaskInd,[1, 1, 1, totalFrames])) = im2double(repmat(refMask(refMaskInd),[1,1,1,totalFrames]));

labelFrames = arrayfun(@(a) round((find(isalmost(bins, a, binSize/1.99),1)-startBin)/numBinsPerFrame), nanmean(allPenetrationSegs), 'UniformOutput', false);
labelFrames(cellfun(@isempty, labelFrames)) = [];
for l = 1:length(labelFrames)-1
    range = labelFrames{l}:labelFrames{l+1}-1;
    for r = 1:length(range)
        allIms(:,:,:,range(r)) = insertText(allIms(:,:,:,range(r)), [600,25],labelNames{l}, 'FontSize', 30, 'AnchorPoint', 'RightTop');
    end
    if (l==length(labelFrames)-1)
        range = labelFrames{l+1}:labelFrames{l+1}+2;
        for r = 1:length(range)
            allIms(:,:,:,range(r)) = insertText(allIms(:,:,:,range(r)), [600,25],labelNames{l+1}, 'FontSize', 30, 'AnchorPoint', 'RightTop');
        end
    end
end
writeVideo(vid,allIms);
close(vid);