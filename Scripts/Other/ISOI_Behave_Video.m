clear all;
close all;

monkey = 'Gilligan';
videoCond = 'ExtraSmallSphere';
% outline intensity values of tTests
outlineType = 'none';
savePlots = false;
% show intensity frames
intensityFrames = true;
fillTtest = false;
clipVals = [-.2 .2]; % values that intensity images are clipped to
smoothISOIKernel = 5; % # of frames
sessionOverlaps = 1:6;
voronoiStripeSize = 6;
frameVideoGap = 100;
BVFontSize = 50;
blankFrames = 15;
monkeyMapPath = ['S:\Lab\ngc14\Working\', monkey,'\'];
vidPath = [monkeyMapPath,'Videos\'];

if(strcmp(monkey, 'Gilligan'))
    pixelToCM = 38;
elseif(strcmp(monkey,'Skipper'))
    pixelToCM = 61;
end

ISOIFR = 10;
behaveFR = 120;
outputFR = 60;
outlineWidth =  4;
LPk = 10;
cropRect = [1, 960, 195, 1308];

numThresholds = length(sessionOverlaps);
thresholdsR = linspace(clipVals(1)/numThresholds,clipVals(1), numThresholds);
thresholdsR = -0.08;
if strcmp(outlineType, 'intensity')
     numThresholds = length(thresholdsR);
end
cm = flipud(jet);
outlineCM = flipud(gray);
outlineCM = interp1([outlineCM(length(outlineCM)/4,:); outlineCM(3*length(outlineCM)/4,:)],...
    (1:1/length(outlineCM):2)');
outlineCM = outlineCM(fix(linspace(1,length(outlineCM),numThresholds)),:);

MM = double(imread(['S:\Lab\ngc14\Figures\Imaging\MM_Simp_',monkey,'-01.png']));
ISOIFramesPath = ['S:\Lab\ngc14\Working\',monkey,'\[',videoCond,']\','TrialAvg.mat'];
tTestFramesPath = [monkeyMapPath,'[',videoCond,']\tTests\'];
behavioralVidPath = ['S:\Lab\ngc14\Figures\Behavior\Behavior_Videos\',...
    monkey,'_',videoCond,'.mp4'];
maskPath = ['S:\Lab\ngc14\Figures\Imaging\M1_Mask_Border_',monkey(1),'-01.png'];
mask = im2double(imread(maskPath));
maskMask = im2double(imread(['S:\Lab\',monkey,'\Mapping\','clean_mask_filled.bmp']))<.70;
maskMask = maskMask | rgb2gray(mask)<.99;
rr = regionprops(rgb2gray(mask)>=.99,'FilledArea','PixelIdxList');
rr = rr(ismember([rr.FilledArea],700:750));
maskMask(vertcat(rr.PixelIdxList)) = true;
maskMask = repmat(maskMask, [1 1 3]);

if(~(exist(vidPath,'dir')))
    mkdir(vidPath);
end
outputVidPath = [vidPath, 'Behave_ISOI_'];

if(strcmp(outlineType, 'intensity'))
    outputVidPath = [outputVidPath,'IntensityOutline_'];
    stringThresholds = arrayfun(@num2str,thresholdsR,'UniformOutput', false);
    outputVidPath = strjoin([outputVidPath,cellfun(@(s) s, stringThresholds,'UniformOutput', false),'_'],'');
elseif(strcmp(outlineType, 'tTest') && fillTtest)
    outputVidPath = [outputVidPath,'tTestFill_'];
    stringSessions = arrayfun(@num2str,sessionOverlaps,'UniformOutput', false);
    outputVidPath = [outputVidPath,stringSessions{1},'to',stringSessions{end},'_'];
elseif(strcmp(outlineType, 'tTest') && ~fillTtest)
    outputVidPath = [outputVidPath,'tTestOutline_'];
    stringSessions = arrayfun(@num2str,sessionOverlaps,'UniformOutput', false);
    outputVidPath = [outputVidPath,stringSessions{1},'to',stringSessions{end},'_'];
end
outputVidPath = [outputVidPath, monkey, '_', videoCond];

BV = VideoReader(behavioralVidPath);
%% get M1 forelimb to quantify overlap
% M1/PMd border is yellow
M1border = bwareafilt(MM(:,:,1)>125 & MM(:,:,2)>125 & MM(:,:,3)<125, 1);
% PMd/PMv border is purple
PMborder = bwareafilt(MM(:,:,1)>125 & MM(:,:,2)<125 & MM(:,:,3)>125, 1);
% M1 mask
[row,col] = find(M1border);
pointsX = unique(row);
pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
PMpointsX = pointsX;
PMpointsY = pointsY;
pointsX(end+1:end+2) = [size(MM,2), 1];
pointsY(end+1:end+2) = [1 1];
M1mask = double(poly2mask(pointsY,pointsX, size(MM,1),size(MM,2)));
% PMd mask
PMpointsX(end+1:end+2) = [size(MM,2), 1];
PMpointsY(end+1:end+2) = [size(MM,2) size(MM,2)];
PMmask = double(poly2mask(PMpointsY,PMpointsX, size(MM,1),size(MM,2)));
[row,col] = find(PMborder);
row(end+1:end+2) = [1, 1];
col(end+1:end+2) = [size(MM,2), 1];
PMdmask = double(poly2mask(col,row,size(MM,1),size(MM,2)));
PMdmask = double(PMmask & PMdmask);
% arm is red
armArea = MM(:,:,1)>125 & MM(:,:,2)<125 & MM(:,:,3)<125;
% hand is green
handArea = MM(:,:,2)>125 & MM (:,:,1)<125 & MM(:,:,3)<125;
LPkMM = voronoiStripeSize+1;
forelimbMM = imfilter(armArea | handArea, fspecial('average', LPkMM-1));
forelimbMM = imfilter(imdilate(forelimbMM,strel('line',LPkMM,90)),ones(LPkMM));
M1mask = logical(M1mask.*forelimbMM);
PMdmask = logical(PMdmask.*forelimbMM);

if(strcmp(monkey, 'Skipper'))
    mask  = mask(:,2:end-1,:);
    maskMask = maskMask(:,2:end-1,:);
    M1mask  = M1mask(:,2:end-1);
    PMdmask  = PMdmask(:,2:end-1);
end

% store total forelimb representation size
M1mask = M1mask.*sum(maskMask,3)==3;
PMdmask = PMdmask.*sum(maskMask,3)==3;
M1total = sum(M1mask(:));
PMdtotal = sum(PMdmask(:));
[M1masks{1:numThresholds}] = deal(M1mask);
[PMdmasks{1:numThresholds}] = deal(PMdmask);
widthRatio = size(mask,2)/size(mask,1);

if(intensityFrames)
    ISOIDir = load(ISOIFramesPath);
    ISOIDir = ISOIDir.frames;
    numISOIFrames = size(ISOIDir,3);
    IMOs = movmean(cat(3,ISOIDir),smoothISOIKernel,3,'omitnan');
end

if strcmp(outlineType,'tTest')
    tTestDir = dir([tTestFramesPath, videoCond,'*.mat']);
    [~,sortInd] = natsort({tTestDir.name});
    tTestDir = tTestDir(sortInd);
    numISOIFrames = length(tTestDir);
end
%%
M1size = [];
PMdsize = [];
IMs = zeros([[BV.Height,round(widthRatio*BV.Height),3],numISOIFrames]);
for d = 1:numISOIFrames
    if(intensityFrames)
        currIm = gray2ind(mat2gray(IMOs(:,:,d),clipVals),length(cm));
    else
        currIm = zeros(BV.Height, round(widthRatio*BV.Height));
    end
    %    if(any(any(checkIm <=cutoffInd & checkIm>0)))
    %       smoothedFrame = imfilter(checkIm, fspecial('gaussian', 25, round(25 /2)));
    %       threshMask = smoothedFrame <= cutoffInd & smoothedFrame > 0 & additionalMask==0;% &...
    %       %(mask(:,:,1)==0 & mask(:,:,2)==0 & mask(:,:,3)==0);
    %       %imfilter(threshMask, fspecial('gaussian', 10, round(10/2)))
    %       threshMask = imfilter(threshMask,fspecial('disk', 5));
    %       maskSmooth = zeros(size(newFrames(:,:,:,d)));
    %       maskSmooth(repmat(contourMask, [1,1,3])) = 1;
    %       contourIdx = contourIdx & maskMask;
    %       contourIdx = (imerode(tTest,ones(8,8)) ~= imdilate(tTest, ones(4,4)));
    %       contourIdx = imdilate(bwperim(imfilter(tTest, fspecial('gaussian', 4, round(4/2)))),ones(4,4));
    %         smallDomains = arrayfun(@(a) a.Area<(size(IMs,1)*size(IMs,2)*.0005), areaMasks);
    %         smallDomains = zeros(1,size(areaMasks,1));
    %         if(any(~smallDomains))
    %             areaMasks = areaMasks(~smallDomains);
    %             contourMask(vertcat(areaMasks.PixelIdxList))= tTest(vertcat(areaMasks.PixelIdxList));
    %         end
    %     end
    if(strcmp(outlineType,'tTest'))
        tTest = load([tTestDir(d).folder, '\',tTestDir(d).name]);
        tTest = tTest.imT;
    end
    for e = 1:numThresholds
        if(strcmp(outlineType, 'tTest'))
            allPerim = imdilate(bwperim(tTest==sessionOverlaps(e)),ones(outlineWidth,outlineWidth));
            M1masks{e}(imfill(allPerim, 'holes')) = 0;
            PMdmasks{e}(imfill(allPerim, 'holes')) = 0;
        elseif(strcmp(outlineType, 'intensity'))
            allPerim = imdilate(bwperim(imfilter(currIm, fspecial('gaussian', ...
                LPk, round(LPk/2)))<=thresholdsR(e)), ones(outlineWidth,outlineWidth));
            M1masks{e}(imfill(allPerim, 'holes')) = 0;
            PMdmasks{e}(imfill(allPerim, 'holes')) = 0;
        end
        if(fillTtest)
            currIm(imfill(allPerim, 'holes')) = e;
        else
            %currIm(allPerim)= e;%clipVals(end) - (e*(diff(clipVals)/evaluateSessions));
        end
        M1size(e,d) = sum(M1masks{e}(:));
        PMdsize(e,d) = sum(PMdmasks{e}(:));
    end
    %
    %currIm(currIm<e) = round((size(cm,1)).*((currIm(currIm<1)-clipVals(1))./range(clipVals)))+numThresholds;
    currIm = im2double(ind2rgb(currIm,cm));%[outlineCM; cm]);
    currIm(maskMask==1) = mask(maskMask==1);
    if(strcmp(monkey, 'Skipper'))
        currIm = currIm(cropRect(1):cropRect(2),cropRect(3):cropRect(4),:);
    end
    IMs(:,:,:,d) = imresize3(currIm, [BV.Height, round(widthRatio*BV.Height),3]);
end
 
if(savePlots)
    figure('Units', 'normalized', 'Position', [0 0 1 1]);
    subplot(1,2,1);
    plot(1:size(M1size,2), 100.*(1-(M1size./M1total)), 'LineWidth', 2);
    ylim([0 100]);
    xticklabels(xticks./10);
    xlabel('Time (s)');
    title('M1 % forelimb activity')
    ylabel('%');
    if(fillTtest)
        l = legend(arrayfun(@(a) num2str(a), sessionOverlaps, 'UniformOutput', false), 'Orientation', 'vertical', 'Location', 'northeast');
        title(l, {'# of tTests', 'overlapped'});
    elseif(intensityFrames)
        l = legend(arrayfun(@(a) num2str(a,'%1.2f'), thresholdsR, 'UniformOutput', false), 'Orientation', 'vertical', 'Location', 'northeast');
        title(l, {'Reflectance change', 'threshold'});
    end
    set(l, 'Units', 'normalized');
    set(l, 'Position', [l.Position(1)+(1.35*l.Position(3)), l.Position(2:4)])
    subplot(1,2,2);
    plot(1:length(PMdsize), 100.*(1-(PMdsize./PMdtotal)), 'LineWidth', 2);
    ylim([0 100]);
    title('PMd % forelimb activity')
    xticklabels(xticks./10);
    xlabel('Time (s)');
    ylabel('%');
    
    if(fillTtest)
        saveas(gcf,[vidPath, 'tTest_Overlap_',monkey,'_',videoCond],'fig');
        saveas(gcf,[vidPath, 'tTest_Overlap_',monkey,'_',videoCond],'png');
        saveas(gcf,[vidPath, 'tTest_Overlap_',monkey,'_',videoCond,'.eps'],'epsc');
    elseif(intensityFrames)
        saveas(gcf,[vidPath, 'Intensity_Overlap_',monkey,'_',videoCond],'fig');
        saveas(gcf,[vidPath, 'Intensity_Overlap_',monkey,'_',videoCond],'png');
        saveas(gcf,[vidPath, 'Intensity_Overlap_',monkey,'_',videoCond,'.eps'],'epsc');
    end
end
clear IMOs ISOIDir currIm MM PMMask currFrame mask maskMask M1masks PMdmasks allPerim armArea forelimbMM handArea M1border M1mask PMborder PMmask col row;
%%
count = 1;
BVFrames = [];
while(BV.hasFrame)
    currFrame = BV.readFrame;
    while(~any(currFrame(:)))
        currFrame= BV.readFrame;
    end
    if(strfind(BV.VideoFormat, 'RGB'))
        BVFrames(:,:,:,count) = im2double(currFrame);
    else
        BVFrames(:,:,:,count) = repmat(im2double(currFrame),[1 1 3]);
    end
    count = count + 1;
end
%%
numTicks = 5;
rulerStroke = BVFontSize /5;
tickHeight = 3*rulerStroke;
rulerFontSz = round(BVFontSize/1.5);
edgeBuffer = rulerStroke/4;
startRulerX = BVFontSize/numTicks;
boundingBoxHeight = BVFontSize+tickHeight+(2*edgeBuffer);
startRulerY = edgeBuffer+tickHeight/2;

FW = VideoWriter([outputVidPath,'_test'], 'Motion JPEG AVI');
FW.FrameRate = outputFR;
FW.Quality = 90;
open(FW);
%%
for d = 1:size(BVFrames,4)
    newFrame = zeros(BV.Height, BV.Width + size(IMs,2) + frameVideoGap,3);
    newFrame(:,1:BV.Width,:) = BVFrames(:,:,:,d);
    newFrame(:,BV.Width+frameVideoGap+1:end,:) = ...
        IMs(:,:,:,min(size(IMs,4),(ISOIFR*+round(d/behaveFR,1))+1));
    % timer bounding box
    newFrame = insertText(newFrame, [BV.Width,BV.Height],...
        [num2str(round(d/behaveFR,1)-1,'%1.1f'),' s'],'AnchorPoint','RightBottom',...
        'TextColor', 'white', 'BoxColor', 'black','Font', 'Lucida Sans Regular',...
        'FontSize',BVFontSize,'BoxOpacity', 0.6);
    % timer
    newFrame = insertText(newFrame,[1,BV.Height],[num2str(outputFR/behaveFR,...
        '%10.2f'),'x'],'AnchorPoint','LeftBottom','TextColor', 'white',...
        'BoxColor', 'black','Font','Lucida Sans Regular','FontSize',...
        BVFontSize,'BoxOpacity', 0.6);
    % Go cue
    if(round(d/behaveFR,2)>=1 && round(d/behaveFR,2)<1.5)
        newFrame = insertText(newFrame, [BV.Width + frameVideoGap/2,1],'Go cue',...
            'AnchorPoint','CenterTop','TextColor','white','BoxColor','black',...
            'Font','Lucida Sans Regular', 'FontSize',50,'BoxOpacity', 0.6);
    end
    % ruler bounding box
    newFrame = insertShape(newFrame,'FilledRectangle',...
        [1,1,pixelToCM*(numTicks+2),BVFontSize+(3*rulerStroke)+(2*edgeBuffer)],...
        'Color', [.6,.6,.6],'Opacity', 1);
    % ruler
    newFrame = insertShape(newFrame,...
        'FilledRectangle', [startRulerX,startRulerY,pixelToCM*6,rulerStroke],...
        'Color', 'white','Opacity', 1);
    % ruler ticks
    for r = 1:numTicks
        newFrame = insertShape(newFrame,'FilledRectangle',...
            [startRulerX+(r*pixelToCM),startRulerY+(rulerStroke/2)-(tickHeight/2),...
            edgeBuffer,rulerFontSz],'Color', 'white','Opacity', 1);
        % ruler numbers
        newFrame = insertText(newFrame,...
            [startRulerX+(r*pixelToCM)+1,edgeBuffer+edgeBuffer+rulerFontSz],...
            num2str(r),'AnchorPoint','CenterTop','TextColor', 'white','Font',...
            'Lucida Sans Regular','FontSize',rulerFontSz,'BoxOpacity', 0);
        if(r==numTicks)
            newFrame = insertText(newFrame,...
                [startRulerX+(numTicks*pixelToCM)+(rulerFontSz/4),...
                edgeBuffer+edgeBuffer+rulerFontSz],'cm','AnchorPoint',...
                'LeftTop','TextColor', 'white','Font', 'Lucida Sans Regular','FontSize',...
                rulerFontSz,'BoxOpacity', 0);
        end
    end
    pause(.01);
    writeVideo(FW,newFrame);
    pause(0.1)
end
close(FW);