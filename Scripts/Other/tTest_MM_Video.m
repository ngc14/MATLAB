clear all;
close all;

monkey = 'Gilligan';
videoCond = 'Photocell';
savePlots = false;
fillTtest = true;
blanking = true;
additive = true;
showBehavior = false;
sessionOverlaps = [1];
voronoiStripeSize = 6;
frameVideoGap = 140;
BVFontSize = 45;
tTestFR = 10;
behaveFR = 120;
outlineWidth = 5;
LPk = 12;
numThresholds = length(sessionOverlaps);
if(strcmp(monkey, 'Skipper'))
    cropWidth = 180;
end
vidPath = ['S:\Lab\ngc14\Working\',monkey,'\Videos\'];
if(~exist(vidPath,'dir'))
    mkdir(vidPath)
end
monkeyMapPath = ['S:\Lab\', monkey,'\Mapping\'];
MM = double(imread([monkeyMapPath,'Motor Maps V3\MM_Simp_Thresh-01.png']));
MMSimp = im2double(imread(['S:\Lab\ngc14\Figures\Imaging\MM_Simp_',monkey,'-01.png']));
mask = im2uint8(imread(['S:\Lab\ngc14\Figures\Imaging\M1_Mask_tTest_',monkey(1),'-01.png']));
maskMask = im2double(imread([monkeyMapPath,'clean_mask_filled.bmp']))<.70;
tTestFramesPath = dir(['S:\Lab\ngc14\Working\',monkey,'\[',videoCond,']\tTests\NaN\tTest_*_p-0.0001.png']);
intensityPath = dir(['S:\Lab\ngc14\Working\',monkey,'\[',videoCond,']\tTests\NaN\Intensity_*.png']);
behavioralVidPath = ['S:\Lab\ngc14\CH1_M1_PMd_Activity_Organization\Up_to_Date_Figures\Behavie_Videos\Behavior\',...
    monkey,'_',videoCond,'.mp4'];
%additionalMaskPath = 'S:\Lab\ngc14\Figures\Imaging\Gililgan_Mask_Activity-01.png';
maskMask = maskMask | rgb2gray(mask)<254;
rr = regionprops(rgb2gray(mask)>=254,'FilledArea','PixelIdxList');
rr = rr(ismember([rr.FilledArea],670:700));
maskMask(vertcat(rr.PixelIdxList)) = true;
maskMask = repmat(maskMask, [1 1 3]);
outputVidPath = [vidPath, 'tTest_MM_',videoCond];
if(strcmp(monkey, 'Gilligan'))
    pixelToCM = 38;
elseif(strcmp(monkey,'Skipper'))
    pixelToCM = 61;
end
if(strcmp(monkey,'Gilligan'))
    if(strcmp(videoCond,'ExtraSmallSphere'))
        blankFrames = 15:35;
    elseif(strcmp(videoCond,'LargeSphere'))
        blankFrames = 15:37;
    elseif(strcmp(videoCond,'Photocell'))
        blankFrames = 15:22;
    end
elseif(strcmp(monkey,'Skipper'))
    if(strcmp(videoCond,'ExtraSmallSphere'))
        blankFrames = 15:36;
    elseif(strcmp(videoCond,'LargeSphere'))
        blankFrames =15:34;
    elseif(strcmp(videoCond,'Photocell'))
        blankFrames = 15:24;
    end
end
if(showBehavior)
    outputFR = behaveFR;
else
    outputFR = tTestFR/3;
end
if(strcmp(monkey, 'Skipper'))
    mask  = mask(:,2:end-1,:);
    MM = MM(:,2:end-1,:);
    MMSimp = MMSimp(:,2:end-1,:);
end
if(fillTtest)
    outputVidPath = [outputVidPath,'_Fill'];
end
if(additive)
    outputVidPath = [outputVidPath,'_Aggregate'];
end
if(showBehavior)
    BV = VideoReader(behavioralVidPath);
end
colors(1,:) = [1 0 0];
colors(2,:) = [0 1 0];
colors(3,:) = [0 0 1];
if(numThresholds>3)
    colors = distinguishable_colors(numThresholds, [0 0 0; 1 1 1;]);
end
colors = rgb2hsv([1 0 0]);
colors(:,2) = 1;
colors(:,3) = .7;
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
% vertical line filter to smooth dual sites
forelimbMM = imfilter(imdilate(forelimbMM,strel('line',LPkMM,90)),ones(LPkMM));
M1mask = logical(M1mask.*forelimbMM);
PMdmask = logical(PMdmask.*forelimbMM);

% store total forelimb representation size
M1total = sum(M1mask(:));
PMdtotal = sum(PMdmask(:));
[M1masks{1:numThresholds}] = deal(M1mask);
[PMdmasks{1:numThresholds}] = deal(PMdmask);
%maskMask = imbinarize(mask,.2) | maskMask==0;
widthRatio = size(mask,2)/size(mask,1);
if(strcmp(monkey, 'Skipper'))
    cropRect = [1, size(mask,1)-cropWidth, cropWidth, size(mask,2)];
end
%% tTest sort
tTestNames = cellfun(@(cc) compose('tTest_%d_p-0.0001.png',cc),...
    string({tTestFramesPath.name}),'UniformOutput',false);
tTestDir = tTestFramesPath(cellfun(@length,tTestNames)<22);
[~,sortInd] = natsort({tTestDir.name});
tTestDir = tTestDir(sortInd);
intensityDir = intensityPath(cellfun(@length,{intensityPath.name})<22);
[~,sortInd] = natsort({intensityDir.name});
intensityDir = intensityDir(sortInd);
% colormap assignments
for e = 1:numThresholds
    startColor = colors(e,:);
        startColor(2) = .5;
        startColor(3) = 1;
    colormaps{e} = im2uint8(hsv2rgb(interp1([colors(e,:); startColor], (1:1/(length(tTestDir)):2)')));
    colormaps{e}(2,:) = 0;
end
%
startColor =  colormaps{e}(1,:);
M1size = M1total.*ones(numThresholds,length(tTestDir));
PMdsize = PMdtotal.*ones(numThresholds,length(tTestDir));
allActive = logical(zeros(size(MMSimp)));
IMs = zeros([[size(mask,1),round(widthRatio*size(mask,1)),3],length(tTestDir)]);
changeMask = zeros(size(mask, [1 2]));
for t = 1:length(tTestDir)
    intensityIm = imread([intensityDir(t).folder,'\',intensityDir(t).name]);
    intensityIm = imresize3d(intensityIm(50:584,115:792,:),[],size(mask),'linear');
    tTest = imread([tTestDir(t).folder, '\',tTestDir(t).name]);
%     tTest = tTest.imT;
    redChannel = intensityIm(:, :, 1);
    greenChannel = intensityIm(:, :, 2);
    blueChannel = intensityIm(:, :, 3);
    if(additive && any(changeMask,'all'))
        outlines = bwperim(imfill(changeMask>0,'holes'));
        addMask =imdilate(outlines,ones(outlineWidth)) - ...
               imerode(outlines,ones(round(outlineWidth/2)));
        addMaskR = im2uint8(addMask).*im2uint8(colormaps{e}(2,1));
        addMaskG = im2uint8(addMask).*im2uint8(colormaps{e}(2,2));
        addMaskB = im2uint8(addMask).*im2uint8(colormaps{e}(2,3));
        redChannel(addMask>0) = addMaskR(addMask>0);
        greenChannel(addMask>0) = addMaskG(addMask>0);
        blueChannel(addMask>0) = addMaskB(addMask>0);
    end
    if(~ismember(t,blankFrames) && blanking)
        for e = 1:numThresholds
            allPerim = false(size(tTest));
            allPerim(imfill(imgaussfilt(double(tTest),round(LPk/2),'FilterSize',...
                LPk+(mod(LPk,2)==0))>=sessionOverlaps(e),'holes')) = true;
            filledPerim = imfill(allPerim, 'holes');
            if(fillTtest)
                changeMask(filledPerim) = t;
            else
                changeMask(allPerim) = t;
            end
            redChannel(changeMask==t) = startColor(1);
            greenChannel(changeMask==t) = startColor(2);
            blueChannel(changeMask==t) = startColor(3);
            M1masks{e}(filledPerim) = 0;
            PMdmasks{e}(filledPerim) = 0;
            M1size(e,t) = sum(M1masks{e}(:));
            PMdsize(e,t) = sum(PMdmasks{e}(:));
        end
    end
    finalIm = cat(3, redChannel, greenChannel, blueChannel);
    finalIm(maskMask==1) = mask(maskMask==1);
%     finalIm(allActive) = 1;
%     if(~ismember(t,blankFrames) && blanking)
%         allActive = allActive | (finalIm(:,:,1) == colormaps{e}(t,1) & ...
%             finalIm(:,:,2) == colormaps{e}(t,2) & finalIm(:,:,3) == colormaps{e}(t,3));
%     end
    if(strcmp(monkey, 'Skipper'))
        finalIm = finalIm(cropRect(1):cropRect(2),cropRect(3):cropRect(4),:);
    end
    IMs(:,:,:,t) = im2double(imresize3(finalIm, [size(mask,1), round(widthRatio*size(mask,1)),3]));
    
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
    l = legend(arrayfun(@(a) num2str(a), sessionOverlaps, 'UniformOutput', false), 'Orientation', 'vertical', 'Location', 'northeast');
    title(l, {'# of tTests', 'overlapped'});
    set(l, 'Units', 'normalized');
    set(l, 'Position', [l.Position(1)+(1.35*l.Position(3)), l.Position(2:4)])
    subplot(1,2,2);
    plot(1:length(PMdsize), 100.*(1-(PMdsize./PMdtotal)), 'LineWidth', 2);
    ylim([0 100]);
    title('PMd % forelimb activity')
    xticklabels(xticks./10);
    xlabel('Time (s)');
    ylabel('%');
    
    saveas(gcf,[vidPath, 'tTest_',monkey,'_',videoCond,num2str(blankFrames)],'fig');
    saveas(gcf,[vidPath, 'tTest_',monkey,'_',videoCond,num2str(blankFrames)],'png');
    saveas(gcf,[vidPath, 'tTest_',monkey,'_',videoCond,num2str(blankFrames),'.eps'],'epsc');
end
%%
numTicks = 5;
rulerStroke = BVFontSize /5;
tickHeight = 3*rulerStroke;
rulerFontSz = round(BVFontSize/1.5);
edgeBuffer = 10;
startRulerX = BVFontSize/numTicks;
boundingBoxHeight = BVFontSize+tickHeight+(2*edgeBuffer);
startRulerY = edgeBuffer+tickHeight/2;
movementPoint = 2*edgeBuffer;
if(showBehavior)
    numFrames = round(BV.FrameRate*BV.Duration);
    newFrames = zeros(BV.Height, BV.Width + size(IMs,2) + frameVideoGap, 3,numFrames);
else
    numFrames = length(tTestDir);
    newFrames = zeros(size(IMs,1) + frameVideoGap, size(IMs,2), 3,numFrames);
    startRulerY = startRulerY + frameVideoGap;
end
FW = VideoWriter(outputVidPath,'Motion JPEG AVI');
FW.FrameRate = outputFR;
FW.Quality = 100;
open(FW);

frameInd = 1;
continueFrames = 1;
incrementMove = (size(newFrames,2)-BVFontSize-(2*edgeBuffer))/numFrames;
moveLine = [2*BVFontSize+movementPoint,  movementPoint+BVFontSize*2.2,...
    2*BVFontSize+movementPoint, movementPoint+BVFontSize*2.2];
onLine = [];
offLine = [];
upLine = [];
downLine = [];
while(continueFrames)
    if(showBehavior)
        % video layout for behavior and tTest results simulatenously
        % plotted
        BVcurrFrame = im2double(BV.readFrame);
        while(~any(BVcurrFrame(:)))
            BVcurrFrame = im2double(BV.readFrame);
        end
        if(size(BVcurrFrame,3)==1)
            BVcurrFrame = repmat(BVcurrFrame,[1 1 3]);
        end
        newFrames(:,1:BV.Width,:,frameInd) = BVcurrFrame;
        newFrames(:,BV.Width+frameVideoGap+1:end,:,frameInd) = ...
            IMs(:,:,:,min(size(IMs,4),(tTestFR*+round(frameInd/behaveFR,1))+1));
        
        % ruler bounding box
        newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
            'FilledRectangle', [1,1,pixelToCM*(numTicks+2)+(4*edgeBuffer),boundingBoxHeight],...
            'Color', [.6,.6,.6],'Opacity', 1);
        
        % ruler
        newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
            'FilledRectangle', [startRulerX,startRulerY,pixelToCM*6,rulerStroke],...
            'Color', 'white','Opacity', 1);
        % ruler ticks
        for r = 1:numTicks
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),'FilledRectangle',...
                [startRulerX+(r*pixelToCM),startRulerY+(rulerStroke/2)-(tickHeight/2),...
                edgeBuffer,rulerFontSz],'Color', 'white','Opacity', 1);
            % ruler numbers
            newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
                [startRulerX+(r*pixelToCM)+1,edgeBuffer+edgeBuffer+rulerFontSz],...
                num2str(r),'AnchorPoint','CenterTop','TextColor', 'white','Font',...
                'Lucida Sans Regular','FontSize',rulerFontSz,'BoxOpacity', 0);
            if(r==numTicks)
                newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
                    [startRulerX+(numTicks*pixelToCM)+(rulerFontSz/4),...
                    edgeBuffer+edgeBuffer+rulerFontSz],'cm','AnchorPoint',...
                    'LeftTop','TextColor', 'white','Font', 'Lucida Sans Regular','FontSize',...
                    rulerFontSz,'BoxOpacity', 0);
            end
        end
        
        % Go cue
        if(round(frameInd/behaveFR,1)>=1 && round(frameInd/behaveFR,2)<1.5)
            newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
                [BV.Width + frameVideoGap/2,1],'Go cue','AnchorPoint','CenterTop',...
                'TextColor', 'white', 'BoxColor', 'black',...
                'Font', 'Lucida Sans Regular', 'FontSize',50,'BoxOpacity', 0.6);
        end
        
        % timer
        newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            [BV.Width,BV.Height],[num2str(round(frameInd/behaveFR,1)-1,'%1.1f'),' s'],...
            'AnchorPoint','RightTop','TextColor', 'white', 'BoxColor', 'black',...
            'Font', 'Lucida Sans Regular', 'FontSize',BVFontSize,'BoxOpacity', 0.6);
    else
        newFrames(frameVideoGap+1:end, :,:,frameInd) = IMs(:,:,:,frameInd);
        %
        %         % behavioral speed
        %         newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
        %             [size(IMs,2),0],[num2str(outputFR/tTestFR,'%10.2f'),'x'],...
        %             'AnchorPoint','RightTop','TextColor', 'white', 'BoxColor', 'black',...
        %             'Font', 'Lucida Sans Regular', 'FontSize',BVFontSize,'BoxOpacity', 0);
        % timer
        newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            [1,size(newFrames,1)],[num2str(round(frameInd/tTestFR,1)-1,'%1.1f'),' s'],...
            'AnchorPoint','LeftBottom','TextColor', 'white', 'BoxColor', 'black',...
            'BoxOpacity',0.6,'Font', 'Lucida Sans Regular', 'FontSize',BVFontSize);
        %30, 170 Skipper
        %19, 120 Gilligan
        newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            [size(newFrames,2)-135,size(newFrames,1)],'0.33x',...
            'AnchorPoint','RightBottom','TextColor', 'white', 'BoxColor', 'black',...
            'BoxOpacity',0.6,'Font', 'Lucida Sans Regular', 'FontSize',BVFontSize);
        
        
        newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            [1,edgeBuffer], 'ON', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
            'BoxColor', 'black','BoxOpacity', 0, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*1));
        newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            [1,(movementPoint+edgeBuffer)+BVFontSize],'OFF', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
            'BoxColor', 'black','BoxOpacity', 0, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*1));
        
        %         newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
        %             [movementPoint, -3*edgeBuffer],'Movement', 'AnchorPoint',...
        %             'LeftTop', 'TextColor', 'white','BoxColor', 'black','BoxOpacity',...
        %             0, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*.6));
        % movement indicator
%         if(frameInd>13 && frameInd<blankFrames)
            % movement ON
            if(frameInd>=14 && frameInd<blankFrames(end))
                if(frameInd==14)
                    upLine = [moveLine(3),moveLine(2), moveLine(3), movementPoint+BVFontSize/2];
                    onLine = [moveLine(3), movementPoint+BVFontSize/2, moveLine(3), movementPoint+BVFontSize/2];
                end
                onLine = onLine + [0, 0, incrementMove, 0];
            % movement OFF
            elseif(frameInd>=blankFrames(end))
                if(frameInd==blankFrames(end))
                    downLine = [onLine(3),onLine(2), onLine(3),moveLine(4)];
                    offLine = [onLine(3), moveLine(4), onLine(3), moveLine(4)];
                end
                offLine = offLine + [0, 0, incrementMove, 0];
            else
               moveLine = moveLine + [0, 0, incrementMove, 0];
            end
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
                'Line',moveLine, 'Color', 'blue', 'LineWidth', 5, 'Opacity', 1);
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
                'Line',upLine,'Color', 'blue', 'LineWidth', 5, 'Opacity', 1);
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
                'Line',onLine,'Color', 'blue', 'LineWidth', 5, 'Opacity', 1);
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
                'Line',offLine,'Color', 'blue', 'LineWidth', 5, 'Opacity', 1);
            newFrames(:,:,:,frameInd) = insertShape(newFrames(:,:,:,frameInd),...
                'Line',downLine,'Color', 'blue', 'LineWidth', 5, 'Opacity', 1);
            %             newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            %                 [movementPoint+edgeBuffer,BVFontSize*.6], 'ON', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
            %                 'BoxColor', 'yellow','BoxOpacity', .6, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*.6));
            %             newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
            %                 [movementPoint+BVFontSize*1.5,BVFontSize*.6],'OFF', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
            %                 'BoxColor', 'black','BoxOpacity', 0, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*.6));
%         else

%             end
%             newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
%                 [movementPoint+edgeBuffer,BVFontSize*.6], 'ON', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
%                 'BoxColor', 'black','BoxOpacity', 0, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*.6));
%             newFrames(:,:,:,frameInd) = insertText(newFrames(:,:,:,frameInd),...
%                 [movementPoint+BVFontSize*1.5,BVFontSize*.6], 'OFF', 'AnchorPoint', 'LeftTop', 'TextColor', 'white',...
%                 'BoxColor', 'yellow','BoxOpacity', .6, 'Font', 'Lucida Sans Regular', 'FontSize', round(BVFontSize*.6));
%         end
    end
    rescaled =  im2uint8(newFrames(:,:,:,frameInd));%rescale(newFrames(:,:,:,frameInd),'InputMin', zeros(size(newFrames,1),size(newFrames,2)),...
        %'InputMax', ones(size(newFrames,1), size(newFrames,2)));
    writeVideo(FW, rescaled);
    frameInd = frameInd + 1;
    if(showBehavior)
        continueFrames = BV.hasFrame;
    else
        continueFrames = frameInd <= length(tTestDir);
    end
end
close(FW);