% function [] = imageAlignment_nsc15()
%% user input
close all;
clearvars;

username = 'ngc14';

animal = 'Gilligan';
date = '11_16_2018';
run = 'run00';
stimSelect = [0 2 3 4];

alignType = 'similarity';   % https://www.mathworks.com/help/images/ref/fitgeotrans.html
% https://www.mathworks.com/examples/computer-vision/mw/vision-ex57604986-image-registration-using-multiple-features


%% variable initialization, find block data
parent = ['\\130.49.229.252\gharbawie\Lab\',username,'\Scripts\'];
path = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
path = ['\\pitt\sni\Gharbawie\Lab\ngc14\MsHowell\Awake_Cutaneous\2014_10_22_MsHowell\run0'];

cd(path);
blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);

anapar = OIHeadRead((blocks(1).name), 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = length(blocks);

timeStart = clock;
dispstat('Loading & Aligning Data... 0.00% Done','init');
% addAttachedFiles(gcp,'C:\Program Files\MATLAB\R2018a\toolbox\vision\vision\');
% addAttachedFiles(gcp,'\\130.49.229.252\gharbawie\Lab\ngc14\Scripts\Imaging\');

%% main alignment loop
for c = 1:length(stimSelect)
for n = 1:size(blocks,1)
    loopStart = clock;
    frames = OIReadStim(blocks(n).name, stimSelect(c), 'v');
    
    % align to first frame of first block
    if n==1
        original = frames(:,:,1);
        original = original - mean(original(:)); % mean-subtract and clip image to increase contrast
        [original, ~, clipLow, clipHigh] = OIClipH2(original, 1, 1.5, ones(768,768));
        frame_original = original;
        
        % detect features in original image
        ptsOriginalBRISK  = detectBRISKFeatures(original, 'MinContrast', 0.9,'MinQuality',0.9);
        ptsOriginalSURF  = detectSURFFeatures(original);
        [featuresOriginalFREAK,  validPtsOriginalBRISK]  = extractFeatures(original,  ptsOriginalBRISK);
        [featuresOriginalSURF,  validPtsOriginalSURF]  = extractFeatures(original,  ptsOriginalSURF);
        outputView = imref2d(size(original));
    end
    
    % align all frames
    for f = 1:NFrames
        % prepare images for alignment
        distorted = frames(:,:,f);
        distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
        [distorted, ~, ~, ~] = OIClipH2(distorted, 9, [clipLow clipHigh], ones(768,768));
        if(any(any(distorted)))
            
            % use BRISK and SURF for image alignment
            ptsDistortedBRISK = detectBRISKFeatures(distorted, 'MinContrast', 0.9,'MinQuality',0.9);
            ptsDistortedSURF = detectSURFFeatures(distorted);
            [featuresDistortedFREAK, validPtsDistortedBRISK] = extractFeatures(distorted, ptsDistortedBRISK);
            [featuresDistortedSURF, validPtsDistortedSURF] = extractFeatures(distorted, ptsDistortedSURF);
            
            
            indexPairsBRISK = matchFeatures(featuresOriginalFREAK, featuresDistortedFREAK, 'MatchThreshold', 40, 'MaxRatio', 0.8);
            indexPairsSURF = matchFeatures(featuresOriginalSURF, featuresDistortedSURF);
            matchedDistortedBRISK = validPtsDistortedBRISK(indexPairsBRISK(:,2));
            matchedDistortedSURF = validPtsDistortedSURF(indexPairsSURF(:,2));
            matchedDistortedXY = [matchedDistortedSURF.Location; matchedDistortedBRISK.Location];
            matchedOriginalBRISK  = validPtsOriginalBRISK(indexPairsBRISK(:,1));
            matchedOriginalSURF  = validPtsOriginalSURF(indexPairsSURF(:,1));
            matchedOriginalXY  = [matchedOriginalSURF.Location; matchedOriginalBRISK.Location];
            
            % initial alignment
            [tformTotal,inlierDistortedXY,inlierOriginalXY] = ...
                estimateGeometricTransform(matchedDistortedXY,matchedOriginalXY,alignType,...
                'maxNumTrials',10000,'confidence',99,'MaxDistance',5);
            
            % get rid of outlier point pairs
            distTemp = sqrt((inlierOriginalXY(:,1)-inlierDistortedXY(:,1)).^2 + ...
                (inlierOriginalXY(:,2)-inlierDistortedXY(:,2)).^2);
            distThreshHigh = median(distTemp) + 3*std(distTemp);
            distThreshLow = median(distTemp) - 3*std(distTemp);
            
            if std(distTemp) > 0
                inlierDistortedXY = inlierDistortedXY((distTemp > distThreshLow) & (distTemp < distThreshHigh),:);
                inlierOriginalXY = inlierOriginalXY((distTemp > distThreshLow) & (distTemp < distThreshHigh),:);
                
                % re-do alignment without outlier points
                [tformTotal,inlierDistortedXY,inlierOriginalXY] = ...
                    estimateGeometricTransform(inlierDistortedXY,inlierOriginalXY,alignType,...
                    'maxNumTrials',10000,'confidence',99,'MaxDistance',5);
                
                distCond{n,f} = sqrt((inlierOriginalXY(:,1)-inlierDistortedXY(:,1)).^2 + ...
                    (inlierOriginalXY(:,2)-inlierDistortedXY(:,2)).^2);
            else
                distCond{n,f} = distTemp;
            end
            
            % x/y displacement between matched point pairs
            displacement_xCond{n,f} = inlierDistortedXY(:,1)-inlierOriginalXY(:,1);
            displacement_yCond{n,f} = inlierDistortedXY(:,2)-inlierOriginalXY(:,2);
            
            %frame_distorted{n,f} = distorted;
            pts_originalCond{n,f} = inlierOriginalXY;
            pts_distortedCond{n,f} = inlierDistortedXY;
            align_XYZRS(:,:,n,f) = tformTotal.T;
        end
    end
    
    % display function progress
    percDone = (c-1)/length(stimSelect) + (1/length(stimSelect) * n/size(blocks,1));
    percDone = round(percDone*10000)/100;
    elapsedTime = etime(clock,timeStart);
    dispstat(['Loading & Aligning Data... ',num2str(percDone),'% Done.']);
end
    dist{c} = distCond;
    displacement_x{c} = displacement_xCond;
    displacement_y{c} = displacement_yCond;
    pts_original{c} = pts_originalCond;
    pts_distorted{c} = pts_distortedCond;
end
totalTime = round(etime(clock,timeStart)*100/60)/100;
dispstat(['Loading & Aligning Data... 100.00% Done. Total time: ',num2str(totalTime),' min.']);

%cd(parent);


%% show movie
if 1==0
    for n = 1:NBlocks
        frames = OIReadStim(blocks(n).name, stimSelect, 'v');
        for f = 1:NFrames
            if n==1 && f==1, figure; end
            
            distorted = frames(:,:,f);
            distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
            [distorted, ~, ~, ~] = OIClipH2(distorted, 9, [clipLow clipHigh], ones(768,768));
            subplot(5,8,[1:4,9:12,17:20,25:28]);
            imagesc(frame_original);
            axis image off; colormap gray;
            hold on;
            plot(pts_original{n,f}(:,1),pts_original{n,f}(:,2),'ro','linewidth',1.5);
            hold off;
            title('Block 1 Frame 1');
            
            subplot(5,8,[5:8,13:16,21:24,29:32]);
            imagesc(distorted);
            axis image off; colormap gray;
            hold on;
            plot(pts_distorted{n,f}(:,1),pts_distorted{n,f}(:,2),'ro','linewidth',1.5);
            hold off;
            title(['Block ',num2str(n),' Frame ',num2str(f)]);
            
            subplot(5,8,33:40);
            histogram(dist{n,f},20,'binlimits',[0 5]);
            xlim([0 5]);
            
            pause(1);
        end
    end
end

%% movie 2
if 1==0
    figure;
    for n = 1:NBlocks
        frames = OIReadStim(blocks(n).name, stimSelect, 'v');
        for f = 1:NFrames
            
            distorted = frames(:,:,f);
            distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
            [distorted, ~, ~, ~] = OIClipH2(distorted, 9, [clipLow clipHigh], ones(768,768));
            imagesc(distorted);
            axis image off; colormap gray;
            title(['Block ',num2str(n),' Frame ',num2str(f)]);
            
            pause(0.1);
        end
    end
end


%% distance box plot
figure;
subplot(3,1,1);
for n = 1:NBlocks
    for f = 1:NFrames
        plot([n+(f-1)/NFrames n+(f-1)/NFrames],[median(dist{n,f})-std(dist{n,f})...
            median(dist{n,f})+std(dist{n,f})],'-','linewidth',0.5,'color',[0.6 0.6 0.6]);
        
        if n==1 && f==1, hold on, end
    end
    
    for f = 1:NFrames
        plot(n+(f-1)/NFrames,median(dist{n,f}),'k.','markersize',20);
    end
end

ylabel('Median Dist Offset (px)');
set(gca,'fontsize',14);

subplot(3,1,2);
for n = 1:NBlocks
    for f = 1:NFrames
        plot([n+(f-1)/NFrames n+(f-1)/NFrames],[median(displacement_x{n,f})-std(displacement_x{n,f})...
            median(displacement_x{n,f})+std(displacement_x{n,f})],'-','linewidth',0.5,'color',[0.8 0.5 0.5]);
        
        if n==1 && f==1, hold on, end
    end
    
    for f = 1:NFrames
        plot(n+(f-1)/NFrames,median(displacement_x{n,f}),'r.','markersize',20);
    end
end

ylabel('Median X Offset (px)');
set(gca,'fontsize',14);

subplot(3,1,3);
for n = 1:NBlocks
    for f = 1:NFrames
        plot([n+(f-1)/NFrames n+(f-1)/NFrames],[median(displacement_y{n,f})-std(displacement_y{n,f})...
            median(displacement_y{n,f})+std(displacement_y{n,f})],'-','linewidth',0.5,'color',[0.5 0.5 0.8]);
        
        if n==1 && f==1, hold on, end
    end
    
    for f = 1:NFrames
        plot(n+(f-1)/NFrames,median(displacement_y{n,f}),'b.','markersize',20);
    end
end

ylabel('Median Y Offset (px)');
set(gca,'fontsize',14);

xlabel('Block #');
save('offsets.mat', 'displacement_x', 'displacement_y', 'dist');

% end