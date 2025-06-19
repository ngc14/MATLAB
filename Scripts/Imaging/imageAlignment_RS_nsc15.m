% function [] = imageAlignment_RS_nsc15()
%% user input
% close all;
clearvars;

username = 'ngc14';

animal = 'Gilligan';
animalType = 'Macaque';
Hemi = 'Right';
date = '11_21_2018';
run = 'run00';

masterBlock = 1;
masterFrame = 1;
FOVx = 135:1120; %400:1100;
FOVy = 50:900; %100:750;

SAVE = true;

alignType = 'similarity';   
    % https://www.mathworks.com/help/images/ref/fitgeotrans.html
    % https://www.mathworks.com/examples/computer-vision/mw/vision-ex57604986-image-registration-using-multiple-features


%% variable initialization, find block data
parent = ['C:\Users\',username,'\Documents\Data'];
path = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];

cd(path);
blocks = dir('*.BLK');

if isempty(blocks) %look on backup drive instead
    disp('Data not found on SSD. Attempting to load from HDD.');
    
    try
        cd(['D:\Data\',animal,'_SqM\',Hemi,'Hemisphere\',date,'\',run]);
    catch error
        error('HDD data path does not exist. Exiting.');
    end
    
    blocks = dir('*.BLK');
    cd(path);
end

if isempty(blocks)
    error('Cannot find data on SSD or HDD. Exiting.');
end

%% sort block files
filenames = blocks;
for f = 1:length(filenames)
    ind1 = strfind(filenames(f).name,'E');
    ind2 = strfind(filenames(f).name,'B'); ind2 = ind2(1);
    ind3 = strfind(filenames(f).name,'.');
    fn_exp_blk(f,1) = str2num(filenames(f).name(ind1+1:ind2-1));
    fn_exp_blk(f,2) = str2num(filenames(f).name(ind2+1:ind3-1));
end

[fn_exp_blk, sortInd] = sortrows(fn_exp_blk,[1 2]);
filenames = filenames(sortInd);
blocks = filenames;
% blocks = blocks(1:10);

%% get block info
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = length(blocks);

timeStart = clock;
dispstat('Loading & Aligning Data... 0.00% Done','init');
addAttachedFiles(gcp,'\\pitt\sni\gharbawie\Lab\ngc14\Scripts\Imaging\');
%% main alignment loop
for n = 1:size(blocks,1)
    loopStart = clock;
    
    % align to first frame of first block
    if n==1
        framesTemp = OIReadStim([blocks(masterBlock).folder,'\',blocks(masterBlock).name], 1, 'v');
        original = framesTemp(:,:,masterFrame); 
%         clear framesTemp;
        original = original(FOVy,FOVx);
        original = original - mean(original(:)); % mean-subtract and clip image to increase contrast
        [original, ~, clipLow, clipHigh] = OIClipH2(original, 1, 1.75, ones(size(original)));
        frame_original = original;

        % detect features in original image
        ptsOriginalBRISK  = detectBRISKFeatures(original, 'MinContrast', 0.9,'MinQuality',0.9);
        ptsOriginalSURF  = detectSURFFeatures(original);
        [featuresOriginalFREAK,  validPtsOriginalBRISK]  = extractFeatures(original,  ptsOriginalBRISK);
        [featuresOriginalSURF,  validPtsOriginalSURF]  = extractFeatures(original,  ptsOriginalSURF);
        outputView = imref2d(size(original));
    end

    frames = OIReadStim([blocks(n).folder,'\',blocks(n).name], 1, 'v');
    % align all frames
    parfor f = 1:NFrames
        % prepare images for alignment
        distorted = frames(:,:,f);
        distorted = distorted(FOVy,FOVx);
        distorted = distorted - mean(distorted(:)); % mean-subtract and clip image to increase contrast
        [distorted, ~, ~, ~] = OIClipH2(distorted, 9, [clipLow clipHigh], ones(size(distorted)));

        % test offset
%             distorted = [zeros(20,768); distorted(1:end-20,:)];

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

            dist{n,f} = sqrt((inlierOriginalXY(:,1)-inlierDistortedXY(:,1)).^2 + ...
                (inlierOriginalXY(:,2)-inlierDistortedXY(:,2)).^2);
        else
            dist{n,f} = distTemp;
        end

        % x/y displacement between matched point pairs
        displacement_x{n,f} = inlierDistortedXY(:,1)-inlierOriginalXY(:,1);
        displacement_y{n,f} = inlierDistortedXY(:,2)-inlierOriginalXY(:,2);

        disp_x_med(n,f) = double(median(displacement_x{n,f}));
        disp_y_med(n,f) = double(median(displacement_y{n,f}));
        
        all_tforms{n,f} = tformTotal;

%         frame_distorted{n,f} = distorted;
%         pts_original{n,f} = inlierOriginalXY;
%         pts_distorted{n,f} = inlierDistortedXY;
%         align_XYZRS(:,:,n,f) = tformTotal.T;
    end

    % display function progress
    percDone = (n-1)/size(blocks,1);
    percDone = round(percDone*10000)/100;
    elapsedTime = etime(clock,timeStart);
    dispstat(['Loading & Aligning Data... ',num2str(percDone),'% Done.']);
end

totalTime = round(etime(clock,timeStart)*100/60)/100;
dispstat(['Loading & Aligning Data... 100.00% Done. Total time: ',num2str(totalTime),' min.']);

if SAVE
    cd(path);
    save('alignment_offset.mat','disp_x_med','disp_y_med','all_tforms');
end
cd(parent);


% %% show movie
% if 1==0
%     for n = 1:NBlocks
%         for f = 1:NFrames
%             if n==1 && f==1, figure; end
% 
%             subplot(5,8,[1:4,9:12,17:20,25:28]);
%             imagesc(frame_original);
%             axis image off; colormap gray;
%             hold on;
%             plot(pts_original{n,f}(:,1),pts_original{n,f}(:,2),'ro','linewidth',1.5);
%             hold off;
%             title('Block 1 Frame 1');
% 
%             subplot(5,8,[5:8,13:16,21:24,29:32]);
%             imagesc(frame_distorted{n,f});
%             axis image off; colormap gray;
%             hold on;
%             plot(pts_distorted{n,f}(:,1),pts_distorted{n,f}(:,2),'ro','linewidth',1.5);
%             hold off;
%             title(['Block ',num2str(n),' Frame ',num2str(f)]);
% 
%             subplot(5,8,33:40);
%             histogram(dist{n,f},20,'binlimits',[0 5]);
%             xlim([0 5]);
% 
%             pause(1);
%         end
%     end
% end
% 
% %% movie 2
% if 1==0
%     figure;
%     for n = 1:NBlocks
%         for f = 1:NFrames
%             imagesc(frame_distorted{n,f});
%             axis image off; colormap gray;
%             title(['Block ',num2str(n),' Frame ',num2str(f)]);
% 
%             pause(0.1);
%         end
%     end
% end


%% offset plot
figure('color',[1 1 1]);
plot(1:(1/NFrames):((NBlocks+1)-1/NFrames), sqrt(reshape(disp_x_med',[NBlocks*NFrames 1]).^2 +...
    reshape(disp_y_med',[NBlocks*NFrames 1]).^2),'.-');
hold on;
box off; grid on;
set(gca,'fontsize',14);
title('XY offset');
xlabel('Block #');
ylabel('Frame #');


% end