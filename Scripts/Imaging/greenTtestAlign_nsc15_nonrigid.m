% function [] = greenTtestAlign_nsc15_nonrigid(varargin)

close all;
clearvars;

%% user input
AUTOMATIC = false;
ADD_POINTS = true;

USE_EXISTING_TFORM = false;
SAVE = true;
PLOT = true;

CLIPMASK = true;

monkey = 'Skipper';

dateMoving = '12_09_2020';
runMoving = 'run00';

PROJECTIVE_THEN_NONRIGID = false; % projective alignment to intermediate step first.
    intStepDate = '07_05_2019';  % projective and nonrigid transforms must already exist.
    intStepRun = 'run00';

CUSTOM_MOVING = false;
CUSTOM_MOVING_PATH = ['C:\Users\nsc15\Documents\Data\Bordeaux_SqM\LeftHemisphere',...
    '\12_11_2018\Green Collage (50-85)\Combined Green Edited.bmp'];

ALIGN_TTESTS = false;
ALIGN_SUNIN_MAT = false;
framesMoving = '8-10';
LPMoving = 5;
HPMoving = 550;
% pValue = 0.001;
clipMoving = .75;

customFixedGreenName = ''; % if empty, load default green name


% fixed green date and run
if strcmpi(monkey,'Gilligan')
        dateFixed = '12_12_2018';
        runFixed = 'run00';
        pathFixed = ['S:\Lab\',monkey,'\Mapping'];
elseif strcmpi(monkey,'Skipper')
    dateFixed = '10_30_2020';
    runFixed = 'run00';
    pathFixed = ['S:\Lab\',monkey,'\All Data\', monkey, '_', dateFixed,'\Imaging\', runFixed];

end


%% display info
if AUTOMATIC
    autoStr = 'Automatically';
else
    autoStr = 'Manually';
end
disp(upper(['%%% ',autoStr,' aligning [',dateMoving(dateMoving~='_'),' ',runMoving,...
    '] to [',dateFixed(dateFixed~='_'),' ',runFixed,'] %%%']));


%% load fixed green

if isempty(customFixedGreenName)
    greenFixed = imread([pathFixed,'\M1_green_edited.bmp']);
else
    greenFixed = imread([pathFixed,'\',customFixedGreenName]);
end

%% load moving green
if ~CUSTOM_MOVING
    pathMoving = ['S:\Lab\',monkey,'\All Data\',monkey,'_',dateMoving,'\Imaging\',runMoving];
    greenMoving = imread([pathMoving,'\green',runMoving(end-1:end),'_edited.bmp']);
%     pathMoving = ['S:\Lab\Gilligan\All Data\Gilligan_12_12_2018\Imaging\run00'];
%     greenMoving = imread([pathMoving,'\green00_edited.bmp']);
    if CLIPMASK
        clipMaskMoving = imread([pathMoving,'\clipMask',runMoving(end-1:end),'.bmp']);
        clipMaskMoving = logical(clipMaskMoving(:,:,1));
    end
else
    pathMoving = CUSTOM_MOVING_PATH;
    greenMoving = imread(CUSTOM_MOVING_PATH);
    
    if CLIPMASK
        clipMaskMoving = imread([pathMoving,'\clipMask',runMoving(end-1:end),'.bmp']);
        clipMaskMoving = logical(clipMaskMoving(:,:,1));
    end
end

%% convert to gray image
if size(greenMoving,3)>1
    greenMoving = rgb2gray(greenMoving);
end

if size(greenFixed,3)>1
    greenFixed = rgb2gray(greenFixed);
end

%% load ttests
if ALIGN_TTESTS
    cd([pathMoving,'\Results\Ttest_nsc15_V4\Frames',framesMoving,...
        '\LP',num2str(LPMoving),'HP',num2str(HPMoving)]);
    ttestFiles = dir('*p-0.*.bmp');
    suninFile = dir('*_submap*.bmp');
    
    for f = 1:length(ttestFiles)
        ttestMoving{f} = imread(ttestFiles(f).name);
    end
    
    for f = 1:length(suninFile)
        suninMoving{f} = imread(suninFile(f).name);
    end
end


%% projective tform to intermediate step, then nonrigid tform to fixed green
if PROJECTIVE_THEN_NONRIGID
    disp(['Projective tform to ',intStepDate,' ',intStepRun]);
    % load projective tform
    projectiveTform = load([pathMoving,'\Results\Aligned_To_',intStepDate(intStepDate~='_'),...
        '_',intStepRun,'\tform.mat']);
    projectiveTform = projectiveTform.tform;
    
    % apply transform to green, clipmask, ttest, sunin
    outputView = imref2d(size(greenMoving));
    greenMoving = imwarp(greenMoving,projectiveTform,'OutputView',outputView);

    if CLIPMASK
        clipMaskMoving = imwarp(clipMaskMoving,projectiveTform,'OutputView',outputView);
    end

    if ALIGN_TTESTS
        for f = 1:length(ttestMoving)
            ttestMoving{f} = imwarp(ttestMoving{f},projectiveTform,'OutputView',outputView);
        end
        for f = 1:length(suninMoving)
            suninMoving{f} = imwarp(suninMoving{f},projectiveTform,'OutputView',outputView);
        end
    end
end


%% adjust size if necessary
if (size(greenMoving,1) < size(greenFixed,1)) || (size(greenMoving,2) < size(greenFixed,2))
    greenMovingTemp = zeros(size(greenFixed));
    greenMovingTemp(1:size(greenMoving,1),1:size(greenMoving,2)) = greenMoving;
    
    if CLIPMASK
        clipMaskMovingTemp = zeros(size(greenFixed));
        clipMaskMovingTemp(1:size(clipMaskMoving,1),1:size(clipMaskMoving,2)) = clipMaskMoving;
    end
    
    if ALIGN_TTESTS
        for i = 1:length(ttestMoving)
            ttestMovingTemp{i} = zeros(size(greenFixed));
            ttestMovingTemp{i}(1:size(ttestMoving{i},1),1:size(ttestMoving{i},2)) = ttestMoving{i};
        end
        
        for i = 1:length(suninMoving)
            suninMovingTemp{i} = zeros(size(greenFixed));
            suninMovingTemp{i}(1:size(suninMoving{i},1),1:size(suninMoving{i},2)) = suninMoving{i};
        end
    end
    
else % no resizing
    greenMovingTemp = greenMoving;
    
    if CLIPMASK
        clipMaskMovingTemp = clipMaskMoving;
    end
    
    if ALIGN_TTESTS
        ttestMovingTemp = ttestMoving;
        suninMovingTemp = suninMoving;
    end
end

parent = 'S:\Lab\ngc14\Scripts\';


%% check if tform exists already
if USE_EXISTING_TFORM
    if ~PROJECTIVE_THEN_NONRIGID
        if ~exist([pathMoving,'/Results/Aligned_To_',dateFixed(dateFixed~='_'),...
                '_',runFixed,'_NONRIGID/tform_nonrigid.mat'])
            USE_EXISTING_TFORM = false;
            disp('No existing tform found.');
        end
    else
        intPathMoving = ['C:\Users\nsc15\Documents\Data\',monkey,'_SqM\',hemi,'Hemisphere\',intStepDate,'\',intStepRun];
        if ~exist([intPathMoving,'/Results/Aligned_To_',dateFixed(dateFixed~='_'),...
                '_',runFixed,'_NONRIGID/tform_nonrigid.mat'])
            USE_EXISTING_TFORM = false;
            disp('No existing tform found.');
        end
    end
end


%% locate / select alignment points
if ~USE_EXISTING_TFORM
    %% MANUAL POINTS
    if ~AUTOMATIC
        [ptFixed, ptMoving] = cpselect(greenFixed,greenMoving,'Wait',true);
        ptMoving_a = cpcorr(ptMoving,ptFixed,greenMoving,greenFixed);
        %pt2_a = pt2;

    else
    %% AUTOMATIC
        original = greenFixed;
        distorted = greenMoving;

%         ptsOriginalBRISK  = detectBRISKFeatures(original, 'MinContrast', 0.01,'MinQuality',0.02);
        ptsOriginalBRISK  = detectBRISKFeatures(original, 'MinContrast', 0.1,'MinQuality',0.05);
        ptsOriginalSURF  = detectSURFFeatures(original);
        [featuresOriginalFREAK,  validPtsOriginalBRISK]  = extractFeatures(original,  ptsOriginalBRISK);
        [featuresOriginalSURF,  validPtsOriginalSURF]  = extractFeatures(original,  ptsOriginalSURF);
        outputView = imref2d(size(original));

%         ptsDistortedBRISK = detectBRISKFeatures(distorted, 'MinContrast', 0.01,'MinQuality',0.02);
        ptsDistortedBRISK = detectBRISKFeatures(distorted, 'MinContrast', 0.1,'MinQuality',0.05);
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

        ptMoving_a = matchedDistortedXY;
        ptFixed = matchedOriginalXY;

        
        if ADD_POINTS
            [~, pt2_inlier, pt_inlier] = estimateGeometricTransform(ptMoving_a,ptFixed,transformType,...
                'maxNumTrials',20000,'confidence',99,'maxDistance',50);

            [ptFixed, ptMoving] = cpselect(greenFixed,greenMoving,double(pt_inlier),double(pt2_inlier),'Wait',true);
            ptMoving_a = cpcorr(ptMoving,ptFixed,greenMoving,greenFixed);
        end
    end
    
    
    %% nonrigid alignment
    % Fit the bspline grid to the corresponding landmarks
    options.Verbose = true;
    options.MaxRef = 9;
    
    [O_trans,Spacing] = point_registration(size(greenMovingTemp),...
        [ptFixed(:,2) ptFixed(:,1)],[ptMoving_a(:,2) ptMoving_a(:,1)],options);
    
else 
    %% load existing tform
    if ~PROJECTIVE_THEN_NONRIGID
        cd([pathMoving,'/Results/Aligned_To_',dateFixed(dateFixed~='_'),'_',runFixed,'_NONRIGID']);
    else
        cd([intPathMoving,'/Results/Aligned_To_',dateFixed(dateFixed~='_'),'_',runFixed,'_NONRIGID']);
    end
    
    nonrigidTform = load('tform_nonrigid');
    O_trans = nonrigidTform.O_trans;
    Spacing = nonrigidTform.Spacing;
    ptFixed = nonrigidTform.ptFixed;
    ptMoving = nonrigidTform.ptMoving;
    ptMoving_a = nonrigidTform.ptMoving_a;
end


%% apply tranformation to images
disp(['Nonrigid tform to ',dateFixed,' ',runFixed]);
% Transform the 2D image  
greenMoving_a = bspline_transform(O_trans,greenMovingTemp,Spacing,3);
greenMoving_a_alpha = bspline_transform(O_trans,ones(size(greenMovingTemp)),Spacing,3);

if CLIPMASK
    clipMaskMoving_a = logical(bspline_transform(O_trans,double(clipMaskMovingTemp),Spacing,3));
end

% create and transform a grid image
gridIm = ones(size(greenMovingTemp));
for x = 1:20:size(gridIm,1)
    gridIm(x,:) = 0;
end
for y = 1:20:size(gridIm,2)
    gridIm(:,y) = 0;
end
gridIm_a = bspline_transform(O_trans,gridIm,Spacing,3);


% apply transformation to ttests
if ALIGN_TTESTS
    for i = 1:length(ttestMovingTemp)
        ttestMoving_a{i} = bspline_transform(O_trans,ttestMovingTemp{i},Spacing,3);
    end
    for i = 1:length(suninMovingTemp)
        suninMoving_a{i} = bspline_transform(O_trans,suninMovingTemp{i},Spacing,3);
    end
end


%% show result
if PLOT
    %% Show the result
    ha = tight_subplot(1,3,0.02,0.03,0.03);
    
%     figure;
    axes(ha(1));
    imagesc(greenFixed); axis image off; colormap gray;
    title('Fixed Image'); 
    hold on;
    for i=1:size(ptFixed,1)
%         plot([ptMoving_a(i,1) ptFixed(i,1)],[ptMoving_a(i,2) ptFixed(i,2)],'y');
%         plot(ptMoving_a(i,1),ptMoving_a(i,2),'ro');
        plot(ptFixed(i,1),ptFixed(i,2),'go');
    end

    axes(ha(2));
    imagesc(greenMoving); axis image off; colormap gray;
    title('Moving Image'); 
    hold on;
    for i=1:size(ptFixed,1)
%         plot([ptMoving_a(i,1) ptFixed(i,1)],[ptMoving_a(i,2) ptFixed(i,2)],'y');
        plot(ptMoving_a(i,1),ptMoving_a(i,2),'ro');
%         plot(ptFixed(i,1),ptFixed(i,2),'go');
    end
    
    axes(ha(3));
    imagesc(greenMoving_a); axis image off; colormap gray;
    title('Registered Image')
    hold on;
    for i=1:size(ptFixed,1)
        plot([ptMoving_a(i,1) ptFixed(i,1)],[ptMoving_a(i,2) ptFixed(i,2)],'y');
        plot(ptMoving_a(i,1),ptMoving_a(i,2),'ro');
        plot(ptFixed(i,1),ptFixed(i,2),'go');
    end
    
    %% show overlap image
    figure; imshowpair(greenFixed,greenMoving_a); title('Overlaid images');
    
    %% Show b-spline grid
    figure; 
    subplot(1,2,1)
    imagesc(greenMoving); axis image off; colormap gray;
    hold on;
    imagesc(gridIm(1:size(greenMoving,1),1:size(greenMoving,2)),'alphadata',...
        gridIm(1:size(greenMoving,1),1:size(greenMoving,2))==0);
    axis image off; colormap gray;
    for i=1:size(ptFixed,1)
        plot(ptMoving_a(i,1),ptMoving_a(i,2),'ro');
    end
    
    subplot(1,2,2);
    imagesc(greenMoving_a); axis image off; colormap gray;
    hold on;
    imagesc(gridIm_a,'alphadata',gridIm_a<0.9);
    axis image off; colormap gray;
    for i=1:size(ptFixed,1)
        plot([ptMoving_a(i,1) ptFixed(i,1)],[ptMoving_a(i,2) ptFixed(i,2)],'y');
        plot(ptMoving_a(i,1),ptMoving_a(i,2),'ro');
        plot(ptFixed(i,1),ptFixed(i,2),'go');
    end
end

%% save tranformation
if SAVE
    cd(pathMoving);
    
    
    if ~USE_EXISTING_TFORM
        save('tform_nonrigid.mat','O_trans','Spacing','greenMovingTemp','ptFixed','ptMoving','ptMoving_a');
    end
        
%     imwrite(ind2rgb(uint8(greenMoving_a),gray(256)),...
%         [dateMoving(dateMoving~='_'),'_green',runMoving(end-1:end),'_Aligned.png'],...
%             'alpha',double(greenMoving_a_alpha>0.99));
    imwrite(norm_to_uint8(greenMoving_a),...
        [dateMoving(dateMoving~='_'),'_green',runMoving(end-1:end),'_Aligned.png'],...
            'alpha',double(greenMoving_a_alpha>0.99));
    imwrite(gridIm_a,[dateMoving(dateMoving~='_'),'_grid',runMoving(end-1:end),'_Aligned.png'],...
            'alpha',double(greenMoving_a_alpha>0.99 & gridIm_a<0.25));
        
    if CLIPMASK
        imwrite(double(clipMaskMoving_a),...
            [dateMoving(dateMoving~='_'),'_clipMask',runMoving(end-1:end),'_Aligned.png'],...
            'alpha',double(greenMoving_a_alpha>0.99 & clipMaskMoving_a>0.99));
    end
    
    if ALIGN_TTESTS        
        for f = 1:length(ttestMoving_a)
            if strfind(ttestFiles(f).name,'_rightTail')
                red(:,:,1) = zeros(size(greenMoving_a,1),size(greenMoving_a,2));
                red(:,:,2) = zeros(size(greenMoving_a,1),size(greenMoving_a,2));
                red(:,:,3) = ones(size(greenMoving_a,1),size(greenMoving_a,2));
            else
                red(:,:,1) = ones(size(greenMoving_a,1),size(greenMoving_a,2));
                red(:,:,2) = zeros(size(greenMoving_a,1),size(greenMoving_a,2));
                red(:,:,3) = zeros(size(greenMoving_a,1),size(greenMoving_a,2));
            end
            
            imwrite(red,['LP',num2str(LPMoving),'HP',num2str(HPMoving),'C',...
                num2str(clipMoving),'_F',framesMoving,'_'...
                ttestFiles(f).name(1:end-4),'_Aligned.png'],...
                'alpha',double(ttestMoving_a{f}));
        end
        
        for f = 1:length(suninMoving_a)
            imwrite(norm_to_uint8(suninMoving_a{f}),['LP',num2str(LPMoving),'HP',num2str(HPMoving),'C',...
                num2str(clipMoving),'_F',framesMoving,'_'...
                suninFile(f).name(1:end-4),'_Aligned.bmp']);
        end
    end
    
    cd(parent);
end


% end