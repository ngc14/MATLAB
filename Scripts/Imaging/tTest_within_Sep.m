function [] = tTest_within_Baseline(varargin)
%% user input
username = 'ngc14';

animal = 'Gilligan';
date = '12_12_2018';
run = 'run00';

ALIGN_OFFSET = true;
stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere', 'LargeSphere', 'Photocell'});
p_values = [0.01 0.001 0.0001];
conditions = [2,4,5];

frameRange = [58:60]; %[16:17];
excludeFrames = []; %[5:12];
style=4;      % 1: simple DC, 2:ffsub (first frame subtraction), 3: subtraction (stim1 -stim2), 4: subtraction after ffsub
fframe=[1];   % first frame (only used for style 2&4)
stim2=[0];    % only used in style 3&4

% Hisashi Filtering
% choose filtering method from the followings
% 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
% 'slowmean' 'disk-like' but better at edge
% 'gaussian' Gaussian filter with half sd
% 'fastmedian' median filter
% 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first,
% apply median filter, then
% expand the image to the original size.
% 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
% 'fft', band-pass filtering using FFT (Fast Fourier Transform).

LPMethod='gaussian';        % Low-pass filtering method, usually 'gaussian'
LPk=10;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=550;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];

numCols = 10;                % Number of columns in big image.
blockSelect = []; % Use only a subset of blocks. Leave empty to use all blocks.


if ~isempty(varargin)
    % animal
    if ~isempty(varargin{1})
        animal = varargin{1};
    end
    
    % hemi
    if ~isempty(varargin{2})
        hemi = varargin{2};
    end
    
    % date
    if ~isempty(varargin{3})
        date = varargin{3};
    end
    
    % run
    if ~isempty(varargin{4})
        run = varargin{4};
    end
    
    % stim #
    if ~isempty(varargin{5})
        stim1 = varargin{5};
    end
    
    % stim name
    if ~isempty(varargin{6})
        stim(stim1) = varargin{6};
    end
    
    % style
    if ~isempty(varargin{7})
        style = varargin{7};
    end
    
    % LPk
    if ~isempty(varargin{8})
        LPk = varargin{8};
    end
    
    % HPk
    if ~isempty(varargin{9})
        HPk = varargin{9};
    end
    
    % clipMethod
    if ~isempty(varargin{10})
        clipMethod = varargin{10};
    end
    
    % clipValue
    if ~isempty(varargin{11})
        clipValue = varargin{11};
    end
    
    % frameRange
    if ~isempty(varargin{12})
        frameRange = varargin{12};
    end
    
    if ~isempty(varargin{13})
        stim2 = varargin{13};
    end
end
%% locate data
dataPath = ['S:\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
greenName = ['green',run(end-1:end),'_edited.bmp'];
parent = ['S:\Lab\',username,'\Scripts\'];
blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);
if isempty(blocks)
    error('Cannot find data on SSD or HDD. Exiting.');
end

if ~isempty(blockSelect)
    disp(['Only using selected blocks (',num2str(length(blockSelect)),' total)']);
    blocks = blocks(blockSelect);
end
try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end

conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
numConds = unique(conds);
badAlign = load('bad_alignments.mat');
align = load('alignment_offset.mat');

% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;

% clip mask
try
    clipMask = imread(clipMaskName);
    clipMask = clipMask(:,:,1)>0;
catch error
    error('Cannot load clip mask.');
end


if size(clipMask,1) == 1082
    clipMask = clipMask(2:end-1,:);
end
if size(clipMask,2) == 1312
    clipMask = clipMask(:,3:end-2);
end
green = imread(greenName);
greenGray = double(rgb2gray(green));
if size(greenGray,1) == 1082
    greenGray = greenGray(2:end-1,:);
end
if size(greenGray,2) == 1312
    greenGray = greenGray(:,3:end-2);
end
clipMask = ~(~clipMask | imgaussfilt(greenGray,1) < median(greenGray(:))-0.5*std(greenGray(:)));
clear green
%% Stim2 average
blockConds2 = blocks(conds==stim2);
alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));
badInds2 = ~ismember([1:size(blockConds2,1)],badAlign.badAlign{alignInd2});
blockConds2 = blockConds2(badInds2);
NBlocks2 = size(blockConds2,1);
if ALIGN_OFFSET
    align.disp_x_med{alignInd2} = align.disp_x_med{alignInd2}(badInds2,:);
    align.disp_y_med{alignInd2} = align.disp_y_med{alignInd2}(badInds2,:);
else
    align = [];
end
maps_stim2 = zeros(H,W,NFrames);
parfor n = 1:NBlocks2
    frames_stim2 = zeros(H,W,NFrames);
    frames_stim2(:,:,:) = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],1,'v');
    if ALIGN_OFFSET
        for f = 1:size(frames_stim2,3)
            frames_stim2(:,:,f) = medfilt2(frames_stim2(:,:,f),[5 5]);
        end
        for f = 1:size(frames_stim2,3)
            frames_stim2(:,:,f) = OIShift(frames_stim2(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
                -1*round(align.disp_y_med{alignInd2}{n,f}));
            frames_stim2(:,:,f) = imageFilter_LPHP_nsc15(frames_stim2(:,:,f), LPk, HPk, clipMask);
        end
    else
        for f = 1:size(frames_stim2,3)
            frames_stim2(:,:,f) = imageFilter_LPHP_nsc15(frames_stim2(:,:,f), LPk, HPk, clipMask);
        end
    end
    maps_stim2 = maps_stim2 + (frames_stim2 - mean(frames_stim2(:,:,fframe),3)) ./...
        mean(frames_stim2(:,:,fframe),3)
    fprintf('!');
end
maps_stim2(isinf(maps_stim2)) = 0;
maps_stim2(isnan(maps_stim2)) = 0;
maps_stim2_avg = maps_stim2 ./ NBlocks2;
avg_stim2 = mean(maps_stim2,3);
clear blockConds2 alignInd2 badInds2 frames_stim2
%% Begin Script
for s=1:length(conditions)
    stim1 = conditions(s);
    disp(['%%%%% BLOCKVIEW_NGC14: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
        num2str(LPk),' HP',num2str(HPk),' f[',num2str(frameRange),'] %%%%%']);
    
    blockConds1 = blocks(conds==stim1);
    alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
    badInds1 = ~ismember([1:size(blockConds1,1)],badAlign.badAlign{alignInd1});
    blockConds1 = blockConds1(badInds1);
    NBlocks = size(blockConds1,1);
    
    if ALIGN_OFFSET
        align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
        align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
    end
    %% load data
    maps_stim1 = zeros(H,W,NFrames,NBlocks);
    % dispstat('','init');
    parfor n = 1:NBlocks
        % reset stim1 and stim2 to zeros
        frames_stim1 = zeros(H,W,NFrames);
        
        % load stim 1
        frames_stim1(:,:,:) = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],1,'v');
        
        if ALIGN_OFFSET
            % median filter to reduce CCD column defects
            for f = 1:size(frames_stim1,3)
                frames_stim1(:,:,f) = medfilt2(frames_stim1(:,:,f),[5 5]);
            end
            for f = 1:size(frames_stim1,3)
                %             disp(['n',num2str(n),' f',num2str(f)]);
                frames_stim1(:,:,f) = OIShift(frames_stim1(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                    -1*round(align.disp_y_med{alignInd1}{n,f}));
                frames_stim1(:,:,f) = imageFilter_LPHP_nsc15(frames_stim1(:,:,f), LPk, HPk, clipMask);
            end
        else
            for f = 1:size(frames_stim1,3)
                frames_stim1(:,:,f) = imageFilter_LPHP_nsc15(frames_stim1(:,:,f), LPk, HPk, clipMask);
            end
        end
        
        % FFS
        frames_stim1 = (frames_stim1 - mean(frames_stim1(:,:,fframe),3)) ./...
            mean(frames_stim1(:,:,fframe),3);
        frames_stim2 = (frames_stim1(:,:,1:10) - mean(frames_stim1(:,:,fframe),3)) ./ ...
            mean(frames_stim1(:,:,fframe),3);
        avg_stim2 = mean(frames_stim2,3);
        maps_stim1(:,:,:,n) = (frames_stim1 - avg_stim2) ./ avg_stim2;
        %dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
        
        fprintf('!');
    end
    % dispstat('Loading Data... 100% complete.');
    maps_stim1(isinf(maps_stim1)) = 0;
    maps_stim1(isnan(maps_stim1)) = 0;
    for f = 1:NFrames
        [h,p(:,:,f)] = ttest2(maps_stim2(:,:,1:7,:), maps_stim1(:,:,f,:), 0.05, 'left', 'equal', 4);
    end
    clear maps_stim1
    %%
    refMask = imread('\\univ.pitt.edu\sni\Gharbawie\Lab\Gilligan\All Data\Mapping\referenceMask', 'bmp')>200;
    refMask = refMask(:,:,1);
    tform = load([dataPath,'\transformation_matrix.mat']);
    tform = inv(tform.tform);
    tform(:,3) = [0 0 1];
    refMask = imwarp(refMask, affine2d(tform),'OutputView',imref2d([W,H]));
    clipMask = ~(~clipMask | ~refMask);
    %% construct big images
    for i = 1:length(p_values)
        pVals = imfill(bwareaopen((p<p_values(i) & clipMask),10),'holes');
        
        numRows = ceil(NFrames/numCols);
        
        big_map_noClip = uint8(zeros(H*(numRows+1),W*(numCols)));
        
        big_map_noClip(numRows*H+1:(numRows+1)*H,1:W) = rgb2gray(imfuse(greenGray,mean(pVals(:,:,frameRange),3)));
        
        k = 1;
        for r = 1:numRows
            for c = 1:numCols
                if k <= NFrames
                    big_map_noClip((r-1)*H+1:r*H,(c-1)*W+1:c*W) = rgb2gray(imfuse(greenGray,pVals(:,:,k)));
                    k = k + 1;
                end
            end
        end
        
        big_map_noClip_small = imresize(big_map_noClip,0.25);
        
        %% save data
        cd(dataPath);
        savePath = ['Results\tTest\',stim(conditions(s)),'\LP',num2str(LPk),'_HP',num2str(HPk),'\p-',num2str(p_values(i)),'\'];
        if ~exist(savePath,'dir')
            mkdir(savePath);
        end
        cd(savePath);
        
        imwrite(big_map_noClip_small,['tTest_p-',num2str(p_values(i)),'.bmp']);
        for f = 1:NFrames
            imwrite(rgb2gray(imfuse(greenGray,pVals(:,:,f))), [num2str(f), '.bmp']);
        end
        saveas(plot(squeeze(sum(sum(pVals)))),['sigPix_p-',num2str(p_values(i)),'.fig']);
    end
    cd(parent);
end
end