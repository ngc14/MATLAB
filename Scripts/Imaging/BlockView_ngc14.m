function [] = BlockView_ngc14_V2(varargin)

%% user input
username = 'ngc14';

animal = 'Gilligan';
date = '10_05_2018';
run = 'run00';

ALIGN_OFFSET = true;
stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere', 'LargeSphere', 'Photocell'});

frameRange = [15:25]; %[16:17];
excludeFrames = []; %[5:12];
style=1;      % 1: simple DC, 2:ffsub (first frame subtraction), 3: subtraction (stim1 -stim2), 4: subtraction after ffsub
fframe=[1];   % first frame (only used for style 2&4)
stim1=[0];    % used in all styles, if more than one stim, will take average
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
LPk=0;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=0;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter
FILTER_SINGLE_CONDITIONS = false; % filter single conditions instead of blank-subtracted images.

% Hisashi's Clipping type
% 0: no clipping;
% 1: clipping at median+-(SD*clipvalue);
% 2: clipping at median+-(SD*clipvalue) with a mask (default.bmp);
% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix [x1, y1; x2, y2])
% 4: clipping at median+-intensity change (intensity change = clipvalue), Ex. 0.0005 or 0.0008;
% 5: clipping at 0+-(SD*clipvalue);
% 6: clipping at 0+-intensity change (intensity change = clipvalue);
clipMethod = 0;
clipValue  = 0;
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];

numCols = 10;                % Number of columns in big image.
avgFrame = true;            % add average frame onto bottom of image?
SAVE_IND_FRAMES = true;     % save individual frames or not?

blockSelect = []; % Use only a subset of blocks. Leave empty to use all blocks.


%% multirun
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
    
    % excludeFrames
    if ~isempty(varargin{13})
        excludeFrames = varargin{13};
    end
    
    clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
end



%% Begin Script
disp(['%%%%% BLOCKVIEW_NGC14: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' f[',num2str(frameRange),'] %%%%%']);


%% locate data
dataPath = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
parent = ['\\pitt\sni\gharbawie\Lab\',username,'\Scripts\'];

try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end

blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);
if isempty(blocks) %look on backup drive instead
    disp('Data not found on SSD. Attempting to load from HDD.');
    
    try
        cd(['D:\Data\',animal,'_SqM\',hemi,'Hemisphere\',date,'\',run]);
    catch error
        error('HDD data path does not exist. Exiting.');
    end
    
    blocks = dir('*.BLK');
    cd(dataPath);
end

if isempty(blocks)
    error('Cannot find data on SSD or HDD. Exiting.');
end

if ~isempty(blockSelect)
    disp(['Only using selected blocks (',num2str(length(blockSelect)),' total)']);
    blocks = blocks(blockSelect);
end

% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = size(blocks,1);


% clip mask
if clipMethod==2
    try
        clipMask = imread(clipMaskName);
    catch error
        error('Cannot load clip mask.');
    end
    
    clipMask = clipMask(:,:,1)>0;
    
    if size(clipMask,1) == 1082
        clipMask = clipMask(2:end-1,:);
    end
    if size(clipMask,2) == 1312
        clipMask = clipMask(:,3:end-2);
    end
    
else
    clipMask = logical(ones(H,W));
end


%% load data
maps_stim1 = zeros(H,W,NFrames);
maps_stim2 = zeros(H,W,NFrames);
maps_avg = zeros(H,W,NFrames);
conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
numConds = unique(conds);

blockConds2 = blocks(conds==stim2);
blockConds1 = blocks(conds==stim1);
if ALIGN_OFFSET
    align = load('alignment_offset.mat');
    badAlign = load('bad_alignments.mat');
    alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
    alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));
    badInds1 = ~ismember([1:size(blockConds1,1)],badAlign.badAlign{alignInd1});
    badInds2 = ~ismember([1:size(blockConds2,1)],badAlign.badAlign{alignInd2});
    blockConds1 = blockConds1(badInds1);
    blockConds2 = blockConds2(badInds2);
    align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
    align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
    if(alignInd1~=alignInd2)
        align.disp_x_med{alignInd2} = align.disp_x_med{alignInd2}(badInds2,:);
        align.disp_y_med{alignInd2} = align.disp_y_med{alignInd2}(badInds2,:);
    end
else
  align = [];  
end

% dispstat('','init');
parfor n = 1:min(size(blockConds2,1), size(blockConds1,1))
    
    % reset stim1 and stim2 to zeros
    frames_stim1 = zeros(H,W,NFrames);
    frames_stim2 = zeros(H,W,NFrames);
    
    % load stim 1
    frames_stim1(:,:,:) = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],1,'v');
    
    
    %if blank sub, load stim 2
    if style==3 || style==4
        frames_stim2(:,:,:) = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],1,'v');
    end
    
    if ALIGN_OFFSET
        % median filter to reduce CCD column defects
        for f = 1:size(frames_stim1,3)
            frames_stim1(:,:,f) = medfilt2(frames_stim1(:,:,f),[5 5]);
            
            if style==3 || style==4
                frames_stim2(:,:,f) = medfilt2(frames_stim2(:,:,f),[5 5]);
            end
        end
        alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
        alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));
        for f = 1:size(frames_stim1,3)
            %             disp(['n',num2str(n),' f',num2str(f)]);
            frames_stim1(:,:,f) = OIShift(frames_stim1(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                -1*round(align.disp_y_med{alignInd1}{n,f}));
            
            if style==3 || style==4
                frames_stim2(:,:,f) = OIShift(frames_stim2(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
                    -1*round(align.disp_y_med{alignInd2}{n,f}));
            end
            
            disp('');
        end
    end
    
    %{
    if FILTER_EACH_TRIAL
        % filter each frame of each stim for this block
        for f = 1:NFrames
            frames_stim1(:,:,f) = OIEasyFilterH2(frames_stim1(:,:,f), LPMethod, LPk, HPMethod, HPk);

            if style==3 || style==4
                frames_stim2(:,:,f) = OIEasyFilterH2(frames_stim2(:,:,f), LPMethod, LPk, HPMethod, HPk);
            end
        end
    end
    %}
    
    switch style
        case 1 % single condition
            maps_stim1 = maps_stim1 + frames_stim1;
            
        case 2 % single condition w/ FFS
            maps_stim1 = maps_stim1 + (frames_stim1 - frames_stim1(:,:,fframe)) ./...
                frames_stim1(:,:,fframe);
            
        case 3 % blank subtraction
            maps_stim1 = maps_stim1 + frames_stim1;
            maps_stim2 = maps_stim2 + frames_stim2;
            
        case 4 % FFS and then blank subtraction
            maps_stim1 = maps_stim1 + (frames_stim1 - frames_stim1(:,:,fframe)) ./...
                frames_stim1(:,:,fframe);
            maps_stim2 = maps_stim2 + (frames_stim2 - frames_stim2(:,:,fframe)) ./...
                frames_stim2(:,:,fframe);
    end
    
    %     dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
end
% dispstat('Loading Data... 100% complete.');


%% average maps and blank subtraction
switch style
    case 1 % single condition
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_avg = maps_stim1_avg;
    case 2 % single condition
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_avg = maps_stim1_avg;
    case 3 % blank subtraction
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_stim2_avg = maps_stim2 ./ NBlocks;
        maps_avg = maps_stim1_avg - maps_stim2_avg;
    case 4 % blank subtraction
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_stim2_avg = maps_stim2 ./ NBlocks;
        maps_avg = maps_stim1_avg - maps_stim2_avg;
end

% subtract average
if style==1
    maps_avg = maps_avg - mean(maps_avg(:));
end

% exclude frames
if ~isempty(excludeFrames)
    maps_avg(:,:,excludeFrames) = 0;
end


%% filtering and clipping
maps_avg_indi = double(zeros(H,W,NFrames));
for f = 1:NFrames
    if (style==2 || style==4)
        if f~=fframe % dont filter fframe if FFS
            maps_avg(:,:,f) = imageFilter_LPHP_nsc15(maps_avg(:,:,f),LPk,HPk,clipMask);
            
            if FILTER_SINGLE_CONDITIONS
                maps_stim1_avg(:,:,f) = imageFilter_LPHP_nsc15(maps_stim1_avg(:,:,f),LPk,HPk,clipMask);
                maps_stim2_avg(:,:,f) = imageFilter_LPHP_nsc15(maps_stim2_avg(:,:,f),LPk,HPk,clipMask);
            end
            
            maps_avg_indi(:,:,f) = OIClipH2(maps_avg(:,:,f), clipMethod, clipValue, clipMask);
        end
        
    else
        maps_avg_indi(:,:,f) = OIClipH2(maps_avg(:,:,f), clipMethod, clipValue, clipMask);
    end
end


% uniform clip
k = 1;
for f = 1:NFrames
    if isempty(find(excludeFrames==f,1))
        temp = maps_avg(:,:,f);
        maps_avg_GF_masked(:,k) = temp(clipMask);
        k = k + 1;
    end
end

if style==2 || style==4
    maps_avg_GF_masked(:,fframe) = []; %remove first frame from consideration if FFS
end

maps_avg_median_uni = median(maps_avg_GF_masked(:));
maps_avg_std_uni = std(maps_avg_GF_masked(:));
uniClipLim = [maps_avg_median_uni-(maps_avg_std_uni*clipValue)...
    maps_avg_median_uni+(maps_avg_std_uni*clipValue)];

parfor f = 1:NFrames
    maps_avg_uni(:,:,f) =  OIClipH2(maps_avg(:,:,f), 9, uniClipLim, []);
end



%% construct big images
numRows = ceil(NFrames/numCols);

if avgFrame
    big_map_uni = double(zeros(H*(numRows+1),W*(numCols)));
    big_map_indi = uint8(zeros(H*(numRows+1),W*(numCols)));
    big_map_noClip = double(zeros(H*(numRows+1),W*(numCols)));
    
    big_map_uni(numRows*H+1:(numRows+1)*H,1:W) = mean(maps_avg_uni(:,:,frameRange),3);
    big_map_indi(numRows*H+1:(numRows+1)*H,1:W) = norm_to_uint8(mean(maps_avg_indi(:,:,frameRange),3));
    big_map_noClip(numRows*H+1:(numRows+1)*H,1:W) = mean(maps_avg(:,:,frameRange),3);
else
    big_map_uni = double(zeros(H*numRows,W*numCols));
    big_map_indi = uint8(zeros(H*numRows,W*numCols));
    big_map_noClip = double(zeros(H*numRows,W*numCols));
end

k = 1;
for r = 1:numRows
    for c = 1:numCols
        if k <= NFrames
            big_map_uni((r-1)*H+1:r*H,(c-1)*W+1:c*W) = maps_avg_uni(:,:,k);
            big_map_indi((r-1)*H+1:r*H,(c-1)*W+1:c*W) = norm_to_uint8(maps_avg_indi(:,:,k));
            big_map_noClip((r-1)*H+1:r*H,(c-1)*W+1:c*W) = maps_avg(:,:,k);
            k = k + 1;
        end
    end
end

big_map_uni_small = imresize(big_map_uni,0.25);
big_map_indi_small = imresize(big_map_indi,0.25);
big_map_noClip_small = imresize(big_map_noClip,0.25);


%% construct structure for .mat file
if style==4
    blockviewStruct.animal = animal;
    blockviewStruct.date = date;
    blockviewStruct.run = run;
    blockviewStruct.filtering = ['LP',num2str(LPk),'HP',num2str(HPk)];
    blockviewStruct.FFrame = fframe;
    blockviewStruct.stim1.frames = maps_stim1_avg;
    blockviewStruct.stim1.name = stim(stim1);
    blockviewStruct.stim2.frames = maps_stim2_avg;
    blockviewStruct.stim2.name = stim(stim2);
end


%% save data
cd(dataPath);

switch style
    case 1
        stimName = ['[',stim(stim1),']'];
        ffName = 'NoFF';
    case 2
        stimName = ['[',stim(stim1),']'];
        ffName = ['FF_',num2str(fframe)];
    case 3
        stimName = ['[',stim(stim1),']-[',stim(stim2),']'];
        ffName = 'NoFF';
    case 4
        stimName = ['[',stim(stim1),']-[',stim(stim2),']'];
        ffName = ['FF_',num2str(fframe)];
end


savePath = ['Results\Blockview_nsc15_V2\',stimName,'\',ffName,'\',...
    'LP',num2str(LPk),'_HP',num2str(HPk),'\'];

if ALIGN_OFFSET
    savePath = [savePath,'Align_Offset\'];
end

if ~exist(savePath,'dir')
    mkdir(savePath);
end
cd(savePath);

% save big images
if ~exist('AllFrames','dir')
    mkdir('AllFrames');
end
imwrite(norm_to_uint8(big_map_uni_small),gray(256),['AllFrames\gray-uni_clip_',num2str(clipValue),'.bmp']);
imwrite(big_map_indi_small,gray(256),['AllFrames\gray-indi_clip_',num2str(clipValue),'.bmp']);
imwrite(mat2gray(big_map_noClip_small),'AllFrames\gray-noClip.bmp');

%% save individual frames
if SAVE_IND_FRAMES
    if ~exist('IndFrames','dir')
        mkdir('IndFrames');
    end
    cd('indFrames');
    if ~exist(['uni_clip_sd_',num2str(clipValue)],'dir')
        mkdir(['uni_clip_sd_',num2str(clipValue)]);
    end
    if ~exist(['indi_clip_sd_',num2str(clipValue)],'dir')
        mkdir(['indi_clip_sd_',num2str(clipValue)]);
    end
    
    cd(['uni_clip_sd_',num2str(clipValue)]);
    uni_minmax = [min(maps_avg_uni(:)) max(maps_avg_uni(:))];
    imwrite(norm_to_uint8b(mean(maps_avg_uni(:,:,frameRange),3),uni_minmax(1),uni_minmax(2)),['avgFrame[',num2str(frameRange),'].bmp']);
    for f = 1:NFrames
        if f<10
            imwrite(norm_to_uint8b(maps_avg_uni(:,:,f),uni_minmax(1),uni_minmax(2)),gray(256),['Frame_0',num2str(f),'_clip_',num2str(clipValue),'.bmp']);
        else
            imwrite(norm_to_uint8b(maps_avg_uni(:,:,f),uni_minmax(1),uni_minmax(2)),gray(256),['Frame_',num2str(f),'_clip_',num2str(clipValue),'.bmp']);
        end
    end
    
    cd('..');
    cd(['indi_clip_sd_',num2str(clipValue)]);
    imwrite(norm_to_uint8(mean(maps_avg_indi(:,:,frameRange),3)),['avgFrame[',num2str(frameRange),'].bmp']);
    for f = 1:NFrames
        if f<10
            imwrite(mat2gray(maps_avg_indi(:,:,f)),['Frame_0',num2str(f),'_clip_',num2str(clipValue),'.bmp']);
        else
            imwrite(mat2gray(maps_avg_indi(:,:,f)),['Frame_',num2str(f),'_clip_',num2str(clipValue),'.bmp']);
        end
    end
    
end

cd(parent);

end