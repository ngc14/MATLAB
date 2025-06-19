function [] = ttest_ngc14_V1(varargin)

%% user input
username = 'ngc14';

animal = 'Gilligan';
date = '12_20_2018';
run = 'run01';

stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere', 'LargeSphere', 'Photocell'});
stim = containers.Map([0 2 3 4 6 7],{'Blank', 'D1_75','pD2', 'P1','D1_50','D1_25'});

ALIGN_OFFSET = true;
FIX_COLUMN_DEFECT = true;

p_values = [0.01 0.001 0.0001];

ALL_FRAMES_IND = true;
frameRange = [10:12]; %[16:17];
style=4;      % 1: simple DC, 2:ffsub (first frame subtraction), 3: subtraction (stim1 -stim2), 4: subtraction after ffsub
fframe=[1];   % first frame (only used for style 2&4)
stim1=[6];    % used in all styles, if more than one stim, will take average
stim2=[3];    % only used in style 3&4

% Filtering
LPMethod='gaussian';        % Low-pass filtering method, usually 'gaussian'
LPk=10;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=550;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter

numCols = 10;
% Hisashi's Clipping type
% 0: no clipping;
% 1: clipping at median+-(SD*clipvalue);
% 2: clipping at median+-(SD*clipvalue) with a mask (default.bmp);
% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix [x1, y1; x2, y2])
% 4: clipping at median+-intensity change (intensity change = clipvalue), Ex. 0.0005 or 0.0008;
% 5: clipping at 0+-(SD*clipvalue);
% 6: clipping at 0+-intensity change (intensity change = clipvalue);
clipMethod = 2;
clipValue  = 1;
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];


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
    
    if ~isempty(varargin{13})
        stim2 = varargin{13};
    end
    
    clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
end

greenName = ['green',run(end-1:end),'_edited.bmp'];

%% Begin Script
timeStart = clock;
disp(['%%%%% TTEST_NGC14_V1: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' f[',num2str(frameRange),'] %%%%%']);


%% locate data
dataPath = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
parent = ['\\pitt\sni\gharbawie\Lab\',username,'\Scripts\'];
%%%%%%%%%%%%%%%%%%%
dataPath = ['S:\Lab\ngc14\MsHowell\Awake_Cutaneous\2014_10_16_MsHowell\run1'];
try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end

blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);

if isempty(blocks)
    error('Cannot find data on SSD or HDD. Exiting.');
end

conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
numConds = unique(conds);

%%%%%%%%%
blockConds2 = blocks;
blockConds1 = blocks;

badAlign = load('bad_alignments.mat');
alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));
badInds1 = ~ismember([1:size(blockConds1,1)],badAlign.badAlign{alignInd1});
badInds2 = ~ismember([1:size(blockConds2,1)],badAlign.badAlign{alignInd2});
blockConds1 = blockConds1(badInds1);
blockConds2 = blockConds2(badInds2);

% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = min(size(blockConds2,1), size(blockConds1,1));


% clip mask
if clipMethod==2
    try
        clipMask = imread(clipMaskName);
        %             green = imread(greenName);
        %     greenGray = double(rgb2gray(green));
        %     if size(greenGray,1) == 1082
        %         greenGray = greenGray(2:end-1,:);
        %     end
        %     if size(greenGray,2) == 1312
        %         greenGray = greenGray(:,3:end-2);
        %     end
        %     clipMask = ~(~clipMask | imgaussfilt(greenGray,1) < median(greenGray(:))-0.5*std(greenGray(:)));
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
    clipMask = logical(ones(W,H));
end


%% load data
frames_sub = zeros(H,W,NFrames,NBlocks);

% dispstat('','init');
for n = 1:NBlocks
    % reset stim1 and stim2 to zeros
    frames_all_stim1 = zeros(H,W,NFrames);
    frames_all_stim2 = zeros(H,W,NFrames);
    
    % load frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frames_all_stim1(:,:,:) = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],alignInd1,'v');
    frames_all_stim2(:,:,:) = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],alignInd2,'v');
    
    % median filter to reduce CCD column defects
    if ALIGN_OFFSET && FIX_COLUMN_DEFECT
        for f = 1:size(frames_all_stim1,3)
            frames_all_stim1(:,:,f) = medfilt2(frames_all_stim1(:,:,f),[5 5]);
            frames_all_stim2(:,:,f) = medfilt2(frames_all_stim2(:,:,f),[5 5]);
        end
    end
    
    % alignment
    if ALIGN_OFFSET
        align = load('alignment_offset.mat');
        
        for f = 1:size(frames_all_stim1,3)
            %             disp(['n',num2str(n),' f',num2str(f)]);
            frames_all_stim1(:,:,f) = OIShift(frames_all_stim1(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                -1*round(align.disp_y_med{alignInd1}{n,f}));
            
            frames_all_stim2(:,:,f) = OIShift(frames_all_stim2(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
                -1*round(align.disp_y_med{alignInd2}{n,f}));
        end
    end
    
    if(n==1)
       greenGray = frames_all_stim1(:,:,fframe); 
    end
    % FFS
    frames_all_stim1 = (frames_all_stim1 - frames_all_stim1(:,:,fframe)) ./ frames_all_stim1(:,:,fframe);
    frames_all_stim2 = (frames_all_stim2 - frames_all_stim2(:,:,fframe)) ./ frames_all_stim2(:,:,fframe);

    
    frames_sub(:,:,:,n) = frames_all_stim1-frames_all_stim2;
    
    parfor f = 1:NFrames
        % filter both stims for each block
        frames_sub(:,:,f,n) = imageFilter_LPHP_nsc15(frames_sub(:,:,f,n), LPk, HPk, clipMask);
    end
    %     dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
    
end
% dispstat('Loading Data... 100% complete.');

for f = 1:NFrames
    [h,p(:,:,f)] = ttest2(frames_sub(:,:,f,:), frames_sub(:,:,2,:), 0.05, 'left', 'equal', 4);
end

%% save data

for i = 1:length(p_values)
    cd(dataPath);
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
    savePath = ['Results\tTest\',stim(stim1),'-',stim(stim2),'\LP',num2str(LPk),'_HP',num2str(HPk),'\p-',num2str(p_values(i)),'\'];
    if ~exist(savePath,'dir')
        mkdir(savePath);
    end
    cd(savePath);
    
    imwrite(big_map_noClip_small,['tTest_p-',num2str(p_values(i)),'.bmp']);
    for f = 1:NFrames
        imwrite(rgb2gray(imfuse(greenGray,pVals(:,:,f))), [num2str(f), '.bmp']);
    end
end

end