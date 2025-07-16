function [] = ttest_ngc14_V1(varargin)

%% user input
username = 'ngc14';

animal = 'Skipper';
date = '12_09_2020';
run = 'run00';

threshold = 10;
blockSelect = [];

ALIGN_OFFSET = true;
FIX_COLUMN_DEFECT = true;

regressVessel = 0;
numVesselMasks = 2;
regressorSmoothKernel = 5;
regressMean = 0;

p_values = [0.01 0.001 0.0001];

ALL_FRAMES_IND = true;
frameRange = [51:56]; %[16:17];
style=4;      % 1: simple DC, 2:ffsub (first frame subtraction), 3: subtraction (stim1 -stim2), 4: subtraction after ffsub
fframe=[1:4];   % first frame (only used for style 2&4)
stim1=[4];    % used in all styles, if more than one stim, will take average
stim2=[0];    % only used in style 3&4
stim = containers.Map([0 2 4 5],{'Blank','Extra Small Sphere','Large Sphere','Photocell'});

% Filtering
LPMethod='gaussian';        % Low-pass filtering method, usually 'gaussian'
LPk=15;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=250;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter

% Hisashi's Clipping type
% 0: no clipping;
% 1: clipping at median+-(SD*clipvalue);
% 2: clipping at median+-(SD*clipvalue) with a mask (default.bmp);
% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix [x1, y1; x2, y2])
% 4: clipping at median+-intensity change (intensity change = clipvalue), Ex. 0.0005 or 0.0008;
% 5: clipping at 0+-(SD*clipvalue);
% 6: clipping at 0+-intensity change (intensity change = clipvalue);
clipMethod = 2;
clipValue  = .75;
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
    
    % style
    if ~isempty(varargin{6})
        style = varargin{6};
    end
    
    % LPk
    if ~isempty(varargin{7})
        LPk = varargin{7};
    end
    
    % HPk
    if ~isempty(varargin{8})
        HPk = varargin{8};
    end
    
    % clipMethod
    if ~isempty(varargin{9})
        clipMethod = varargin{9};
    end
    
    % clipValue
    if ~isempty(varargin{10})
        clipValue = varargin{10};
    end
    
    % frameRange
    if ~isempty(varargin{11})
        frameRange = varargin{11};
    end
    
    if ~isempty(varargin{12})
        stim2 = varargin{12};
    end
    
    if ~isempty(varargin{13})
        stim = varargin{13};
    end
    
    clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
end

greenName = ['green',run(end-1:end),'_edited.bmp'];

if(strcmp(date, '11_09_2020') || strcmp(date, '12_09_2020'))
    stim = containers.Map([0 4 5],{'Blank', 'LargeSphere', 'Photocell'});
elseif(strcmp(date, '10_30_2020'))
    stim = containers.Map([0 2],{'Blank', 'Extra Small Sphere'});
end

%% Begin Script
timeStart = clock;
disp(['%%%%% TTEST_NGC14_V1: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' f[',num2str(frameRange),'] %%%%%']);


%% locate data
dataPath = ['S:\Lab\',animal,'\All Data\',animal,'_',date, '\Imaging\',run];
parent = ['S:\Lab\',username,'\Scripts\'];

try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end
if regressVessel
    vesselMaskRR = load('regMaskVessel_SS.mat');
    vesselMaskRR = vesselMaskRR.vesselMaskRR;
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

blockConds1 = blocks(conds==stim1);
alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
%blockConds2 = blocks(conds==stim2);
%alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));

badAlign = load('offsets.mat');
if ~isempty(blockSelect)
    blockConds1 = blockConds1(blockSelect);
    badAlign.dist{alignInd1} = badAlign.dist{alignInd1}(blockSelect,:);
    %blockConds2 = blockConds2(blockSelect);
    %badAlign.dist{alignInd2} = badAlign.dist{alignInd2}(blockSelect,:);
end

[row,col] =find(cellfun(@(a) nanmean(a)>threshold, badAlign.dist{alignInd1}));
row = row(~(ismember(col,11:20)));
badInds = unique(row);
badInds1 = ~ismember([1:size(blockConds1,1)], badInds);
% [row,col] =find(cellfun(@(a) nanmean(a)>threshold, badAlign.dist{alignInd2}));
% row = row(~(ismember(col,11:20)));
% badInds = unique(row);
% badInds2 = ~ismember([1:size(blockConds2,1)], badInds);

blockConds1 = blockConds1(badInds1,:);
%blockConds2 = blockConds2(badInds2,:);
align = load('alignment_offset.mat');
if ~isempty(blockSelect)
    align.disp_x_med{alignInd1} =  align.disp_x_med{alignInd1}(blockSelect,:);
    align.disp_y_med{alignInd1} =  align.disp_y_med{alignInd1}(blockSelect,:);
    %     align.disp_x_med{alignInd2} =  align.disp_x_med{alignInd2}(blockSelect,:);
    %     align.disp_y_med{alignInd2} =  align.disp_y_med{alignInd2}(blockSelect,:);
end

align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
% align.disp_x_med{alignInd2} = align.disp_x_med{alignInd2}(badInds2,:);
% align.disp_y_med{alignInd2} = align.disp_y_med{alignInd2}(badInds2,:);
shiftX1 = align.disp_x_med{alignInd1};
shiftY1 = align.disp_y_med{alignInd1};

% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = size(blockConds1,1);

% clip mask
if clipMethod==2
    try
        clipMask = imread(clipMaskName);
        clipMask = clipMask(:,:,1)>0;
        green = imread(greenName);
        greenGray = double(rgb2gray(green));
        if size(greenGray,1) == 1082
            greenGray = greenGray(2:end-1,:);
            clipMask = clipMask(2:end-1,:);
        end
        if size(greenGray,2) == 1312
            greenGray = greenGray(:,3:end-2);
            clipMask = clipMask(:,3:end-2);
        end
        clipMask = ~(~clipMask | imgaussfilt(greenGray,1) < median(greenGray(:))-0.5*std(greenGray(:)));
    catch error
        error('Cannot load clip mask.');
    end
    if size(clipMask,1) == 1082
        clipMask = clipMask(2:end-1,:);
    end
    if size(clipMask,2) == 1312
        clipMask = clipMask(:,3:end-2);
    end
else
    clipMask = logical(ones(W,H));
end
clear anapar blocks badAlign badInds1 clipMaskName conds greenGray green greenName indx align
%% load data
frames_all_stim1N = zeros(H,W,NFrames,NBlocks);
% dispstat('','init');
for n = 1:NBlocks
    % reset stim1 and stim2 to zeros
    %   frames_all_stim1 = zeros(H,W,NFrames);
    %    frames_all_stim2 = zeros(H,W,NFrames);
    
    % load frames
    frames_all_stim1 = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],0,'v');
    %    frames_all_stim2 = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],0,'v');
    shiftX = [shiftX1{n,:}];
    shiftY = [shiftY1{n,:}];
    % median filter to reduce CCD column defects
    parfor f = 1:NFrames
        if  FIX_COLUMN_DEFECT
            frames_all_stim1(:,:,f) = medfilt2(frames_all_stim1(:,:,f),[5 5]);
            %            frames_all_stim2(:,:,f) = medfilt2(frames_all_stim2(:,:,f),[5 5]);
        end
        % alignment
        if ALIGN_OFFSET
            frames_all_stim1(:,:,f) = OIShift(frames_all_stim1(:,:,f),-1*round(shiftX(f)),...
                -1*round(shiftY(f)));
            %frames_all_stim2(:,:,f) = OIShift(frames_all_stim2(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
            %    -1*round(align.disp_y_med{alignInd2}{n,f}));
        end
    end
    
    %frames_all_stim1(:,:,f) = bspline_transform(tform.O_trans,frames_all_stim1(:,:,f),tform.Spacing,3);
    %frames_all_stim2(:,:,f) = bspline_transform(tform.O_trans,frames_all_stim2(:,:,f),tform.Spacing,3);
    %             if regressVessel
    %                 tempFrame = frames_all_stim1(:,:,f);
    %                 for j = 1:numVesselMasks
    %                      vesselSig(f,j) = mean(tempFrame(vesselMaskRR(:,:,j)));
    %                 end
    %             end
    
    % FFS
    frames_all_stim1N(:,:,:,n) = (frames_all_stim1 - mean(frames_all_stim1(:,:,fframe),3)) ./ mean(frames_all_stim1(:,:,fframe),3);
end


    %frames_all_stim2 = (frames_all_stim2 - frames_all_stim2(:,:,fframe)) ./ frames_all_stim2(:,:,fframe);
    
    %     if regressVessel || regressMean
    %         if regressVessel, vesselSig = movmedian(vesselSig,regressorSmoothKernel); end
    %
    %         % reshape 3d image matricies into 2d ones
    %         frames_reshape = reshape(frames_all_stim1,[size(frames_all_stim1,1)*size(frames_all_stim1,2) size(frames_all_stim1,3)])';
    %
    %         % construct regression matrix X
    %         X = [ones(size(frames_reshape,1),1)];
    %
    %         if regressVessel
    %             for j = 1:numVesselMasks
    %                 X = [X, vesselSig(:,j)];
    %             end
    %         end
    %
    %         if regressMean
    %             frames_reshape(find(isinf(frames_reshape))) = NaN;
    %             X = [X, nanmean(frames_reshape,2)];
    %         end
    %
    %         %     perform regression
    %         parfor i = 1:size(frames_reshape,2)
    %             b(i,:) = regress(frames_reshape(:,i),X);
    %             maps_avg_reshape_reg(:,i) = frames_reshape(:,i)-X*b(i,:)';
    %         end
    %
    %         % reshape back into 3d image
    %         frames_all_stim1 = reshape(maps_avg_reshape_reg',[size(frames_all_stim1,1) size(frames_all_stim1,2) size(frames_all_stim1,3)]);
    %     end
    
    %         if(stim2==0)
    %             frames_avg_stim2(:,:,n) = nanmean(frames_all_stim2(:,:,1:NFrames),3);
    %         else
    %             frames_avg_stim2(:,:,n) = nanmean(frames_all_stim2(:,:,frameRange),3);
    %         end
    
frames_avg_stim2 = zeros(H,W,NBlocks);
clear frames_all_stim1 blockConds align;
%%
for frameRanges = 1:1:70
    frameRange = max(1,frameRanges-2):min(70,frameRanges+2);
    frames_avg_stim1(:,:,1:NBlocks) = squeeze(nanmean(frames_all_stim1N(:,:,frameRange,:),3));
    if(frameRanges==1)
        frames_avg_stim2(:,:,1:NBlocks) = squeeze(nanmean(frames_all_stim1N(:,:,setdiff(1:10,fframe),:),3));
    end
    
    parfor n = 1:NBlocks
        % filter both stims for each block
        frames_avg_stim1(:,:,n) = imageFilter_LPHP_nsc15(frames_avg_stim1(:,:,n), LPk, HPk, []);
        if(frameRanges==1)
            frames_avg_stim2(:,:,n) = imageFilter_LPHP_nsc15(frames_avg_stim2(:,:,n), LPk, HPk, []);
        end
        %     dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
    end
    
    % dispstat('Loading Data... 100% complete.');
    
    subMap = nanmean(frames_avg_stim1,3)-nanmean(frames_avg_stim2,3);
    subMapC = norm_to_uint8(OIClipH2(subMap,clipMethod,clipValue,clipMask));
    
    [h, p] = ttest2(frames_avg_stim1,frames_avg_stim2,0.05,'left','equal',3);
    
    
    %% save data
    cd(dataPath);
    
    stimName = ['[',stim(stim1),']-[BASELINE]'];
    
    savePath = ['Results\Ttest_nsc15_V3\Frames',num2str(frameRange(1)),...
        '-',num2str(frameRange(end)),'\LP',num2str(LPk),'HP',num2str(HPk),'\'];
    
    if ALIGN_OFFSET
        savePath = [savePath,'Align_Offset\FF_1\'];
    end
    
    if regressVessel
        savePath = [savePath,'Regress\'];
    end
    if regressMean
        savePath = [savePath,'Regress_Mean\'];
    end
    if ~isempty(blockSelect)
        savePath = [savePath,'Blocks_',num2str(blockSelect(1)),'_',num2str(blockSelect(end)),'\'];
    end
    if ~exist(savePath,'dir')
        mkdir(savePath);
    end
    cd(savePath);
    
    save([stimName,'_pMap.mat'],'p');
    
    for i = 1:length(p_values)
        cleanedP = imfill(bwareaopen((p<p_values(i) & clipMask),10),'holes');
        imwrite(cleanedP,[stimName,'_p-',num2str(p_values(i)),'.bmp']);
        %imwrite(imwarp(cleanedP,affine2d(tform),'OutputView',imref2d(size(cleanedP))),...
    end
    
    imwrite(subMapC,gray(256),[stimName,'_subMap_C',num2str(clipValue),'.bmp']);
    imwrite(norm_to_uint8(subMap),gray(256),[stimName,'_subMap_C0.bmp']);
end
cd(parent);

timeTotal = round(etime(clock,timeStart)/60*100)/100;
disp(['Script took ',num2str(timeTotal),' minutes total.']);

end