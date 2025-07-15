function [] = ttest_ngc14_V1(varargin)

%% user input
username = 'ngc14';

animal = 'Gilligan';
date = '12_20_2018';
run = 'run01';

stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere', 'LargeSphere', 'Photocell'});

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
LPk=15;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=550;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter

% Hisashi's Clipping type
% 0: no clipping;
% 1: clipping at median+-(SD*clipvalue);
% 2: clipping at median+-(SD*clipvalue) with a mask (default.bmp);
% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix [x1, y1; x2, y2])
% 4: clipping at median+-intensity change (intensity change = clipvalue), Ex. 0.0005 or 0.0008;
% 5: clipping at 0+-(SD*clipvalue);
% 6: clipping at 0+-intensity change (intensity change = clipvalue);
clipMethod = 0;
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
dataPath = ['S:\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
parent = ['S:\Lab\',username,'\Scripts\'];
%%%%%%%%%%%%%%%%%%%
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
blockConds1 = blocks;

badAlign = load('bad_alignments.mat');
alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
badInds1 = ~ismember([1:size(blockConds1,1)],badAlign.badAlign{alignInd1});
blockConds1 = blockConds1(badInds1);

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
maps_stim1 = zeros(H,W,1);
maps_stim2 = zeros(H,W,1);
maps_avg = zeros(H,W,1);

frames_avg_stim1 = zeros(H,W,NBlocks);
frames_avg_baseline = zeros(H,W,NBlocks);

% dispstat('','init');
parfor n = 1:NBlocks
    % reset stim1 and stim2 to zeros
    frames_all_stim1 = zeros(H,W,NFrames);
    frames_all_baseline = zeros(H,W,1:10);
    
    % load frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frames_all_stim1(:,:,:) = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],alignInd1,'v');
    
    % median filter to reduce CCD column defects
    if ALIGN_OFFSET && FIX_COLUMN_DEFECT
        for f = 1:size(frames_all_stim1,3)
            frames_all_stim1(:,:,f) = medfilt2(frames_all_stim1(:,:,f),[5 5]);
            frames_all_baseline(:,:,f) = medfilt2(frames_all_baseline(:,:,f),[5 5]);
        end
    end
    
    % alignment
    if ALIGN_OFFSET
        align = load('alignment_offset.mat');
        
        for f = 1:size(frames_all_stim1,3)
            %             disp(['n',num2str(n),' f',num2str(f)]);
            frames_all_stim1(:,:,f) = OIShift(frames_all_stim1(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                -1*round(align.disp_y_med{alignInd1}{n,f}));
        end
    end
    
    % FFS
    frames_all_stim1 = (frames_all_stim1 - frames_all_stim1(:,:,fframe)) ./ frames_all_stim1(:,:,fframe);
    
    % average framerange
    frames_avg_stim1(:,:,n) = mean(frames_all_stim1(:,:,frameRange),3);
    frames_avg_baseline(:,:,n) = mean(frames_all_baseline(:,:,frameRange),3);
    
    % filter both stims for each block
    frames_avg_stim1(:,:,n) = imageFilter_LPHP_nsc15(frames_avg_stim1(:,:,n), LPk, HPk, clipMask);
    frames_avg_baseline(:,:,n) = imageFilter_LPHP_nsc15(frames_avg_baseline(:,:,n), LPk, HPk, clipMask);
    %     dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
    
end
% dispstat('Loading Data... 100% complete.');

subMap = mean(frames_avg_stim1,3)-mean(frames_avg_baseline,3);
subMapC = norm_to_uint8(OIClipH2(subMap,clipMethod,clipValue,clipMask));

[h, p] = ttest2(frames_avg_stim1,frames_avg_baseline,0.05,'left','equal',3);


%% save data
cd(dataPath);

stimName = ['[',stim(stim1),']-[',stim(stim2),']'];

savePath = ['Results\Ttest_nsc15_V3\Frames',num2str(frameRange(1)),...
    '-',num2str(frameRange(end)),'\LP',num2str(LPk),'HP',num2str(HPk),'\'];

if ALIGN_OFFSET
    savePath = [savePath,'Align_Offset\'];
end

if ~exist(savePath,'dir')
    mkdir(savePath);
end
cd(savePath);

save([stimName,'_pMap.mat'],'p');
% tform = load([dataPath,'\transformation_matrix.mat']);
% tform = tform.tform;
% if ~exist('Warped\', 'dir')
%     mkdir('Warped\');
% end
for i = 1:length(p_values)
    cleanedP = imfill(bwareaopen((p<p_values(i) & clipMask),10),'holes');
    imwrite(cleanedP,[stimName,'_p-',num2str(p_values(i)),'.bmp']);
    imwrite(ones(W,H,3).*cat(3,1,0,0),[stimName,'_p-',num2str(p_values(i)),'.png'],...
        'alpha',double(cleanedP));
end

imwrite(subMapC,gray(256),[stimName,'_subMap_C',num2str(clipValue),'.bmp']);
imwrite(norm_to_uint8(subMap),gray(256),[stimName,'_subMap_C0.bmp']);

cd(parent);

timeTotal = round(etime(clock,timeStart)/60*100)/100;
disp(['Script took ',num2str(timeTotal),' minutes total.']);

end