function [] = BlockView_ngc14_V1(varargin)
%% user input
animal = 'Skipper';
date = '11_24_2020';
run = 'run00';
threshold = 10;
stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere' ,'LargeSphere', 'Photocell'});
frameRange = [50:55]; %[16:17];
fframe=[1:10];   % first frame (only used for style 2&4)
stim1=[5];    % used in all styles, if more than one stim, will take average
LPk=5;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPk=250;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter

%% multirun
if ~isempty(varargin)
    % animal
    if ~isempty(varargin{1})
        animal = varargin{1};
    end
    % date
    if ~isempty(varargin{2})
        date = varargin{2};
    end
    % run
    if ~isempty(varargin{3})
        run = varargin{3};
    end
    % stim #
    if ~isempty(varargin{4})
        stim1 = varargin{4};
    end
    % stim name
    if ~isempty(varargin{5})
        stim = varargin{5};
    end
    % LPk
    if ~isempty(varargin{6})
        LPk = varargin{6};
    end
    % HPk
    if ~isempty(varargin{7})
        HPk = varargin{7};
    end
    
end
dataPath = ['S:\Lab\',animal,'\All Data\',animal,'_',date, '\Imaging\',run];
savePath = [dataPath,'\Results\Blockview\','[',stim(stim1),']\FF_1\','LP',...
    num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
if(~exist(savePath,'dir'))
    mkdir(savePath)
end
%% Begin Script
disp(['%%%%% BLOCKVIEW_NGC14: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
    num2str(LPk),' HP',num2str(HPk),' f[',num2str(frameRange),'] %%%%%']);
%% locate data
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
blockConds1 = blocks(conds==stim1);
alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
badAlign = load('offsets.mat');

[row,col] =find(cellfun(@(a) nanmean(a)>threshold, badAlign.dist{alignInd1}));
row = row(~(ismember(col,11:20)));
badInds = unique(row);
badInds1 = ~ismember([1:size(blockConds1,1)], badInds);
blockConds1 = blockConds1(badInds1);
% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NBlocks = size(blockConds1,1);
%% load data
align = load('alignment_offset.mat');
align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
alignX = align.disp_x_med{alignInd1};
alignY = align.disp_y_med{alignInd1};

% allTrials = repmat({zeros(H,W,NFrames)},1,NBlocks);
fileName2 = matfile([savePath,'trialFrames.mat'],'Writable', true);

for n = 1:NBlocks
    % reset stim1 and stim2 to zeros
    %frames_stim1 = zeros(H,W,NFrames);
    %framesFilt = zeros(H,W,NFrames);
    
    % load stim 1
    frames_stim1(:,:,n) = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],0,'v',ff);
    % median filter to reduce CCD column defects
    for f = 1:size(frames_stim1,3)
        frames_stim1(:,:,f) = medfilt2(frames_stim1(:,:,f),[5 5]);
        frames_stim1(:,:,f) = OIShift(frames_stim1(:,:,f),-1*round(alignX{n,f}),...
            -1*round(alignY{n,f}));
    end
    frames_stim1 = 100*(frames_stim1 - mean(frames_stim1(:,:,fframe),3)) ./...
        mean(frames_stim1(:,:,fframe),3);
    %     parfor f = 1:size(framesFFS,3)
    %         framesFilt(:,:,f) = imageFilter_LPHP_nsc15(framesFFS(:,:,f),LPk,HPk,[]);
    %     end
    frames_stim1(isinf(frames_stim1)) = 0;
    frames_stim1(isnan(frames_stim1)) = 0;
    
    subsref(fileName2.allTrials(1,n)= frames_stim1(:,:,ff);
    disp(['Trial ', num2str(n), ' of ', num2str(NBlocks)]);
end
disp('Par:Done');
% fileName = [savePath,'AllTrials.mat'];
% save(fileName, 'allTrials', '-v7.3');
% allTrials = num2cell(cat(4,allTrials{:}),[1 2 4]);
% save(fileName2, 'allTrials', '-v7.3');
end
