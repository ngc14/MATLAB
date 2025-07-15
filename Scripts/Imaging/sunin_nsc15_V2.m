function [] = sunin_nsc15_V2(varargin)

%% user input
username = 'nsc15';                 % windows username

animal = 'Merele';                  % animal name
hemi = 'Left';
date = '03_21_2017';                % experiment date
run = 'run02';

stim(1,1)={'Blank'};
% stim(2,1)={'1500ms_Blue_LED'};
stim(2,1)={'150BiPulses_40microAmps'};
stim(3,1)={'150BiPulses_60microAmps'};
stim(4,1)={'150BiPulses_80microAmps'};

frameRange = [7:9];
style=4;      % 1: simple DC, 2:ffsub (first frame subtraction), 3: subtraction (stim1 -stim2), 4: subtraction after ffsub
fframe=[1];   % first frame (only used for style 2&4)
stim1=[3];    % used in all styles, if more than one stim, will take average
stim2=[1];    % only used in style 3&4

ALIGN_OFFSET = true;

% Filtering
LPMethod='gaussian';        % Low-pass filtering method, usually 'gaussian'
LPk=0;                     % Low-pass kernel size (usually small) set as 0 if no LP filter
HPMethod='fastmedian';      % High-pass filtering method, usually 'fastmedian'
HPk=0;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter
FILTER_EACH_TRIAL = false;  % this will take forever.

% Hisashi's Clipping type
    % 0: no clipping;
    % 1: clipping at median+-(SD*clipvalue);
    % 2: clipping at median+-(SD*clipvalue) with a mask (default.bmp);
    % 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix [x1, y1; x2, y2])
    % 4: clipping at median+-intensity change (intensity change = clipvalue), Ex. 0.0005 or 0.0008;
    % 5: clipping at 0+-(SD*clipvalue);
    % 6: clipping at 0+-intensity change (intensity change = clipvalue);
clipMethod = 2;
clipValue  = 2;
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
        stim{stim1} = varargin{6};
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
    
    clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
end



%% Begin Script
disp(['%%%%% SUNIN_NSC15_V2: ',animal,' ',date(date~='_'),' ',run,' [',stim{stim1},'] LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' f[',num2str(frameRange),'] %%%%%']);


%% locate data
dataPath = ['C:\Users\',username,'\Documents\Data\',animal,'_SqM\',hemi,'Hemisphere\',date,'\',run];
parent = ['C:\Users\',username,'\Documents\Data'];

try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end

blocks = dir('*.BLK');

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
    clipMask = logical(ones(W,H));
end

%% load data
maps_stim1 = zeros(H,W,1);
maps_stim2 = zeros(H,W,1);
maps_avg = zeros(H,W,1);

% dispstat('','init');
parfor n = 1:NBlocks
    if ALIGN_OFFSET
        align = load('alignment_offset.mat');
    end
    
    % reset stim1 and stim2 to zeros
    frames_stim1 = zeros(H,W,NFrames);
    frames_stim2 = zeros(H,W,NFrames);
    
    % load stim 1
    frames_stim1(:,:,:) = OIReadStim([blocks(n).folder,'\',blocks(n).name],stim1,'v');
    
    %if blank sub, load stim 2
    if style==3 || style==4 
        frames_stim2(:,:,:) = OIReadStim([blocks(n).folder,'\',blocks(n).name],stim2,'v');
    end
    
%     % median filter to reduce CCD column defects
%     for f = 1:size(frames_stim1,3)
%         frames_stim1(:,:,f) = medfilt2(frames_stim1(:,:,f),[3 9]);
%         
%         if style==3 || style==4 
%             frames_stim2(:,:,f) = medfilt2(frames_stim2(:,:,f),[3 9]);
%         end
%     end

    if ALIGN_OFFSET
        for f = 1:size(frames_stim1,3)
%             disp(['n',num2str(n),' f',num2str(f)]);
            frames_stim1(:,:,f) = OIShift(frames_stim1(:,:,f),-1*round(align.disp_x_med{stim1}(n,f)),...
                -1*round(align.disp_y_med{stim1}(n,f)));
            
            if style==3 || style==4 
                frames_stim2(:,:,f) = OIShift(frames_stim2(:,:,f),-1*round(align.disp_x_med{stim2}(n,f)),...
                    -1*round(align.disp_y_med{stim2}(n,f)));
            end
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
            maps_stim1 = maps_stim1 + mean(frames_stim1(:,:,frameRange),3);
            
        case 2 % single condition w/ FFS
            maps_stim1 = maps_stim1 + (mean(frames_stim1(:,:,frameRange),3) -...
                frames_stim1(:,:,fframe)) ./ frames_stim1(:,:,fframe);
            
        case 3 % blank subtraction
            maps_stim1 = maps_stim1 + mean(frames_stim1(:,:,frameRange),3);
            maps_stim2 = maps_stim2 + mean(frames_stim2(:,:,frameRange),3);
            
        case 4 % FFS and then blank subtraction
            maps_stim1 = maps_stim1 + (mean(frames_stim1(:,:,frameRange),3) -...
                frames_stim1(:,:,fframe)) ./ frames_stim1(:,:,fframe);
            maps_stim2 = maps_stim2 + (mean(frames_stim2(:,:,frameRange),3) -...
                frames_stim2(:,:,fframe)) ./ frames_stim2(:,:,fframe);
    end
    
%     dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);

end
% dispstat('Loading Data... 100% complete.');


%% average maps and blank subtraction
switch style
    case {1, 2} % single condition
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_avg = maps_stim1_avg;

    case {3, 4} % blank subtraction
        maps_stim1_avg = maps_stim1 ./ NBlocks;
        maps_stim2_avg = maps_stim2 ./ NBlocks;
        maps_avg = maps_stim1_avg - maps_stim2_avg;
end


%% filtering and clipping

if (style==2 || style==4)
    maps_avg_F = imageFilter_LPHP_nsc15(maps_avg,LPk,HPk,clipMask);
    maps_avg_FC = OIClipH(maps_avg_F, clipMethod, clipValue, clipMask);
    
    % just some testing
    %{
    maps_avg_C = OIClipH(maps_avg, clipMethod, clipValue, clipMask);
    
    maps_avg_F2 = imageFilter_LPHP_nsc15(maps_avg,LPk,0,clipMask);
    maps_avg_FC2 = OIClipH(maps_avg_F2, clipMethod, clipValue, clipMask);
    %}
        
else
    maps_avg_FC = OIClipH(maps_avg_F, clipMethod, clipValue, clipMask);
end


% testing distance correction
%{ 
seed = [310 440];
clear distMat;
for x = 1:size(maps_avg_F,1)
    for y = 1:size(maps_avg_F,2)
        distMat(x,y) = sqrt((seed(2)-x)^2 + (seed(1)-y)^2) / (74);
    end
end

maps_avg_DC = maps_avg + maps_avg .* distMat;
maps_avg_DC(maps_avg>0) = maps_avg(maps_avg>0);
maps_avg_DC_F = imageFilter_LPHP_nsc15(maps_avg_DC,LPk,HPk,clipMask);
maps_avg_DC_FC = OIClipH(maps_avg_DC_F, clipMethod, clipValue, clipMask);
%}


%% save data
cd(dataPath);

switch style
    case 1
        stimName = ['[',stim{stim1},']'];
        ffName = 'NoFF';
    case 2
        stimName = ['[',stim{stim1},']'];
        ffName = ['FF_',num2str(fframe)];
    case 3
        stimName = ['[',stim{stim1},']-[',stim{stim2},']'];
        ffName = 'NoFF';
    case 4
        stimName = ['[',stim{stim1},']-[',stim{stim2},']'];
        ffName = ['FF_',num2str(fframe)];
end

savePath = ['Results\Sunin_nsc15_V2\',stimName,'\',ffName,'\',...
    'LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue),'\'];

if ALIGN_OFFSET
    savePath = [savePath,'Align_Offset\'];
end

if ~exist(savePath,'dir')
    mkdir(savePath);
end
cd(savePath);

% uint8 images
imwrite(norm_to_uint8(maps_avg_FC),gray(256),['sunin_F',num2str(frameRange(1)),'-',...
    num2str(frameRange(end)),'.bmp']);
imwrite(norm_to_uint8(maps_avg),gray(256),['sunin_F',num2str(frameRange(1)),'-',...
    num2str(frameRange(end)),'_unFilt.bmp']);

% double .mat files
save(['sunin_F',num2str(frameRange(1)),'-',num2str(frameRange(end)),'.mat'],'maps_avg_FC');
save(['sunin_F',num2str(frameRange(1)),'-',num2str(frameRange(end)),'_unFilt.mat'],'maps_avg');


cd(parent);

end