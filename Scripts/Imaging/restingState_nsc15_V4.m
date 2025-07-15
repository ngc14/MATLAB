% function [] = restingState_nsc15_V4()

% close all;
clearvars;

%% user input
username = 'ngc14';                 % windows username

datasetNum = 5;

% dataset details
animal = {'Merele','Merele','Merele','Bordeaux','Gilligan'};                  % animal name
animalType = {'SqM','SqM','SqM','SqM','Macaque'};
hemi = {'Left','Left','Left','Left','Right'};
date = {'11_28_2017','01_22_2018','05_07_2018','09_24_2018','11_20_2018'};             % experiment date
run = {'run01','run01','run16','run00','run01'};                     % run name
fps = {10,10,2,10,10};    %frames per second
blockSelect = {[1:22],[1:51],[37:75],[1:40],[1:50]};

SPATIAL_BIN = 2;      %factor to scale resolution down by
TEMPORAL_BIN = 5;     %factor to scale time down by

if datasetNum==1
    seed = [144 212]; %dataset 1
    
elseif datasetNum==2
    seed = [162 196]; %dataset 2
    
elseif datasetNum==3
    seed = [206 218]; %dataset 3
    
elseif datasetNum==4
    seed = [200 200];
    
elseif datasetNum==5
    seed = [328 235];
end

bandPass = [0.03 0.2]; % butterworth filter for signal timecourse
buttOrder = 3;

ALIGN_OFFSET = true;

corrType = 'pearson'; %'pearson','partial'

HPk = 0;                          % High pass filter kernel
HPMethod = 'fastmedian';            % High pass filter method
LPk = 3;                           % Low pass filter kernel
LPMethod = 'gaussian';              % Low pass filter method

%% begin
disp(['%%%%%  RESTING STATE V4: ',animal{datasetNum},' ',date{datasetNum}(date{datasetNum}~='_'),...
    ' ',run{datasetNum},'  %%%%%']);
timeStart = clock;

dataPath = ['\\pitt\sni\Gharbawie\Lab\',animal{datasetNum},'\All Data\',animal{datasetNum},'_',date{datasetNum},'\Imaging\',run{datasetNum}];
parent = ['C:\Users\',username,'\Documents\Data\'];


%% load data info
cd(dataPath);
green = imread(['green',run{datasetNum}(end-1:end),'_edited.bmp']);
red = imread(['red',run{datasetNum}(end-1:end),'.bmp']);
clipMask = imread(['clipMask',run{datasetNum}(end-1:end),'.bmp']);
clipMask = clipMask(:,:,1)>0;

greenGray = double(rgb2gray(green));
clipMask = ~(~clipMask | imgaussfilt(greenGray,1) < median(greenGray(:))-0.5*std(greenGray(:)));

%% sort files
filenames = dir('Data*.BLK');
for f = 1:length(filenames)
    ind1 = strfind(filenames(f).name,'E');
    ind2 = strfind(filenames(f).name,'B'); ind2 = ind2(1);
    ind3 = strfind(filenames(f).name,'.');
    fn_exp_blk(f,1) = str2num(filenames(f).name(ind1+1:ind2-1));
    fn_exp_blk(f,2) = str2num(filenames(f).name(ind2+1:ind3-1));
end

[fn_exp_blk, sortInd] = sortrows(fn_exp_blk,[1 2]);
filenames = filenames(sortInd);

filenames = filenames(blockSelect{datasetNum});
anpar = OIHeadRead(filenames(1).name,'v');


%% load blk files, downsample to 2Hz
% frames = zeros(anpar.FrameHeight,anpar.FrameWidth,anpar.FramesPerStim*length(blockSelect{datasetNum}));
dispstat('','init');
k = 1;
for n = 1:length(filenames)
    framesTemp = OIReadStim([filenames(n).folder,'\',filenames(n).name], 1, 'v');
%     frames(:,:,((n-1)*anpar.FramesPerStim+1):(n*anpar.FramesPerStim)) = framesTemp;

    % alignment
    if ALIGN_OFFSET
        if n==1
            align = load('alignment_offset.mat');
            outputView = imref2d(size(framesTemp(:,:,1)));
        end
        
        for f = 1:size(framesTemp,3)
%             disp(['n',num2str(n),' f',num2str(f)]);
%             framesTemp(:,:,f) = OIShift(framesTemp(:,:,f),-1*round(align.disp_x_med(n,f)),...
%                 -1*round(align.disp_y_med(n,f)));
            
            framesTemp(:,:,f) = imwarp(framesTemp(:,:,f),align.all_tforms{n,f},'outputView',outputView);
        end
    end
    
    for i = 1:(size(framesTemp,3)/TEMPORAL_BIN)
        frames_2Hz_temp = mean(framesTemp(:,:,((i-1)*TEMPORAL_BIN+1):(i*TEMPORAL_BIN)),3);
        frames_2Hz(:,:,k) = imresize(frames_2Hz_temp,1/SPATIAL_BIN);
        k = k + 1;
    end

    dispstat(['LOADING DATA... ',num2str(round(n/length(filenames)*100*100)/100),'% Done']);
end
cd(parent);
dispstat('LOADING DATA... 100% Done');
clear framesTemp;


%% image filtering
if HPk>0 || LPk>0
    dispstat('','init');
    for f = 1:size(frames_2Hz,3)
        frames_2Hz(:,:,f) = OIEasyFilterH2(frames_2Hz(:,:,f), LPMethod, LPk, HPMethod, HPk);

        if rem(f-1,10)==0
            dispstat(['FILTERING IMAGES... ',num2str(round(f/size(frames_2Hz,3)*100*100)/100),'% Done']);
        end
    end
    dispstat('FILTERING IMAGES... 100% Done');
end


%% timecourse filtering
[b1,a1] = butter(buttOrder, [bandPass(1) bandPass(2)]/(2/2), 'bandpass');
frames_2Hz_filt = zeros(size(frames_2Hz));
% frames_filt = zeros(size(frames));

dispstat('','init');
for x = 1:size(frames_2Hz,1)
    for y = 1:size(frames_2Hz,2)
        frames_2Hz_filt(x,y,:) = filtfilt(b1,a1,frames_2Hz(x,y,:));
%         frames_filt(x,y,:) = filtfilt(b1,a1,frames(x,y,:));
    end

    if rem(x-1,10)==0
        dispstat(['BANDPASS FILTERING DATA... ',num2str(round(x/size(frames_2Hz,1)*100*100)/100),'% Done']);
    end
end
dispstat('BANDPASS FILTERING DATA... 100% Done.');


%% calculate seed signal (gaussian average of a radius around chosen point)
% rows are observations, columns are variables
seed = round([546 556] ./ SPATIAL_BIN);
seedRad = 3; % 6 = 250um radis, 12 = 500um radius, 24 = 1000um radius
gaussianDisk = fspecial('gaussian',2*seedRad+1,seedRad);
green_seed = green(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad,:);
clipMask_seed = clipMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
gaussianDisk = gaussianDisk .* double(clipMask_seed);
gaussianDisk = gaussianDisk ./ sum(gaussianDisk(:));

seedSig = zeros(size(frames_2Hz_filt,3),1);
seedSig_unfilt = zeros(size(frames_2Hz,3),1);
k = 1;
for x = 1:2*seedRad+1
    for y = 1:2*seedRad+1    
        seedSig = seedSig + gaussianDisk(x,y) * squeeze(frames_2Hz_filt(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));
        seedSig_unfilt = seedSig_unfilt + gaussianDisk(x,y) * squeeze(frames_2Hz(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));
    end
end

%% calculate pearson correlation
frames_filt_reshape = reshape(frames_2Hz_filt,...
    [size(frames_2Hz_filt,1)*size(frames_2Hz_filt,2) size(frames_2Hz_filt,3)])';

c_filt = corr(seedSig,frames_filt_reshape);
c_filt = reshape(c_filt,[size(frames_2Hz_filt,1) size(frames_2Hz_filt,2)]);

%% plot correlation map
figure; imagesc(c_filt); 
hold on;
plot(seed(1),seed(2),'k.','markersize',30);
plot(seed(1),seed(2),'w.','markersize',23);
axis image off; 
colormap jet; 
colorbar; 
caxis([0 1]);
title([animal{datasetNum},' ',date{datasetNum}(date{datasetNum}~='_'),...
    ' ',run{datasetNum},' - Pearson corr map - ALIGN:',num2str(ALIGN_OFFSET)]);

%% movies
if 1==0
    % unfiltered movie
    figure('color',[1 1 1]);
    for f = 1:1:size(frames_2Hz,3)
        imagesc(frames_2Hz(:,:,f));
        axis image off; colormap gray;
        caxis([0 4E4]);
        title([num2str(0.1*f),' seconds']);
        pause(0.02);
    end
    
    % filtered movie
    figure('color',[1 1 1]);
    for f = 1:1:size(frames_2Hz_filt,3)
        imagesc(frames_2Hz_filt(:,:,f));
        axis image off; colormap gray;
        caxis([-750 750]);
        title([num2str(0.1*f),' seconds']);
        pause(0.05);
    end
    
     % both
    figure('color',[1 1 1]);
    for f = 1:1:size(frames_2Hz_filt,3)
        subplot(1,2,1);
        imagesc(frames_2Hz(:,:,f));
        axis image off; colormap gray;
        caxis([0 4E4]);
        title([num2str(0.1*f),' seconds']);
        
        subplot(1,2,2);
        imagesc(frames_2Hz_filt(:,:,f));
        axis image off; colormap gray;
        caxis([-1000 1000]);
        title([num2str(0.1*f),' seconds']);
        pause(0.1);
    end
end
