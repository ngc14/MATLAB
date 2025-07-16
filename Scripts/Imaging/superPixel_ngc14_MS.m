 %function [] = superPixel_nsc15()
clearvars;
% close all;

%% user input
username = 'ngc14';

monkey = 'Gilligan';                   % animal name
hemi = 'Left';
date = '04_18_2018';                % experiment date
run = 'run01';                      % run name
blockSelect = [];                     % leave empty for all blocks
excludeFrames = []; %[5:12];

ImgHz = 10;

FILTER = true;                % HP/LP filter and clip?
HPk = 425;                          % High pass filter kernel
HPMethod = 'fastmedian';            % High pass filter method
LPk = 4;                           % Low pass filter kernel
LPMethod = 'gaussian';              % Low pass filter method

CLIP = false;
clipMethod = 2;                     % Clip Method
clipValue = 1;                         % Clip value (SDs)
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
greenName = ['green',run(end-1:end),'_edited.bmp'];

frameLimit = 100;
blockLimit = 100;

stimSelect = [1];                   % stim #s to do superpixel for
FFS = true;
fframe = 1;
SP_INDI_BLKS = true;               % do superpixel for each block (takes a lot longer), or the averaged blocks?
USE_MASKS = true;                  % if true, uses masks instead of SP at each pixel. masks go in run/SPMasks folder
GRID_MASKS = false;
GRID_MASKS_INVERT = false;
LPFILTER = false;
NORMALIZE = false;                   % if true, normalize masks to have a peak at -1

%stims = containers.Map([0 2 3 4],{'Blank', 'D2', 'D3', 'D4'});
%stims = containers.Map([0 2 3 6 7],{'Blank', 'DThenar_75', 'pD3', 'DThenar_50', 'DThenar_25'});
stims = containers.Map([2 1],{'Rest','Large Sphere'});
SAVE = false;                        % Save data?
EMAIL_UPON_COMPLETION = false;


%% begin script
c0 = clock;

if ~FILTER
    LPk = 0;
    HPk = 0;
end

if ~CLIP
    clipValue = 0;
end

% display data that will be analyzed
disp(['%%%%% SUPERPIXEL: ',monkey,' ',date(date~='_'),' ',run,' LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' %%%%%']);

% data path
path = ['\\pitt\sni\gharbawie\Lab\',monkey,'\All Data\', monkey, '_',date, '\Imaging\',run];
path = ['\\univ.pitt.edu\sni\Gharbawie\Lab\ngc14\MsHowell\Left_Chamber\2013_04_23_MsHowell\run1'];
parent = ['\\pitt\sni\gharbawie\Lab\',username,'\Scripts\'];
% load(['C:\Users\nsc15\Documents\Data\',animal,'_SqM\',hemi,'Hemisphere\dateRunSite_',animal,'.mat']);

% get additional data properties
cd(path);

% block identification
blocks = dir('*.BLK');


if isempty(blocks) %look on backup drive instead
    disp('Data not found on SSD. Attempting to load from HDD.');
    
    try
        cd(['D:\Data\',monkey,'_SqM\',hemi,'Hemisphere\',date,'\',run]);
    catch error
        error('HDD data path does not exist. Exiting.');
    end
    
    blocks = dir('*.BLK');
    cd(path);
end
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);
conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.NStim];  %12_20_2019 Omar changed anapar.Cond to anapar.NStim
end
if ~isempty(blockSelect)
    blocks = blocks(blockSelect);
end
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;


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

if NFrames > frameLimit
    NFrames = frameLimit;
end

xtime = ((1:NFrames)-1)*(1/ImgHz);

% green = imread(greenName);
% greenGray = double(rgb2gray(green));
% clipMask = ~(~clipMask | imgaussfilt(greenGray,1) < median(greenGray(:))-0.5*std(greenGray(:)));


% load masks if necessary
if USE_MASKS
%     tform = load('transformation_matrix.mat');
%     tform = inv(tform.tform);
%     tform(:,3) = [0 0 1];
    if ~GRID_MASKS
        %cd('\\pitt\sni\Gharbawie\Lab\Gilligan\Mapping\SPMasks');
        cd('spmask');
    else
        if ~GRID_MASKS_INVERT
            cd([path,'\SPMasks\gridMasks']);
        else
            cd([path,'\SPMasks\gridMasks_InvertTtest']);
        end
    end
    
    maskFiles = dir('mask*.bmp');
    fileNames = {maskFiles.name};
    undInds = cellfun(@(a) regexp(a, '_'), fileNames, 'UniformOutput', false);
    domains = 1;
    for i = 1:length(maskFiles)
        %mask{i} = imwarp(imread(maskFiles(i).name), affine2d(tform),'OutputView',imref2d([W,H]));
        mask{i} = imread(maskFiles(i).name);
        mask{i} = mask{i}(:,:,1)>0;
    end
    cd(path);
    
    if GRID_MASKS
        DRS_ind = find(dateRunSite(:,1) == str2num(date(date~='_')) &...
            dateRunSite(:,2) == str2num(run(end-1:end)));
        DRS_xy = [dateRunSite(DRS_ind,4), 768-dateRunSite(DRS_ind,5)];
        
        for m = 1:length(mask)
            [xm, ym] = find(mask{m}>0);
            mask_xy_avg(m,:) = [mean(ym) 768-mean(xm)];
            maskDist(m,1) = sqrt((mask_xy_avg(m,1) - DRS_xy(1))^2 + (mask_xy_avg(m,2) - DRS_xy(2))^2);
        end
        
        [maskDist, maskSort] = sort(maskDist);
        mask = mask(maskSort);
        mask_xy_avg = mask_xy_avg(maskSort,:);
    end
end

% colors = hsv(length(mask));
% figure;
% for i=1:length(mask)
%     plot(mask_xy_avg(i,1), mask_xy_avg(i,2), 'o', 'color', colors(i,:));
%     hold on;
% end
% plot(DRS_xy(1),DRS_xy(2),'ko')
% ylim([0 768])
% xlim([0 768])
% axis square;
%
% figure;
% for i = 1:length(mask)
%     if i==1
%         imagesc(mask{i});
%         hold on;
%     else
%         imagesc(mask{i},'alphadata',mask{i});
%     end
% end
% axis image off;
% colormap(gray);
badAlign = load('bad_alignments.mat');
align = load('alignment_offset.mat');
%% begin image loading/filtering loop
for s = 1:length(stimSelect)
%     blankBlocks = blocks(conds==0);
%     blocksCond = blocks(conds==stimSelect(s));
    blankBlocks = blocks;
    blocksCond = blocks;
    alignInd2 = find(cellfun(@(a) a==stimSelect(s), stims.keys));
    badInds2 = ~ismember([1:size(blocksCond,1)],badAlign.badAlign{alignInd2});
    blocksCond = blocksCond(badInds2);
    align.disp_x_med{alignInd2} = align.disp_x_med{alignInd2}(badInds2,:);
    align.disp_y_med{alignInd2} = align.disp_y_med{alignInd2}(badInds2,:);
    
    if(s==1)
        alignInd1 = find(cellfun(@(a) a==2, stims.keys));
        badInds1 = ~ismember([1:size(blankBlocks,1)],badAlign.badAlign{alignInd1});
        blankBlocks = blankBlocks(badInds1);
        align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
        align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
    end
    
    NBlocks = min(size(blocksCond,1), size(blankBlocks,1));
    if NBlocks > blockLimit
        NBlocks = blockLimit;
    end
    %% preallocating
    if SP_INDI_BLKS
        framesFFS_BS_FILT = cell(length(stimSelect),1);
        
        framesFFS_BS_FILT{s} = cell(NBlocks,1);
        
        framesFFS_B_FILT = cell(length(stimSelect),1);
        
        framesFFS_B_FILT{s} = cell(NBlocks,1);
        for n = 1:NBlocks
            framesFFS_BS_FILT{s}{n} = zeros(H,W,NFrames);
            framesFFS_B_FILT{s}{n} = zeros(H,W,NFrames);
        end
    end
    
    framesFFS_AVG_BS_FILT{s} = zeros(H,W,NFrames);
    %% preallocate for FFS
    if s==1
        framesFFS_blank = cell(NFrames,1);
        for f = 1:NFrames
            framesFFS_blank{f} = zeros(H,W,NBlocks);
        end
    end
    
    framesFFS = cell(NFrames,1);
    for f = 1:NFrames
        framesFFS{f} = zeros(H,W,NBlocks);
    end
    
    
    %% begin loop
    c1 = clock;
    for n = 1:NBlocks
        
        %% load data
        if s==1
            frames_blank(:,:,:) = OIReadStim([blankBlocks(n).folder,'\',blankBlocks(n).name], alignInd1, 'v');
            frames_blank(:,:,excludeFrames) = 0;
        end
        frames(:,:,:) = OIReadStim([blocksCond(n).folder,'\',blocksCond(n).name], alignInd2, 'v');
        frames(:,:,excludeFrames) = 0;
        
        %% FFS
        if SP_INDI_BLKS
            framesFFS_BS = zeros(H,W,NFrames);
            framesFFS_B = zeros(H,W,NFrames);
        end
        
        for f = 1:NFrames
            
              frames(:,:,f) = OIShift(frames(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
                  -1*round(align.disp_y_med{alignInd2}{n,f}));
            if FFS
                if s==1
                      frames_blank(:,:,f) = OIShift(frames_blank(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                          -1*round(align.disp_y_med{alignInd1}{n,f}));
                    framesFFS_blank{f}(:,:,n) = (frames_blank(:,:,f) - frames_blank(:,:,fframe)) ./ frames_blank(:,:,fframe) * 100;
                end
                framesFFS{f}(:,:,n) = (frames(:,:,f) - frames(:,:,fframe)) ./ frames(:,:,fframe) * 100;
                
            else % no FFS
                if s==1
                    frames_blank(:,:,f) = OIShift(frames_blank(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                        -1*round(align.disp_y_med{alignInd1}{n,f}));
                    framesFFS_blank{f}(:,:,n) = frames_blank(:,:,f);
                end
                framesFFS{f}(:,:,n) = frames(:,:,f);
            end
            
            if SP_INDI_BLKS  % filter individual blocks for superpixel
                framesFFS_BS(:,:,f) = (framesFFS{f}(:,:,n));
                framesFFS_B(:,:,f) = framesFFS_blank{f}(:,:,n);
                if FILTER
                    framesFFS_BS_FILT{s}{n}(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_BS(:,:,f),LPk,HPk,clipMask);
                    framesFFS_B_FILT{s}{n}(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_B(:,:,f),LPk,HPk,clipMask);
                else
                    framesFFS_BS_FILT{s}{n}(:,:,f) = framesFFS_BS(:,:,f);
                end
                
                if CLIP
                    [framesFFS_BS_FILT{s}{n}(:,:,f), ~, ~, ~] = OIClipH(framesFFS_BS_FILT{s}{n}(:,:,f),...
                        clipMethod, clipValue, clipMask);
                end
            end
        end
        
        clear frames frames_blank framesFFS_BS;  % clear up some memory
        
        
        %progress and time estimation
        if n==1
            time_loop = etime(clock,c1);
            time_est = round(time_loop*100/60)/100 * NBlocks;
            disp(['Loading stim #',num2str(stimSelect(s)),' is estimated to take ',num2str(time_est),' minutes.']);
            dispstat(['Loading stim #',num2str(stimSelect(s)),'... 0.00% Done'],'init');
        end
        
        percDone = (n-1)/NBlocks;
        percDone = round(percDone*10000)/100;
        dispstat(['Loading stim #',num2str(stimSelect(s)),'... ',num2str(percDone),'% Done']);
    end
    
    time_loop = etime(clock,c1);
    time_est = round(time_loop*100/60)/100;
    dispstat(['Loading stim #',num2str(stimSelect(s)),'... 100% Done. Took ',num2str(time_est),' minutes.']);
    
    
    %% average across blocks
    dispstat('','init');
    dispstat('Averaging across blocks... 0.0% Done.');
    if s==1
        framesFFS_blank_AVG = zeros(H,W,NFrames);
    end
    framesFFS_AVG{s} = zeros(H,W,NFrames);
    
    for f = 1:NFrames
        if s==1
            framesFFS_blank_AVG(:,:,f) = mean(framesFFS_blank{f},3);
        end
        framesFFS_AVG{s}(:,:,f) = mean(framesFFS{f},3);
        
        dispstat(['Averaging across blocks... ',num2str(round(f/NFrames*100*100)/100),'% Done.']);
    end
    
    dispstat('Averaging across blocks... 100% Done.');
    
    
    %% blank subtraction and filtering
    framesFFS_AVG_BS{s} = zeros(H,W,NFrames);
    
    c1 = clock;
    for f = 1:NFrames
        framesFFS_AVG_BS{s} = framesFFS_AVG{s}(:,:,f) - framesFFS_blank_AVG(:,:,f);
        
        if FILTER
            framesFFS_AVG{s}(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_AVG{s}(:,:,f),LPk,HPk,clipMask);
            framesFFS_blank_AVG(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_blank_AVG(:,:,f),LPk,HPk,clipMask);
            framesFFS_AVG_BS_FILT{s}(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_AVG_BS{s},LPk,HPk,clipMask);
        else
            framesFFS_AVG_BS_FILT{s}(:,:,f) = framesFFS_AVG_BS{s};
        end
        
        if CLIP
            framesFFS_AVG_BS_FILT{s}(:,:,f) = OIClipH(framesFFS_AVG_BS_FILT{s}(:,:,f), 2, clipValue, clipMask);
        end
        
        if f==1
            time_loop = etime(clock,c1);
            time_est = round(time_loop*100/60)/100 * NFrames;
            disp(['Filtering stim #',num2str(stimSelect(s)),' is estimated to take ',num2str(time_est),' minutes.']);
            dispstat(['Filtering stim #',num2str(stimSelect(s)),'... 0.00% Done'],'init');
        end
        
        percDone = (f-1)/NFrames;
        percDone = round(percDone*10000)/100;
        dispstat(['Filtering stim #',num2str(stimSelect(s)),'... ',num2str(percDone),'% Done']);
    end
    
    clear framesFFS_blank framesFFS;  % clear out memory
    
    
    time_loop = etime(clock,c1);
    time_est = round(time_loop*100/60)/100;
    dispstat(['Filtering stim #',num2str(stimSelect(s)),'... 100% Done. Took ',num2str(time_est),' minutes.']);
    %disp(' ');
end
cd(parent);


%% sp calculation
if SP_INDI_BLKS  % superpixel for each block
    if ~USE_MASKS
        spPeak = cell(length(stimSelect),1);
        spPeakInd = cell(length(stimSelect),1);
        spPeakAvg = cell(length(stimSelect),1);
        spPeakIndAvg = cell(length(stimSelect),1);
        spPeakAvg2 = cell(length(stimSelect),1);
        spPeakIndAvg2 = cell(length(stimSelect),1);
    end
    
    for s = 1:length(stimSelect)
        if ~USE_MASKS
            spPeak{s} = zeros(H,W,NBlocks);
            spPeakInd{s} = zeros(H,W,NBlocks);
            spPeakAvg{s} = zeros(H,W);
            spPeakIndAvg{s} = zeros(H,W);
            spPeakAvg2{s} = zeros(H,W);
            spPeakIndAvg2{s} = zeros(H,W);
        end
        
        c1 = clock;
        for n = 1:NBlocks
            
            if ~USE_MASKS  % if no masks, do superpixel at every pixel
                for x = 1:H
                    for y = 1:W
                        [spPeak{s}(x,y,n), spPeakInd{s}(x,y,n)] = min(squeeze(framesFFS_BS_FILT{s}{n}(x,y,:)));
                    end
                end
                
            else  % use masks
                for m = 1:length(mask)
                    for f = 1:NFrames
                        tempData = framesFFS_BS_FILT{s}{n}(:,:,f);
                        maskMean_BS{s}(m,f,n) = mean(tempData(mask{m}>0));
                        
                        tempData2 = framesFFS_B_FILT{s}{n}(:,:,f);
                        maskMean_B{s}(m,f,n) = mean(tempData2(mask{m}>0));
                    end
                end
            end
            
            if n==1
                time_loop = etime(clock,c1);
                time_est = round(time_loop*100/60)/100 * NBlocks;
                disp(['Computing superpixel for stim #',num2str(stimSelect(s)),' is estimated to take ',num2str(time_est),' minutes.']);
                dispstat('','init');
                dispstat(['Computing superpixel for stim #',num2str(stimSelect(s)),'... 0.00% Done']);
            end
            
            percDone = (n-1)/NBlocks;
            percDone = round(percDone*10000)/100;
            dispstat(['Computing superpixel for stim #',num2str(stimSelect(s)),'... ',num2str(percDone),'% Done']);
        end
        
        time_loop = etime(clock,c1);
        time_est = round(time_loop*100/60)/100;
        dispstat(['Computing superpixel for stim #',num2str(stimSelect(s)),'... 100% Done. Took ',num2str(time_est),' minutes.']);
        disp(' ');
        
        
        if USE_MASKS  % average mask superpixel over blocks
            maskMean_BS_Avg{s} = mean(maskMean_BS{s},3);
            maskMean_BS_Std{s} = std(maskMean_BS{s},0,3)/sqrt(size(maskMean_BS{s},3));
            
            if LPFILTER
                for m = 1:size(maskMean_BS_Avg{s},1)
                    if m==1
                        %                         if ImgHz == 5
                        %                             fc = ImgHz/15;
                        %                         elseif ImgHz == 20
                        %                             fc = ImgHz/6;
                        %                         end
                        %                         %fc = ImgHz/15; %ImgHz/6;
                        %fc = 0.35;
                        %fc = 1.2;
                        %fs = ImgHz;
                        %[b1, a1] = butter(2,fc/(fs/2),'low');
                        [b1, a1] = butter(2,0.14,'low');
                    end
                    maskMean_BS_Avg{s}(m,:) = filtfilt(b1,a1,maskMean_BS_Avg{s}(m,:));
                    maskMean_BS_Std{s}(m,:) = filtfilt(b1,a1,maskMean_BS_Std{s}(m,:));
                    
                    INTERP_FACTOR = 4;
                    maskMeanAvg_interp{s}(m,:) = interp(maskMean_BS_Avg{s}(m,:),INTERP_FACTOR);
                    xtime_interp = linspace(xtime(1), xtime(end), size(maskMeanAvg_interp{s},2));
                    [maskMinVal(m), maskMinLoc(m)] = min(maskMeanAvg_interp{s}(m,:));
                    %maskMinLoc(m) = ceil(maskMinLoc(m)/INTERP_FACTOR);
                    maskMinLoc(m) = xtime_interp(maskMinLoc(m));
                end
            end
            
            if NORMALIZE
                for m = 1:size(maskMean_BS_Avg{s},1)
                    maskMean_BS_Std{s}(m,:) = maskMean_BS_Std{s}(m,:)/abs(min(maskMean_BS_Avg{s}(m,:)));
                    maskMean_BS_Avg{s}(m,:) = maskMean_BS_Avg{s}(m,:)/abs(min(maskMean_BS_Avg{s}(m,:)));
                end
            end
        else
            for x = 1:H
                for y = 1:W
                    spPeakAvg{s}(x,y) = mean(spPeak{s}(x,y,:));
                    spPeakIndAvg{s}(x,y) = mean(spPeakInd{s}(x,y,:));
                end
            end
        end
    end
end

%% superpixel for averaged blocks
for s = 1:length(stimSelect)
    if ~USE_MASKS
        for x = 1:H
            for y = 1:W
                [spPeakAvg2{s}(x,y), spPeakIndAvg2{s}(x,y)] = min(squeeze(framesFFS_AVG_BS_FILT{s}(x,y,:)));
            end
        end
    else
        for m = 1:length(mask)
            for f = 1:NFrames
                tempData2 = framesFFS_AVG{s}(:,:,f);
                maskMean_stim_Avg2{s}(m,f) = mean(tempData2(mask{m}>0));
                
                tempData2 = framesFFS_blank_AVG(:,:,f);
                maskMean_blank_Avg2{s}(m,f) = mean(tempData2(mask{m}>0));
                
                tempData2 = framesFFS_AVG_BS_FILT{s}(:,:,f);
                maskMean_BS_Avg2{s}(m,f) = mean(tempData2(mask{m}>0));
            end
            
            if LPFILTER
                maskMean_BS_Avg2{s}(m,:) = filtfilt(b1,a1,maskMean_BS_Avg2{s}(m,:));
                maskMean_BS_Avg2_interp{s}(m,:) = interp(maskMean_BS_Avg2{s}(m,:),INTERP_FACTOR);
                [maskMinVal2(m), maskMinLoc2(m)] = min(maskMean_BS_Avg2_interp{s}(m,:));
                %maskMinLoc2(m) = ceil(maskMinLoc2(m)/INTERP_FACTOR);
                maskMinLoc2(m) = xtime_interp(maskMinLoc2(m));
            end
            
            if NORMALIZE
                maskMean_BS_Avg2{s}(m,:) = maskMean_BS_Avg2{s}(m,:)/abs(min(maskMean_BS_Avg2{s}(m,:)));
            end
        end
    end
end


%% total loop time
time_loop2 = etime(clock,c0);
time_est = round(time_loop2*100/60)/100;
disp(['Script took ',num2str(time_est),' minutes total.']);
disp(' ');
cd('S:\Lab\ngc14\MsHowell\Left_Chamber\2013_04_23_MsHowell\run1\Results\Superpixel_Figure');
save('trialMask_NEW', 'maskMean_BS', 'maskMean_B');
%% plot
if 1==1
    if ~USE_MASKS  % no masks
        figure;
        ha = tight_subplot(1,2,[0.06 0.01],[0.01 0.04],0.01);
        
        axes(ha(1));
        imagesc(spPeakAvg{1});
        axis image; colorbar;
        title('Superpixel Peak of invididual blocks, averaged');
        %caxis([-255 0]);
        
        axes(ha(2));
        imagesc(spPeakIndAvg{1});
        axis image; colorbar;
        title('Superpixel Peak Index of invididual blocks, averaged');
        caxis([0 20]);
        
    else  % masks
        
        for i = 0:size(maskMean_BS_Avg2{s},1)-1
            lgnd{i+1} = ['Mask #',num2str(i)];
        end
        
        figure('color',[1 1 1]);
        
        for s = 1:length(stimSelect)
            if SP_INDI_BLKS
                if ~GRID_MASKS
                    colors = distinguishable_colors(size(maskMean_BS_Avg{s},1));
                else
                    colors = jet(size(maskMean_BS_Avg{s},1));
                end
                
                %% indi plot
                subplot(length(stimSelect),2,s*2);
                
                for i = 1:size(maskMean_BS_Avg{s},1)
                    ax = shadedErrorBar(xtime,maskMean_BS_Avg{s}(i,:),maskMean_BS_Std{s}(i,:));
                    ax.edge(1).Visible = 'off';
                    ax.edge(2).Visible = 'off';
                    
                    if i==1, hold on, end;
                end
                
                xlabel('Time (s)')
                ylabel('Reflectance Change');
                title([stims(stimSelect(s)),'  -  SP Indi']);
                xlim([min(xtime) max(xtime)]);
                %legend(lgnd);
                set(gca,'fontsize',16);
                ylim([-0.4 0.2]);
                
                if LPFILTER
                    for i = 1:size(maskMean_BS_Avg{s},1)
                        plot([maskMinLoc(i) maskMinLoc(i)],ylim,'--','color',colors(i,:),'linewidth',1);
                        plot(maskMinLoc(i),maskMinVal(i),'o','color',colors(i,:),'linewidth',1.5);
                    end
                end
                
                
                %% avg plot
                subplot(length(stimSelect),2,s*2-1);
                for i = 1:size(maskMean_BS_Avg2{s},1)
                    plot(xtime,maskMean_BS_Avg2{s}(i,:),'linewidth',2,'color',colors(i,:));
                    if i==1, hold on, end;
                end
                
                xlabel('Time (s)')
                ylabel('Reflectance Change');
                title([stims(stimSelect(s)),'  -  SP Avg']);
                xlim([min(xtime) max(xtime)]);
                legend(lgnd);
                set(gca,'fontsize',16);
                ylim([-0.4 0.2]);
                
                if LPFILTER
                    for i = 1:size(maskMean_BS_Avg2{s},1)
                        plot([maskMinLoc2(i) maskMinLoc2(i)],ylim,'--','color',colors(i,:),'linewidth',1);
                        plot(maskMinLoc2(i),maskMinVal2(i),'o','color',colors(i,:),'linewidth',1.5);
                    end
                end
                
                if GRID_MASKS
                    colormap jet;
                    colorbar;
                    caxis([min(maskDist)*0.013 max(maskDist)*0.013]);
                end
                
            else % only indi plot
                subplot(1,length(stimSelect),s);
                plot(xtime,maskMean_BS_Avg2{s},'linewidth',2);
                
                xlabel('Time (s)');
                ylabel('Reflectance Change');
                title([stims(stimSelect(s)),'  -  SP Avg']);
                xlim([min(xtime) max(xtime)]);
                legend(lgnd);
                set(gca,'fontsize',16);
            end
        end
        
        suptitle([monkey,' ',date(date~='_'),' ',run,'  -  HP',num2str(HPk),...
            ' LP',num2str(LPk),' C',num2str(clipValue)]);
        
        
        if GRID_MASKS
            figure('color',[1 1 1]);
            subplot(3,1,1);
            plot(maskDist*0.013,maskMinVal,'.','markersize',12);
            %ylim([-0.5 0]);
            xlabel('Distance Between Mask and Electrode (mm)');
            ylabel('Minimum R-Change (%)');
            box off;
            set(gca,'fontsize',12);
            
            subplot(3,1,2);
            plot(maskDist*0.013,maskMinLoc,'.','markersize',12);
            %ylim([1 2.2]);
            xlabel('Distance Between Mask and Electrode (mm)');
            ylabel('Time of Peak R-Change (s)');
            box off;
            set(gca,'fontsize',12);
            
            for m = 1:length(mask)
                maskMeanAvg_diff{1}(m,:) = diff(maskMean_BS_Avg{1}(m,:));
                maskMeanAvg_diffMin(m) = min(maskMeanAvg_diff{1}(m,:));
            end
            
            subplot(3,1,3);
            plot(maskDist*0.013,maskMeanAvg_diffMin,'.','markersize',12);
            %ylim([-0.025 0]);
            xlabel('Distance Between Mask and Electrode (mm)');
            ylabel('Max slope of descent');
            box off;
            set(gca,'fontsize',12);
            
            suptitle([monkey,' ',date(date~='_'),' ',run,'  -  HP',num2str(HPk),...
                ' LP',num2str(LPk),' C',num2str(clipValue)]);
        end
    end
end

%% individual condition maps
if ~exist([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)])
    mkdir([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)]);
end
cd([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)]);
for i = 0:size(maskMean_BS_Avg2{s},1)-1
    lgnd{i+1} = ['Mask #',num2str(i)];
end
for s =1:length(stimSelect)
    fig1 = figure('color',[1 1 1]);
    subplot(2,2,1);
    p = plot(maskMean_blank_Avg2{s}','linewidth',1.25);
    ylabel('R/dR %');
    set(gca,'fontsize',14);
    if FFS, ylim([-0.2 0.1]); end
    title(['[',stims(2),']']);
    
    colors = distinguishable_colors(size(maskMean_blank_Avg2{s},1));
    for i = 1:length(p)
        set(p(i),'color',colors(i,:));
    end
    
    subplot(2,2,3);
    p = plot(maskMean_stim_Avg2{s}','linewidth',1.25);
    hold on;
    baseline_stim = mean(maskMean_stim_Avg2{s}(:,1:10),2);
    vars_stim = std(maskMean_stim_Avg2{s}(:,1:10),0,2);
    sigInds_stim = maskMean_stim_Avg2{s} < baseline_stim-3.5*vars_stim;
    sigInds_stim(:,1:25) = 0;
    [~, sigFrame_stim] = min(maskMean_stim_Avg2{s}, [],2);
    scatter(sigFrame_stim,maskMean_stim_Avg2{s}...
        (sub2ind(size(maskMean_stim_Avg2{s}), 1:length(sigFrame_stim),sigFrame_stim')),50*ones(size(sigFrame_stim)),colors);
    xlabel('Frame #');
    ylabel('R/dR %');
    set(gca,'fontsize',14);
    if FFS, ylim([-0.2 0.1]); end
    legend(lgnd);
    title(['[',stims(stimSelect(s)),']']);
    
    for i = 1:length(p)
        set(p(i),'color',colors(i,:));
    end
    
    subplot(1,2,2);
    % plot((maskMean_stim_Avg2{1}-maskMean_blank_Avg2{1})','linewidth',1.25);
    p = plot((maskMean_BS_Avg2{s})','linewidth',1.25);
    hold on;
    baseline_BS = mean(maskMean_BS_Avg2{s}(:,1:10),2);
    vars_BS = std(maskMean_BS_Avg2{s}(:,1:10),0,2);
    sigInds_BS = maskMean_BS_Avg2{s} < baseline_BS-3.5*vars_BS;
    sigInds_BS(:,1:25) = 0;
    [~, sigFrame_BS] = max(sigInds_BS, [],2);
    scatter(sigFrame_BS,maskMean_BS_Avg2{s}...
        (sub2ind(size(maskMean_BS_Avg2{s}), 1:length(sigFrame_BS),sigFrame_BS')),50*ones(size(sigFrame_BS)),colors);
    ylabel('R/dR %');
    set(gca,'fontsize',14);
    if FFS, ylim([-0.2 0.1]); end
    title(['[',stims(stimSelect(s)),'] - [',stims(2),']']);
    
    for i = 1:length(p)
        set(p(i),'color',colors(i,:));
    end
    
    suptitle([monkey,' ',date(date~='_'),' ',run,'  -  HP',num2str(HPk),' LP',num2str(LPk),' C',num2str(clipValue)]);
    saveas(fig1, [stims(stimSelect(s)), '_masksByConds'], 'fig');
    %%
    fig2 = figure('color',[1 1 1]);
    subplotSize = numSubplots(size(maskMean_BS_Avg2{s},1));
    for m = 1:size(maskMean_BS_Avg2{s},1)
        subplot(subplotSize(1), subplotSize(2),m);
        hold on;
        plot(smooth(maskMean_blank_Avg2{s}(m,:)),'r','linewidth',1.25);
        plot(smooth(maskMean_stim_Avg2{s}(m,:)),'g','linewidth',1.25);
        plot(smooth(maskMean_BS_Avg2{s}(m,:)),'b','linewidth',1.25);
        if m==1
            legend([{stims(2)}, {stims(stimSelect(s))}, {[stims(stimSelect(s)), '-', stims(2)]}]);
        end
%         scatter(sigFrame_stim(m),maskMean_stim_Avg2{s}(m,sigFrame_stim(m)),35,'g');
%         scatter(sigFrame_BS(m),maskMean_BS_Avg2{s}(m,sigFrame_BS(m)),35,'b');
        title(['Mask #',num2str(m-1)]);
        xlabel('Frame #');
        ylabel('R/dR %');
        set(gca,'fontsize',14);
    end
    for m = 1:length(mask)
        [dm_rows,dm_cols] = find(mask{m});
        mask_stim_frames = squeeze(mean(mean(framesFFS_AVG{s}(dm_rows,dm_cols,:))));
        mask_BS_frames = squeeze(mean(mean(framesFFS_AVG_BS_FILT{s}(dm_rows,dm_cols,:))));
        
%         mask_stim_var = std2(framesFFS_AVG{s}(dm_rows,dm_cols,1:10));
%         mask_BS_var = std2(framesFFS_AVG_BS_FILT{s}(dm_rows,dm_cols,1:10));
%         mask_stim_t = find(mask_stim_frames<mean(mask_stim_frames(1:10))-2.5*mask_stim_var);
%         mask_BS_t = find(mask_BS_frames<mean(mask_BS_frames(1:10))-2.5*mask_BS_var);
%         mask_stim{m} = mask_stim_t(find(mask_stim_t>25,1));
%         mask_BS{m} = mask_BS_t(find(mask_BS_t>25,1));
        
                max_stim = max([inf(1,25) -inf(1,length(mask_stim_frames)-25)]',mask_stim_frames);
                max_BS = max([inf(1,25) -inf(1,length(mask_BS_frames)-25)]',mask_BS_frames);
                max_stim(isinf(max_stim)) = NaN;
                max_BS(isinf(max_BS)) = NaN;
                [~,mask_stim{m}] = nanmin(max_stim);
                [~,mask_BS{m}] = nanmin(max_BS);
    end
%     for d = 1:length(domains)
%         domain = domains{d};
%         maskInds = cellfun(@(a) strcmp(a, domain), maskNames);
%         domain_stim = mask_stim(maskInds);
%         domain_BS = mask_BS(maskInds);
%         save([stims(stimSelect(s)), '_',domain], 'domain_stim', 'domain_BS');
%     end
    saveas(fig2, [stims(stimSelect(s)), '_individualMasks'], 'fig');
end
%% save
if SAVE && USE_MASKS
    if ~exist([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)])
        mkdir([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)]);
    end
    cd([path,'\Results\Superpixel\LP',num2str(LPk),'HP',num2str(HPk),'C',num2str(clipValue)]);
    
    SPData.Parameters.animal = monkey;
    SPData.Parameters.date = date;
    SPData.Parameters.run = run;
    SPData.Parameters.stim = stimSelect;
    SPData.Parameters.NBlocks = NBlocks;
    SPData.Parameters.NFrames = NFrames;
    SPData.Parameters.ImgHz = ImgHz;
    SPData.Parameters.LPk = LPk;
    SPData.Parameters.HPk = HPk;
    SPData.Parameters.clip = clipValue;
    SPData.Parameters.LPFILTER = LPFILTER;
    SPData.Parameters.NORMALIZE = NORMALIZE;
    SPData.Parameters.GRID_MASKS = GRID_MASKS;
    SPData.Parameters.GRID_MASKS_INVERT = GRID_MASKS_INVERT;
    
    SPData.maskMeanAvg = maskMean_BS_Avg;
    SPData.maskMeanStd = maskMean_BS_Std;
    SPData.maskMean_BS_Avg2 = maskMean_BS_Avg2;
    SPData.mask = mask;
    
    if GRID_MASKS
        if GRID_MASKS_INVERT
            save('SPData_GridMasks_Invert.mat','SPData');
        else
            save('SPData_GridMasks.mat','SPData');
        end
    else
        save('SPData.mat','SPData');
    end
    cd(parent);
end


%% email
if EMAIL_UPON_COMPLETION
    mailSubject = ['Superpixel for ',monkey,' ',date(date~='_'),' ',run,' is complete!'];
    mailBody = ['Superpixel for ',monkey,' ',date(date~='_'),' ',run,' is complete and took '...
        num2str(time_est),' minutes to run.'];
    
    sendEmail('nsc15@pitt.edu',mailSubject,mailBody);
end

% clearvars;

%end