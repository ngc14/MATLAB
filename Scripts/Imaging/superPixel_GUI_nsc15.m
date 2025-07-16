function [] = superPixel_GUI_nsc15()
clearvars;
% close all;

%% user input
username = 'ngc14';

animal = 'Gilligan';                   % animal name
hemi = 'Right';
date = '01_07_2019';                % experiment date
run = 'run00';                      % run name
blockSelect = [];                     % leave empty for all blocks
excludeFrames = []; %[5:12];

ImgHz = 5;

FILTER = true;                % HP/LP filter and clip?
HPk = 0;                          % High pass filter kernel
HPMethod = 'fastmedian';            % High pass filter method
LPk = 0;                           % Low pass filter kernel
LPMethod = 'gaussian';              % Low pass filter method

CLIP = false;
clipMethod = 2;                     % Clip Method
clipValue = 2;                         % Clip value (SDs)
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];

frameLimit = 100;
blockLimit = 100;

stimSelect = [2];                   % stim #s to do superpixel for
FFS = true;
fframe = 1;
SP_INDI_BLKS = false;               % do superpixel for each block (takes a lot longer), or the averaged blocks? 
LPFILTER = false;
NORMALIZE = false;                   % if true, normalize masks to have a peak at -1

stims{1} = 'Blank';                 % Stim condition names
stims{2} = '500ms_Blue_LED';
% stims{2} = '150BiPulses_60microAmps';
% stims{3} = '150BiPulses_60microAmps';
% stims{4} = '150BiPulses_80microAmps';

diskRad = 5;
imXLim = [1 1308]; %[1 768]; %[350 550]; %[1 768];
imYLim = [1 1080]; %[1 768]; %[150 350]; %[1 768];

TTEST_OVERLAY = false;
ttestLP = '0';
ttestHP = '0';
ttestFR = '8-10';
ttestP = '0.001';
ttestStim = '150BiPulses_40microAmps';


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
disp(['%%%%% SUPERPIXEL: ',animal,' ',date(date~='_'),' ',run,' LP',...
    num2str(LPk),' HP',num2str(HPk),' C',num2str(clipValue),' %%%%%']);

% data path
path = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
parent = ['\\pitt\sni\gharbawie\Lab\',username,'\Scripts\'];
% load(['C:\Users\nsc15\Documents\Data\',animal,'_SqM\',hemi,'Hemisphere\dateRunSite_',animal,'.mat']);

% get additional data properties
cd(path);

% block identification
blocks = dir('*.BLK');

if isempty(blocks) %look on backup drive instead
    disp('Data not found on SSD. Attempting to load from HDD.');
    
    try
        cd(['D:\Data\',animal,'_SqM\',hemi,'Hemisphere\',date,'\',run]);
    catch error
        error('HDD data path does not exist. Exiting.');
    end
    
    blocks = dir('*.BLK');
    cd(path);
end

if ~isempty(blockSelect)
    blocks = blocks(blockSelect);
end

anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = size(blocks,1);


%% green
green = imread(['green',run(end-1:end),'_edited.bmp']);

if TTEST_OVERLAY
    ttestPath = ['Results\Ttest_nsc15_V4\Frames',ttestFR,'\LP',ttestLP,'HP',ttestHP,...
        '\[',ttestStim,']-[Blank]_p-',ttestP,'.png'];
    [~,~,ttestAlpha] = imread(ttestPath);
    ttestColor(:,:,1) = ones(size(ttestAlpha,1),size(ttestAlpha,2))*1;
    ttestColor(:,:,2) = ones(size(ttestAlpha,1),size(ttestAlpha,2))*1;
    ttestColor(:,:,3) = ones(size(ttestAlpha,1),size(ttestAlpha,2))*0;
end

%% clip mask
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
if NBlocks > blockLimit
    NBlocks = blockLimit;
end

xTime = ((1:NFrames)-1)*(1/ImgHz);


%% preallocating
framesFFS_blank_AVG = zeros(H,W,NFrames);
framesFFS_AVG = zeros(H,W,NFrames);

framesFFS_AVG_BS_FILT = cell(length(stimSelect),1);
for s = 1:length(stimSelect)
    framesFFS_AVG_BS_FILT{s} = zeros(H,W,NFrames);
end

%% begin image loading/filtering loop
for s = 1:length(stimSelect)
    
    %% begin loop
    c1 = clock;
    for n = 1:NBlocks
        
        %% load data
        if s==1
            frames_blank(:,:,:) = OIReadStim([blocks(n).folder,'\',blocks(n).name], 1, 'v');
            frames_blank(:,:,excludeFrames) = 0;
        end
        frames(:,:,:) = OIReadStim([blocks(n).folder,'\',blocks(n).name], 0, 'v');
        frames(:,:,excludeFrames) = 0;
        
        %% FFS       
        for f = 1:NFrames
            
            if FFS
                if s==1
                    framesFFS_blank_AVG(:,:,f) = framesFFS_blank_AVG(:,:,f) +...
                        ((frames_blank(:,:,f) - frames_blank(:,:,fframe)) ./ frames_blank(:,:,fframe) * 100);
                end
                framesFFS_AVG(:,:,f) = framesFFS_AVG(:,:,f) +...
                    ((frames(:,:,f) - frames(:,:,fframe)) ./ frames(:,:,fframe) * 100);
                disp('');
            else % no FFS
                if s==1
                    framesFFS_blank_AVG{f} = framesFFS_blank_AVG(:,:,f) + frames_blank(:,:,f);
                end
                framesFFS_AVG(:,:,f) = framesFFS_AVG(:,:,f) + frames(:,:,f);
            end
        end      
        
        clear frames frames_blank;
        
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
    
    framesFFS_blank_AVG = framesFFS_blank_AVG / NBlocks;
    framesFFS_AVG = framesFFS_AVG / NBlocks;
    
    time_loop = etime(clock,c1);
    time_est = round(time_loop*100/60)/100;
    dispstat(['Loading stim #',num2str(stimSelect(s)),'... 100% Done. Took ',num2str(time_est),' minutes.']);
    
    
    %% blank subtraction and filtering
    framesFFS_AVG_BS = zeros(H,W,NFrames);
    dispstat('','init');
    c1 = clock;
    
    for f = 1:NFrames
        framesFFS_AVG_BS(:,:,f) = framesFFS_AVG(:,:,f) - framesFFS_blank_AVG(:,:,f);
        
        if FILTER
            framesFFS_AVG(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_AVG(:,:,f),LPk,HPk,clipMask);
            framesFFS_blank_AVG(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_blank_AVG(:,:,f),LPk,HPk,clipMask);
            framesFFS_AVG_BS_FILT{s}(:,:,f) = imageFilter_LPHP_nsc15(framesFFS_AVG_BS(:,:,f),LPk,HPk,clipMask);
        else
            framesFFS_AVG_BS_FILT{s}(:,:,f) = framesFFS_AVG_BS(:,:,f);
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
    
    time_loop = etime(clock,c1);
    time_est = round(time_loop*100/60)/100;
    dispstat(['Filtering stim #',num2str(stimSelect(s)),'... 100% Done. Took ',num2str(time_est),' minutes.']);
    %disp(' ');
end
cd(parent);

%% movie
if 1==0
    figure;
    for f = 1:NFrames
        imagesc(framesFFS_blank_AVG(:,:,f));
        axis image off;
        colormap gray;
        title(num2str(f));
        pause(0.2);
    end
end 

%% GUI
GUIDims = [1800 950];
imDims = [H W];

if imDims(2)>GUIDims(1) || imDims(1)>GUIDims(2)
    imScale = min([GUIDims(1)/imDims(2) GUIDims(2)/imDims(1)]);
else
    imScale(1) = GUIDims(1)/imDims(1);
    imScale(2) = (GUIDims(2)+100)/imDims(2);
    imScale = min(imScale);
end

S.fh = figure('units','pixels',...
    'position',[50 50 GUIDims(1)+50 GUIDims(2)+50],...
    'menubar','none',...
    'numbertitle','off',...
    'name','Superpixel GUI',...
    'resize','on');

S.ax = axes('units','pixels',...
    'position',[20 20 round(imScale*imDims(2))-100 round(imScale*imDims(1))-100]);
imagesc(green);
hold on;
if TTEST_OVERLAY, imagesc(ttestColor,'alphadata',ttestAlpha*0.3); end
axis image off;
title([animal,' ',date(date~='_'),' ',run]);
set(gca,'fontsize',12);
ylim(imYLim)
xlim(imXLim)


S.ax_blank = axes('units','pixels',...
    'position',[1040 625 770 225]);
set(gca,'fontsize',12);

S.ax_stim = axes('units','pixels',...
    'position',[1040 350 770 225]);
set(gca,'fontsize',12);

S.ax_BS = axes('units','pixels',...
    'position',[1040 75 770 225]);
xlabel('Time (s)');
set(gca,'fontsize',12);


pickPoint = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[1040 900 90 40],...
    'string','Pick Point',...
    'fontsize',12);
set(pickPoint,'call',{@pointPickFcn, S});

pickPointCorr = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[1540 900 90 40],...
    'string','Correlation',...
    'fontsize',12);
set(pickPointCorr,'call',{@pointPickCorrFcn, S});

clearPlots = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[1715 900 90 40],...
    'string','Clear Plots',...
    'fontsize',12);
set(clearPlots,'call',{@clearPlotFcn, S});


diskRadTxt1 = uicontrol('style','text',...
    'unit','pix',...
    'position',[1150 890 150 40],...
    'string','Point Radius:',...
    'fontsize',12);

diskRadBox1 = uicontrol('style','edit',...
    'unit','pix',...
    'position',[1280 900 50 40],...
    'string',num2str(diskRad),...
    'fontsize',12);

set(diskRadBox1,'call',{@diskRadFcn1, S});



%% point pick function
function pointPickFcn(pickPoint,event,S)
    [x, y] = getpts(S.ax);
    x = round(x);
    y = round(y);
    
    cla(S.ax);
    cla(S.ax_blank);
    cla(S.ax_stim);
    cla(S.ax_BS);
    
    for i = 1:length(x)
        xt = x(i);
        yt = y(i);

        disk = fspecial('disk',diskRad);
        diskImage = zeros(size(green,1),size(green,2));
        diskImage(yt-diskRad:yt+diskRad,xt-diskRad:xt+diskRad) = disk;

        for f = 1:NFrames
            blankIm = framesFFS_blank_AVG(:,:,f);
            blankSig(i,f) = mean(blankIm(diskImage>0));

            stimIm = framesFFS_AVG(:,:,f);
            stimSig(i,f) = mean(stimIm(diskImage>0));

            BS_Im = framesFFS_AVG_BS_FILT{1}(:,:,f);
            BS_Sig(i,f) = mean(BS_Im(diskImage>0));
        end

        ylimTemp = [-0.3 0.15];

        numLines = length(get(S.ax_blank,'Children'));
        colors = lines(numLines+1);
%         colors = cmapMaker([0 0 0.8],[0 1 1],length(x));
%         colors = cmapMaker([1 0 0],[0.3 0 0.3],length(x));
        
        colorIm(:,:,1) = ones(size(green,1),size(green,2))*colors(numLines+1,1);
        colorIm(:,:,2) = ones(size(green,1),size(green,2))*colors(numLines+1,2);
        colorIm(:,:,3) = ones(size(green,1),size(green,2))*colors(numLines+1,3);

        % plot
        axes(S.ax);
%         cla(S.ax);
        if i==1, imagesc(green); end
        hold on;
        imagesc(colorIm,'alphadata',diskImage>0);
        axis image off;
        title([animal,' ',date(date~='_'),' ',run]);
        set(gca,'fontsize',12);
        ylim(imYLim)
        xlim(imXLim)

        axes(S.ax_blank);
%         cla(S.ax_blank);
        plot(xTime,blankSig(i,:),'color',colors(numLines+1,:),'linewidth',1);
        hold on;
        ylim(ylimTemp);
        set(gca,'fontsize',12);
        title('Blank');

        axes(S.ax_stim);
%         cla(S.ax_stim);
        plot(xTime,stimSig(i,:),'color',colors(numLines+1,:),'linewidth',1);
        hold on;
        ylim(ylimTemp);
        set(gca,'fontsize',12);
        title('Stim');

        axes(S.ax_BS);
%         cla(S.ax_BS);
        plot(xTime,BS_Sig(i,:),'color',colors(numLines+1,:),'linewidth',1);
        hold on;
        ylim(ylimTemp);
        set(gca,'fontsize',12);
        title('Stim - Blank');
        xlabel('Time (s)');
    end
    disp('');
end


%% point pick correlation function
function pointPickCorrFcn(pickPointCorr,event,S)
    [x, y] = getpts(S.ax);
    x = round(x);
    y = round(y);
    
    x = x(1);
    y = y(1);

    disk = fspecial('disk',diskRad);
    diskImage = zeros(size(green,1),size(green,2));
    diskImage(y-diskRad:y+diskRad,x-diskRad:x+diskRad) = disk;
    
    for f = 1:NFrames
        blankIm = framesFFS_blank_AVG(:,:,f);
        blankSig(f) = mean(blankIm(diskImage>0));
        
        stimIm = framesFFS_AVG(:,:,f);
        stimSig(f) = mean(stimIm(diskImage>0));
        
        BS_Im = framesFFS_AVG_BS_FILT{1}(:,:,f);
        BS_Sig(f) = mean(BS_Im(diskImage>0));
    end
    
    c = zeros(size(green,1),size(green,2));
    dispstat('','init');
    for x = imXLim(1):imXLim(2)
        for y = imYLim(1):imYLim(2)
            c(y,x) = corr(BS_Sig',squeeze(framesFFS_AVG_BS_FILT{1}(y,x,:)));
        end
        
        if rem(x,5)==1
            dispstat(['Calculating Correlations... ',...
                num2str(round((x-imXLim(1))/(imXLim(2)-imXLim(1))*10000)/100),'% Done.']);
        end
    end
    dispstat(['Calculating Correlations... 100% Done.']);
    
    axes(S.ax);
    cla(S.ax);
    imagesc(c);
    hold on;
    imagesc(ones(size(green)),'alphadata',diskImage>0);
    axis image off;
    colormap jet;
    colorbar;
    ylim(imYLim)
    xlim(imXLim)
    
    disp('');
end


%% plot clear function
function clearPlotFcn(clearPlots,event,S)
    axes(S.ax);
    cla(S.ax);
    imagesc(green);
    hold on;
    if TTEST_OVERLAY, imagesc(ttestColor,'alphadata',ttestAlpha*0.3); end
    axis image off;
    title([animal,' ',date(date~='_'),' ',run]);
    set(gca,'fontsize',12);
    ylim(imYLim)
    xlim(imXLim)
    
    cla(S.ax_blank);
    cla(S.ax_stim);
    cla(S.ax_BS);
end

%% disk rad function
function diskRadFcn1(diskRadBox1,event,S)
    diskRad = str2double(diskRadBox1.String);
end

end