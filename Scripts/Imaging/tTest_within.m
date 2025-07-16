%% user input
username = 'ngc14';

animal = 'Gilligan';
date = '01_07_2019';
run = 'run00';

ALIGN_OFFSET = true;
stim = containers.Map([0 2 4 5],{'Rest', 'ExtraSmallSphere', 'LargeSphere', 'Photocell'});
p_values = [0.01 0.001 0.0001];
conditions = [2];

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
HPk=700;                    % High-pass kernel size (usually greater than 250) set as 0 if no HP filter

numCols = 10;                % Number of columns in big image.
blockSelect = []; % Use only a subset of blocks. Leave empty to use all blocks.

%% locate data
maps_avg = [];
dataPath = ['\\univ.pitt.edu\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run];
try
    cd(dataPath);
catch error
    error('SSD data path does not exist. Exiting.');
end
clipMaskName = ['clipMask',run(end-1:end),'.bmp'];
greenName = ['green',run(end-1:end),'_edited.bmp'];
parent = ['\\univ.pitt.edu\sni\gharbawie\Lab\',username,'\Scripts\'];
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

conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
numConds = unique(conds);
badAlign = load('bad_alignments.mat');

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

%% Begin Script
for s=1:length(conditions)
    cd(dataPath);
    stim1 = conditions(s);
    disp(['%%%%% BLOCKVIEW_NGC14: ',animal,' ',date(date~='_'),' ',run,' [',stim(stim1),'] LP',...
        num2str(LPk),' HP',num2str(HPk),' f[',num2str(frameRange),'] %%%%%']);
    
    blockConds1 = blocks(conds==stim1);
    alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
    badInds1 = ~ismember([1:size(blockConds1,1)],badAlign.badAlign{alignInd1});
    blockConds1 = blockConds1(badInds1);
    blockConds2 = blocks(conds==stim2);
    alignInd2 = find(cellfun(@(a) a==stim2, stim.keys));
    badInds2 = ~ismember([1:size(blockConds2,1)],badAlign.badAlign{alignInd2});
    blockConds2 = blockConds2(badInds2);
    NBlocks = size(blockConds1,1);
    
    %% load data
    maps_stim1 = single(zeros(H,W,NFrames,NBlocks));
        if(s==1)
            maps_stim2 = single(zeros(H,W,NBlocks));
        end
    
    if ALIGN_OFFSET
        align = load('alignment_offset.mat');
        align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
        align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
                if(alignInd1~=alignInd2)
                    align.disp_x_med{alignInd2} = align.disp_x_med{alignInd2}(badInds2,:);
                    align.disp_y_med{alignInd2} = align.disp_y_med{alignInd2}(badInds2,:);
                end
    else
        align = [];
    end
    hbar = parfor_progressbar(NBlocks,'Loading and Aligning Data...'); %create the progress bar
    % dispstat('','init');
    for n = 1:NBlocks
        
        % reset stim1 and stim2 to zeros
        frames_stim1 = zeros(H,W,NFrames);
                if(s==1)
                    frames_stim2 = zeros(H,W,NFrames);
                end
        
        % load stim 1
        frames_stim1(:,:,:) = single(OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],1,'v'));
        
        
        %if blank sub, load stim 2
                if style==3 || style==4 && s==1
                    frames_stim2(:,:,:) = single(OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],1,'v'));
                end
        
        if ALIGN_OFFSET
            % median filter to reduce CCD column defects
            for f = 1:size(frames_stim1,3)
                frames_stim1(:,:,f) = medfilt2(frames_stim1(:,:,f),[5 5]);
                
                                if style==3 || style==4 && s==1
                                    frames_stim2(:,:,f) = medfilt2(frames_stim2(:,:,f),[5 5]);
                                end
            end
            for f = 1:size(frames_stim1,3)
                %             disp(['n',num2str(n),' f',num2str(f)]);
                frames_stim1(:,:,f) = OIShift(frames_stim1(:,:,f),-1*round(align.disp_x_med{alignInd1}{n,f}),...
                    -1*round(align.disp_y_med{alignInd1}{n,f}));
                
                                if style==3 || style==4 && s==1
                                    frames_stim2(:,:,f) = OIShift(frames_stim2(:,:,f),-1*round(align.disp_x_med{alignInd2}{n,f}),...
                                        -1*round(align.disp_y_med{alignInd2}{n,f}));
                                end
            end
        end
        
        switch style
            %             case 1 % single condition
            %                 maps_stim1 = maps_stim1 + frames_stim1;
            %
            case 2 % single condition w/ FFS
                maps_stim1(:,:,:,n) = single((frames_stim1 - mean(frames_stim1(:,:,fframe),3)) ./...
                    mean(frames_stim1(:,:,fframe),3));
                %
                %             case 3 % blank subtraction
                %                 maps_stim1 = maps_stim1 + frames_stim1;
                %                 if(s==1)
                %                     maps_stim2 = maps_stim2 + frames_stim2;
                %                 end
                
            case 4 % FFS
                maps_stim1(:,:,:,n) = single((frames_stim1 - mean(frames_stim1(:,:,fframe),3)) ./...
                    mean(frames_stim1(:,:,fframe),3));
                if(s==1)
                    maps_stim2(:,:,n) = mean(single((frames_stim2 - mean(frames_stim2(:,:,fframe),3)) ./...
                        mean(frames_stim2(:,:,fframe),3)),3);
                end
                maps_stim1(:,:,:,n) = single((maps_stim1(:,:,:,n) - maps_stim2(:,:,n)));
                
        end
        %dispstat(['Loading Data... ',num2str(round(n/NBlocks*10000)/100),'% complete.']);
        
        hbar.iterate(1); % update progress by one iteration
    end
    close(hbar); % close the progress bar
    
    clear frames_stim1 frames_stim2 blockConds2
    % dispstat('Loading Data... 100% complete.');
    maps_stim1(isinf(maps_stim1)) = 0;
    maps_stim1(isnan(maps_stim1)) = 0;
        maps_stim2(isinf(maps_stim2)) = 0;
        maps_stim2(isnan(maps_stim2)) = 0;
    %     %% average maps and blank subtraction
    %     switch style
    %         case 1 % single condition
    %             maps_stim1_avg = maps_stim1 ./ NBlocks;
    %             maps_avg = maps_stim1_avg;
    %         case 2 % single condition
    %             maps_stim1_avg = maps_stim1 ./ NBlocks;
    %             maps_avg = maps_stim1_avg;
    %         case 3 % blank subtraction
    %             maps_stim1_avg = maps_stim1 ./ NBlocks;
    %             if(s==1)
    %                 maps_stim2_avg = maps_stim2 ./ NBlocks;
    %             end
    %             maps_avg = maps_stim1_avg - maps_stim2_avg;
    %         case 4 % blank subtraction
    %             maps_stim1_avg = maps_stim1 ./ NBlocks;
    %             if(s==1)
    %                 maps_stim2_avg = maps_stim2 ./ NBlocks;
    %                 maps_stim2_avg = mean(maps_stim2_avg,3);
    %             end
    %             maps_avg = maps_stim1_avg - maps_stim2_avg;
    %     end
    %     % subtract average
    %     if style==1
    %         maps_avg = maps_avg - mean(maps_avg(:));
    %     end
    %
    %     % exclude frames
    %     if ~isempty(excludeFrames)
    %         maps_avg(:,:,excludeFrames) = 0;
    %     end
    %
    %     cd(dataPath);
    %     %% filtering
    %
    %%

        maps_stim1_avg = maps_stim1;
    clear maps_stim2 maps_stim1 blockConds1;
    maps_stim1_avg_FILTERED = single(zeros(H,W,NFrames,NBlocks));
    hbar = parfor_progressbar(NBlocks,'Filtering Data...'); %create the progress bar
    parfor n =1:NBlocks
        for f = 1:NFrames
            %if ~any(ismember(fframe,f))% dont filter fframe if FFS
            maps_stim1_avg_FILTERED(:,:,f,n) = imageFilter_LPHP_nsc15(maps_stim1_avg(:,:,f,n),LPk,HPk,clipMask);
            %end
        end
        hbar.iterate(1); % update progress by one iteration
    end
    close(hbar); % close the progress bar
    clear maps_stim1_avg
    baseline = squeeze(mean(maps_stim1_avg_FILTERED(:,:,2:10,:),3));
    pl = zeros(H,W,NFrames);
    pr = zeros(H,W,NFrames);
    for f = 1:NFrames
        [~, pl(:,:,f)] = ttest2(baseline,squeeze(maps_stim1_avg_FILTERED(:,:,f,:)),0.05,'left','equal',3);
        [~, pr(:,:,f)] = ttest2(baseline,squeeze(maps_stim1_avg_FILTERED(:,:,f,:)),0.05,'right','equal',3);
    end
    %%
    refMask = imread('\\univ.pitt.edu\sni\Gharbawie\Lab\Gilligan\All Data\Mapping\referenceMask', 'bmp')>200;
    refMask = refMask(:,:,1);
    tform = load([dataPath,'\transformation_matrix.mat']);
    tform = inv(tform.tform);
    tform(:,3) = [0 0 1];
    refMask = imwarp(refMask, affine2d(tform),'OutputView',imref2d([W,H]));
    clipMask = ~(~clipMask | ~refMask);

    clear big_map_green
    %% construct big images
    for i = 1:length(p_values)
        pVals = imfill(bwareaopen((pl<p_values(i) & clipMask),10),'holes');
        pVals = pVals.* 2;
        pVals = pVals + imfill(bwareaopen((pr<p_values(i) & clipMask),10), 'holes');
        
        numRows = ceil(NFrames/numCols);
        
        big_map_noClip = uint8(zeros(H*(numRows+1),W*(numCols)));
        big_map_green = uint8(zeros(H*(numRows+1),W*(numCols)));
        big_map_noClip(numRows*H+1:(numRows+1)*H,1:W) = mean(pVals(:,:,frameRange),3);
        big_map_green(numRows*H+1:(numRows+1)*H,1:W) = greenGray;
        k = 1;
        for r = 1:numRows
            for c = 1:numCols
                if k <= NFrames
                    big_map_noClip((r-1)*H+1:r*H,(c-1)*W+1:c*W) = pVals(:,:,k);
                    big_map_green((r-1)*H+1:r*H,(c-1)*W+1:c*W) = greenGray;
                    k = k + 1;
                end
            end
        end
        
        big_map_noClip_small = imresize(big_map_noClip,0.25);
    big_map_green_small = imresize(big_map_green,0.25);

        %% save data
        cd(dataPath);
        savePath = ['Results\tTest\',stim(conditions(s)),'\LP',num2str(LPk),'_HP',num2str(HPk),'\p-',num2str(p_values(i)),'\'];
        if ~exist(savePath,'dir')
            mkdir(savePath);
        end
        cd(savePath);
        
        figure();
        ax_g = axes;
        imagesc(big_map_green_small);
        colormap(ax_g, 'gray');
        ax_m = axes;
        imagesc(ax_m, big_map_noClip_small,'AlphaData',big_map_noClip_small>0);
        colormap(ax_m,hsv);
        caxis(ax_m,[1 3]);
        ax_m.Visible = 'off';
        linkprop([ax_g ax_m],'Position');
        
        saveas(gcf,['tTest_p-',num2str(p_values(i)),'.bmp']);
        pause(.5);
        close all;
        for f = 1:NFrames
            ax1 = axes;
            imagesc(greenGray)
            colormap(ax1, 'gray');
            ax2 = axes;
            imagesc(ax2, pVals(:,:,f),'AlphaData',pVals(:,:,f)>0);
            ax2.Visible = 'off';
            linkprop([ax1 ax2],'Position');
            saveas(gcf, [num2str(f), '.bmp']);
        end
        saveas(plot(squeeze(sum(sum(pVals)))),['sigPix_p-',num2str(p_values(i)),'.fig']);
    end
    cd(parent);
end