function [] = onlineImagingAnalysis_ngc14_V1()
    % computes blockview, sunin, and ttest analysis online
%     close all;
%     clearvars;
    
    dataPath = 'C:\Data\Skipper_152-17\11_06_2020\run00\';
    parent = 'C:\Users\OmarLab\Desktop\Online Imaging Code';
    
    stim1 = 2;
    stim2 = 0;
    frameRange = 51:56;
    fframe = 1;
    pValue = 0.001;
    LP = 5;
    HP = 550;
    clip = 1.75;
    
    excludeFrames = [];
    blockSelect = [];%[1:20,23:25,27:31,33:35,38:72];
    
    CONTINUOUS_MODE = true;
    
    
    %% identify starter block files
    cd(dataPath);
    blkFiles = dir('*.BLK');
    conds = [];
    % sort filenames
    for n = 1:length(blkFiles)
        ind1 = strfind(blkFiles(n).name,'_E');
        ind2 = strfind(blkFiles(n).name,'B');
        ind3 = strfind(blkFiles(n).name,'.BLK');
        
        blockNum(n,1) = str2num(blkFiles(n).name((ind1+2):(ind2(1)-1)));
        blockNum(n,2) = str2num(blkFiles(n).name((ind2(1)+1):(ind3(1)-1)));
        anapar = OIHeadRead((blkFiles(n).name), 'v');
        conds = [conds, anapar.Cond];
    end
    
    [blockNum, sortInd] = sortrows(blockNum,[1 2]);
    blkFiles = blkFiles(sortInd);
    if ~isempty(blockSelect)
        blkFiles = blkFiles(blockSelect);
        conds = conds(blockSelect);
    end
    disp(['Found ',num2str(length(blkFiles)),' BLK files to start.']);
    
    blockConds1 = blkFiles(conds==stim1);
    blockConds2 = blkFiles(conds==stim2);
    %% data properties
    anapar = OIHeadRead([blkFiles(1).folder,'\',blkFiles(1).name], 'v');
    W = anapar.FrameWidth;
    H = anapar.FrameHeight;             % anapar contains information about the block files
    NFrames = anapar.FramesPerStim;
    NConds = anapar.NStim;
    NBlocks = min(length(blockConds1),length(blockConds2));
    
    %% starter variables
    frames_stim1_temp = zeros(H,W,NFrames);
    frames_stim2_temp = zeros(H,W,NFrames);
    
    frames_stim1_avg = zeros(H,W,NFrames);
    frames_stim2_avg = zeros(H,W,NFrames);
    
    frames_stim1_ttest = zeros(H/2,W/2,NBlocks);
    frames_stim2_ttest = zeros(H/2,W/2,NBlocks);

    %% load starter data
    for n = 1:NBlocks
        frames_stim1_temp = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],0,'v');
        frames_stim2_temp = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],0,'v');
        
        if ~isempty(excludeFrames)
            frames_stim1_temp(:,:,excludeFrames) = zeros(H,W,length(excludeFrames));
            frames_stim2_temp(:,:,excludeFrames) = zeros(H,W,length(excludeFrames));
        end
        
        frames_stim1_temp = (frames_stim1_temp - frames_stim1_temp(:,:,fframe)) ./ frames_stim1_temp(:,:,fframe) * 100;
        frames_stim2_temp = (frames_stim2_temp - frames_stim2_temp(:,:,fframe)) ./ frames_stim2_temp(:,:,fframe) * 100;
        
        frames_stim1_avg = frames_stim1_avg + frames_stim1_temp;
        frames_stim2_avg = frames_stim2_avg + frames_stim2_temp;
        
        frames_stim1_ttest(:,:,n) = imageFilter_LPHP_nsc15(imresize(mean(frames_stim1_temp(:,:,frameRange),3),0.5),...
            round(LP/2),round(HP/2),[]);
        frames_stim2_ttest(:,:,n) = imageFilter_LPHP_nsc15(imresize(mean(frames_stim2_temp(:,:,frameRange),3),0.5),...
            round(LP/2),round(HP/2),[]);
    end

    
    %% generate figures and axes
    fig_all = figure('color',[1 1 1],'position',[150 100 1500 850]);
    ax_ttest = axes('position',[0.1 0.1 0.3 0.4]);
    ax_sunin = axes('position',[0.1 0.55 0.3 0.4]);
    ax_blockview = axes('position',[0.45 0.1 0.50 0.85]);
    
    suninImage = imaging_ttest(frames_stim1_ttest,frames_stim2_ttest,pValue);
    clipLims = imaging_blockview(frames_stim1_avg,frames_stim2_avg,LP,HP,clip);
%     suptitle(['Stim',num2str(stim1),'-Stim',num2str(stim2),'     LP',num2str(LP),...
%         ' HP',num2str(HP),' C',num2str(clip),'     Blocks:',num2str(length(blkFiles))]);
    
    pickPoint = uicontrol('style','pushbutton',...
        'unit','pix',...
        'position',[300 870 120 40],...
        'string','Pick SP Points',...
        'fontsize',12);
    set(pickPoint,'call',{@pointPickFcn});
    
    
    %% sumview figure
    nRows = 10;
    nCols = ceil(length(blockNum)/nRows);
    sumviewImage = zeros(nCols*size(frames_stim1_ttest,1),nRows*size(frames_stim1_ttest,2));
    
    k2 = 1;
    for x = 1:nCols
        for y = 1:nRows
            if k2 <= length(NBlocks)
                sumviewImage(((x-1)*size(frames_stim1_ttest,1)+1):(x*size(frames_stim1_ttest,1)),...
                    ((y-1)*size(frames_stim1_ttest,2)+1):(y*size(frames_stim1_ttest,2))) =...
                    frames_stim1_ttest(:,:,k2) - frames_stim2_ttest(:,:,k2);
                k2 = k2 + 1;
            end
        end
    end
    
    clipMed = median(sumviewImage(:));
    clipStd = std(sumviewImage(:));
    clipLims(1) = clipMed - clip*clipStd;
    clipLims(2) = clipMed + clip*clipStd;
    
    sumviewImage = OIClipH2(imresize(sumviewImage,1/4), 9, clipLims, []);
    fig_sumview = figure('color',[1 1 1],'position',[250 100 1500 850]);
    imagesc(sumviewImage); 
    axis image off;
    colormap gray;
    title('Sumview Image');
    disp('');
    
    %% update with new data
    if CONTINUOUS_MODE
        pause(5);
        
        while 1==1
            
            NBlocksOld = NBlocks;
            cd(dataPath);
            blkFiles = dir('*.BLK');
            for n = 1:length(blkFiles)
                anapar = OIHeadRead((blkFiles(n).name), 'v');
                conds = [conds, anapar.Cond];
            end
            blockConds1 = blkFiles(conds==stim1);
            blockConds2 = blkFiles(conds==stim2);
            NBlocks = min(length(blockConds1),length(blockConds2));
            
            if NBlocks > NBlocksOld
                disp(['Updating plots with ',num2str(NBlocks-NBlocksOld),' new blocks...']);

                for n = (NBlocksOld+1):NBlocks
                    frames_stim1_temp = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],0,'v');
                    frames_stim2_temp = OIReadStim([blockConds2(n).folder,'\',blockConds2(n).name],0,'v');

                    frames_stim1_temp = (frames_stim1_temp - frames_stim1_temp(:,:,fframe)) ./...
                        frames_stim1_temp(:,:,fframe) * 100;
                    frames_stim2_temp = (frames_stim2_temp - frames_stim2_temp(:,:,fframe)) ./...
                        frames_stim2_temp(:,:,fframe) * 100;

                    frames_stim1_avg = frames_stim1_avg + frames_stim1_temp;
                    frames_stim2_avg = frames_stim2_avg + frames_stim2_temp;

                    frames_stim1_ttest(:,:,n) = imageFilter_LPHP_nsc15(imresize(mean(frames_stim1_temp(:,:,frameRange),3),0.5),...
                        round(LP/2),round(HP/2),[]);
                    frames_stim2_ttest(:,:,n) = imageFilter_LPHP_nsc15(imresize(mean(frames_stim2_temp(:,:,frameRange),3),0.5),...
                        round(LP/2),round(HP/2),[]);
                end

                imaging_ttest(frames_stim1_ttest,frames_stim2_ttest,pValue);
                imaging_blockview(frames_stim1_avg,frames_stim2_avg,LP,HP,clip);
                suptitle(['Stim',num2str(stim1),'-Stim',num2str(stim2),'     LP',num2str(LP),...
                    ' HP',num2str(HP),' C',num2str(clip),'     Blocks:',num2str(length(blkFiles))]);

                disp('Plots updated.');
            end

            disp('Waiting for new blocks...');
            pause(15); %wait 15 seconds and recheck
        end

    end
    
    
    cd(parent);
    
    %% ttest function
    function [suninImage] = imaging_ttest(ttest_stim1,ttest_stim2,pValue)
        [~,p] = ttest2(ttest_stim1,ttest_stim2,0.05,'left','equal',3);
        
        suninImage = mean(ttest_stim1,3)-mean(ttest_stim2,3);
        suninImage = OIClipH(suninImage, 1, 2, []);
        
        axes(ax_sunin);
        hold off;
        imagesc(suninImage);
        axis image off;
        colormap gray;
        title(['stim#',num2str(stim1),' - stim#',num2str(stim2),', F',...
            num2str(frameRange(1)),'-',num2str(frameRange(end)),', ',...
            num2str(length(blockSelect)),' Blocks']);
        
        axes(ax_ttest);
        hold off;
        imagesc(suninImage);
        hold on;
        imagesc(ones(size(p,1),size(p,2),3).*cat(3,1,0,0),'alphadata',p<pValue);
        axis image off;
        colormap gray;
%         title(['p < ',num2str(pValue)]);
    end
    
    %% blockview function
    function [uniClipLim] = imaging_blockview(blockview_stim1,blockview_stim2,LP,HP,clip)
        numCols = 5;
        numRows = ceil(NFrames/numCols);
        
        big_map_uni = double(zeros(size(blockview_stim1,1)*numRows,size(blockview_stim1,2)*numCols));
        
        k = 1;
        for r = 1:numRows
            for c = 1:numCols
                if k <= NFrames
                    big_map_uni((r-1)*H+1:r*H,(c-1)*W+1:c*W) = blockview_stim1(:,:,k)-blockview_stim2(:,:,k);
                    k = k + 1;
                end
            end
        end
        
        big_map_uni = imresize(big_map_uni,0.25);
        big_map_uni = imageFilter_LPHP_nsc15(big_map_uni,round(LP/4),round(HP/4),big_map_uni~=0);
        
        clipMed = median(big_map_uni(:));
        clipStd = std(big_map_uni(:));
        uniClipLim = [clipMed-clipStd*clip clipMed+clipStd*clip];
        big_map_uni = OIClipH2(big_map_uni, 9, uniClipLim, []);
        
        axes(ax_blockview);
        hold off;
        imagesc(big_map_uni);
        axis image off;
        colormap gray;
        
        warning off;
        title(dataPath);
        warning on;
    end


    function [fig_SP] = pointPickFcn(pickPoint,event,S)
        [x, y] = getpts(ax_sunin);
        x = round(x);
        y = round(y);
        
        if ~exist('fig_SP')
            fig_SP = figure('color',[1 1 1]);
            ax_SP = axes;
        end
        
        cla(ax_SP);
        axes(ax_SP);
        colors = distinguishable_colors(length(x));
        
        for i = 1:length(x)
            subplot(3,1,1);
            plot(squeeze(frames_stim2_avg(x(i),y(i),:))/NBlocks,'color',colors(i,:));
            ylim([-0.3 0.15]);
            hold on;
            title('stim2');
            
            subplot(3,1,2);
            plot(squeeze(frames_stim1_avg(x(i),y(i),:))/NBlocks,'color',colors(i,:));
            ylim([-0.3 0.15]);
            hold on;
            title('stim1');
            
            subplot(3,1,3);
            plot((squeeze(frames_stim1_avg(x(i),y(i),:)) - ...
                squeeze(frames_stim2_avg(x(i),y(i),:)))/NBlocks,'color',colors(i,:));
            ylim([-0.3 0.15]);
            hold on;
            title('stim1 - stim2');
            xlabel('Frame #');
        end
        
        cla(ax_sunin);
        axes(ax_sunin);
        imagesc(suninImage);
        axis image off;
        colormap gray;
        hold on;
        for i = 1:length(x)
            plot(x(i),y(i),'.','markersize',20,'color',colors(i,:));
        end
    end
    
end