%function [] = normalizeAverages(varargin)
% close all;
% clear all;
monkey = 'Skipper';                   % animal name
if(strcmp(monkey,'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
        HPk = 550;

elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
        HPk = 250;
end
saveVar = true;
tTest = false;

pvs = [0.01, 0.001, 0.0001];
clipMethod = 9;
clipValue = 0.3;
LPk = 5;
LPt = 5;
cm = [flipud(jet(254));.5 .5 .5];
conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
saveStartPath = strcat("S:\Lab\",monkey,"\Mapping\Intensity Maps\");

refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'],'bmp')>180;
refMask = refMask(:,:,1);
W = size(refMask,2);
H = size(refMask,1);
saveStartPath = strcat(saveStartPath,"TrialAVGS\");
%%
for cn = 1:length(conds)
    ims = {};
%     if(tTest)
%         if(cn==1)
%             if(strcmp(monkey,'Gilligan'))
%                 loadFrames = {1:10,49:53};
%             else
%                 loadFrames = {1:10,46:50};
%             end
%         elseif(cn==2)
%             if(strcmp(monkey,'Gilligan'))
%                 loadFrames = {1:10,56:60};
%             else
%                 loadFrames = {1:10,58:62};
%             end
%         elseif(cn==3)
%             if(strcmp(monkey,'Gilligan'))
%                 loadFrames = {1:10,46:50};
%             else
%                 loadFrames = {1:10,58:62};
%             end
%         elseif(cn==4)
%             if(strcmp(monkey,'Gilligan'))
%                 loadFrames = {1:10,50:54};
%             else
%                 loadFrames = {1:10,54:58};
%             end
%         end
        if(strcmp(monkey,'Gilligan'))
            loadFrames = num2cell(1:70);
        else
            loadFrames =  num2cell(1:70)
        end
%     else
%         loadFrames = num2cell(1:70);
%     end
    currCond = conds{cn};
    for f = 1:length(loadFrames)
        imT = {};
        currF = loadFrames{f};
        startFrame = clock;
        parfor d = 1:length(dates)
            framePointers = [];
            sessionPath = ['S:\Lab\',monkey,'\All Data\', monkey, '_',dates{d}, '\Imaging\'];
            runFolder = [sessionPath,'run0', num2str(runs{d})];
            savePath = [runFolder,'\Results\Blockview\',currCond,'\FF_1\','LP',...
                num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
            if(exist([savePath,'filtered_HP', num2str(HPk),'_LP',num2str(LPk),'.mat'],'file'))
                tform = matfile([runFolder,'\tform_nonrigid.mat']);
                m = matfile([savePath,'filtered_HP', num2str(HPk),'_LP',num2str(LPk),'.mat']);
                framePointers = squeeze(cell2mat(m.allTrials(1,currF)));
                for t = 1:size(framePointers,3)
                    framePointers(:,:,t) = bspline_transform(tform.O_trans,...
                        framePointers(:,:,t), tform.Spacing,3);
                end
                imT{d} = framePointers;
            end
        end
        clear framePointers
        disp([num2str(currF),': ',num2str(etime(clock,startFrame))]);
        %         if(trialsAvg)
        ims{f} = tall(squeeze(cat(3,imT{:})));
        clear imT
        %         else
        %             ims{f} = cellfun(@(t) nanmean(cat(3,t{:}),3), num2cell([imT{:}],2), 'UniformOutput',false);
        %         end
    end
        mapreducer(0);
    saveFolder = strcat("S:\Lab\ngc14\Meeting\",monkey,"\",string(conds{cn}),'_All\');
    if(~exist(saveFolder,'dir'))
        mkdir(saveFolder);
    end
    %%
    if(tTest)
        %%
        bFrames = ims(1:10);
        baselineFrames = gather(nanmean(cat(4,bFrames{:}),4));
        for f = 1:length(ims)
            for v = 1:length(pvs)
                p_value = pvs(v);
                [~, p] = ttest2(baselineFrames,gather(ims{f}),'alpha',0.05,'tail','right','dim',3);
                cleanedP = imgaussfilt(im2double(bwareaopen(imfill((p<p_value & refMask),'holes'),10)),round(LPt/2),'FilterSize', LPt+(mod(LPt,2)==0));
                if(saveVar)
                    imwrite(cleanedP & refMask,strjoin([saveFolder,'tTest_',num2str(f),'_p-',num2str(p_value),'.bmp'],''));
                    imwrite(cleanedP & refMask,strjoin([saveFolder,'tTest_',num2str(f),'_p-',num2str(p_value),'.png'],''));
                    save(strjoin([saveFolder,'tTest',num2str(f),'_p-',num2str(p_value),'.mat'],''),'cleanedP');
                end
            end
            disp([num2str(f),': ',num2str(gather(size(ims{f})))]);
        end
    else
        parfor f = 1:length(loadFrames)
            filtIms{f} = imageFilter_LPHP_nsc15(nanmean(cat(3,ims{f}{:}),3),LPk,HPk,[]);
        end
        filtIms = cat(3,filtIms{:});
        filtIms(isnan(filtIms)) = 0;
        filtIms(isinf(filtIms)) = 0;
        filtMasked = [];
        
        for k = 1:size(filtIms,3)
            temp = filtIms(:,:,k);
            filtMasked(:,k) = temp(refMask);
        end
        numCols = 10;
        numRows = ceil(size(filtIms,3)/numCols);
        
        maps_avg_median_uni = nanmedian(filtMasked(:));
        maps_avg_std_uni = nanstd(filtMasked(:));
        uniLowClip = maps_avg_median_uni-(maps_avg_std_uni*clipValue);
        uniHighClip = maps_avg_median_uni+(maps_avg_std_uni*clipValue);
        big_map_indi = uint8(zeros(H*numRows,W*numCols));
        k=1;
        for r = 1:numRows
            for c = 1:numCols
                if k <= size(filtIms,3)
                    imClip = rescale(filtIms(:,:,k),1,length(cm)-1,'InputMin',-clipValue ,'InputMax',clipValue);
                    imClip(refMask==0) = length(cm);
                    imClip = im2uint8(imClip,'indexed');
                    big_map_indi((r-1)*H+1:r*H,(c-1)*W+1:c*W) = imClip;
                    if(saveVar)
                        if(~exist(strjoin([saveFolder,'Frames\'],''),'dir'))
                            mkdir(strjoin([saveFolder,'Frames\'],''));
                        end
                        imwrite(imClip,cm,strcat(saveFolder,'Frames\', num2str(k),'_Clip_',num2str(clipValue,'%0.2f'),'.png'));
                        saveIm = filtIms(:,:,k);
                        save(strcat(saveFolder,'Frames\', num2str(k),'.mat'),'saveIm');
                    end
                    k = k + 1;
                end
            end
        end
        big_map_indi_small = imresize(big_map_indi,0.25);
        
        if(saveVar)
            imwrite(big_map_indi_small,[gray(254); .5 .5 .5],strcat(saveFolder,'AllFrames_Clip_',num2str(uniHighClip,'%0.2f'),'.bmp'));
            imwrite(big_map_indi_small,cm,strcat(saveFolder,'AllFrames_Color_Clip_',num2str(uniHighClip,'%0.2f'),'.bmp'));
        end
    end
    
end