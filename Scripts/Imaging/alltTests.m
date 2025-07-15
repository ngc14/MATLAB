monkey = 'Gilligan';
if(strcmp(monkey,'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    mm = MotorMapping(42);

elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    mm = MotorMapping(56);
end
pvs = [0.0001, .000000001] ;
LPk = 5;
HPk = 250;
conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
baselineFrames = 1:10;
lastFrame = {32:70};%num2cell(1:70);
saveDir = ['S:\Lab\ngc14\Working\',monkey,'\'];
refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'],'bmp')>180;
clipLim = [-0.3 0.3];
cMap = flipud(jet(255));
cleanArea = 10;
%%
for cn =1:length(conds)
    currCond = conds{cn};
    saveFolder = strcat(saveDir, string(currCond),'\tTests\NaN\HP', num2str(HPk),'\');
    if(~exist(saveFolder,'dir'))
        mkdir(saveFolder);
    end
    ci = datetime('now');
    baseline = {};
    for d = 1:length(dates)
        bPath = ['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
            '\Imaging\run0',num2str(runs{d}),'\Results\Blockview\',currCond,...
            '\FF_1\LP',num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
        if(exist([bPath, 'NaN_filtered_HP',num2str(HPk),'_LP',num2str(LPk),'.h5'],'file'))
            baseline{d} = squeeze(mean(h5read([bPath,'NaN_filtered_HP',...
                num2str(HPk),'_LP',num2str(LPk),'.h5'],'/allTrials',...
                [1 1 baselineFrames(1) 1],[size(refMask,1),size(refMask,2),...
                length(baselineFrames),Inf]),3,'omitnan'));
            tform = matfile(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},...
                '\Imaging\run0',num2str(runs{d}),'\tform_nonrigid.mat']);
            Otrans{d} = tform.O_trans;
            spacing{d} = tform.Spacing;
            tform = [];
            for t = 1:size(baseline{d},3)
                baseline{d}(:,:,t) = bspline_transform(Otrans{d},baseline{d}(:,:,t), spacing{d},3);
            end
        end
    end
    baseline = tall(cat(3,baseline{:}));
    for f = 1:length(lastFrame)
        ims = {};
        startFrame = datetime('now');
        currFrame = lastFrame{f};
        for (d = 1:length(dates))
            runFolder = ['S:\Lab\',monkey,'\All Data\',monkey, '_',dates{d},...
                '\Imaging\run0', num2str(runs{d})];
            condPath = [runFolder,'\Results\Blockview\',currCond,'\FF_1\','LP',...
                num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\'];
            if(exist([condPath,'NaN_filtered_HP',num2str(HPk),'_LP',num2str(LPk),'.h5'],'file'))
                imR = squeeze(mean(h5read([condPath,'NaN_filtered_HP',...
                    num2str(HPk),'_LP',num2str(LPk),'.h5'],'/allTrials',...
                    [1 1 currFrame(1) 1],[size(refMask,1),size(refMask,2),...
                    length(currFrame),Inf]),3,'omitnan'));
                for t = 1:size(imR,3)
                    ims{d}(:,:,t) = bspline_transform(Otrans{d},imR(:,:,t), spacing{d},3);
                end
            end
        end
        ims = cat(3,ims{:});
        disp(['Frame ', num2str(f), ' loaded (', num2str(time2num(datetime('now')-startFrame,'seconds')),' s)'])
        [~, p] = ttest2(gather(baseline),ims,'alpha',0.05,'tail','right','dim',3);
        figure();
        imagesc(p)
        colormap(flipud(colormap("jet")));
        caxis([0 max(pvs)]);
        colorbar;
        saveas(gcf,strjoin([saveFolder,'tTest_Nomask',num2str(lastFrame{f}),'_p.bmp'],''));

        for v = 1:length(pvs)
            p_value = pvs(v);
            cleanedP = im2double(bwareaopen(imfill((p<p_value),'holes'),cleanArea));
%                 frameSA(f) = sum(cleanedP & forelimbMask,'all')/sum(forelimbMask,'all');
                imwrite(im2double(bwareaopen(imfill(p<p_value,'holes'),cleanArea)),...
                    strjoin([saveFolder,'tTest_Nomask',num2str(lastFrame{f}),...
                    '_p-',num2str(p_value),'.bmp'],''));
        end
        p = [];
        figure();
        mIm = mean(ims,3,'omitnan');
        clippedIntensity = OIClipH2(mIm,9,clipLim, []);
        imagesc(clippedIntensity);
        colormap(cMap);
        imwrite(ind2rgb(im2uint8(mat2gray(clippedIntensity,clipLim)),cMap),cMap,strjoin([saveFolder,'Intensity_',num2str(lastFrame{f}),'.png'],''));
        save(strjoin([saveFolder,'Intensity_',num2str(lastFrame{f}),'.mat'],''),'mIm','-mat');
        close all;
    end
    disp([currCond,': ',num2str(time2num(datetime('now')-ci,'minutes'))]);
end