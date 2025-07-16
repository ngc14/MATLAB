allConds = {'[ExtraSmallSphere]', '[LargeSphere]', '[Photocell]', '[Rest]'};
monkey = 'Gilligan';
dirPath = ['S:\Lab\',monkey,'\Mapping\Intensity Maps\'];
refMask = imread(['S:\Lab\', monkey,'\Mapping\clean_mask_filled'], 'bmp')<250;
refMask = refMask(:,:,1);
savePath = ['S:\Lab\ngc14\Working\',monkey,'\'];
framesSmooth = 3;
histBins = 50;
threshPixels = '';
colorCond = {'r','g','b','k'};
numFrames=70;
blankFrames = 1:23;
HPk = 250;
LPk=5;
H = size(refMask,1);
W = size(refMask,2);
if(strcmp(monkey,"Skipper"))
    mm = MotorMapping(56);
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
else
    mm = MotorMapping(42);
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
end
[datesTable, masksR, ~] = getMonkeyInfo('S:\Lab\',monkey,["M1", "PMd"],false);
refMask = masks{1};
siteMask = getVoronoiMask(datesTable,mm,refMask,["Arm","Hand"]);
%%
close all
delete(gcp('nocreate'))
parpool(length(dates));
for cn = 1:length(allConds)
    saveCond = [savePath,allConds{cn},'\Trial_Peak\'];
    if(~exist(saveCond,'dir'))
        mkdir(saveCond);
    end
    for d = 1:length(dates)
        runFolder = ['S:\Lab\',monkey,'\All Data\', monkey, '_',dates{d}, '\Imaging\run0', num2str(runs{d})];
        tform = matfile([runFolder,'\tform_nonrigid.mat']);
        tMat = tform.O_trans;
        spacing = tform.Spacing;
        tform = [];
        frameFile = [runFolder,'\Results\Blockview\',allConds{cn},'\FF_1\','LP',...
            num2str(LPk),'_HP',num2str(HPk),'\Align_Offset\noClip\',...
            'NaN_filtered_HP', num2str(HPk),'_LP',num2str(LPk),'.h5'];
        if(exist(frameFile,'file'))
            hFinfo = h5info(frameFile);
            sessionSz = hFinfo.Datasets.Dataspace.Size;
            hSize{d} = sessionSz;
            numTrials = sessionSz(end);
            sessionDS = {};
            for t = 1:numTrials
                sessionDS{t} = transform(fileDatastore(frameFile,"ReadFcn",@(x) ...
                    squeeze(h5read(x,'/allTrials',[1,1,1,t],[sessionSz(1:3),1]))),...
                    @(x) arrayfun(@(n) bspline_transform(tMat,x(:,:,n),spacing,3),1:size(x,3), 'UniformOutput', false));
            end
            sessionTrials{d} = combine(sessionDS{:});
        end
    end
    sessionTrials = sessionTrials(~cellfun(@isempty,sessionTrials));
    for s = 1:length(sessionTrials)
        df = datetime('now');
        mvs = {};
        mts = {};
        currSession = sessionTrials{s}.UnderlyingDatastores;
        parfor t = 1:size(currSession,2)
            smoothFrames = movmean(cell2mat(reshape(readall(currSession{t}),1,1,[])),framesSmooth,3,'omitnan','EndPoints','fill');
            smoothFrames(:,:,blankFrames) = NaN;
            [minVal,minFrame] = min(smoothFrames,[],3);
            mvs{t} = minVal;
            mfs{t} = minFrame;
        end
        sessionFrames{s} = tall(mfs);
        sessionVals{s} = tall(mvs);
        disp(['Session ', num2str(s), ' loaded (',num2str(seconds(duration(datetime('now')-df))),')']);
    end
    sessionFrames = cellfun(@gather, sessionFrames, 'UniformOutput', false);
    sessionVals = cellfun(@gather, sessionVals, 'UniformOutput', false);
    minFrame = mean(cell2mat(reshape([sessionFrames{:}],1,1,[])),3);
    minVal = mean(cell2mat(reshape([sessionVals{:}],1,1,[])),3);

    figure();
    histogram(minFrame(:),1:numFrames+1);
    saveas(gcf,[saveCond,'Hist_Frame',threshPixels,'.png']);

    figure();
    imagesc(minFrame);
    hold on;
    colormap([colormap(jet(numFrames));0 0 0;.65 .65 .65;]);
    colorbar
    [rI,cI] = ind2sub(size(refMask),find(imdilate(bwperim(~refMask & siteMask),ones(3,3))));
    scatter(cI,rI,'g.');
    caxis([blankFrames(end)+1 numFrames+2])
    minFrame(isnan(minFrame)|refMask|~siteMask|isinf(minFrame)) = NaN;
    title(strcat("Median Frame = ", num2str(median(minFrame(:),'omitnan')), "; Mean Frame = ", num2str(round(mean(minFrame(:),'omitnan')))));
    saveas(gcf,[saveCond,'Min_Frame',threshPixels,'.png']);

    figure();
    histogram(minVal(:));
    saveas(gcf,[saveCond,'Hist_Val',threshPixels,'.png']);
    xlim([-0.2 0.05]);

    figure();
    minVal(isnan(minVal))  = -Inf;
    minVal(refMask) = Inf;
    im = imagesc(minVal);
    cm = colormap(hot(histBins));
    hold on;
    colormap([0 0 0; cm(round(length(cm)/5):end,:); .65 .65 .65;]);
    caxis([-0.15 0.05]);
    colorbar
    scatter(cI,rI,'g.');
    minVal(isnan(minVal)|refMask|~siteMask|isinf(minVal)) = NaN;
    title(strcat("Median RChange % = ", num2str(median(minVal(:),'omitnan')), "; Mean RChange % = ", num2str(mean(minVal(:),'omitnan'))));
    saveas(gcf,[saveCond,'Min_Val',threshPixels,'.png']);
end