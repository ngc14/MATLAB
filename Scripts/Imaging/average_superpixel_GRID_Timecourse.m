close all;
clear all;
monkey = 'Gilligan';
areaNames = {'M1','PMd'};
MMNames = {'Arm', 'Hand','Face', 'Trunk'};
if strcmp(monkey, 'Gilligan')
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    
    MM = double(imread(['S:\Lab\Gilligan\Mapping\', 'Motor Maps V3\', 'MM_Simp_RGB-01.png']));
    activationFilenames = {'ESS_HSV.png', 'LS_HSV.png', 'PC_HSV.png', 'Rest_HSV.png'};

    %border is yellow
    M1border = MM(:,:,1)>200 & MM(:,:,2) > 200 & MM(:,:,3) < 100;
    PMborder = MM(:,:,1)>100 & MM(:,:,2) < 100 & MM(:,:,3) > 150;
    
    [row,col] = find(M1border);
    pointsX = unique(row);
    pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
    PMpointsX = pointsX;
    PMpointsY = pointsY;
    
    pointsX(end+1:end+2) = [768, 1];
    pointsY(end+1:end+2) = [1 1];
    M1mask = double(poly2mask(pointsY,pointsX, 768,768));
    M1mask(M1mask==0) = NaN;
    M1mask = repmat(M1mask,[1,1,3]);
    
    PMpointsX(end+1:end+2) = [768, 1];
    PMpointsY(end+1:end+2) = [768 768];
    PMmask = double(poly2mask(PMpointsY,PMpointsX, 768,768));
    
    [row,col] = find(PMborder);
    row(end+1:end+2) = [1, 1];
    col(end+1:end+2) = [768, min(col)];
    PMdmask = double(poly2mask(col,row,768,768));
    PMdmask = double(PMmask & PMdmask);
    PMdmask(PMdmask==0) = NaN;
    PMdmask = repmat(PMdmask,[1,1,3]);
    
    areaMasks = {M1mask, PMdmask};
    
    %arm is red
    armArea = MM(:,:,1) > 100 & MM(:,:,2) < 100 & MM(:,:,3) < 100;
    %hand is green
    handArea = MM(:,:,2) > 100 & MM (:,:,1) < 100 & MM(:,:,3) < 100;
    %face is blue
    faceArea = MM(:,:,3) > 100 & MM (:,:,1) < 100 & MM(:,:,2) < 100;
    %trunk is black
    trunkArea = MM(:,:,1) < 100 & MM (:,:,2) < 100 & MM(:,:,3) < 100;
    
    allMMs = {armArea, handArea, faceArea, trunkArea};
    allMMs = cellfun(@(a) conv2(a, ones(5)/5^2, 'same')>.5, allMMs, 'UniformOutput', false);
    allMMs = cellfun(@(a) imfilter(a, ones(15)), allMMs, 'UniformOutput', false);
    
    armArea = allMMs{1};
    handArea = allMMs{2};
    faceArea = allMMs{3};
    trunkArea = allMMs{4};
else
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0]};
    %runs = {[0], [1]};
end
stimNames = {'ExtraSmallSphere','LargeSphere', 'Photocell'};

refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'], 'bmp')>200;
consecutiveFrames = 7;
windowSmoothSize = 5;
stdThresh = 0;
areaThresh = 0.25;

colors = {'r', 'g', 'b', 'k'};
for c = 1:length(stimNames)+1
    toPlot = {};
    for d = 1:length(dates)
        if(strcmp(monkey,'Gilligan'))
            sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP550C0\SPData_Best.mat']);
        else
            sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP250C0\SPData_Best.mat']);
        end
        sp = sp.SPData;
        if(strcmp(monkey, 'Skipper') && c<length(stimNames)+1 && ((d==1 && c>1) || (d==2)))
            if(d==1)
                toPlot{d} = NaN(size(sp.maskMean_stim_Avg2{1}));
            elseif(d==2)
                if(c==1)
                    toPlot{d} = NaN(size(sp.maskMean_stim_Avg2{1}));
                else
                    toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c-1},2,'movmean',windowSmoothSize);
                end
            end
        elseif(c==length(stimNames)+1)
            toPlot{d} = smoothdata(sp.maskMean_blank_Avg2{1},2,'movmean',windowSmoothSize);
        else
            %toPlot{d} = sp.maskMeanAll{c};
            toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c},2,'movmean',windowSmoothSize);
        end
        if(strcmp(monkey, 'Gilligan') && d == length(dates))
            toPlot{d} = toPlot{d}(:,1:70);
        end
    end
    allMasks = reshape(cell2mat(toPlot),[size(toPlot{1},1), size(toPlot{1},2), length(toPlot)]);
    subplotSize = numSubplots(length(areaMasks)*length(allMMs));
    masks = sp.mask;
    for m = 1:length(areaMasks)*length(allMMs)
        currMask = areaMasks{floor(m/(length(allMMs)+1))+1};
        currMask(isnan(currMask)) = 0;
        
        currInd = m-floor(m/(length(allMMs)+1))*length(allMMs);
        MMmap = allMMs{currInd};
        others = allMMs(setdiff(1:end,currInd));
        exclusiveMaps = sum(cat(3,others{:}),3)>0;
        exclusiveMaps = zeros(size(MMmap));
         activation = imresize(imread(['S:\Lab\', monkey,'\Mapping\tTests\',activationFilenames{c}]),[768,768]);
        %activation with 50% overlap is red
        activation = activation(:,:,1) > 100 & activation(:,:,2) < 100 & activation (:,:,3) < 100;
        mask = currMask(:,:,1) & MMmap & ~activation;%~exclusiveMaps;
        
        overlapInds = cellfun(@(a) sum(sum(a&mask))>sum(sum(masks{1}))*.5, masks);
        currMasks = allMasks(overlapInds,:,:);
        timecourse = nanmean(nanmean(currMasks,3),1);
        subplot(subplotSize(1), subplotSize(2),m);
        hold on;
        plot(timecourse, colors{c});
        if(c==3)
            title([areaNames{floor(m/(length(allMMs)+1))+1},': ',MMNames{currInd}]);
        end
    end
    
end