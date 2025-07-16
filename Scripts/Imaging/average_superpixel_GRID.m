close all;
clear all;
monkey = 'Gilligan';
stimNames = {'ExtraSmallSphere','LargeSphere', 'Photocell'};
parentPath = ['S:\Lab\', monkey,'\Mapping\'];
activationFilenames = {'ESS_HSV.png', 'LS_HSV.png', 'PC_HSV.png'};

windowSmoothSize = 5;
maskSpacing = 15;
saveVar = 1;
negative = 0;
paramVals  = [1, 1.1, 1.2, 1.3, 1.4, 1.5];
errPlot = [];

if strcmp(monkey, 'Gilligan')
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    %dates = {'12_07_2018', '12_12_2018', '12_13_2018', '12_14_2018'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    %runs = {[0], [0], [0], [1]};
    
    MM = double(imread(['S:\Lab\Gilligan\Mapping\', 'Motor Maps V3\', 'MM_Simp_RGB.png']));
    
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
    
    allMMsLabs = ["Arm", "Hand", "Face", "Trunk"];
    allMMs = {armArea, handArea, faceArea, trunkArea};
    allMMs = cellfun(@(a) conv2(a, ones(5)/5^2, 'same')>.5, allMMs, 'UniformOutput', false);
    allMMs = cellfun(@(a) imfilter(a, ones(15)), allMMs, 'UniformOutput', false);
    clear armArea handArea faceArea trunkArea PMdmask M1mask M1border PMborder MM PMmask points* PMpoints*
else
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0]};
end

refMask = imread(['S:\Lab\', monkey, '\Mapping\clean_mask_filled.bmp']);
refMask = refMask(:,:,1)>50;
%refMask = (refMask & (armArea | handArea));

refMask(:,1) = 0;
refMask(:,end) = 0;
refMask(1,:) = 0;
refMask(end-4:end,:) = 0;
count = 1;
H = size(refMask,1);
W = size(refMask,2);


% activation = imresize(imread([parentPath,'tTests\',activationFilenames{c}]),[H,W]);
% activation = activation(:,:,1) > 100 & activation(:,:,2) < 100 & activation (:,:,3) < 100;
%activation = imread([parentPath, 'M1_Forelimb_Estimate.bmp']);
%refMask = activation(:,:,1) > 100 & refMask;

for m = 1:ceil(H/maskSpacing)
    for n = 1:ceil(W/maskSpacing)
        masks{count} = zeros(H,W);
        masks{count}(((m-1) * maskSpacing) + 1:min((m * maskSpacing) + 1,H) , ((n-1) * maskSpacing) + 1:min((n * maskSpacing) + 1 ,W)) = 1;
        masks{count}(refMask==0) = 0;
        count = count + 1;
    end
end
badMasks = cellfun(@(a) any(any(a)), masks);
badMasks = ~badMasks;
masks(badMasks) = [];
if(strcmp(monkey, 'Skipper'))
%masks(end-373:end) = [];
end

%%
for c = 1:length(stimNames)+1
    toPlot = {};
    for d = 1:length(dates)
        if(strcmp(monkey,'Gilligan'))
            if(maskSpacing==15)
                sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP550C0\SPData_GridMasks_SMALL.mat']);
            else
                sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP550C0\SPData_GridMasks.mat']);
            end
            
        else
            sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP250C0\SPData_GridMasks.mat']);
        end
        sp = sp.SPData;
        
        
        sp.maskMean_stim_Avg2 = cellfun(@(a) a(~badMasks,:), sp.maskMean_stim_Avg2, 'UniformOutput', false);
        sp.maskMean_blank_Avg2 = cellfun(@(a) a(~badMasks,:), sp.maskMean_blank_Avg2, 'UniformOutput', false);
        
%         allMasksMeans = vertcat(sp.maskMean_stim_Avg2{:});
%         allMasksMeans = [allMasksMeans; sp.maskMean_blank_Avg2{1}];
%         meanVals = nanmean(allMasksMeans(:));
%         stdVals = nanstd(allMasksMeans(:));
%         median_uni = nanmedian(allMasksMeans(:));
%         std_uni = nanstd(allMasksMeans(:));
%         
%        for z = 1:length(sp.maskMean_stim_Avg2)
%            sp.maskMean_stim_Avg2{z} = (sp.maskMean_stim_Avg2{z} - meanVals) ./ stdVals;
%            sp.maskMean_blank_Avg2{z} = (sp.maskMean_blank_Avg2{z} - meanVals) ./ stdVals;
%        end
       
%        allMaskMeans = [];
%        allMasksMeans = vertcat(sp.maskMean_stim_Avg2{:});
%        allMasksMeans = [allMasksMeans; sp.maskMean_blank_Avg2{1}];
%        median_uni = nanmedian(allMasksMeans(:));
%        std_uni = nanstd(allMasksMeans(:));
%         for z = 1:length(sp.maskMean_stim_Avg2)
%             for f = 1:size(sp.maskMean_stim_Avg2{z},1)
%                 sp.maskMean_stim_Avg2{z}(f,:) = OIClipH2(sp.maskMean_stim_Avg2{z}(f,:),...
%                     9, [median_uni-(2*std_uni), median_uni + (2*std_uni)], []);
%                 sp.maskMean_blank_Avg2{z}(f,:) = OIClipH2(sp.maskMean_blank_Avg2{z}(f,:),...
%                     9, [median_uni-(2*std_uni), median_uni + (2*std_uni)], []);
%             end
%         end
        
        if(strcmp(monkey, 'Skipper') && c<length(stimNames)+1 && ((d==1 && c>1) || (d==2)))
            if(d==1)
                toPlot{d} = NaN(size(sp.maskMean_blank_Avg2{1}));
            elseif(d==2)
                if(c==1)
                    toPlot{d} = NaN(size(sp.maskMean_blank_Avg2{1}));
                else
                    toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c-1},2,'movmean',windowSmoothSize);
                end
            end
        elseif(c==length(stimNames)+1)
            toPlot{d} = smoothdata(sp.maskMean_blank_Avg2{1},2,'movmean',windowSmoothSize);
        else
            toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c},2,'movmean',windowSmoothSize);
        end
        %masks = sp.mask{1};
        clear sp;
    end
    
    if(strcmp(monkey,'Gilligan'))
        toPlot{end} = toPlot{end}(:,1:70);
    end
    badMasksR = cellfun(@(a) double(any(isnan(a),2)), toPlot, 'UniformOutput', false);
    if(strcmp(monkey, 'Skipper'))
        if (c==1)
            badMasksR = badMasksR{1} + badMasksR{3} + badMasksR{4} + badMasksR{5} + badMasksR{6} + badMasksR{7} + badMasksR{8} + badMasksR{9}>1;
        elseif (c==4)
            badMasksR = badMasksR{1} + badMasksR{2} + badMasksR{3} + badMasksR{4} + badMasksR{5} + badMasksR{6} + badMasksR{7} + badMasksR{8} + badMasksR{9}>1;
        else
            badMasksR = badMasksR{2} + badMasksR{3} + badMasksR{4} + badMasksR{5} + badMasksR{6} + badMasksR{7} + badMasksR{8} + badMasksR{9}>1;
        end
    else
        badMasksR = badMasksR{1} + badMasksR{2} + badMasksR{3} + ...
            badMasksR{4} + badMasksR{5} + badMasksR{6} + badMasksR{7} + badMasksR{8} >1;
    end
    toPlot = cellfun(@(a) a(~badMasksR,:), toPlot, 'UniformOutput', false);
    
    allMasks = reshape(cell2mat(toPlot),[size(toPlot{1},1), size(toPlot{1},2), length(toPlot)]);
    meanMasks = nanmedian(allMasks,3);
    stdMasks = nanstd(allMasks, [], 3);
    
    DG = [];
    %for d = 1:length(dates)
    %     for z = 1:length(meanMasks)
    %         for f = 1:length(meanMasks)
    %             xCorrVals = xcorr(meanMasks(z,:),meanMasks(f,:));
    %             %xCorrVals = xcorr(allMasks(z,:,d),allMasks(f,:,d));
    %             [~, maxMagInd] = max(abs(xCorrVals));
    %             if(negative)
    %                 %DG(z,f,d) = xCorrVals(maxMagInd);
    %                 DG(z,f) = xCorrVals(maxMagInd);
    %             else
    %                 maxVal = max(xcorr(meanMasks(z,:),meanMasks(f,:))) ;
    %                 %maxVal = max(xcorr(allMasks(z,:,d),allMasks(f,:,d)));
    %                 %DG(z,f,d) = max(0, maxVal);
    %                 DG(z,f) = max(0,maxVal);
    %             end
    %         end
    %     end
    %end
    DG = nanmean(DG,3);
    
    %     for m = 1:size(meanMasks,1)
    %         curve = allMasks(m,:,:);
    %         if(any(~isnan(meanMasks(m,:))))
    %          parfor n = 1:50
    % %             for d = 1:5
    %             [~,r] = polyfit(1:70,curve,n);
    %             err(n) = r.normr;
    % %             [nC{n,d},r] = fit([[1:70]', meanMasks(m,:)'],stdMasks(m,:)',['poly',num2str(n),num2str(d)],'Exclude', [13:23]);
    % %             err(n,d) = r.rmse;
    % %             end
    %          end
    % %         [bestfitRVals,bestfitR] = nanmin(err);
    % %         [~,bestfitC] = nanmin(bestfitRVals);
    % %             bestfit(m,1) = bestfitR(bestfitC);
    % %             bestfit(m,2) = bestfitC;
    % %             numCs(m) = numcoeffs(nC{bestfit(m,1),bestfit(m,2)});
    %         if(~any(~isnan(err)) || ~any(err<1))
    %             bestfit(m) = 1;
    %         else
    %             bestfit(m) = find(err<1, 1,'first');
    %         end
    %         disp(m);
    %         else
    %            %bestfit(m,:) = [NaN, NaN];
    %         end
    %     end
    %     bestfit(badMasks) = 1;
    %     largePoly = max(bestfit);
    %     eqFit = zeros(size(meanMasks,1),largePoly);
    %     for m = 1:size(eqFit,1)
    %         if(any(~isnan(meanMasks(m,:))))
    %             %polyCoefs = fit([[1:70]', meanMasks(m,:)'],stdMasks(m,:)',['poly',num2str(bestfit(m,1)),num2str(bestfit(m,2))],'Exclude', [13:23]);
    %             polyCoefs =  polyfit(1:70,meanMasks(m,:),5);
    % %             names = coeffnames(polyCoefs);
    % %             for n = 1:length(names)
    % %                 eqFit(m,n) = polyCoefs.(names{n});
    % %             end
    %         end
    %         eqFit(m,end-length(polyCoefs)+1:end) = polyCoefs;
    %     end
    
    meanMasksF = [];
    
    %             meanMasksM = meanMasks;
    %             %meanMasksM(:,13:23) = NaN;
    %             meanMasksF(:,end+1) = nanmean(meanMasksM');
    %             meanMasksF(:,end+1) = nanstd(meanMasksM');
    %             meanMasksF(:,end+1) = nanmin(meanMasksM');
    %             meanMasksF(:,end+1) = nanmax(meanMasksM');
    %             [~, meanMasksF(:,end+1)]  = nanmin(meanMasksM');
    %             [~, meanMasksF(:,end+1)]  = nanmax(meanMasksM');
    %             meanMasksF(:,end+1) = sum(meanMasksM'<0 & ~isnan(meanMasksM'));
    %             meanMasksF(:,end+1) = sum(meanMasksM'>0 & ~isnan(meanMasksM'));
    %             meanMasksF(:,end+1) = nansum(abs(meanMasksM'.*(meanMasksM'<0)));
    %             meanMasksF(:,end+1) = nansum(abs(meanMasksM'.*(meanMasksM'>0)));
    
    %             meanMasksF(:,end+1:end+size(eqFit,2)) = eqFit;
    %             meanMasksF(:,end+1) = sum(diff(sign(meanMasksM),1,2)~=0,2);
    
    meanMasksF = (meanMasks-repmat(nanmean(meanMasks),[length(allMasks),1]))./repmat(nanstd(meanMasks')',[1,size(meanMasks,2)]);
    meanMasksF = smoothdata(meanMasks,2,'movmean',3);
    
    %      meanMasksF = (allMasks-repmat(nanmean(allMasks,1),[size(allMasks,1),1]))./repmat(nanstd(allMasks,0,1),[size(allMasks,1),1]);
    %      meanMasksF = reshape(meanMasksF,[size(meanMasksF,1),size(meanMasksF,2)*size(meanMasksF,3)]);
    %      meanMasksF(:,end+1:end+70) = (meanMasks-repmat(nanmean(meanMasks),[length(allMasks),1]))./repmat(nanstd(meanMasks')',[1,size(meanMasks,2)]);
    %      [wcoeff,score,latent,~,explained] = pca((meanMasks-repmat(nanmean(meanMasks),[length(allMasks),1]))./repmat(nanstd(meanMasks')',[1,size(meanMasks,2)]));
    %%
    sumD = {};
    inter = {};
    information = {};
    clus = {};
    sol = [];
    parfor cl=1:10
        [idxC, clu,sumM,interM,~,infor]  = kmedoids(meanMasksF,cl,'Distance',@dtwf);
        %[idxC, clu,sumM,interM,]  = kmeans(meanMasksF,cl);
        clus(cl) = {clu};
        sumD(cl) = {sumM};
        inter(cl)= {interM};
        %information(cl) = {infor};
        sol(cl,:) = idxC;
    end
    sumW = cellfun(@(a) sum(a),sumD);
    %     tots = cellfun(@(a) mean(a(:)),inter);
    %     betw = tots-sumW;
    %     exVar = betw./tots;
    %     exVar = exVar.*((size(X,1)-[1:10])./([1:10]-1));
    %     sumD = cellfun(@sum,sumD);
    %     variance=sumD(1:end-1)-sumD(2:end);
    %     distortion_percent=cumsum(variance)/(sumD(1)-sumD(end));
    %     r=find(exVar>0.9,1);
    %     optimalK(z) = r; %find(variance < variance(1)/10)+1
    
    optimalK = triangle_threshold(sumW,'right', 0);
    %%
    clusters = kmedoids(meanMasksF,optimalK,'Distance',@dtwf);
    %     for p = 1:length(paramVals)
    %         if(negative)
    %             clusters = community_louvain(DG,paramVals(p),[], 'negative_asym');
    %         else
    %             clusters = community_louvain(DG,paramVals(p),[], 'modularity');
    %         end
    %         currK = length(unique(clusters));
    %         errMetric = 0;
    %         for m = 1:currK
    %             currGrids = nanmean(allMasks(clusters==m,:,:),3);
    %             errs = nanstd(currGrids,0,1);
    %
    %             if(~any(errs))
    %                 errs = ones(size(errs));
    %             end
    %
    %             errMetric = errMetric + log(sum(errs)) + log(sum(clusters==m)/(length(clusters)/(currK*maskSpacing)));
    %         end
    %         errPlot(p) = errMetric;
    %     end
    %     optimalP = triangle_threshold(errPlot,'right', 0);
    %     if(negative)
    %         clusters = community_louvain(DG,paramVals(optimalP),[], 'negative_asym');
    %     else
    %         clusters = community_louvain(DG,paramVals(optimalP),[], 'modularity');
    %     end
    %     optimalK = length(unique(clusters));
    %[~, minFrames] = min(meanMasks');
    %%
    finalMask = zeros(size(refMask,1),size(refMask,2));

    for m = 1:length(masks)
        %currMask(find(currMask)) = minFrames(m);
        %minVal = find(meanMasks(m,:)<0);
        %         if(any(meanMasks(m,:)<mean(mean(meanMasks))-stdThresh*std(std(meanMasks))))
        %             minVal = diff(meanMasks(m,:));
        %         else
        %             minVal = [];
        %         end
        %         if(isempty(minVal))
        %             currMask(find(currMask))= NaN;
        %         else
        %             diffVals = find(diff(meanMasks(m,:))<0)+1;
        %             areas = cumsum(abs(meanMasks(m,diffVals)));
        %             areas = [1 2 3];
        %             minValInd = strfind(minVal<0, ones(1, consecutiveFrames));
        %             %minValInd = [1 2 3];
        %             if(isempty(minValInd) || isempty(areas))
        %                 currMask(find(currMask)) = NaN;
        %             else
        %                 %minValInd = diffVals(areas > areaThresh*areas(end));
        %                 currMask(find(currMask)) = minValInd(1)+1;
        %             end
        %         end
        if(~ismember(m,find(badMasks)))
            [rs,cs] = find(masks{m}>0);
            finalMask(rs,cs) = clusters(m);
        end
    end
    spIm{c}=finalMask;
    
    figure();
    h = imshow(uint8(spIm{c}));
    hold on;
    set(h,'alphadata',~isnan(spIm{c}))
    cmap = colormap(flipud(jet(optimalK)));
    colormap([1 1 1 ;cmap]);
    caxis([0,optimalK]);
    h = imagesc(refMask);
    set(h, 'alphadata', refMask(:,:,1)==0);
    savePath = ['S:\Lab\',monkey,'\Mapping\Clustering\15px_Grid\XCorr_Elbow\Biased\No_Norm_No_Clip\'];
     if(saveVar)
        if(c==length(stimNames)+1)
            saveas(gcf, [savePath, 'Rest_Grid.fig']);
            ts = getframe(gca);
            imwrite(ts.cdata,[savePath, 'Rest_Grid.png']);
        else
            title(stimNames{c});
            saveas(gcf, [savePath, stimNames{c},'_Grid.fig']);
            ts = getframe(gca);
            imwrite(ts.cdata,[savePath,stimNames{c},'_Grid.png']);
        end
        for z = 1:optimalK
            figure();
            imshow(spIm{c}==z);
            if (c > length(stimNames))
                export_fig(gcf,[savePath,'Rest_',num2str(z),'.bmp']);
            else
                export_fig(gcf,[savePath,stimNames{c},'_',num2str(z),'.bmp']);
            end
            pause(1);
        end
    end
    figure();
    subplotSize = numSubplots(optimalK);
    colorP = colormap(flipud(jet(optimalK)));
    %     fID = fopen([savePath,'\errors.txt'],'a');
    %     fprintf( fID, 'Parameter Value %f \n', paramVals(p));
    for m = 1:optimalK
        currGrids = nanmean(allMasks(clusters==m,:,:),3);
        errs = nanstd(currGrids,0,1);
        timecourse = nanmean(currGrids,1);
        
        if(~any(errs))
            errs = ones(size(errs));
        end
        
        %      errMetric = errMetric + log(sum(errs)) + log(sum(clusters==m)/(length(clusters)/(optimalK*maskSpacing)));
        
        subplot(subplotSize(1), subplotSize(2),m);
        hold on;
        shadedErrorBar(1:length(errs),timecourse,errs,'lineprops',{'Color', [colorP(m,:)]});
    end
    
    if(saveVar)
        if(c==length(stimNames)+1)
            saveas(gcf, [savePath, 'Rest_Curves.fig']);
        else
            title(stimNames{c});
            saveas(gcf, [savePath, stimNames{c},'_Curves.fig']);
        end
    end
    
    %     if(c==length(stimNames)+1)
    %         fprintf( fID, '%s : %f\n', 'Rest', errMetric);
    %     else
    %         fprintf( fID, '%s : %f\n', stimNames{c}, errMetric);
    %     end
    
    %     figure();
    %     plot(paramVals,errPlot);
    %     if(c==length(stimNames)+1)
    %         saveas(gcf, [savePath, 'Error_Rest_Curve.fig']);
    %     else
    %         saveas(gcf, [savePath, 'Error_',stimNames{c},'_Curve.fig']);
    %
    %     end
end
for c = 1:length(stimNames)
    vals = spIm{c};
    figure();
    percs = [];
    for m = 1:length(areaMasks)
        currMask = areaMasks{m};
        currMask(isnan(currMask)) = 0;
        
        armPerc = nanmean(vals(armArea&~handArea&currMask(:,:,1)));
        handPerc = nanmean(vals(~armArea&handArea&currMask(:,:,1)));
        facePerc = nanmean(vals(faceArea&currMask(:,:,1)));
        trunkPerc = nanmean(vals(trunkArea&currMask(:,:,1)));
        percs(m,:) = [armPerc, handPerc,facePerc,trunkPerc];
        %         percs(m,:) = [(armPerc+handPerc)/2,facePerc,trunkPerc];
        errors(m,:) = [armPerc/nanstd(vals(armArea&~handArea&currMask(:,:,1))),...
            handPerc/nanstd(vals(~armArea&handArea&currMask(:,:,1))),...
            facePerc/nanstd(vals(faceArea&currMask(:,:,1))),trunkPerc/nanstd(vals(trunkArea&currMask(:,:,1)))];
        %          errors(m,:) = [(armPerc/nanstd(vals(armArea&~handArea&currMask(:,:,1)))+...
        %              handPerc/nanstd(vals(~armArea&handArea&currMask(:,:,1))))/2,facePerc/nanstd(vals(faceArea&currMask(:,:,1))),...
        %              trunkPerc/nanstd(vals(trunkArea&currMask(:,:,1)))];
    end
    
    bar(percs);
    hbar = bar(percs);
    for p = 1:size(percs,2)
        ctr(p,:) = bsxfun(@plus, hbar(1).XData, hbar(p).XOffset);
    end
    hold on
    errorbar(ctr', percs, errors, '.k')
    title(stimNames{c});
    xticklabels({'M1', 'PMd'});
    legend({'Arm','Hand', 'Face', 'Trunk'});
end
% function clus = kmedoidsF(X,k)
%     clus = kmedoids(X,k,'Distance',@dtwf);
%     disp(k);
% end
function dist = dtwf(x,y)
% n = numel(x);
m2 = size(y,1);
dist = zeros(m2,1);
for i=1:m2
    %     coeffs = [];
    %     for s = 1:8
    %         str = (s-1)*70;
    %         if(any([x(:,str+1:str+70)==-1000, y(i,str+1:str+70)==-1000]))
    %             coeffs(s) = NaN;
    %         else
    %           coeffs(s) = max(xc orr(x(:,str+1:str+70),y(i,str+1:str+70), 10,'coeff'));
    %           coeffs(s) = 100*max(0,1-coeffs(s));
    %          %   coeffs(s) = dtw(x(:,str+1:str+70),y(i,str+1:str+70));
    %         end
    %     end
    %      dist(i) = nanmedian(coeffs);
    
    %dist(i) = dtw(x,y(i,:));
    dist(i) = max(0,max(xcorr(x,y(i,:),'biased')));
    %     h = corrcoef(x,y(i,:));
    %     dist(i) = max(0,h(1,2));
end
end 