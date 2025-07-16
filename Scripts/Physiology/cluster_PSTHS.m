function clusterPSTHS(saveFigs, monkey, singleOrAll, sessionOrUnitPSTHS,...
    jointName, jointClass, PSTHS, activation)
vMap = load(['S:\Lab\', monkey, '\Mapping\Encoding Maps\Voronoi\Recording_Voronoi.mat']);
refMask = imread(['S:\Lab\', monkey, '\Mapping\clean_mask_filled'], 'bmp')>200;
refMask = double(refMask(:,:,1));
%MMHand = imdilate(bwperim(imdilate(bwareafilt(MM(:,:,1)<125 & MM(:,:,2)>125 & MM(:,:,3)<125, 1),strel('disk',24))),ones(24,24));
refMask = refMask;% & ~MMHand;
%%
taskLength = cellfun(@length,PSTHS);
PSTHSNaN = NaN(length(PSTHS), max(taskLength));
minSeg = min(cellfun(@(a) a(1), PSTHSegs))-1;
for i = 1:length(PSTHS)
    alignSegs = PSTHSegs{i}(1)-minSeg;
    PSTHSNaN(i,alignSegs:alignSegs+length(PSTHS{i})-1) = PSTHS{i};
    %PSTHS{i}(end+1:meanTaskLength) = eps;
    %PSTHS{i} =  interp1(1:length(PSTHS{i}), PSTHS{i}, 1:meanTaskLength);
end
allPSTHS = PSTHSNaN;
featureMatrix = allPSTHS;

% allPSTHS = reshape(cell2mat(PSTHSNaN'), length(PSTHSNaN), meanTaskLength);
% allPSTHSW = reshape(cell2mat(PSTHSW'), length(PSTHS), length(PSTHSW{1}));
% allPSTHS = [allPSTHS; allPSTHSW];
% condInds = [ones(length(PSTHS),1); 2*ones(length(PSTHS),1)];
% featureMatrixW = allPSTHSW;
% featureMatrix(isnan(featureMatrix))=eps;
% [score, explained, V] = svd(featureMatrix');
% [scoreW,explainedW,V] = svd(featureMatrixW');
% explained = diag(explained);
% [~,score,~,~,explained,~] = pca(featureMatrix, 'Economy', false);
% [~,scoreW,~,~,explainedW,~] = pca(featureMatrixW, 'Economy', false);
% dimInd = find(cumsum(explained)>95,1);
% dimIndW = find(cumsum(explainedW)>95,1);
% toCluster = score(:,1:dimInd);
% toClusterW = scoreW(:,1:dimIndW);
% toCluster = (score*allPSTHS)';
% toCluster = toCluster(:,1:dimInd);
% toClusterW = (scoreW'*allPSTHSW')';
% toClusterW = toClusterW(:,1:dimInd);
% toCluster = haart(uint8(featureMatrix(:,1:end-1)'), 5,'noninteger')';
%%
sumA = [];
clus = {};
toCluster = featureMatrix;
%toClusterW = featureMatrixW;
%sumW = [];
parfor cl=1:10
    [~, clu,sumM,~,~,~]  = kmedoids(toCluster,cl,'Distance',@dtwf);
    clus(cl) = {clu};
    sumA(cl) = sum(sumM);
    %[~, ~,sumN,~,~,~]  = kmedoids(toClusterW,cl,'Distance',@dtwf);
    %sumW(cl) = sum(sumN);
end
% kFunc = @(X,K)(kmedoids(X,K,'Distance',@dtwf));
% eva = evalclusters(toCluster,kFunc,'DaviesBouldin','KList',[1:5])
% optimalK = eva.OptimalK;
% optimalKW = triangle_threshold(sumW,'right',0);
% clustersW  = kmedoids(toClusterW,optimalKW,'Distance',@dtwf);
optimalK = triangle_threshold(sumA,'right', 1);
clusters  = kmedoids(toCluster,optimalK,'Distance',@dtwf);
%%
saveFigPath = ['S:\Lab\', monkey,'\Mapping\Encoding Maps\Clustering\'];
if(~exist(saveFigPath,'dir'))
    mkdir(saveFigPath);
end
saveFigPath = [saveFigPath, singleOrAll, '_', sessionOrUnitPSTHS, '_'];
cMap = [0 0 0;  hsv(optimalK); .9 .75 .2;];
% figure();
% subplot(1,2,1);
% scatter(score(:,1), score(:,2), 36, cMap(clusters+1,:), 'filled');
% xlabel('PCA 1');
% ylabel('PCA 2');
% title(['95% variance captured in first ',  num2str(dimInd), ' dimensions']);
% subplot(1,2,2);
% scatter3(score(:,1), score(:,2), score(:,3), 36, cMap(clusters+1,:), 'filled');
% xlabel('PCA 1');
% ylabel('PCA 2');
% zlabel('PCA 3');
% saveas(gcf,[saveFigPath, 'PCA_Graphs.fig']);
%for c = 1:optimalK
clusteredImage = NaN(size(vMap.siteMask,1),size(vMap.siteMask,2));
for i = 1:length(clusters)
    for x = 1:size(vMap.siteMask,1)
        for y = 1:size(vMap.siteMask,2)
            if vMap.siteMask(x,y,i) == 1
                [in, on] = inpolygon(x,y,vMap.v(vMap.c{i},2),vMap.v(vMap.c{i},1));
                if (in || on)
                    clusteredImage(x,y,1) = clusters(i)+1;
                end
            end
        end
    end
end

figure();
im = imshow(clusteredImage,cMap);
hold on;
set(im, 'AlphaData', ~isnan(clusteredImage));
h = imagesc(refMask);
set(h, 'AlphaData', refMask(:,:,1)==0);
h2 = imagesc(activation.*(optimalK+1));
set(h2, 'AlphaData', (activation==1).*.7);
% if(c==1)
%   title('ESS');
%   saveas(gcf,[saveFigPath, 'ESS_Spatial_Clustering.fig']);
% else
%   title('Rest');
%   saveas(gcf,[saveFigPath, 'Rest_Spatial_Clustering.fig']);
% end
% title(['Cluster ', num2str(c)]);
% saveas(gcf,[saveFigPath, 'Cluster_', num2str(c), '.fig']);
% hold off;
% end
% figure();
% im = imshow(clusteredImage,cMap);
% hold on;
% set(im, 'AlphaData', ~isnan(clusteredImage));
% h = imagesc(refMask);
% set(h, 'AlphaData', refMask(:,:,1)==0);
% title('ESS');
% hold off;
if(saveFigs)
    saveas(gcf,[saveFigPath, 'ESS_Spatial_Clustering.fig']);
end

figure();
hold on;
subplotSize = numSubplots(optimalK);
for c = 1:optimalK
    subplot(subplotSize(1), subplotSize(2),c);
    plot(allPSTHS(clusters==c,:)', 'Color', cMap(c+1,:));
    hold on;
    plot(nanmean(allPSTHS(clusters==c,:))', 'Color', 'k', 'LineWidth', 2);
    title(sprintf('%d sites of %d', sum(clusters==c), length(clusters)));
end
hold off;
if(saveFigs)
    saveas(gcf,[saveFigPath, 'Clusters_All.fig']);
end

figure();
hold on;
subplotSize = numSubplots(2*optimalK);
for a =0:1
    if (a==0)
        currActiveInds = ~activationInds;
        activeName = 'Inactive sites';
    elseif (a==1)
        currActiveInds = activationInds;
        activeName = 'Active sites';
    end
    for c = 1:optimalK
        subplot(subplotSize(1), subplotSize(2),(a*optimalK)+c);
        plot(PSTHS(currActiveInds & clusters==c,:)', 'Color', cMap(c+1,:));
        hold on;
        plot(nanmean(PSTHS(currActiveInds & clusters==c,:))', 'Color', 'k', 'LineWidth', 2);
        title(sprintf('%d sites of %d %s', sum(currActiveInds & clusters==c), ...
            sum(currActiveInds), activeName));
    end
end
if(saveFigs)
    saveas(gcf,[saveFigPath, 'Clusters_Activity.fig']);
end

figure();
hold on;
subplotSize = numSubplots(length(jointToName));
for j = 1:length(jointName)
    subplot(subplotSize(1), subplotSize(2),j);
    plot(PSTHS(cellfun(@(a) strcmp(a, jointName{j}), jointClass),:)', 'Color', cMap(j+1,:));
    hold on;
    plot(nanmean(PSTHS(cellfun(@(a) strcmp(a, jointName{j}), jointClass),:))', 'Color', 'k', 'LineWidth', 2);
    title(sprintf('%s: %d sites of %d %s', jointName{j}, sum(cellfun(@(a) strcmp(a, jointName{j}), jointClass)), ...
        length(jointClass)));
end
if(saveFigs)
    saveas(gcf,[saveFigPath, 'Clusters_Representation.fig']);
end

% for c = 1:2
%     figure();
%     hold on;
%     subplotSize = numSubplots(optimalK);
%     for k = 1:optimalK
%         subplot(subplotSize(1), subplotSize(2),k);
%         plot(allPSTHS(clusters==k & condInds==c,:)', 'Color', cMap(k+1,:));
%         hold on;
%         plot(nanmean(allPSTHS(clusters==k & condInds==c,:))', 'Color', 'k', 'LineWidth', 2);
%         title(sprintf('%d sites of %d', sum(clusters==k & condInds==c), numSites));
%     end
%     hold off;
%     if(c==1)
%         saveas(gcf,[saveFigPath, 'Clusters_ESS.fig']);
%     else
%         saveas(gcf,[saveFigPath, 'Clusters_Rest.fig']);
%     end
% end

% cMapW = [0 0 0; hsv(optimalKW)];
% figure();
% subplot(1,2,1);
% scatter(scoreW(:,1), scoreW(:,2), 36, cMapW(clustersW+1,:), 'filled');
% xlabel('PCA 1');
% ylabel('PCA 2');
% title(['Rest; 95% variance captured in first ',  num2str(dimIndW), ' dimensions']);
% subplot(1,2,2);
% scatter3(scoreW(:,1), scoreW(:,2), scoreW(:,3), 36, cMapW(clustersW+1,:), 'filled');
% xlabel('PCA 1');
% ylabel('PCA 2');
% zlabel('PCA 3');
% saveas(gcf,[saveFigPath, 'Rest_PCA_Graphs.fig']);

% clusteredImageW = NaN(size(vMap.siteMask,1),size(vMap.siteMask,2));
% for i = 1:numSites
%     for x = 1:size(vMap.siteMask,1)
%         for y = 1:size(vMap.siteMask,2)
%             if vMap.siteMask(x,y,i) == 1
%                 [in, on] = inpolygon(x,y,vMap.v(vMap.c{i},2),vMap.v(vMap.c{i},1));
%                 if (in || on)
%                     clusteredImageW(x,y,1) = clustersW(i)+1;
%                 end
%             end
%         end
%     end
%
% end
% figure();
% im2 = imshow(clusteredImageW,cMapW);
% hold on;
% set(im2, 'AlphaData', ~isnan(clusteredImageW));
% h2 = imagesc(refMask);
% set(h2, 'AlphaData', refMask(:,:,1)==0);
% title('Rest');
% hold off;
% saveas(gcf,[saveFigPath, 'Rest_Spatial_Clustering.fig']);

% figure();
% hold on;
% subplotSize = numSubplots(optimalKW);
% for c = 1:optimalKW
%     subplot(subplotSize(1), subplotSize(2),c);
%     plot(allPSTHSW(clustersW==c,:)', 'Color', cMapW(c+1,:));
%     hold on;
%     plot(nanmean(allPSTHSW(clustersW==c,:))', 'Color', 'k', 'LineWidth', 2);
%     title(sprintf('%d sites of %d', sum(clustersW==c), length(clustersW)));
% end
% saveas(gcf,[saveFigPath, 'Rest_Clusters.fig']);
end


function dist = dtwf(x,y)
m2 = size(y,1);
dist = zeros(m2,1);
for i=1:m2
    currY = y(i,:);
    dist(i) = dtw(x(~isnan(x)),currY(~isnan(currY)));
    %dist(i) = dtw((x-min(x))/(max(x)-min(x)),(currY-min(currY))/(max(currY)-min(currY)));
    %r = corrcoef(x,currY, 'Rows', 'pairwise');
    %dist(i) = 100-10*r(2);
    %lagRs = xcorr(x,currY);
    %dist(i) = max(0,max(lagRs(round(length(lagRs)/2)-2:round(length(lagRs)/2)+2)));
    %dist(i) = corr(x',currY','Type','Kendall');
    %[~,mx1] = max(x);
    %[~,mx2] = max(currY);
    %dist(i) = abs(mx1-mx2);
    %dist(i) = sum(abs(x-currY));
    %dist(i) = norm(x-currY);
end
end