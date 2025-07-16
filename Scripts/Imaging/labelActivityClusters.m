clear all;
close all;

monkey = 'Gilligan';
parentPath = ['S:\Lab\', monkey,'\Mapping\'];
MMFilename = 'MM_RGB_Simp-01.png';
activationFilenames = {'ESS_HSV.png', 'LS_HSV.png', 'PC_HSV.png'};
clusterFilenames = {'ExtraSmallSphere_Grid.fig', 'LargeSphere_Grid.fig', 'Photocell_Grid.fig'};
clusterPath = [parentPath,'\Clustering\15px_Grid\XCorr_Elbow\Best\'];
refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled.bmp'])>200;
refMask(:,1) = 0;
refMask(:,end) = 0;
refMask(1,:) = 0;
refMask(end-4:end,:) = 0;
refMaskD = refMask;
refMask = logical(refMask(:,:,1));

parentPath = ['S:\Lab\', monkey,'\Mapping\'];

MM = double(imread([parentPath, 'Motor Maps V3\', 'MM_Simp_RGB_Thresh-01.png']));

%border is yellow
M1border = MM(:,:,1)>200 & MM(:,:,2) > 200 & MM(:,:,3) < 100;
PMborder = MM(:,:,1)>100 & MM(:,:,2) < 100 & MM(:,:,3) > 150;

MM = double(imread([parentPath, 'Motor Maps V3\', MMFilename]));

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

areaMasks = {logical(rgb2gray(M1mask)), logical(rgb2gray(PMdmask))};

%arm is red
armArea = MM(:,:,1) > 100 & MM(:,:,2) < 100 & MM(:,:,3) < 100;
%hand is green
handArea = MM(:,:,2) > 100 & MM (:,:,1) < 100 & MM(:,:,3) < 100;
%face is blue
faceArea = MM(:,:,3) > 100 & MM (:,:,1) < 100 & MM(:,:,2) < 100;
%trunk is black
trunkArea = MM(:,:,1) < 100 & MM (:,:,2) < 100 & MM(:,:,3) < 100;

allMMs = {armArea, handArea, faceArea, trunkArea};
for a = 1:length(allMMs)
    otherInds = 1:length(allMMs);
    otherInds(a) = [];
    allMMsEx{a} = allMMs{a} > .5 & ...
        ~logical(sum(reshape(cell2mat(allMMs(otherInds)),[size(allMMs{a},1), size(allMMs{a},2),length(otherInds)]),3));
end
allMMsEx = cellfun(@(a) imfilter(a, ones(5)), allMMsEx, 'UniformOutput', false);


repNames = {'Arm', 'Hand', 'Face', 'Trunk'};
areaNames = {'M1', 'PMd'};
%%
for a = 1:length(activationFilenames)
    activation = imresize(imread([parentPath,'tTests\',activationFilenames{a}]),[size(refMask)]);
    %activation with 50% overlap is red
    activation = activation(:,:,1) > 100 & activation(:,:,2) < 100 & activation (:,:,3) < 100;
    activation = activation & refMask;
    clusters = openfig([clusterPath, clusterFilenames{a}],'invisible');
    clusters = gca(clusters);
    imIndex = find(cellfun(@length, arrayfun(@(a) size(a.CData), ...
        clusters.Children, 'UniformOutput', false))==2);
    clusters = clusters.Children(2).CData;
    numClusters = unique(max(clusters(:)));
    
    combined = zeros(size(activation));
    clusterSize = [];
    for c = 1:numClusters
        combined(activation==1 & clusters==c & refMask) = c;
        clusterSize(c) = sum(sum(clusters==c & refMask));
    end
    
    [~,minCluster] = min(clusterSize);
    goodClusters = 1:numClusters;
    goodClusters(minCluster) = [];
    totalClusters = ismember(clusters,goodClusters);
    contained = 100*sum(sum(activation & totalClusters))/sum(sum(activation));
    
    figure();
    %subplot(2,2,[1 3]);
    
    imagesc(combined);
    hold on;
    cmap = colormap(flipud(jet(numClusters)));
    colormap([1 1 1; cmap]);
    colorbar;
    caxis([0,numClusters]);
    h = imagesc(refMaskD);
    set(h, 'AlphaData',refMaskD(:,:,1)==0);
end
for a = 1:length(activationFilenames)
    subplot(2,2,2);
    hold on;
    percOfCluster = [];
    xlabs = {};
    for m1 = 1:length(allMMsEx)
        rep = allMMsEx{m1};
        for m2 = 1:length(areaMasks)
            area = areaMasks{m2};
            for c = 1:numClusters
                clus = (logical(clusters==c & refMask) & rep & area);
                percOfCluster((length(allMMsEx)*(m2-1)) + m1,c) = 100*sum(sum(clus))/clusterSize(c);
                xlabs(length(allMMsEx)*(m2-1) + m1) = {[areaNames{m2} ' ',repNames{m1}]};
            end
        end
    end
    
    title(sprintf('Cluster by representation'));
    b =bar(percOfCluster);
    for c = 1:numClusters
        b(c).FaceColor = cmap(c,:);
    end
    xticks([1:length(xlabs)]);
    xticklabels(xlabs);
    ylabel('%');
    
    
    subplot(2,2,4);
    hold on;
    percOfActivation = [];
    xlabs = {};
    for m1 = 1:length(allMMsEx)
        rep = allMMsEx{m1};
        for m2 = 1:length(areaMasks)
            area = areaMasks{m2};
            for c = 1:numClusters
                clus = (logical(clusters==c & refMask) & rep & area & activation);
                percOfActivation((length(allMMsEx)*(m2-1)) + m1,c) = 100*sum(sum(clus))/sum(sum(activation));
                xlabs(length(allMMsEx)*(m2-1) + m1) = {[areaNames{m2} ' ',repNames{m1}]};
            end
        end
    end
    
    title(sprintf('Activation by representation\n %.2f%% of activation in largest group(s)',contained));
    b =bar(percOfActivation);
    for c = 1:numClusters
        b(c).FaceColor = cmap(c,:);
    end
    xticks([1:length(xlabs)]);
    xticklabels(xlabs);
    ylabel('% ');
    saveas(gcf,['S:\Lab\Gilligan\Mapping\Clustering\15px_Grid\XCorr_Elbow\tTests\' ,activationFilenames{a}(1:end-8) ,'_Representation.fig']);
    
    
    %     for z = 1:numClusters
    %         figure();
    %         if(~exist([clusterPath,'tTests\'], 'dir'))
    %             mkdir([clusterPath,'tTests\']);
    %         end
    %         imshow(combined==z);
    %         export_fig(gcf,[clusterPath,'tTests\',clusterFilenames{a}(1:end-8),num2str(z),'.bmp']);
    %         pause(1);
    %     end
    
end