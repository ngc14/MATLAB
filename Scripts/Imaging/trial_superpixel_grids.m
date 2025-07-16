close all;
clear all;
monkey = 'Gilligan';
stimNames = ["[ExtraSmallSphere]","[LargeSphere]", "[Photocell]","[Rest]"];
meanMed = 'Mean';
maskSpacing = 5;
medsOrMeans = 'medoids';
maxClusters = 7;
nCoeff = 1:7;
LPk = 5;
savePath = ["S:\Lab\ngc14\Working\",monkey,"\Clustering\",num2str(maskSpacing),"px\Trials\DFunc\STD\",string(meanMed),"\"];
refMask = uint8(imbinarize(imread(['S:\Lab\', monkey, '\Mapping\clean_mask_filled.bmp'])));
if(strcmp(monkey,"Gilligan"))
    refMask(end-85:end,:) = 0;
else
    refMask(1:75,:) = 0;
    refMask = rgb2gray(refMask);
end
[refMask(:,1),refMask(:,end),refMask(1,:),refMask(end,:)] = deal(0);
refMask = mat2gray(refMask);
H = size(refMask,1);
W = size(refMask,2);
if strcmp(monkey,'Gilligan')
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    HPk=550;
    tFrames = 51:55;
else
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    HPk = 250;
    tFrames = 49:53;
end

if(isempty(gcp('nocreate')))
    parpool();
end
optionSet = statset('UseParallel',true,'UseSubstreams',true,'Streams',RandStream('mlfg6331_64'),'Display','iter');
if(strcmp(medsOrMeans,'medoids'))
    clusterOpts= struct('Distance','spearman',...%@(xi,xj)dFunc(xi,xj,5),...
        'Algorithm','large','NumSamples',[],'PercentNeighbors',0.001,...
        'Replicates',3,'OnlinePhase','on','Options',optionSet);
elseif(strcmp(medsOrMeans,'means'))
    clusterOpts= struct('Distance','cosine',...
        'Display','iter','EmptyAction','Singleton','MaxIter',1000,...
        'Replicates',3,'OnlinePhase','on','Options',optionSet);
    savePath = [savePath,"kmeans\"];
end
savePath = strjoin(savePath,'');
if(~exist(savePath,'dir'))
    mkdir(savePath);
    mkdir(strjoin([savePath, 'Grids\'],''));
    mkdir(strjoin([savePath, 'Curves\'],''));
end

mm = MotorMapping(42);
count = 1;
linRef = find(refMask==0);
for m = 1:ceil(H/maskSpacing)
    rows =  ((m-1) * maskSpacing) + 1:min((m * maskSpacing) + 1,H);
    for n = 1:ceil(W/maskSpacing)
        cols = ((n-1) * maskSpacing) + 1:min((n * maskSpacing) + 1 ,W);
        [rV,cV] = ndgrid(rows,cols);
        lin = sub2ind(size(refMask),rV(:), cV(:));
        [rs,cs] = ind2sub(size(refMask),lin(~ismember(lin,linRef)));
        masks{count} = [rs,cs];
        masksLin{count} = lin(~ismember(lin,linRef));
        count = count + 1;
    end
end
badMasks = cellfun(@isempty, masks);
masks(badMasks) = [];
allMasks = masks;
%%
[cSpIm,spIm] = deal(cell(1,length(stimNames)));
for c = 1:length(stimNames)
    %%
    tTestIm = imread(['S:\Lab\ngc14\Working\',monkey,'\',char(stimNames(c)),'\tTests\NaN\tTest_Nomask',num2str(tFrames),'_p-0.0001.bmp']) & refMask;
    tTestMasks = cellfun(@(m) mode(tTestIm(m(:,1),m(:,2)),'all'),masks)';
%     masks = allMasks(tTestMasks);
    condFrames = loadCondFrames([meanMed,stimNames(c)],monkey);
    maskFrames = cellfun(@(ms) squeeze(mean(condFrames(ms(:,1),ms(:,2),:),...%[1:15, 26:70]),...
        [1 2],'omitnan')), masks,'UniformOutput',false);
    meanMasks = smoothdata(cell2mat(maskFrames),1,'movmean',3)';

    shiftR = randi(size(meanMasks,2)+1,1,size(meanMasks,1))-1;
    for cm = 1:size(meanMasks,1)
        shiftedSigs(cm,:)= circshift(meanMasks(cm,:),shiftR(cm));%randperm(size(meanMasks,2)));
%       shiftedSigs(cm,:) = meanMasks(cm,randperm(size(meanMasks,2)));
    end
%   fa = fft(meanMasks,[],2); % Fourier decomposition
%   fa(:,(nCoeff(end)+1:end)) = 0;
%   clusterSig = real(ifft(fa(:,nCoeff),[],2));
    clusterSig = meanMasks;
    % Condition clustering map and time series
    % perturb time series signal to control generated map
    startSample = randsample(length(masks),maxClusters*clusterOpts.Replicates);
    startMasks = permute(reshape(clusterSig(startSample,:)',[],maxClusters,clusterOpts.Replicates),[2 1 3]);
    startShifts = permute(reshape(shiftedSigs(startSample,:)',[],maxClusters,clusterOpts.Replicates),[2 1 3]);
    % determine optimal # of k clusters for perturbed signal
    [sumD,idxC,clk,sumDS,idxS,cls] = deal(cell(1,maxClusters));
    for cl=1:maxClusters
        clusterOpts.Start = startMasks(1:cl,:,:);
        clk{cl} = namedargs2cell(clusterOpts);
        clusterOpts.Start = startShifts(1:cl,:,:);
        cls{cl} = namedargs2cell(clusterOpts);
    end
    %%
    parfor cl=1:maxClusters
        if(strcmp(medsOrMeans,'medoids'))
            [idxC{cl},clu,sumD{cl},interM,midx,info]= kmedoids(clusterSig,cl,clk{cl}{:});
            [idxS{cl},~,sumDS{cl},~,~,~]= kmedoids(shiftedSigs,cl,cls{cl}{:});
        elseif(strcmp(medsOrMeans,'means'))
            [idxC{cl},clu,sumD{cl},interM]= kmeans(clusterSig,cl,clk{cl}{:});
            [idxS{cl},~,sumDS{cl},~]= kmeans(shiftedSigs,cl,cls{cl}{:});
        end
        disp(cl)
    end
    %%
    spIm{c} = plotClusters(maxClusters,idxC,idxS,clusterSig,shiftedSigs,masks,refMask,savePath,stimNames(c));
    %%
    figure();
%   optimalK = triangle_threshold(cellfun(@(a) sum(a),sumD),'right', 1);
    clusterP = evalclusters(clusterSig,@clusterFunc,'gap','KList',[1:10]);%cell2mat(idxC),'gap');
%     clusterP = [];
%     groupP = {};
%     for cl = 1:maxClusters
%         clusterCoeffs = {};
%         allGComb = {};
%         for cl2=1:cl
%             allGComb{cl2} = nchoosek(1:cl,cl2);
%             allGComb{cl2} = [1:cl]';
%             for g = 1:size(allGComb{cl2},1)
%                 ag = allGComb{cl2}(g,:);
%                 ig = setdiff(1:cl,ag);
%                 gMask = NaN(size(tTestMasks));
%                 gMask(ismember(idxC{cl},ag)) = 1;% & tTestMasks==1) = 1;
%                 gMask(ismember(idxC{cl},ig)) = 0;% & tTestMasks==0) = 0;
%                 if(cl==1 && g==1)
%                     gMask(1) = 0;
%                 end
%                 outHold = corrcoef(double(gMask),double(tTestMasks));
%                 clusterCoeffs{cl2}(g) = outHold(2);
%             end
%         end
%         [intermM,mIndInterim] = cellfun(@max,clusterCoeffs);
%         [clusterP(cl),mInd] = max(intermM);
%         groupP{cl} = allGComb{mInd}(mIndInterim(mInd),:);
% %       [percs,groupN] = sort(histcounts(idxC{cl}(tTestMasks==1),'Normalization','probability'),'descend');
% %       activeGroups = percs>.20;
% %       ag = groupN(activeGroups);
% %       ig = groupN(~activeGroups);
% %       (sum(ismember(idxC{cl}(tTestMasks==1),ag)) - sum(ismember(idxC{cl}(tTestMasks==1),ig)))/sum(tTestMasks==1) + ...
% %       (sum(ismember(idxC{cl}(tTestMasks==0),ig)) - sum(ismember(idxC{cl}(tTestMasks==1),ag)))/sum(tTestMasks==0);
%     end
     plot(clusterP);
     ylim([0 1]);
%     xNames = cellfun(@(xs,a) [num2str(a),': ', num2str(xs,'%d,')], groupP, num2cell(1:maxClusters), 'UniformOutput',false);
%     xticklabels(cellfun(@(s)s(1:end-1), xNames,'UniformOutput',false))
    saveas(gcf,[savePath+"Match"+stimNames(c)],'png')
end

function spIm = plotClusters(maxClusters,clu,cluS,clusterSigs,shiftedSigs,masks,refMask,savePath,condName)
saveVar = true;
for cl =1:maxClusters
    [currK, currKC] = deal(cl);
    clusters = clu{currK};
    controlClusters = cluS{currKC};
    % spatial map fig
    spIm= zeros(size(refMask,1),size(refMask,2));
    cSpIm= zeros(size(refMask,1),size(refMask,2));
    for m = 1:length(masks)
        spIm(masks{m}(:,1),masks{m}(:,2)) = clusters(m);
        cSpIm(masks{m}(:,1),masks{m}(:,2)) = controlClusters(m);
    end

    g = figure();
    h = imagesc(spIm);
    hold on;
    set(h,'alphadata',spIm.*spIm>0)
    cmap = colormap(flipud(jet(currK)));
    colormap([.6 .6 .6 ;cmap]);
    clim([0,currK]);
    h = imagesc(refMask);
    set(h, 'alphadata', refMask==0);
    title(condName+" "+num2str(cl));

    % cluster time course avgs fig
    ts = figure();
    subplotSize = numSubplots(currK);
    colorP = colormap(flipud(jet(currK)));

    for m = 1:currK
        currGrids = mean(clusterSigs(clusters==m,:,:),3,'omitnan');
        errs = std(currGrids,0,1,'omitnan');
        timecourse = mean(currGrids,1,'omitnan');

        subplot(subplotSize(2), subplotSize(1),m);
        hold on;
        title(strjoin([condName, num2str(m)]));
        shadedErrorBar(1:length(errs),timecourse,errs,'lineprops',{'Color', [colorP(m,:)]});
        ylim([-.25 .25])
    end
    if(saveVar)
        fd = getframe(gca(g));
        fileTypeName = '.png';
        imwrite(fd.cdata,strjoin([savePath,'Grids\',condName,num2str(cl),fileTypeName],''));

        tsc = get(ts,'Children');
        linkaxes(tsc(1:end))
        saveas(gcf, strjoin([savePath,'Curves\',condName,num2str(cl),fileTypeName],''));
    end
    %% Control mask
    cg = figure();
    h = imagesc(cSpIm);
    hold on;
    set(h,'alphadata',cSpIm.*cSpIm>0)
    cmap = colormap(flipud(jet(currKC)));
    colormap([.6 .6 .6 ;cmap]);
    clim([0,currKC]);
    h = imagesc(refMask);
    set(h, 'alphadata', refMask==0);
    title(strjoin([condName,"Control",num2str(cl)]));

    % time series of control clusters
    ct = figure();
    subplotSize = numSubplots(2*currKC);
    colorP = colormap(flipud(jet(currKC)));
    for m = 1:currKC
        currGrids = mean(clusterSigs(controlClusters==m,:,:),3,'omitnan');
        errs = std(currGrids,0,1,'omitnan');
        timecourse = mean(currGrids,1,'omitnan');
        subplot(subplotSize(2), subplotSize(1),(2*m)-1);
        hold on;
        shadedErrorBar(1:length(errs),timecourse,errs,'lineprops',{'Color', [colorP(m,:)]});
        title(strjoin([condName, "C"+num2str(m)]));

        subplot(subplotSize(2), subplotSize(1),2*m);
        hold on;
        shiftedT = mean(shiftedSigs(controlClusters==m,:,:),3,'omitnan');
        title(strjoin([condName,"shuffled", num2str(m)]));
        shadedErrorBar(1:length(errs), mean(shiftedT,1,'omitnan'),...
            std(shiftedT,0,1,'omitnan'),'lineprops',{'Color', [colorP(m,:)],...
            'Linestyle','--'},'patchSaturation',0.1);
        ylim([-1 1])
    end

    if(saveVar)
        fileTypeName = '.png';
        cf = getframe(gca(cg));
        imwrite(cf.cdata,strjoin([savePath,'Grids\',condName,num2str(cl),'Shuffled',fileTypeName],''));
        tcc= get(ct,'Children');
        linkaxes(tcc(1:2:end))
        linkaxes(tcc(2:2:end))
        saveas(ct, strjoin([savePath,'Curves\',condName,num2str(cl),'Shuffled',fileTypeName],''));
    end
end
end

function condFrames = loadCondFrames(condName,monkey)
LPk = 5;
condFrames = load(['S:\Lab\ngc14\Working\',monkey,'\',char(condName(end)),'\Trial',char(condName(1)),'.mat']);
condFrames = condFrames.condFrames;
stdFrames = load(['S:\Lab\ngc14\Working\',monkey,'\',char(condName(end)),'\STDFrames.mat']);
stdFrames = imgaussfilt3(cat(3,stdFrames.stdCFrames{:}),round(LPk/2),'FilterSize',...
    [LPk, LPk, 1], 'Padding', 0,'FilterDomain', 'spatial');
condFrames = condFrames ./ stdFrames;
 
% stdPix = arrayfun(@(r) repmat(r,size(stdFrames,[1 2])),...
%     squeeze(std(stdFrames,0,[1 2],'omitnan')), 'UniformOutput', false);
% meanPix = arrayfun(@(r) repmat(r,size(stdFrames,[1 2])),...
%     squeeze(mean(condFrames,[1 2],'omitnan')), 'UniformOutput', false);
% condFrames = condFrames-cat(3, meanPix{:})./ cat(3, stdPix{:});
end

