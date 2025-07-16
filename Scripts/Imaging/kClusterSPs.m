 close all;
clear all;
monkey = 'Gilligan';
stimNames = ["[ExtraSmallSphere]","[LargeSphere]", "[Photocell]","[Rest]"];
meanMed = 'Mean';
maskSpacing = 5;
medsOrMeans = 'medoids';
maxClusters = 10;
savePath = ["S:\Lab\ngc14\Working\",monkey,"\Clustering\",num2str(maskSpacing),"px\Trials\XCorr\STD\Gap\",string(meanMed),"\"];
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
    mm = MotorMapping(42);
else
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    HPk = 250;
    tFrames = 49:53;
    mm = MotorMapping(56);
end
if(strcmp(medsOrMeans,'means'))
    savePath = [savePath,"kmeans\"];
end
savePath = strjoin(savePath,'');
if(~exist(savePath,'dir'))
    mkdir(savePath);
    mkdir(strjoin([savePath, 'Grids\'],''));
    mkdir(strjoin([savePath, 'Curves\'],''));
end

if(isempty(gcp('nocreate')))
    parpool('threads');
end

[datesTable, masksR, ~] = getMonkeyInfo('S:\Lab\',monkey,["M1", "PMd"],false);
repThreshInd = cellfun(@(r)  ~(strcmp(r,"Face") | strcmp(r,"Trunk")), datesTable.SiteRep,'UniformOutput',false);
datesTable.SiteRep = cellfun(@(s,r) s(r), datesTable.SiteRep, repThreshInd, 'UniformOutput',false);
datesTable.Thresh = cellfun(@(t,r) t(r), datesTable.Thresh, repThreshInd, 'UniformOutput',false);
datesTable = datesTable(~cellfun(@isempty,datesTable.SiteRep),:);
siteLocation = [datesTable.x,datesTable.y];
[verticies, vCells] = voronoin(fliplr([siteLocation; [0 size(refMask,2); ...
    size(refMask,1) 0; 0 0; size(refMask,1) size(refMask,2)]]));
simpRep = string(cellfun(@(s,t) s(find(t==min(t),1)),  datesTable.SiteRep,...
    datesTable.Thresh,'UniformOutput',false));
siteMasks = [];
for i = 1:length(siteLocation)
    currSite = siteLocation(i,:);
    tempCircle = zeros(size(refMask) + 2*mm.tileBuffer);
    tempCircle((currSite(2)-mm.siteRadius+mm.tileBuffer):...
        (currSite(2)+mm.siteRadius+mm.tileBuffer),...
        (currSite(1)-mm.siteRadius+mm.tileBuffer):...
        (currSite(1)+mm.siteRadius+mm.tileBuffer)) = mm.poolCircle;
    tempCircle = tempCircle(mm.tileBuffer:end-(mm.tileBuffer+1),...
        mm.tileBuffer:end-(mm.tileBuffer+1));
    siteMasks(:,:,i) = tempCircle & poly2mask(verticies(vCells{i},2),...
        verticies(vCells{i},1),size(tempCircle,1),size(tempCircle,2));
end
forelimbMask = any(siteMasks,3);

count = 1;
linRef = find(refMask==0)% | forelimbMask==0);
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
    tTestMasks = cellfun(@(m) mode(tTestIm(m(:,1),m(:,2)),'all'),allMasks)';
%     masks = allMasks(tTestMasks);
    condFrames = loadCondFrames([meanMed,stimNames(c)],monkey);
    maskFrames = cellfun(@(ms) squeeze(mean(condFrames(ms(:,1),ms(:,2),:),...%[1:15, 26:70]),...
        [1 2],'omitnan')), allMasks,'UniformOutput',false);
    meanMasks = smoothdata(cell2mat(maskFrames),1,'movmean',3)';

    shiftR = randi(size(meanMasks,2)+1,1,size(meanMasks,1))-1;
    shiftedSig = NaN(size(meanMasks));
    for cm = 1:size(meanMasks,1)
        shiftedSig(cm,:)= circshift(meanMasks(cm,:),shiftR(cm));%randperm(size(meanMasks,2)));
    end
    clusterSig = meanMasks;
    %
    clearvars -except clusterSig maxClusters shiftedSig tTestMasks stimNames savePath c monkey tFrames refMask masks meanMed
    clusterP = evalclusters(clusterSig,@clusterFunc,'gap','KList',1:maxClusters,'B',3,...
        'ReferenceDistribution','PCA','Distance','correlation','SearchMethod','firstMaxSE');
    for s = 1:length(clusterP.InspectedK)
        sCluster{s} = clusterFunc(shiftedSig,s);
    end
    %
    clusterG = [];
    groupP = {};
    for cl = 1:maxClusters
        clusterCoeffs = {};
        allGComb = {};
        for cl2=1:cl
            allGComb{cl2} = nchoosek(1:cl,cl2);
            for g = 1:size(allGComb{cl2},1)
                ag = allGComb{cl2}(g,:);
                ig = setdiff(1:cl,ag);
                gMask = NaN(size(tTestMasks));
                gMask(ismember(clusterP.Y{cl},ag)) = 1;
                gMask(ismember(clusterP.Y{cl},ig)) = 0;
                if(cl==1 && g==1)
                    gMask(1) = 0;
                end
                outHold = corrcoef(double(gMask),double(tTestMasks));
                clusterCoeffs{cl2}(g) = outHold(2);
            end
        end
        groupP{cl} = allGComb;
        clusterG{cl}= clusterCoeffs;    
    end
    [maxPerCluster,maxGroup] =  cellfun(@(g) cellfun(@max,g(1),'UniformOutput',true),clusterG(2:end),'UniformOutput',false);
    maxPerCluster(2:end+1) = maxPerCluster(1:end);
    maxPerCluster(1) = clusterG{1};
    [intermM,mIndInterim] = cellfun(@max,maxPerCluster);
    [maxCorr,mInd] = max(intermM);
    gro = groupP{mInd}{1}(mIndInterim(mInd),:);
    figure(); plot(clusterP); hold on;
    plot(intermM); 

    corrChold = corrcoef(double(sCluster{2}==1),double(tTestMasks));
    corrC = corrChold(2);
    corrChold = corrcoef(double(sCluster{2}~=1),double(tTestMasks));
    corrC(end+1) = corrChold(2);
    title("Optimal: "+ num2str(clusterP.OptimalK)...
         +". Corr (" +num2str(mInd) + ": " +num2str(mIndInterim(mInd)) +" (" + num2str(maxGroup{mInd}) + ")):  "...
         +num2str(maxCorr, '%.3f') + " (Control: " + num2str(max(corrC),'%.3f')+")");

    saveas(gcf,savePath+"Match"+stimNames(c),'png');
    saveas(gcf,savePath+"Match"+stimNames(c),'eps');
    save(savePath+"Cluster"+stimNames(c),'clusterP','-mat');
%
plotClusters(clusterP.Y,sCluster,clusterSig,shiftedSig,masks,refMask,savePath,stimNames(c));
end

function spIm = plotClusters(clu,cluS,clusterSigs,shiftedSigs,masks,refMask,savePath,condName)
saveVar = true;

for cl =1:length(clu)
    NClusters = cl;
    clusters = clu{NClusters};
    controlClusters = cluS{NClusters};
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
    cmap = colormap(flipud(jet(NClusters)));
    colormap([.6 .6 .6 ;cmap]);
    clim([0,NClusters]);
    h = imagesc(refMask);
    set(h, 'alphadata', refMask==0);
    title(condName+" "+num2str(NClusters));

    % cluster time course avgs fig
    ts = figure();
    subplotSize = numSubplots(NClusters);
    colorP = colormap(flipud(jet(NClusters)));

    for m = 1:NClusters
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
        imwrite(fd.cdata,strjoin([savePath,'Grids\',condName,num2str(NClusters),fileTypeName],''));

        tsc = get(ts,'Children');
        linkaxes(tsc(1:end))
        saveas(gcf, strjoin([savePath,'Curves\',condName,num2str(NClusters),fileTypeName],''));
    end
    %% Control mask
    cg = figure();
    h = imagesc(cSpIm);
    hold on;
    set(h,'alphadata',cSpIm.*cSpIm>0)
    cmap = colormap(flipud(jet(NClusters)));
    colormap([.6 .6 .6 ;cmap]);
    clim([0,NClusters]);
    h = imagesc(refMask);
    set(h, 'alphadata', refMask==0);
    title(strjoin([condName,"Control",num2str(NClusters)]));

    % time series of control clusters
    ct = figure();
    subplotSize = numSubplots(2*NClusters);
    colorP = colormap(flipud(jet(NClusters)));
    for m = 1:NClusters
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
        imwrite(cf.cdata,strjoin([savePath,'Grids\',condName,num2str(NClusters),'Shuffled',fileTypeName],''));
        tcc= get(ct,'Children');
        linkaxes(tcc(1:2:end))
        linkaxes(tcc(2:2:end))
        saveas(ct, strjoin([savePath,'Curves\',condName,num2str(NClusters),'Shuffled',fileTypeName],''));
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