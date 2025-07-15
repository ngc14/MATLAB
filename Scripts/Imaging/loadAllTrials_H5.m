function im = loadAllTrials_H5(dirPath,monkey,dates,runs,LPk,HPk,cond,frames)
mDir = strjoin([dirPath,monkey,"All Data"],'\\');
dPaths = arrayfun(@(d,r) strjoin([mDir,monkey+"_"+d,"Imaging\run0"+...
    num2str(r)],'\\'), dates, runs, 'UniformOutput',false);

frameFiles = cellfun(@(m) strjoin([m,"Results\Blockview",cond,"FF_1",...
    strjoin(["LP",num2str(LPk),"_HP",num2str(HPk)],''),"Align_Offset\noClip",...
    "NaN_filtered_HP"+num2str(HPk)+"_LP"+num2str(LPk)+".h5"],'\\'),dPaths,'UniformOutput', false);
validDates = cellfun(@(e) isfile(e),frameFiles);
frameFiles = frameFiles(validDates);
dPaths = dPaths(validDates);
%%
BFS = matlab.io.datastore.BlockedFileSet(cellstr(frameFiles));
im = cell(1,BFS.NumBlocks);
for d = 1:BFS.NumBlocks
    bIM = blockedImage(BFS.BlockInfo(d).Filename,Adapter=customH5Adapter('allTrials'));
    FRAMEIND = find(bIM.BlockSize ==1);
    offsetL = ones(1,numel(bIM.BlockSize));
    offsetL(:,FRAMEIND) = length(frames);
    bIM.BlockSize = bIM.BlockSize .*offsetL;
    offsetL(:,FRAMEIND) = min(frames);
    bLs = blockLocationSet(1,offsetL,bIM.BlockSize);
    bIs = blockedImageDatastore(bIM,BlockLocationSet=bLs,PadPartialBlocks=false,PadMethod=NaN);
    sessionTrials = readall(bIs, 'UseParallel', false);
    avgFrameTrials = squeeze(mean(cat(3,sessionTrials{:}),3,'omitnan'));
    
    ts = load(strjoin([dPaths(d),"tform_nonrigid.mat"],'\\'));
    rpt = cellfun(@(b) cat(3,ts.O_trans,b),arrayfun(@(a) a.*ones(size(...
        ts.O_trans,[1 2])),1:size(avgFrameTrials,FRAMEIND),'UniformOutput',false),'UniformOutput', false);
    nTransform = permute(cat(4,rpt{:}),[1 2 4 3]);
    im{d} = bspline_transform(nTransform,avgFrameTrials,[ts.Spacing,mode(ts.Spacing)],3);
end
end