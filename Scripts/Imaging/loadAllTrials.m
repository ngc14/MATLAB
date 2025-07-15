function im = loadAllTrials(dirPath,monkey,dates,runs,LPk,HPk,cond,frames,FRAMEIND)
dPaths = arrayfun(@(d,r) strjoin([dirPath,monkey,"All Data",monkey+"_"+d,...
    "Imaging\run0"+num2str(r)],'\\'), dates, runs, 'UniformOutput',false);
frameFiles = cellfun(@(m) strjoin([m,"Results\Blockview",cond,"FF_1",...
    strjoin(["LP",num2str(LPk),"_HP",num2str(HPk)],''),"Align_Offset\noClip",...
    "NaN_filtered_HP"+num2str(HPk)+"_LP"+num2str(LPk)+".h5"],'\\'),dPaths,'UniformOutput', false);
validDates = cellfun(@(e) isfile(e),frameFiles);
frameFiles = frameFiles(validDates);
dPaths = dPaths(validDates);

if(~exist('FRAMEIND','var'))
    FRAMEIND = 3;
end
nFiles = sum(validDates);
im = cell(1,nFiles);
for d = 1:nFiles
    fileInfo = h5info(frameFiles{d},'/allTrials');
    offset = ones(1,numel(fileInfo.ChunkSize));
    offset(FRAMEIND) = min(frames);
    blockSz = Inf(1,numel(fileInfo.ChunkSize));
    blockSz(FRAMEIND) = length(frames);
    avgFrameTrials = squeeze(mean(h5read(frameFiles{d},'/allTrials',offset,...
        blockSz),FRAMEIND,'omitnan'));

    NDIM = ndims(avgFrameTrials);
    ts = load(strjoin([dPaths(d),"tform_nonrigid.mat"],'\\'));
    spacing = ts.Spacing;
    ts = ts.O_trans;
    parfor t = 1:size(avgFrameTrials,3)
        avgFrameTrials(:,:,t) = bspline_transform(ts,avgFrameTrials(:,:,t),spacing,NDIM);
    end
    im{d} = avgFrameTrials;
end
end
% rpt = cellfun(@(b) cat(NDIM,ts.O_trans,b),arrayfun(@(a) a.*ones(size(...
% ts.O_trans,1:NDIM-1)),1:size(avgFrameTrials,NDIM),'UniformOutput',false),'UniformOutput', false);
% nTransform = permute(cat(NDIM+1,rpt{:}),[1:NDIM-1,NDIM+1,NDIM]);