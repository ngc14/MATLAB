function [] = BlockView_ngc14_V1(varargin)
%% user input
animal = 'Skipper';
date = '11_30_2020';
run = 'run00';
threshold = 10;
stim = containers.Map([0 2 4 5],{'Rest','ExtraSmallSphere','LargeSphere', 'Photocell'});
fframe= 1:10;
stim1= 0;
LPk= 5; % Low-pass kernel size set as 0 if no LP filter
HPk= 250; % High-pass kernel size (usually >= 250) set as 0 if no HP filter
%% multirun
if ~isempty(varargin)
    % animal
    if ~isempty(varargin{1})
        animal = varargin{1};
    end
    % date
    if ~isempty(varargin{2})
        date = varargin{2};
    end
    % run
    if ~isempty(varargin{3})
        run = varargin{3};
    end
    % stim #
    if ~isempty(varargin{4})
        stim1 = varargin{4};
    end
    % stim name
    if ~isempty(varargin{5})
        stim = varargin{5};
    end
    % LPk
    if ~isempty(varargin{6})
        LPk = varargin{6};
    end
    % HPk
    if ~isempty(varargin{7})
        HPk = varargin{7};
    end
end
%% Begin Script
disp(['%%%%% BLOCKVIEW_NGC14: ',animal,' ',date(date~='_'),' ',run,...
    ' [',stim(stim1),'] LP',num2str(LPk),' HP',num2str(HPk),' %%%%%']);
dataPath = ['S:\Lab\',animal,'\All Data\',animal,'_',date, '\Imaging\',run];
savePath = ['C:\Users\ngc14',dataPath(3:end),'\Results\Blockview\','[',stim(stim1),']\'];
    if(~exist(savePath,'dir'))
    mkdir(savePath)
end
%% locate data
blocks = dir([dataPath,'\*.BLK']);
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);
if isempty(blocks)
    error('Cannot find data on SSD. Exiting.');
end
conds = NaN(size(blocks,1),1);
for n =1:size(blocks,1)
    anapar = OIHeadRead([blocks(n).folder,'\',blocks(n).name], 'v');
    conds(n)= anapar.Cond;
end
blockConds1 = blocks(conds==stim1);
% alignment
alignInd1 = find(cellfun(@(a) a==stim1, stim.keys));
badAlign = load([dataPath,'\offsets.mat']);
[row,col] =find(cellfun(@(a) nanmean(a)>threshold, badAlign.dist{alignInd1}));
row = row(~(ismember(col,11:20)));
badInds = unique(row);
badInds1 = ~ismember(1:size(blockConds1,1), badInds);
blockConds1 = blockConds1(badInds1);
% data properties
anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
W = anapar.FrameWidth;
H = anapar.FrameHeight;
NFrames = anapar.FramesPerStim;
NBlocks = size(blockConds1,1);
%% load data
align = load([dataPath,'\alignment_offset.mat']);
align.disp_x_med{alignInd1} = align.disp_x_med{alignInd1}(badInds1,:);
align.disp_y_med{alignInd1} = align.disp_y_med{alignInd1}(badInds1,:);
alignX = round(cell2mat(align.disp_x_med{alignInd1}));
alignY = round(cell2mat(align.disp_y_med{alignInd1}));
clear badAlign align
hFile = ['NaN_filtered_HP', num2str(HPk),'_LP',num2str(LPk),'.h5'];
if(exist([savePath,hFile],'file'))
    hinfo = h5info([savePath,hFile]);
    if(~contains([hinfo.Datasets.Name],'allTrials'))
        h5create([savePath,hFile], '/allTrials', [H W NFrames NBlocks],'Datatype','double','Chunksize',[H W 1 NBlocks]);
    end
else
    h5create([savePath,hFile], '/allTrials', [H W NFrames NBlocks],'Datatype','double','Chunksize',[H W 1 NBlocks]);
    fapl = H5P.create('H5P_FILE_ACCESS');
    H5P.set_fapl_sec2(fapl);
    H5P.set_fclose_degree(fapl,'H5F_CLOSE_STRONG');
    H5P.set_libver_bounds(fapl,'H5F_LIBVER_LATEST','H5F_LIBVER_LATEST')

    fileID = H5F.create([savePath,hFile],'H5F_ACC_TRUNC','H5P_DEFAULT',fapl);
    datatypeID = H5T.copy('H5T_NATIVE_DOUBLE');
    dataspaceID = H5S.create_simple(4,fliplr( [H W NFrames NBlocks]),[]);
    dcpl = H5P.create('H5P_DATASET_CREATE');
    H5P.set_chunk(dcpl,fliplr([H W NFrames 1]))
    H5P.set_alloc_time(dcpl,'H5D_ALLOC_TIME_INCR');
    datasetID = H5D.create(fileID,"allTrials",datatypeID,dataspaceID,dcpl);
    H5D.close(datasetID);
    H5P.close(dcpl);
    H5S.close(dataspaceID);
    H5T.close(datatypeID);
    H5F.close(fileID);
end
fileID = H5F.open([savePath,hFile],"H5F_ACC_RDWR","H5P_DEFAULT");
dsID = H5D.open(fileID,"/allTrials");
memspaceID = H5S.create_simple(4,fliplr([H W NFrames 1]),[]);

%%
for n = 1:NBlocks
    ni = datetime('now');
    frames_stim1 = OIReadStim([blockConds1(n).folder,'\',blockConds1(n).name],0,'v');
    % median filter to reduce CCD column defects
    for f = 1:size(frames_stim1,3)
        frames_stim1(:,:,f) = OIShift(medfilt2(frames_stim1(:,:,f),[5,5]),-1*alignX(n,f),-1*alignY(n,f));
    end
    %%
    frames_stim1 = imageFilter_LPHP_ngc14_NaN(100*(frames_stim1 -mean(frames_stim1(:,:,fframe),3))...
        ./mean(frames_stim1(:,:,fframe),3),LPk,HPk,[]);
    frames_stim1(isinf(frames_stim1) | isnan(frames_stim1)) = 0;
    dspaceID = H5D.get_space(dsID);
    H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr([0 0 0 n-1]),[],[],fliplr([H W NFrames 1]));
    H5D.write(dsID,'H5ML_DEFAULT',memspaceID,dspaceID,"H5P_DEFAULT",frames_stim1);
    H5F.flush(dsID,'H5F_SCOPE_GLOBAL');

    %h5write([savePath,hFile],'/allTrials',frames_stim1,[1 1 1 n],[H W NFrames 1]);
    disp(['Trial ', num2str(n), ' of ', num2str(NBlocks), ' (',num2str(time2num(datetime('now')-ni,'seconds')),')']);
end
H5S.close(dspaceID);
H5D.close(dsID);
H5F.close(fileID);
end