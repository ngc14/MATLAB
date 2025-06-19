%% Sunincore	Optical Imaging Data Processing 	- HDL
%%
%% See 'Sunin_readme.txt' for detail useage. 

%%		(060502) Add trial by trial filtering (flagtrialfilter)
%%      Temporally delete 'mean.txt' and 'mapinfo.txt' output, will put back later
%%      Add 'datasource' into suninuser so can process previous saved ivf files besides block file. 
%%      Add 'operation' into suninuser
%% V3.3	Add high-pass, low-pass filtering feature for average maps	(051128-060123)		
%		Also add 'tDCFrames=tDCFrames+1;' around line 520 for solving the problem that some pixels are 0 in dye imaging
%% V3.2v add vector analysis, 051025: correct one error in sp value calculating: previous one is scaled by number of blocks (relative shape doesn't change)
%% V3.2 (050519-20) Changed defination of 'flagmap' to better control output maps. 1 -> 4, less maps -> more maps
%%      Add 'blockselect' for control which block to process, delete 'startblk' and 'startstim'
%%      Add 'blockfilenum' to separate from 'blocknum', the latter one is the number in 'blockselect'
%%		Add 'DataType' to 'anapar', and modified all blk reading function to recognize different data type (12, 13, 14)
%%      Add 'flagsaveblk' for save averaged block file, 'flagquantify' to control fun.txt output
%%      Corrected 'goodstim' usage, was shaped wrong before.
%%      Add function that will search under '\expfolder\masks' for *.bmp or *.txt for masking (if name is not provided)
%% V3.1 Add randomnize stim sequence option 'flagrandom', will read '_stimseq.txt' from data folder.
%%      Correction: Modified 'operation' part: sum --> mean
%%		Modify OIAlign2, to include method 4 (normxcorr2.m)
%%		Modified filenames for accummap and tempmap output 
%% V3.0 Add 'goodstim' matrix, to select only good conditions for averaging. 
%%      add 'stimsumnum' to keep track of stim number been summed (useful in 'flagrandom' and '..' cases
%%      unfinished: need change the way to calculate domainavg and domainstdev

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check and create directories
FrameRate=4;    %hz
resultdriver=datadriver;
resultfolder=[resultdriver, 'expresult\', expname, runname];
if ~isdir([resultdriver, 'expresult\'])
    mkdir(resultdriver, 'expresult\');    
end
if ~isdir([resultdriver, 'expresult\', expname])
    mkdir([resultdriver, 'expresult\'], expname);    
end
if ~isdir(resultfolder)
    mkdir([resultdriver, 'expresult\', expname], runname);    
end
if ~isdir([resultfolder, 'blkavg\'])
    mkdir(resultfolder, 'blkavg\');    
end
if ~isdir([resultfolder, 'blkavg\singlecond\'])
    mkdir([resultfolder, 'blkavg\'], 'singlecond\');    
end
if ~isdir([resultfolder, 'singleblk\'])
    mkdir(resultfolder, 'singleblk\');    
end
if ~isdir([resultfolder, 'singleblk\singlecond\'])
    mkdir([resultfolder, 'singleblk'], '\singlecond\');    
end
if tempsaveblocks~=0 & ~isdir([resultfolder, 'tempmap\'])
    mkdir(resultfolder, 'tempmap\');    
end
if saveaccum~=0 & ~isdir([resultfolder, 'accummap\'])
    mkdir(resultfolder, 'accummap\');    
end
fprintf('Results will be saved in   "%s"\n', resultfolder);
blockfolder = strcat(datadriver, datafolder, expname, runname); 
if flagalign~=0
    shiftlog = strcat(resultfolder, 'shiftlog.txt');     %for output log file during processing
    fidshiftlog = fopen(shiftlog, 'w');  %for shiftlog output
end


% default value for vector analysis, 
vectclipsd=0;	% clip image befor calculating, to remove some extreams, may not necessory, put 0 for no clipping
vectmask=ones(504); % mask used for Bosking normalization
%vectmask=imread('masks\v1mask.bmp');
polarclipsd=1;  % angle*map map is usually too dark, use lower clipsd to make it brighter, used 1 for 040113GarRun2

% If this program is called only for creating vector maps based on previous
% calculation then do vector map only and exit.
% check if vector map (single condition map) exist
if flagvector==2
    tempfilename=struct2cell(dir([resultfolder, 'vector\vect*', ext, '.ivf']));
    filename=sort(tempfilename(1,:)');
    Nvect=size(filename, 1);      % how many conditions
    if Nvect==0
        fprintf('Could not find any single condition map for vector analysis, you need run map analysis first (set flagvector=1 in suninuser file)\r');
        return;
    else    
        fprintf('Use following files for vector analysis:\r');
        for i=1:size(filename,1)
            fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
        end
    end
	t1=clock;
	for k = 1:Nvect
        if k<10
            filename=strcat('vect0', num2str(k), ext, '.ivf');  
        else
            filename=strcat('vect', num2str(k), ext, '.ivf');  
        end
       	filename=strcat(resultfolder, 'vector\', filename);
		b=OIReadIVF(filename);
        if (k==1)
		    [height, width]=size(b);
            singlemap = zeros(height, width, Nvect); 
        end
        singlemap(:, :, k) = b; %norm_to_01(b);
    end
    vectmask=ones(height, width);
	lut=textread('bwfix.lut');  % this color table should be in sunin folder
	Polarmap=OIPolar(singlemap, lut, vectmask, vectclipsd, lowpass, highpass, filtermethod);
%	figure; image(norm_to_uint8(Polarmap.ang));  axis equal; axis off; colormap(lut);  colorbar; %%load billcolorfix; 
%	figure; imagesc(norm_to_uint8(Polarmap.mag)); axis equal; axis off; colormap(gray(256)); colorbar; 
	filename1 = strcat('Angle', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename1 = strcat(resultfolder, 'vector\', filename1);
	imwrite(norm_to_uint8(Polarmap.ang), lut, filename1, 'tiff');
%	OIWriteIVF(Polarmap.ang, [filename1(1:end-4), ext, '.ivf']);
	filename2 = strcat('Mag', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename2 = strcat(resultfolder, 'vector\', filename2);
	imwrite(norm_to_uint8(OIClip(Polarmap.mag, clipmethod, clipvalue)), filename2, 'tiff');
%	OIWriteIVF(Polarmap.mag, [filename2(1:end-4), ext, '.ivf']);
	filename3 = strcat('Polar', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename3 = strcat(resultfolder, 'vector\', filename3);
    mag=norm_to_01(Polarmap.mag);
    mag=OIClip(mag, 1, polarclipsd);   % to adjust map darkness    
    ang=double(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
    polarmap=ang;
    polarmap(:,:,1)=ang(:,:,1).*mag;
    polarmap(:,:,2)=ang(:,:,2).*mag;
    polarmap(:,:,3)=ang(:,:,3).*mag;
    imwrite(norm_to_uint8(polarmap), lut, filename3, 'tiff');
%	OIWriteIVF(polarmap, [filename3(1:end-4), ext, '.ivf']);
	filename4 = strcat('Sum', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename4 = strcat(resultfolder, 'vector\', filename4);
	imwrite(norm_to_uint8(OIClip(Polarmap.sum, clipmethod, clipvalue)), filename4, 'tiff');
%	OIWriteIVF(Polarmap.sum, [filename4(1:end-4), ext, '.ivf']);
	fprintf('\rDone (time: %f minutes)\r', etime(clock, t1)/60);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determine data source & check availability
clipmask = strcat(resultdriver, 'expresult\', expname, 'masks\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2
if (clipmethod==2)      %use masking (to exclude blood vessels or only include certain area)
    immask = imread (clipmask, 'bmp');
    %        immask = immask (cropsize+1:end-cropsize, cropsize+1:end-cropsize);
else 
    immask = 0;
end

if datasource==1
    %check;
    if isempty(filename) % Search for block files if filenames are not provided
        if system=='v'
            tempfilename=struct2cell(dir([blockfolder, '*.blk']));
        elseif system=='r'
            fprintf('Note: you may need delete non-block "*.da" files in data folder\n');
            tempfilename=struct2cell(dir([blockfolder, '*.da']));
        end
        filename=sort(tempfilename(1,:)');
        for i=1:size(filename,1)
            fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
        end
        fprintf('\nFound %d blk files(sorted, check sequence).\n', size(filename,1));
    end
    blockfilenum=size(filename, 1);      % how many blocks
elseif datasource==2    % if source data is 'ivf' files.
    smapnum = size(Smap, 1);       % Subtraction map number
    % check 'ivf' files
    ivffolder=[resultfolder, 'blkavg\singlecond\'];
    tempfilename=struct2cell(dir([ivffolder, '*.ivf']));
    filename=sort(tempfilename(1,:)');
    fprintf('Found following ''ivf'' files:\r');
    for i=1:size(filename,1)
        fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
    end
    NStim=size(filename,1);
    fprintf('Subtraction maps will be based on above %d ivf files(sorted, check sequence).\r', NStim);
    fprintf('reading ivf: ');
	for i=1:NStim
		maps(:,:,i)=OIReadIVF([ivffolder, getfield(cell2struct(filename(i), 'junk'), 'junk')]);
        fprintf('%d ', i);
    end
    fprintf('\r');
    % following info has to be defined before process
    flagsavedata=0;
    FrameHeight=size(maps(:,:,1), 1);
    FrameWidth=size(maps(:,:,1),2);
    for i=1:NStim
        if ~isempty(stim)
            mapnames(i)=stim(i);
        else
            if NStim~=size(stim, 1)
                fprintf('Error, condition number not match!');
                exit;
            end
            mapnames(i)={num2str(i)};
        end
    end
    for j=1:size(Smap,1)
        mapnames(i+j)=Smap(j, 1);
    end
    smapnum = size(Smap, 1);       % Subtraction map number
    mapnum = NStim + smapnum;   % total map number

    % Process
    SNCCalculateSMap;
    if flaglpfilter | flaghpfilter>0
        SNCFilterSMap;
    end
    if flagmap>0
        SNCSaveSMap;
    end
    fprintf('\rProcess finished\r');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read  info
anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(1), 'junk'), 'junk')), system)
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;
DataType=anapar.DataType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check & define parameters
if isempty(fframe)
    fprintf('note: you specified "fframe=[]", which means no ffsub (equal to opeation=0)\n');
    operation=0;
end
if isempty(superpixel)
    flagsp=0;
end
for i=1:NStim
    if ~isempty(stim)
        mapnames(i)=stim(i);
    else
        mapnames(i)={num2str(i)};
    end
end
for j=1:size(Smap,1)
    mapnames(i+j)=Smap(j, 1);
end
smapnum = size(Smap, 1);       % Subtraction map number
mapnum = NStim + smapnum;   % total map number

if isempty(blockselect)
    blockselect=1:blockfilenum;
end
blocknum=size(blockselect,2);

subtraction = 1;    % 0: division, 1: subtraction   (division to be implemented)


% === shifting ===
flagdarksub = 0;    % whether or not subtract dark frame
%cropsize = 0;   % how many pixel crop off along 4 edges, should be used only when doing image alignment (which causes 'blank' margins), (currently not used), may couse problem in clipping

% === data/map saving ===
flagsavedata=0;    % 1: save data file inaddition to the maps, note: slow
SDtype=0;          % 0: SD is across blocks, 1: SD is across different dots (only for only for 'flagmasking=2', i.e. dots from coordinates)
flagafunout=0;     % controls fun1-3 output: 0: quantify only single-condition map, 1: quantify both single-condition and analysis maps
flagsavedetailfun = 0;  % save detailed averaged (across different spots) domain values for each stim and block.  0: no detail save
dotdetail=[             % first colum is mapnum (single cond), second column is masknum
      ];     
if flagquantify~=0
    if flagmasking==0
        fprintf('Error: Need masks for quantification (flagquantify~=0 and flagmasking==0)!\n');
    end
    % open fun file here
end


% === masking/Quantifying ===
savemaskmap=1;     % 1: save maskmap based on coordinates (only for 'flagmasking=2')
if flagmasking>0 
    maskfolder=strcat(resultdriver, 'expresult\', expname, 'masks\');
    switch flagmasking
    case 1      % for map-loading, Mask bmp files are 0-1 monochrome bmp files, desired areas are 1's
        if isempty(maskname)        
            tempfilename=struct2cell(dir([maskfolder, '*.bmp']));
            maskfilename=sort(tempfilename(1,:)');
            masknum=size(maskfilename, 1);
            for i=1:masknum
                masknametemp(i,:)=getfield(cell2struct(maskfilename(i), 'junk'), 'junk');
            end
            maskname=char(masknametemp(:,1:end-4));         % somehow this line can't be put together with the previous i loop
        else
            masknum=size(maskname,1);
        end
        for i=1:masknum
            mask(:,:,i) = imread (strcat(maskfolder, maskname(i,:)), 'bmp');
        end
        mask=uint8(mask);   % need check if mask is binary.
        fprintf('\nFound %d bmp masks.\n', masknum);
    case 2     % masks are from coordinates, create masks in memory here
        if isempty(domainradius)
            domainradius=zeros(masknum)+5;
        end
        if isempty(maskname)
            tempfilename=struct2cell(dir([maskfolder, '*x.txt']));  % only look for x.txt
            maskfilename=sort(tempfilename(1,:)')
            masknum=size(maskfilename, 1)
            dotsinmap=200;       % assume maximum 100 dots in each domain map, only for 'flagmask==2'
            xx=zeros(dotsinmap, masknum);       % store the dot center coordinates          
            yy=zeros(dotsinmap, masknum);
            for i=1:masknum
                masknametemp(i,:)=getfield(cell2struct(maskfilename(i), 'junk'), 'junk');
            end
            maskname=char(masknametemp(:,1:end-5));
        else
            masknum=size(maskname,1);
        end
        for i=1:masknum
            domainfidx=fopen (strcat(resultdriver, 'expresult\', expname, 'masks\', maskname(i,:), 'x.txt'), 'r');
            domainfidy=fopen (strcat(resultdriver, 'expresult\', expname, 'masks\', maskname(i,:), 'y.txt'), 'r');
            xtemp=fscanf(domainfidx, '%f');
            ytemp=fscanf(domainfidy, '%f');
            domaindotnum(i)=size(xtemp, 1);
            for j=1:domaindotnum(i)
                xx(j, i)=xtemp(j);      % j is dots, i is mask type
                yy(j, i)=ytemp(j);
            end
            fclose(domainfidx);
            fclose(domainfidy);        
            mask(:,:,i)=OIMaskMake(0, xx(:,i), yy(:,i), domainradius(i), domaindotnum(i), FrameWidth, FrameHeight); % maskmake(manual, xs, ys, pixradius, dotnum, xdim, ydim)
            if (savemaskmap==1)     % output memory masks for evaluating if selected 
                maskmap=mask(:,:,i)*255;
                imwrite(uint8(maskmap), strcat(resultdriver, 'expresult\', expname, 'masks', maskname(i,:), '_1.bmp'), 'bmp');
            end
        end
        clear maskmap;
    end
end
mapinfo=zeros(3, blocknum, mapnum);      %save (median, lowclip, highclip) for each maps, blknum, mapnum 

% === pixel distribution histograms ===
histout=0;         % gray value histogram, 0: no histogram output, 1: histogram from entire map, 2: histogram from masked map
histindex=[8, 8];  % histindex(1, 1) is the mask number used for output histograms across all maps, -1 is used for non-blood-vessel masks, 0 for whole map, output to 'hist1.txt'
                   % histindex(1, 2) is the map number used for output histograms across all masks, output to 'hist2.txt'
hist3index=[        % pair of map-mask numbers to export histgram, output to 'hist3.txt'
    1 2           % e.g. '5 10' means using map #5 and mask #10
];
binnum=50;        % bin number for histogram
binrange=[-0.0025, 0.0025];   % histogram bin range (typical single condition map: [-0.006, 0.002], typical difference map: [-0.0005, 0.0005])
                            % also can change program to set bin range as the maximum and min of all the maps (around line 500).
nobvmask = strcat(resultdriver, 'expresult\', expname, 'masks\nobvlarge.bmp');     %used when histindex(1,1)=0, which means calculate histograms on nobv rigion, if no such bmp exist, calculate the whole map and give warning by the end.
% need improve, delete nobvmask

if (flagrandom==1);
    stimseqfile=[blockfolder, '_stimseq.txt'];
    fidtemp1=fopen(stimseqfile, 'r');
    if fidtemp1==-1
        fprintf('Error: couldn''t open file ''_stimseq.txt''\n');
    end
    stimseq=fscanf(fidtemp1, '%d');
	fclose(fidtemp1);
	if size(stimseq,1)~=blockfilenum*NStim
        fprintf('Error: number of stim in ''_stimseq.txt'' is %d\n', size(stimseq));
    end
    stimseq=reshape(stimseq, [NStim, blockfilenum])';
    % check if it's a valid stim sequence file
    if max(max(stimseq))>NStim | min(min(stimseq))<1
        fprintf('Error: stimid in ''_stimseq.txt'' must >0 and <%d\n', NStim);
    end
    for i=1:blockfilenum
        checksum=zeros(1, NStim);
        checksum(stimseq(i,:))=1;
        if sum(checksum)~=NStim
            fprintf('Error: checksum of ''_stimseq.txt'' is not equal to stim number (%d) at line %d\n', NStim, i);
        end
    end
end
    
if (flaggoodstim==1) % load good stim condition file
	fidgoodstim=fopen(strcat(resultfolder, 'goodstim.txt'), 'r');
	goodstim=fscanf(fidgoodstim, '%d');
	fclose(fidgoodstim);
	if size(goodstim, 1)~=NStim*blockfilenum
		fprintf('ERROR: check goodstim.txt, make sure it has %d blocks and %d conditions\r', blockfilenum, NStim);
		return;
	end
	if sum(sum(goodstim))~=sum(sum(goodstim.*goodstim))	% only 0*0=0 and 1*1=1
		fprintf('Error: ''goodstim.txt'' should contain only "1" and "0"'\n);
	end
	goodstim=reshape(goodstim, NStim, blockfilenum);    %was wrong before 050520: goodstim=reshape(goodstim, blockfilenum, NStim);
	goodstim=goodstim';		%goodstim(x, y) is block x, stim y
    if prod(sum(goodstim, 2))==0
        fprintf('Worning: too few good stims, one or more condition has no data!\n');
    end
	fprintf('\r Good stim parameters loaded\r');
    if flagmap>=3   % user reqires single block maps?
        fprintf('Warning: ''goodstim=1'' and ''flagmap>=3'', you may not get all single block maps, anykey to continue...');
        pause;
    end
else
    goodstim=ones(blockfilenum, NStim);
end	


% check
if (flagmasking==1 & SDtype==1)
	fprintf('ERROR: check flagmasking and SDtype\r');
end
if flaggoodstim==1 & flagtrialfilter==1
	fprintf('Can not perform trial-by-trial filtering if using goodstim, flagtrialfilter disabled\r')
	flagtrialfilter=0;
end
finalnotice='Process Finished: ';    % Programmer: add notice into this string, it will be displayed at the end of the calculation, use <finalnotice=strcat(finalnotice, 'addhere')>
timer0=clock;       % used for calculate how much time this program runing

% Log the parameters for this process secession.
fidmasterlog = fopen(strcat(resultfolder, 'masterlog.txt'), 'a');  
%fprintf(fidmasterlog, '%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\tstart%d-%d-%d %d:%d.%d\r',ext, blockfolder, blocknum, operation, flagalign, shiftrange(1), shiftrange(2), shiftrange(3), shiftrange(4), shiftstep1, clipmethod, clipvalue, sframe, fframe, bframe, flagmasking, fix(clock));

% Open & Initiate log files
%fprintf(fidshiftlog, '__________________Start Date: %d %d %d    Time: %d %d %d\r', fix(clock));

%Read Header info from one block file

% Other variables used
%        DCframes[y,x,j]: Stores DC frames in one trial
%        fmShifts[j, i, k]: shifts for best alignment
%

% Print out important information from the block header: width, height in terms of pixels; frame number per stim, stim numbers
fprintf('\nFrameWidth=%d FrameHeight=%d NFramesPerStim=%d NStimuli=%d\n',FrameWidth, FrameHeight, FramesPerStim, NStim);

% if do shift alignment, readin the very first frame
if flagalign==1
    frameone=OIReadFrame(strcat(blockfolder, firstframe), 'v', 1, 1);
%    shiftmask = imread(strcat(resultdriver, 'expresult\', expname, 'masks\shiftmask1.bmp'), 'bmp');     % bv mask for shifting method3
end
if flagalign==11    %% Load shift parameters
    [sfname, sfblock, sfstim, sfframe, sfx, sfy, sfcoor]=textread(strcat(resultfolder, 'shiftinput.txt'), '%s %d %d %d %f %f %f');
    if size(sfname,1)~=blocknum*NStim*FramesPerStim
        fprintf('\r Wrong shift parameternum');
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blocknum
        for i=1:NStim
            for j=1:FramesPerStim
                if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim*FramesPerStim+(i-1)*FramePerStim+j).sfname
                    fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname((k-1)*NStim*FramesPerStim+(i-1)*FramePerStim+j).sfname);
                end
            end
        end
    end
end
if flagalign==41  % Load shift parameters generated by flagalign==4 process (include '.1' file for all first frame shift and '.3' file for all rest frames)
    [sfname, sfblock, sfstim, sfframe, sfx1, sfy1, sfcoor]=textread(strcat(resultfolder, 'shiftinput.2.txt'), '%s %d %d %d %f %f %f');
    if size(sfname,1)~=blockfilenum*NStim
        fprintf('Error: flagalign===41, .2 file doesnot contain blockfilenum*NStim number of entries!\n');
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blocknum
        for i=1:NStim
            if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim+i).sfname
                fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstruct(k*NStim*FramesPerStim).sfname);
            end
        end
    end
%     [sfname, sfblock, sfstim, sfframe, sfx3, sfy3, sfcoor]=textread(strcat(resultfolder, 'shiftinput.3.txt'), '%s %d %d %d %f %f %f');
%     if size(sfname,1)~=blockfilenum*NStim*(FramesPerStim-1)
%         fprintf('Error: flagalign===41, .3 file doesnot contain blockfilenum*NStim*(FramePerStim-1) number of entries!\n');
%     end
%     shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
%     for k=1:blocknum
%         for i=1:NStim
%             for j=1:FramesPerStim
%                 if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramePerStim-1)+j).sfname
%                     fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramePerStim-1)+j).sfname);
%                 end
%             end
%         end
%     end
end
    


%*********************************************************************************************************************
%********************************        Start Process     ***********************************************************
%*********************************************************************************************************************

fprintf('Wait...');
if system=='v'
    headerlength = 1716;    % byte
    switch DataType
    case 12
        stimlength=FrameWidth*FrameHeight*FramesPerStim*2;
    case 13
        stimlength=FrameWidth*FrameHeight*FramesPerStim*4;
    case 14
        stimlength=FrameWidth*FrameHeight*FramesPerStim*4;
    otherwise
        fprintf('Error: data type must be 12, 13 or 14!\n');
        return;
    end        
elseif system=='r'
    headerlength = 5120;
    stimlength=FrameWidth*FrameHeight*FramesPerStim*2+8*FramesPerStim*2+FrameWidth*FrameHeight*2+8*2;
elseif system=='d'	% for Dan Ts's old data 
    headerlength = 0;
    stimlength=FrameWidth*FrameHeight*FramesPerStim*2;
end
timertrigger=1;         %just for estimate time to run
if (flagmasking>0) 
    domainvalues=zeros(masknum, mapnum, blocknum);    %store average value for each mask, each stim, and each block, can expand to include analysis map values
end
maps=zeros(FrameHeight, FrameWidth, mapnum);     %for average maps
if tempsaveblocks~=0
    tempmaps=zeros(FrameHeight, FrameWidth, mapnum);     %for tempory maps
end
NewSum=zeros(FrameHeight,FrameWidth,NStim);		  % a new way to calcultate average trials (can be different number for different stim)
stimsumnum=zeros(NStim,1);	%keep how many stim was summed. (in randomized & goodstim cases, it may not same for each condition).
if flagsaveblk
    avgblock=zeros(FrameHeight, FrameWidth, FramesPerStim, NStim);   % note: memory waste
end
timertrigger=1;
ctime=clock;
ctime0=clock;
% Start read/process block by block
fprintf('total block number = %d\r\n', blocknum);
for k=blockselect  % start block by block processing
	if blocknum>1
		if k==blockselect(2)    % just for estimate process time
			blockproctime=etime(clock, ctime0);
		    fprintf('\rTime for one block: %7.4f secs', blockproctime);
    		fprintf('\rTime for all (%d) blocks: %5.2f hours\r', blocknum, blockproctime*blocknum/3600);
		end
	end
    SumFrame=zeros([FrameHeight, FrameWidth, NStim]);
    fprintf('\rblock#%d(%s): ',k, getfield(cell2struct(filename(k), 'junk'), 'junk'));
    fid=fopen(strcat(blockfolder, getfield(cell2struct(filename(k), 'junk'), 'junk')),'rb','ieee-le'); % open block file
    fseek(fid, headerlength,'bof');    % skip header
    % Here start read/process each stim in block 'k'
    for i = 1:NStim
      if goodstim(k, i)==0    % skip this stim
        fseek(fid, stimlength,'cof');    
      else
        if flagrandom   % if stim seq in data is random, use i for reading, and use 'currentstimid' for averageing
            currentstimid=stimseq(k, i);
        else
            currentstimid=i;
        end
        fprintf(' stim%d ', currentstimid);
        switch system
        case 'v'    % for VDAQ data 
		    switch DataType
    		case 12
	            tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>double');   % read a long vector of integers from block file, ('t' is for temp, may need shifting)
    		case 13
	            tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint32=>double');   % read a long vector of integers from block file, ('t' is for temp, may need shifting)
    		case 14
	            tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'float32=>double');   % read a long vector of integers from block file, ('t' is for temp, may need shifting)
    		otherwise
    			fprintf('Error: data type must be 12, 13 or 14!\n');
                return;
    		end
            tDCFrames=reshape(tDCFrames, FrameWidth, FrameHeight, FramesPerStim);     % reshape the long vector into 3-D matrix
            tDCFrames=permute(tDCFrames, [2,1,3]); % Since matlab read column first. Same purpose as the old method (flip&rotate). Note, the width and height changed
        case 'r'    % for RedShirt data
            tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>double');  % Redshirt use 'uint16' data type
            tBNCFrames=fread(fid,[8*FramesPerStim], 'uint16');   %NBNC=8 Read 8 BNCs
            tDarkFrames=fread(fid,[FrameWidth*FrameHeight],'uint16');    % read dark frames
            tBNCDark=fread(fid, [8], 'uint16');  % read BNC dark 
            tDCFrames=reshape(tDCFrames, FramesPerStim, FrameWidth, FrameHeight); % reshape into 3-D matrix
            tDCFrames=permute(tDCFrames, [3,2,1]);    % so the dimension order is (height, width, frame)
        case 'd'
            tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>double');   % read a long vector of integers from block file, ('t' is for temp, may need shifting)
            tDCFrames=reshape(tDCFrames, FrameWidth, FrameHeight, FramesPerStim);     % reshape the long vector into 3-D matrix
            tDCFrames=permute(tDCFrames, [2,1,3]); % Since matlab read column first. Same purpose as the old method (flip&rotate). Note, the width and height changed
        end
        % add a small value here to prevent 'Devide by Zero' error for some 'zero pixel' case in dye imaging (060123)
        tDCFrames=tDCFrames+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test dark frame sub
        if flagdarksub==1
            AvgFrame2=zeros(504, 504);
%% test: subtract dark frame
%             tempfid=fopen(strcat(resultfolder, 'darkframe.txt'), 'r');
%             AvgFrame2= fscanf(tempfid, '%f\t');
%             AvgFrame2= permute(reshape(AvgFrame2, 504, 504), [2,1]);
%% test: subtract dark frame from dark line
            tempfid=fopen(strcat(resultfolder, 'darkline.txt'), 'r');
            line= fscanf(tempfid, '%f\t');
            line=line';
            size(line);
            for ii=1:504
                AvgFrame2(ii, :)=line;
            end
            fclose(tempfid);

            for eachf=1:FramesPerStim
                tDCFrames(:,:,eachf)=tDCFrames(:,:,eachf)-AvgFrame2;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image alignment 
        switch flagalign
        case 0 %if (flagalign==0)   % no shifting
            DCFrames = tDCFrames;   % how much time will it cost for this operation?
        case 1 %elseif (flagalign==1)        % shifting 1: all aligned with 'frameone'
            pixoff=[0 0]; 
            maxcor=0;
            % note: fmshifted frame is the same size as before shifting, may contain zeros due to the shifting
            for j = 1:FramesPerStim            
                if timertrigger==1|timertrigger==2  % estimate time
                    [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*FramesPerStim);
                end
                [tDCFrames(:,:,j), pixoff, maxcor] = OIAlign2(frameone, tDCFrames(:,:,j), shiftrange, shiftstep1, shiftstep2, shiftmethod); % find the best alignment
                fmShifts(:, j, i, k) = pixoff;       % didn't use
                fprintf(fidshiftlog, '%s\t%d\t%d\t%d\t%3.1f\t%3.1f\t%6.6f\n', getfield(cell2struct(filename(k), 'junk'), 'junk'), k, i, j, pixoff(1), pixoff(2), maxcor);
         %       DCFrames(:,:,j) = OIShift(tDCFrames(:,:,j), pixoff(1), pixoff(2));
            end
            % Unsolved problem here: shifting cause 'blank' margins, crop will affect masking, no crop at this time, use median to fill the margin in 'OIShift.m'
            DCFrames = tDCFrames;
%            DCFrames=tDCFrames(cropsize+1:end-cropsize, cropsize+1:end-cropsize,:);  %crop off blank margins caused by shifting, img size changed from here on...
%            FrameWidth=FrameWidth-2*cropsize;
%            FrameHeight=FrameHeight-2*cropsize;
        case 2 %elseif (flagalign==2)        % shifting within stim frames (align all other frames with the frame#1)
            pixoff=[0 0]; 
            maxcor=0;
            for j = 1:FramesPerStim      
                if (j<sumrange(1))|(j>sumrange(2))   % do not shift those frames not summing to save time
                    pixoff=[0 0]; 
                    maxcor=1;
                else
                    if timertrigger==1|timertrigger==2  % estimate time
                        [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*FramesPerStim);
                    end
                    [tDCFrames(:,:,j), pixoff, maxcor] = OIAlign2(tDCFrames(:,:,1), tDCFrames(:,:,j), shiftrange, shiftstep1, shiftstep2, shiftmethod); % find the best alignment
                end
                fprintf(fidshiftlog, '%s\t%d\t%d\t%d\t%3.1f\t%3.1f\t%6.6f\n', getfield(cell2struct(filename(k), 'junk'), 'junk'), k, i, j, pixoff(1), pixoff(2), maxcor);
            end
            DCFrames = tDCFrames;
        case 11 %elseif (flagalign==11) % loading shift parameters and do shifting
            for j = 1:FramesPerStim 
                 DCFrames(:,:,j) = OIShift(tDCFrames(:,:,j), sfx((k-1)*NStim*FramesPerStim+(i-1)*FramesPerStim+j), sfy((k-1)*NStim*FramesPerStim+(i-1)*FramesPerStim+j));
                 fprintf(fidshiftlog, '%s\t%d\t%d\t%d\t%3.1f\t%3.1f\t%6.6f\n', getfield(cell2struct(filename(k), 'junk'), 'junk'), k, i, j, sfx((k-1)*NStim*FramesPerStim+(i-1)*FramesPerStim+j), sfy((k-1)*NStim*FramesPerStim+(i-1)*FramesPerStim+j), k);
            end
        case 41 %elseif (flagalign==41) % loading 2 shift parameters and do shifting
            DCFrames = tDCFrames;   % how much time will it cost for this operation?
            % wait until later to do shifting
        end
        clear tDCFrames;    
        tSumFrame=mean(DCFrames(:,:,sumrange), 3);   % tSumFrame stores frame-summed maps.
        tSumFrame=reshape(tSumFrame, FrameHeight, FrameWidth);     % don't know why this is necessory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Operations to compute the single condition map and quantification measurement
        switch operation    
        case 0      % 0: Rx: no operation (original average frame values)
            SumFrame(:,:,currentstimid)=tSumFrame;           
        case 1      % 1: dR: (Rx-R0), Rx is average frames, and R0 is the first frame of each stim condition
            thefframe=mean(DCFrames(:,:,fframe),3);
            SumFrame(:,:,currentstimid)=(tSumFrame-thefframe);
        case 2      % 2: dR/R0: percentage change, R0 is the first frame of each stim condition
            thefframe=mean(DCFrames(:,:,fframe),3);
            SumFrame(:,:,currentstimid)=(tSumFrame-thefframe)./thefframe;
        case 3      % 3: dR/Rb: (Rx-Rb)/Rb: Rb is the averaged blank frames     ---- currently not implemented, since blank is read last
           % SumFrame(:,:,i)=(SumFrame-SumFrame(:,:,NStim))./SumFrame(:,:,NStim);      % Note: here assume the last frame is blank
           fprintf ('\r\roperation=3 doesn''t work yet!');
           return;
        case 4      % 4: dR/R0end: Devide by the average of first and last frames, see Chen-Bee__Frostig 1996
            SumFrame(:,:,currentstimid)=(tSumFrame-mean(DCFrames(:,:,fframe),3))./((mean(DCFrames(:,:,fframe),3)+DCFrames(:,:,FramesPerStim))/2);    
        otherwise
            fprintf('\n  Error: need specify "operation" type!\n');
            return
        end
        if flagalign==41
            SumFrame(:,:,currentstimid)=OIShift(SumFrame(:,:,currentstimid), sfx1((k-1)*NStim+i), sfy1((k-1)*NStim+i));
        end
        if flagsaveblk==1
            avgblock(:,:,:,currentstimid)=avgblock(:,:,:,currentstimid)+DCFrames;       % Note, simple summaton, not dRR, need compare with sum of dRR
        end
        NewSum(:,:,currentstimid)=NewSum(:,:,currentstimid)+SumFrame(:,:,currentstimid);    % Add this sum because of use of 'goodstim'
        stimsumnum(currentstimid)=stimsumnum(currentstimid)+1;     % record how many repeats for each stim
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Pixel time course recording
        if (flagsp>0)
            for j=1:masknum
                for f=1:FramesPerStim
                    maskmap= double(DCFrames(:,:,f)).*double(mask(:, :, j));
                    spvalues(f,j,currentstimid,stimsumnum(currentstimid))=sum(sum(maskmap))/sum(sum(double(mask(:,:,j))));
                end
            end
        end
      end     %% if goodstim(k, i)==0
    end       %%i = 1:NStim    
    status = fclose(fid);
    if flagalign~=0
        fclose(fidshiftlog);
        fidshiftlog = fopen(shiftlog, 'a'); %close and reopen to save current results
    end
    
    singlemap(:,:,1:NStim)=SumFrame(:,:,1:NStim);   % use singlemap to save both single condition and subtraction map from this block
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Generate subtraction maps    
    for j=1:size(Smap, 1)						% note: not differentiate goodstim=1 or 0, get single block map anyway
        Sum1=zeros(FrameHeight, FrameWidth);
        Sum2=zeros(FrameHeight, FrameWidth);
        condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
        for n=condlist1
            Sum1=Sum1+SumFrame(:,:,n);
        end
        Sum1=Sum1/size(condlist1,2);
        condlist2=getfield(cell2struct(Smap(j, 3), 'junk'), 'junk');
        for n=condlist2
            Sum2=Sum2+SumFrame(:,:,n);
        end
        if size(condlist2,2)~=0
            Sum2=Sum2/size(condlist2,2);
        end
        if flagtrialfilter & (flaghpfilter|flaglpfilter)
	        singlemap(:,:,i+j)=OIEasyFilter(Sum1-Sum2, lpfmethod, lpkernelsize, hpfmethod, hpkernelsize);
	    else
		    singlemap(:,:,i+j)=Sum1-Sum2;    
		end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sum over blocks
%     if (flaggoodstim==1)    % do not sum subtraction map (will re-calculate later), only sum single condition map
% 	    for i = 1:NStim
%             maps(:,:,i) =  maps(:,:,i) + singlemap(:,:,i).*goodstim(k, i);  % only good stim are sumed (need work on avg map.
%         end
%     else
	    for i = 1:mapnum   
    	    maps(:,:,i) =  maps(:,:,i) + singlemap(:,:,i); 
            if tempsaveblocks~=0
                tempmaps(:,:,i) =  tempmaps(:,:,i) + singlemap(:,:,i); 
            end
        end
%	end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  save single block maps for each blocks
    if flagmap==3
        for i=NStim+1:mapnum
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            [singlemap(:,:,i), mapinfo(1, k, i), mapinfo(2, k, i), mapinfo(3, k, i)]=OIClip(singlemap(:,:,i), clipmethod, clipvalue, immask);       % mapinfo will be write into 'mapinfo.txt'
            map=norm_to_uint8(singlemap(:,:,i));             
            if (flaggoodstim==1) & (i<=NStim) & (goodstim(k, i)==0)
                map(5:20, 5:20)=1;  % to indicate this is a bad stim
            end
            if (i<=9)   % just for nice looking of the file names.
                if (k<=10)
                    savefilename=strcat('0', num2str(i), '_', mapname, '_k0', num2str(k-1), ext, '.bmp');
                else
                    savefilename=strcat('0', num2str(i), '_', mapname, '_k', num2str(k-1), ext, '.bmp');
                end
            else 
                if (k<=10)
                    savefilename=strcat(num2str(i), '_', mapname, '_k0', num2str(k-1), ext, '.bmp');
                else
                    savefilename=strcat(num2str(i), '_', mapname, '_k', num2str(k-1), ext, '.bmp');
                end                
            end
            savefilename=strcat(resultfolder, 'singleblk\', savefilename);   %for saveing single-block-analysis map
            imwrite(map, savefilename, 'bmp');
        end
    end
    if flagmap==4
        for i=1:mapnum
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            [singlemap(:,:,i), mapinfo(1, k, i), mapinfo(2, k, i), mapinfo(3, k, i)]=OIClip(singlemap(:,:,i), clipmethod, clipvalue, immask);       % mapinfo will be write into 'mapinfo.txt'
            map=norm_to_uint8(singlemap(:,:,i));             
            if (flaggoodstim==1) & (i<=NStim) & (goodstim(k, i)==0)
                map(5:20, 5:20)=1;  % to indicate this is a bad stim
            end
            if (i<=9)   % just for nice looking of the file names.
                if (k<=10)
                    savefilename=strcat('0', num2str(i), '_', mapname, '_k0', num2str(k-1), ext, '.bmp');
                else
                    savefilename=strcat('0', num2str(i), '_', mapname, '_k', num2str(k-1), ext, '.bmp');
                end
            else 
                if (k<=10)
                    savefilename=strcat(num2str(i), '_', mapname, '_k0', num2str(k-1), ext, '.bmp');
                else
                    savefilename=strcat(num2str(i), '_', mapname, '_k', num2str(k-1), ext, '.bmp');
                end                
            end
            if (i<=NStim)
                savefilename=strcat(resultfolder, 'singleblk\singlecond\', savefilename);    %for saving single-block-single-condition map
            else
                savefilename=strcat(resultfolder, 'singleblk\', savefilename);   %for saveing single-block-analysis map
            end
            imwrite(map, savefilename, 'bmp');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantify (get avg value from masked maps), need these block value in 'f1detail.txt' output
    if (flagquantify~=0)
        for i=1:mapnum %NStim       % currently only output single-condition values, change to 'mapnum' here to include analysis maps
            for j=1:masknum
                maskmap= singlemap(:,:,i).*double(mask(:, :, j));
                domainvalues(j,i,k)=sum(sum(maskmap))/sum(sum(double(mask(:,:,j)))); % avg value for block 'k', stim 'i', and mask 'j'.
            end
        end
        clear maskmap;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save temp sum maps & data (for every x blocks, then clear matrix), useful for animal status changing)
    if (mod(k, tempsaveblocks)==0)  % set tempsaveblocks=0 or 999 to avoid save tempmaps
        for i=1:mapnum    
           mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
	       if (flagsavedata) 
               tempname = strcat(resultfolder, 'tempdata\', num2str(i), '_k', num2str(k), '.txt');
               dlmwrite(tempname, tempmaps(:,:,i), '\t'); 
           end
           [tempmaps(:,:,i), framemedian, lowClip, highClip]=OIClip(tempmaps(:,:,i), clipmethod, clipvalue, immask);
           map=norm_to_uint8(tempmaps(:,:,i));  
           if (k<=10)
               savefilename=strcat(mapname, '_tk00', num2str(k-1), ext, '.bmp');
           elseif k>10&k<=100
               savefilename=strcat(mapname, '_tk0', num2str(k-1), ext, '.bmp');
           else	%101-999
               savefilename=strcat(mapname, '_tk', num2str(k-1), ext, '.bmp');
           end 
           savefilename=strcat(resultfolder, 'tempmap\', savefilename);
           imwrite(map, savefilename, 'bmp');
        end
        tempmaps=zeros(FrameHeight, FrameWidth, mapnum);
    end %  if (mod(k, tempsaveblocks)==0)

    % to save temp accumulate map
    if (saveaccum)  % save the temporary accumulate maps
        for i=accummapnum+NStim    % use 'mapnum' to save all maps (a lot)
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            [tempaccu, framemedian, lowClip, highClip]=OIClip(maps(:,:,i), clipmethod, clipvalue, immask);
            map=norm_to_uint8(tempaccu);
            if (k<=10)
                savefilename=strcat(mapname, '_ak00', num2str(k-1), ext, '.bmp');
            elseif k>10&k<=100
                savefilename=strcat(mapname, '_ak0', num2str(k-1), ext, '.bmp');
            else	%101-999
                savefilename=strcat(mapname, '_ak', num2str(k-1), ext, '.bmp');
            end 
            savefilename=strcat(resultfolder, 'accummap\', savefilename);
            imwrite(map, savefilename, 'bmp');
        end
    end  
end % for k1:blocknum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Done with file reading and initial processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tempmaps;
clear DCFrames;
if flagalign~=0
    fclose(fidshiftlog);
end

% average 
if (flagsaveblk)    % save averaged block file 
    for i=1:NStim
        avgblock(:,:,:,i)=avgblock(:,:,:,i)/stimsumnum(i);
    end
    fidlastblock=fopen(strcat(blockfolder,getfield(cell2struct(filename(end), 'junk'), 'junk')), 'r');   % simply copy an existing file header since they are the same, except data type
    header=fread(fidlastblock, 1716, 'uint8');
    fclose(fidlastblock);
    header(29)=uint8(14);   % force to float type
    fidnewblock=fopen(strcat(blockfolder, '_avgblk', ext, '.aBLK'), 'wb');
    fwrite(fidnewblock, header, 'uint8');
    avgblock=permute(avgblock, [2,1,3,4]);
    fwrite(fidnewblock, avgblock, 'float32');
    fclose(fidnewblock);    
end
    
if (flaggoodstim==1)        % if use goodstim, avg map is from averaged single condition map,  if goodstim=0, use averaged subtraction map, should be the same. 
	for i=1:NStim
		maps(:,:,i)=NewSum(:,:,i)/stimsumnum(i); %maps(:,:,i)/sum(goodstim(k,:));
	end
	% for 'flaggoodstim==1', need calculate subtraction map based on selected stims
    SNCCalculateSMap;
else 	%         
	for i=1:mapnum
	    maps(:,:,i)=maps(:,:,i)/blocknum;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output mapinfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output histograms (pixel gray value distributions) for each map
if (histout==2) % need implement histout==1 (or delete)
    % 1. output pixel gray value distributions for each map (using one mask)
% modify following 6 lines to change bin range style.
%    allmapmax=max(max(max(maps)));
%    allmapmin=min(min(min(maps)));
%    binmax=allmapmax;      
%    binmin=allmapmin;
    binmin=binrange(1,1);
    binmax=binrange(1,2);
    binsize=abs(binmax-binmin)/binnum;
    bins=(binmin:binsize:binmax);
    k=histindex(1,1);   % mask number
    if k==-1            % use non-blood-vessel region
        histomask=imread(nobvmask, 'bmp');
    elseif k==0         % use the whole map
        histomask=ones(FrameHeight,FrameWidth);
    else                % use specific domains
        histomask=mask(:,:,k);    
    end
    for i=1:mapnum     % get the values from the masked maps
        pixnum(i)=sum(sum(double(histomask)));              % following 6 lines put masked pixels into 'maskary', any better way?
        maskmap= maps(:,:,i) .* double(histomask);
        maskmap= maskmap + (double(histomask)-1)*(-999);    %% NOTE: -999 need change if the elements of the map can be that value or higher
        maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
        maskmap= sort(maskmap, 1);
        maskary= maskmap(1:pixnum(i), 1);
        maphist(:,i)=(hist(maskary, bins))';
    end
    % output
    histfid1= fopen(strcat(resultfolder, strcat('hist1', '.txt')), 'w');  % histo type 1: across all maps but using one mask
    fprintf(histfid1, '#\tbin\t');  
    for i=1:mapnum
       fprintf(histfid1, '%s\t', getfield(cell2struct(mapnames(i), 'junk'), 'junk'));
    end
%    fprintf(histfid1, '\r\t\t');
%    for i=1:mapnum
%        fprintf(histfid1, '%d\t', pixnum(i));
%    end
    fprintf(histfid1, '\r');
    for j=1:binnum
        fprintf(histfid1, '%d\t%f\t', j, bins(j));
        for i=1:mapnum
            fprintf(histfid1, '%f\t', maphist(j, i)/pixnum(i));
        end
        fprintf(histfid1, '\r');
    end
    clear maphist;
    clear pixnum;
    fclose(histfid1);
    
    % 2. hist type 2: output pixel value distrubution from one map but for each masks
    i=histindex(1,2);
    for j=1:masknum
        pixnum(j)=sum(sum(double(mask(:,:,j))));         % following 6 lines put masked pixels into 'maskary', better way?
        maskmap= maps(:,:,i) .* double(mask(:, :, j));
        maskmap= maskmap + (double(mask(:, :, j))-1)*(-999);    
        maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
        maskmap= sort(maskmap, 1);
        maskary= maskmap(1:pixnum(j), 1);
        maphist(:,j)=(hist(maskary, bins))';
    end
    % output
    histfid2= fopen(strcat(resultfolder, strcat('hist2', '.txt')), 'w');  
    fprintf(histfid2, '#\tbin\t');
    for i=1:masknum
        fprintf(histfid2, '%s\t',maskname(i,:));
    end
%    fprintf(histfid2, '\r\t\t');
%    for i=1:masknum
%        fprintf(histfid2, '%d\t', pixnum(i));
%    end
    fprintf(histfid2, '\r');
    for j=1:binnum 
        fprintf(histfid2, '%d\t%f\t', j, bins(j));
        for i=1:masknum
            fprintf(histfid2, '%f\t', maphist(j, i)/pixnum(i));
        end
        fprintf(histfid2, '\r');
    end
    clear maphist;
    clear pixnum;
    fclose(histfid2);
    
    % 3. hist type 3: output gray value histograms from provided map-mask pairs
    pairs=size(hist3index, 1);  % number of pairs
    for i=1:pairs
        imap =hist3index(i, 1); %index of map
        imask=hist3index(i, 2);
        pixnum(i)=sum(sum(double(mask(:,:,imask))));         % following 6 lines put masked pixels into 'maskary', better way?
        maskmap= maps(:,:,imap) .* double(mask(:, :, imask));
        maskmap= maskmap + (double(mask(:, :, imask))-1)*(-999);    
        maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
        maskmap= sort(maskmap, 1);
        maskary= maskmap(1:pixnum(i), 1);
        maphist(:,i)=(hist(maskary, bins))';
    end
    % output
    histfid3= fopen(strcat(resultfolder, strcat(getfield(cell2struct(filename(1), 'junk'), 'junk'), 'hist3', '.txt')), 'w');  %here also used file names
    fprintf(histfid3, '\t\t');
    for i=1:pairs
        mapname=getfield(cell2struct(mapnames(hist3index(i,1)), 'junk'), 'junk');
        fprintf(histfid3, '%s\t',mapname);
    end
    fprintf(histfid3, '\r');
    fprintf(histfid3, '#\tbin\t');
    for i=1:pairs
        fprintf(histfid3, '%s\t',maskname(hist3index(i,2),:));
    end
%    fprintf(histfid2, '\r\t\t');
%    for i=1:masknum
%        fprintf(histfid2, '%d\t', pixnum(i));
%    end
    fprintf(histfid3, '\r');
    for j=1:binnum 
        fprintf(histfid3, '%d\t%f\t', j, bins(j));
        for i=1:pairs
            fprintf(histfid3, '%f\t', maphist(j, i)/pixnum(i));
        end
        fprintf(histfid3, '\r');
    end
    clear maphist;
    clear pixnum;
    fclose(histfid3);
end

%save spvalues.mat spvalues;			
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Pixel time course output, only for 'flagsp>0' and 'flagmasking>0'
if (flagsp>0)
    xlabel=(1:FramesPerStim)/FrameRate;     % in unit of second
    spavg=zeros(FramesPerStim, masknum, NStim);
    for i=1:NStim  
        spavg(:,:,i)=sum(spvalues(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
    end
    spavg=reshape(spavg, [FramesPerStim, masknum, NStim]);
	spstd=zeros([FramesPerStim, masknum, NStim]);	
	for i=1:NStim		
		spstd(:,:,i)=std(spvalues(:,:,i,1:stimsumnum(i)),0,4);
	end
    spstd=reshape(spstd, [FramesPerStim, masknum, NStim]);
    superfid1= fopen(strcat(resultfolder, strcat('superpixel1', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    superfid2= fopen(strcat(resultfolder, strcat('superpixel2', '.txt')), 'w');
    fprintf(superfid1, 'original value\r');
    fprintf(superfid1, '%d\t', 1:FramesPerStim);
    fprintf(superfid1, '\r', f);
    for (j=1:masknum)
        fprintf(superfid1, '%s\t', maskname(j, :));
        if j==1
            fprintf(superfid1, '%f\t', xlabel);
        end
        fprintf(superfid1, '\r');
        for (i=1:NStim)
            fprintf(superfid1, 'stim%d\t', i);
            for (f=1:FramesPerStim)
                fprintf(superfid1, '%6.5e\t', spavg(f, j, i));			% original value
            end
            fprintf(superfid1, '\r');
        end
    end
    fprintf(superfid1, '\r\r_______________SD_________\r\t');
    for (f=1:FramesPerStim)
        fprintf(superfid1, 'f%d\t', f);
    end
    for (j=1:masknum)
        fprintf(superfid1, '%s\r', maskname(j, :));
        for (i=1:NStim)
            fprintf(superfid1, 'stim%d\t', i);
            for (f=1:FramesPerStim)
                fprintf(superfid1, '%6.5e\t', spstd(f, j, i));
            end
            fprintf(superfid1, '\r');
        end
    end

    % dRR 	
    superfid1a= fopen(strcat(resultfolder, strcat('superpixel1dRR', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    superfid2a= fopen(strcat(resultfolder, strcat('superpixel2dRR', '.txt')), 'w');
    superfid3a= fopen(strcat(resultfolder, strcat('superpixel3dRR', '.txt')), 'w');  
    spvaluestemp=spvalues;
   	for i=1:NStim
	    for k=1:stimsumnum(i)
    		for j=1:masknum
    			spvaluestemp(:, j, i, k)=(spvalues(:,j,i,k)-spvalues(1,j,i,k))/spvalues(1,j,i,k);
    		end
    	end
    end
    for i=1:NStim  
        spavg(:,:,i)=sum(spvaluestemp(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
    end
    spavg=reshape(spavg, [FramesPerStim, masknum, NStim]);
%	spstd=zeros([FramesPerStim, masknum, NStim]);
% 	for i=1:NStim		
% 		spstd(:,:,i)=std(spvalues(:,:,i,1:stimsumnum(i)),0,4);
% 	end
 	for i=1:NStim		
 		spstd(:,:,i)=std(spvaluestemp(:,:,i,1:stimsumnum(i)),0,4);
 	end
    spstd=reshape(spstd, [FramesPerStim, masknum, NStim]);

    % superpix1 (same mask are grouped together)
    fprintf(superfid1a, 'dR/R value\t');
    fprintf(superfid1a, 'f%d\t', 1:FramesPerStim);
    fprintf(superfid1a, '\r');
    for (j=1:masknum)
        fprintf(superfid1a, '%s\t', maskname(j, :));
        if j==1
            fprintf(superfid1a, '%f\t', xlabel);
        end
        fprintf(superfid1a, '\r');
        for (i=1:NStim)
            fprintf(superfid1a, 'stim%d\t', i);
            for (f=1:FramesPerStim)
                fprintf(superfid1a, '%4.3f%%\t', spavg(f, j, i)*100);			% original value
            end
            fprintf(superfid1a, '\r');
        end
    end
    fprintf(superfid1a, '\r\r_______________SD_________\r\t');
    fprintf(superfid1a, 'f%d\t', 1:FramesPerStim);
    fprintf(superfid1a, '\r');
    for (j=1:masknum)
        fprintf(superfid1a, '%s\r', maskname(j, :));
        for (i=1:NStim)
            fprintf(superfid1a, 'stim%d\t', i);
            for (f=1:FramesPerStim)
                if blocknum>1
	                fprintf(superfid1a, '%6.5e\t', spstd(f, j, i));
	            end
            end
            fprintf(superfid1a, '\r');
        end
    end
    
	% superpix2 is just another format for sp1 (same stim are grouped together)
    fprintf(superfid2a, 'dR/R value\t');
    fprintf(superfid2a, 'f%d\t', 1:FramesPerStim);
    fprintf(superfid2a, '\r');
    for (i=1:NStim)
        fprintf(superfid2a, 'stim%d\t', i);
        if i==1
            fprintf(superfid2a, '%f\t', xlabel);
        end
        fprintf(superfid2a, '\r');
	    for (j=1:masknum)
	        fprintf(superfid2a, '%s\t', maskname(j, :));
            for (f=1:FramesPerStim)
                fprintf(superfid2a, '%6.5e\t', spavg(f, j, i));			% original value
            end
            fprintf(superfid2a, '\r');
        end
    end
    fprintf(superfid2a, '\r\r_______________SD_________\r\t');
    for (f=1:FramesPerStim)
        fprintf(superfid2a, 'f%d\t', f);
    end
    fprintf(superfid2a, '\r');
    for (i=1:NStim)
        fprintf(superfid2a, 'stim%d\r', i);
	    for (j=1:masknum)
	        fprintf(superfid2a, '%s\t', maskname(j, :));
            for (f=1:FramesPerStim)
                if blocknum>1
	                fprintf(superfid2a, '%6.5e\t', spstd(f, j, i));
	            end
            end
            fprintf(superfid2a, '\r');
        end
    end

    % superfid3a is another version of superfid1a, it plot frames in vertical instead 050928 per Kristof
    fprintf(superfid3a, 'dR/R value\t');
    for (j=1:masknum)
        fprintf(superfid3a, 'mask: %s\t', maskname(j, :));
        for (i=1:NStim)
        	fprintf(superfid3a, '\t');
        end
	end
	fprintf(superfid3a, '\r');
    for (j=1:masknum)
        fprintf(superfid3a, '\t');
        for (i=1:NStim)
            fprintf(superfid3a, 'stim%d\t', i);
        end
	end
	fprintf(superfid3a, '\r');	
    for (f=1:FramesPerStim)
        fprintf(superfid3a, 'f%d\t', f);
	    for (j=1:masknum)
    	    for (i=1:NStim)
                fprintf(superfid3a, '%6.5e\t', spavg(f, j, i));			% original value
            end
            fprintf(superfid3a, '\t');
        end
        fprintf(superfid3a, '\r');
    end
    clear spavg;
    clear spstd;
    fclose(superfid1);
    fclose(superfid2);
    fclose(superfid1a);
    fclose(superfid2a);
    fclose(superfid3a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save average values across domains, only for 'flagmasking>0' 
if (flagquantify>0)
	% output detailed (including every block) averaged masked values for each domain type (masks) and stimulus. 
	if (flagquantify==2)
        fundfid1=fopen(strcat(resultfolder, 'fun1detail.txt'), 'w');
        for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            fprintf(fundfid1, '%s\t', mapname);
            for k=1:blocknum
                fprintf(fundfid1,'blk%d\t', k);
            end
            fprintf(fundfid1, '\r');
            for j=1:masknum
                fprintf(fundfid1, '%s\t', maskname(j, :));
                for k=1:blocknum
                    fprintf(fundfid1, '%f\t', domainvalues(j, i, k));
                end
                fprintf(fundfid1, '\r');
            end
            fprintf(fundfid1, '\r');
        end
        fclose(fundfid1);
        fundfid2=fopen(strcat(resultfolder, 'fun1detail2.txt'), 'w');
        for j=1:masknum
            fprintf(fundfid1, '%s\t', maskname(j, :));
            for k=1:blocknum
                fprintf(fundfid2,'%d\t', k-1);
            end
            fprintf(fundfid2, '\r');
            for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map
		        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
                fprintf(fundfid2, '%s\t', mapname);
                for k=1:blocknum
                    fprintf(fundfid2, '%f\t', domainvalues(j, i, k));
                end
                fprintf(fundfid2, '\r');
            end
            fprintf(fundfid2, '\r');
        end
        fclose(fundfid2);
	end
	
	% save averaged map value covered by the mask. calculate SD 
	if (flagmasking>0 & SDtype~=1)  % calculate mean and SD across blocks
        domainavg=zeros(masknum, mapnum);
        domainstd=zeros(masknum, mapnum);
        if (flaggoodstim==1)    % note: for selected stim presentation analysis, no analysis map quantatification since there is no such map for every block
            for j=1:masknum
                newdomainvalues=reshape(domainvalues(j,1:NStim,:), NStim, blocknum);
                domainavg(j,1:NStim)=(sum(newdomainvalues.*goodstim', 2)./sum(goodstim', 2))';
                % sd?
            end
        else
            domainavg=mean(domainvalues, 3); % mean across blocks (3rd dimension)
            domainstd=std(domainvalues, 0, 3); 
        end
	%    [0, p, 95 sig level] = 
        %[hypothesis, pvalue, percentile, siglevel]=ttest2(domainvalues(1,7,:),domainvalues(2,7,:), 0.05, 0);
        %pvalue
	elseif (flagmasking==2&SDtype==1) % calculate mean & SD across dots (use 'maps' matrix instead 'domainvalues')
        for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map     
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            for j=1:masknum
                dotvalue=zeros(domaindotnum(j),1);
                for k=1:domaindotnum(j)
                    L=1;
                    pixelvalue=0;
                    for x=round(xx(k, j)-domainradius(j)):round(xx(k, j)+domainradius(j))
                        for y=round(yy(k, j)-domainradius(j)):round(yy(k, j)+domainradius(j))
                            if ((x-xx(k,j))^2+(y-yy(k,j))^2)<=domainradius(j)^2
                                pixelvalue=pixelvalue+maps(y,x,i);
	%                            if (i==1)  % for testing
	%                                testmap(y,x,j)=255;
	%                            end
                                L=L+1;
                            end
                        end
                    end
                    dotvalue(k)=pixelvalue/(L-1);
                end
                domainavg(j, i)=mean(dotvalue);     % mask in column and condition in row
                domainstd(j, i)=std(dotvalue);
                if (flagsavedetailfun==1)
                    dotindextemp=dotdetail(:,1)*100+dotdetail(:,2);     % just for easy finding
                    if ~isempty(find(dotindextemp == (i*100+j)))    % in dotdetail, first colum is mapnum, second column is masknum
                        dotfid=fopen(strcat(resultfolder, mapname, '-', maskname(j,:), '.txt'),'w');
                        for ii=1:domaindotnum(j)
                            fprintf(dotfid, '%f\r',dotvalue(ii));
                        end
                        fclose(dotfid);
                    end
                end
                clear dotvalue;    
            end
        end
	else
	%    fprintf ('\r!!!Error001, check flag "flagmasking" & "SDtype"\r');
        %return;
	end
	% output
	
	
	funfid1= fopen(strcat(resultfolder, strcat('fun1', '.txt')), 'a');  % for domain avg & stdev output, 1, 2 & 3 are three different formats
	funfid2= fopen(strcat(resultfolder, strcat('fun2', '.txt')), 'a');
	funfid3= fopen(strcat(resultfolder, strcat('fun3', '.txt')), 'a');
	stdfid1 =fopen(strcat(resultfolder, strcat('std1', '.txt')), 'a');
	stdfid2 =fopen(strcat(resultfolder, strcat('std2', '.txt')), 'a');
	stdfid3 =fopen(strcat(resultfolder, strcat('std3', '.txt')), 'a');
	if flagafunout==1
        outnum=mapnum;  %output all map values
	else
        outnum=NStim;   %output only stim map values
	end
	if (flagmasking>0)
        fprintf(funfid1, '\t');
        fprintf(stdfid1, '\t');
        for i=1:outnum      %output to fun1, group by domains after all blocks
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            fprintf(funfid1, '\t%s', mapname);
            fprintf(stdfid1, '\t%s', mapname);
        end
        fprintf(funfid1, '\r');
        fprintf(stdfid1, '\r');
        for i=1:masknum       
            fprintf(funfid1, '\t%s\t', maskname(i,:));
            fprintf(stdfid1, '\t%s\t', maskname(i,:));
            for j=1:outnum
                fprintf(funfid1, '%10.7f\t', domainavg(i, j));
                fprintf(stdfid1, '%10.7f\t', domainstd(i, j));
            end
	%        fprintf(funfid1, '\r');    %put these two line on for single block output
	%        fprintf(stdfid1, '\r');
        end
        fprintf(funfid1, '\r');
        fprintf(stdfid1, '\r');
        
        fprintf(funfid2, '\t');
        fprintf(stdfid2, '\t');
        for i=1:mapnum           %output to fun2, group by SF
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            fprintf(funfid2, '\t%s', mapname);
            fprintf(stdfid2, '\t%s', mapname);
        end
        fprintf(funfid2, '\r');
        fprintf(stdfid2, '\r');
        for i=1:masknum     
            fprintf(funfid2, '\t%s\t', maskname(i,:));
            fprintf(stdfid2, '\t%s\t', maskname(i,:));
            for j=1:outnum
                fprintf(funfid2, '%10.7f\t', domainavg(i, j));
                fprintf(stdfid2, '%10.7f\t', domainstd(i, j));
            end
            fprintf(funfid2, '\r');
            fprintf(stdfid2, '\r');
        end
        fprintf(funfid2, '\r');
        fprintf(stdfid2, '\r');
	
	%    for i=1:masknum           %output to fun3, group by contrast
	%        fprintf(funfid3, '\t%s', maskname(i,:));
	%        fprintf(stdfid3, '\t%s', maskname(i,:));
	%    end
        fprintf(funfid3, '\r');
        fprintf(stdfid3, '\r');
        for i=1:outnum     
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            fprintf(funfid3, '\t%s\t', mapname);
            fprintf(stdfid3, '\t%s\t', mapname);
            for j=1:masknum
                fprintf(funfid3, '%10.7f\t', domainavg(j, i));
                fprintf(stdfid3, '%10.7f\t', domainstd(j, i));
            end
        end
        fprintf(funfid3, '\r');
        fprintf(stdfid3, '\r');
	end
    fclose(funfid1);
	fclose(stdfid1);
	fclose(funfid2);
	fclose(stdfid2);
	fclose(funfid3);
	fclose(stdfid3);
end     %if (flagquantify>0)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Filtering
if flaglpfilter | flaghpfilter>0
    SNCFilterSMap;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save sum maps
if flagmap>0
    SNCSaveSMap;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vector analysis
if flagvector==1
    Nvect=size(vect, 1);
	% calculate vector map
	for j=1:Nvect
        Sum1=zeros(FrameHeight, FrameWidth);
        Sum2=zeros(FrameHeight, FrameWidth);
        condlist1=getfield(cell2struct(vect(j, 2), 'junk'), 'junk');
        for n=condlist1
            Sum1=Sum1+maps(:,:,n);
        end
        Sum1=Sum1/size(condlist1,2);
        condlist2=getfield(cell2struct(vect(j, 3), 'junk'), 'junk');
        for n=condlist2
            Sum2=Sum2+maps(:,:,n);
        end
        if size(condlist2,2)~=0
            Sum2=Sum2/size(condlist2,2);
        end
        vmaps(:,:,j)=Sum1-Sum2;
	end
    if ~isdir([resultfolder, 'vector\'])
        mkdir(resultfolder, 'vector\');    
    end
	for i=1:Nvect
        [maptemp, framemedian, lowClip, highClip] = OIClip(vmaps(:,:,i), clipmethod, clipvalue, immask);
        maptemp = norm_to_uint8(maptemp);         
        if i<10
            savefilename=strcat('vect0', num2str(i), ext, '.bmp');  
        else
            savefilename=strcat('vect', num2str(i), ext, '.bmp');  
        end
    	savefilename=strcat(resultfolder, 'vector\', savefilename);
    	imwrite(maptemp, savefilename, 'bmp');
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(vmaps(:,:,i), savefilename);	
    end

    % Generate vector map
	t1=clock;
	for k = 1:Nvect
        if k<10
            filename=strcat('vect0', num2str(k), ext, '.ivf');  
        else
            filename=strcat('vect', num2str(k), ext, '.ivf');  
        end
       	filename=strcat(resultfolder, 'vector\', filename);
		b=OIReadIVF(filename);
        if (k==1)
		    [height, width]=size(b);
            singlemap = zeros(height, width, Nvect); 
        end
        singlemap(:, :, k) = b; %norm_to_01(b);
	end
	lut=textread('bwfix.lut');
	Polarmap=OIPolar(singlemap, lut, vectmask, vectclipsd, lowpass, highpass, filtermethod);
%	figure; image(norm_to_uint8(Polarmap.ang));  axis equal; axis off; colormap(lut);  colorbar; %%load billcolorfix; 
%	figure; imagesc(norm_to_uint8(Polarmap.mag)); axis equal; axis off; colormap(gray(256)); colorbar; 
	filename1 = strcat('Angle', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename1 = strcat(resultfolder, 'vector\', filename1);
	imwrite(norm_to_uint8(Polarmap.ang), lut, filename1, 'tiff');
%	OIWriteIVF(Polarmap.ang, [filename1(1:end-4), ext, '.ivf']);
	filename2 = strcat('Mag', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename2 = strcat(resultfolder, 'vector\', filename2);
	imwrite(norm_to_uint8(OIClip(Polarmap.mag, clipmethod, clipvalue)), filename2, 'tiff');
	OIWriteIVF(Polarmap.mag, [filename2(1:end-4), ext, '.ivf']);
	filename3 = strcat('Polar', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename3 = strcat(resultfolder, 'vector\', filename3);
    mag=norm_to_01(Polarmap.mag);
    mag=OIClip(mag, 1, polarclipsd);   % to adjust map darkness    
    ang=double(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
    polarmap=ang;
    polarmap(:,:,1)=ang(:,:,1).*mag;
    polarmap(:,:,2)=ang(:,:,2).*mag;
    polarmap(:,:,3)=ang(:,:,3).*mag;
    imwrite(norm_to_uint8(polarmap), lut, filename3, 'tiff');
%	OIWriteIVF(polarmap, [filename3(1:end-4), ext, '.ivf']);
	filename4 = strcat('Sum', num2str(lowpass), '-', num2str(highpass), filtermethod, ext, '.tif');
    filename4 = strcat(resultfolder, 'vector\', filename4);
	imwrite(norm_to_uint8(OIClip(Polarmap.sum, clipmethod, clipvalue)), filename4, 'tiff');
%	OIWriteIVF(Polarmap.sum, [filename4(1:end-4), ext, '.ivf']);
    
    % output a color table
    colorbarmap=zeros(60, 256, 3);
    for i=1:30
        colorbarmap(i, :, :)=256.*lut(:,:);
    end
    for i=1:Nvect
        x=(i-1)*floor(256/Nvect)+1;
        colorbarmap(36:60, x:x+9, 1)=256*ones(25, 10)*lut(x,1);
        colorbarmap(36:60, x:x+9, 2)=256*ones(25, 10)*lut(x,2);
        colorbarmap(36:60, x:x+9, 3)=256*ones(25, 10)*lut(x,3);
    end
    imwrite(uint8(colorbarmap), strcat(resultfolder, 'vector\colortable.tif'), 'tiff'); 
    
	fprintf('\rDone (time: %f minutes)\r', etime(clock, t1)/60);
end

% finishing up
fprintf(fidmasterlog, '\tend%d-%d-%d %d:%d.%d\n', fix(clock));
fclose(fidmasterlog);
totaltime=etime(clock, timer0);
fprintf('\r%s', finalnotice);
fprintf('\n\n\nTotal time used: %8.4f secs (%8.4f hours).\r\r', totaltime, totaltime/3600);

return;

% ToDo
%% number every masking dot for computing SD (not yet, difficult)
%% resize maps (binning)
%% high/low pass filtering
% Cliping: 20% around mean (or median)
% Polar map
% use matrix for fun.txt output
% Notes:
% 1) didn't initialize most matrixs, need do that if memory size becomes a problem, also need clean temporary vaiables
% 2) 