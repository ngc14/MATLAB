%% Sunincore	Optical Imaging Data Processing 	- HDL
%%
%%      See 'Sunin_readme.txt' for detail useage. 

%%      (090807) H3.0 SNCSaveMapH was fixed. It now clips images correctly at center value +- range
%%                    Added 'circular' and 'replicate' for PADOPT in medfilt2 (medfilt2H). The default is 'symmetric'. 091011 
%%      (090207  H2.9 Parameter 'stimselect' was introduced to select stims. This will reduce data reading time if you don't need to read all stim data.
%%        -090805)    Parameter 'gsfilename' was introduced to specify goodstim file.
%%                    Parameter 'outputdriver' was introduced to determine the output folder flexibly.
%%      (090129) H2.8 The parameters, 'flagsaveaverageblk' and 'datasource=3'  were omitted. They had been introduced since Ver H2.1 
%%                    Several post-processes (ex. histout, flagsp, flagsavespdataall, flagspttest, SNCQuantify, SNCVectorMap)
%%                      were removed from the core program to sub-programs.
%%                    Parameter 'spoperation' was introduced to specify calculation type for super pixel analysis
%%                    Parameter 'flagspcond' was introduced.
%%                    SNCFilterSMapH3.m was introduced.
%%      (081014) H2.7 'flagsavespdataall' was introduced for save trial-by-trial raw data in a text file     
%%                                  'flagsaveframeave' was introduced
%%      (080502) H2.5 'flagvectormaskfilter' was introduced.
%%                                  SNCFilterSMapH2.m was introduced.
%%      (080303) H2.3 Superpixel analysis part was updated and t-test function was introduced                  
%%      (080303) H2.2 The default input and output folder name, 'Expt' and 'expt1' was changed to 'expt0' and 'expt1', respectively.
%%      (080228) H2.1 You can averaged block file for analysis which you saved previously by 'flagsaveaverageblk=1'. 
%%                          The extention of this type of file is '.aBLK'. Set 'datasource=3'.
%%      (080130) H2.0 You can choose precision type (precisiontype='single' or 'double') . Using 'single' precision may reduce memory usage.
%%      (080115) H1.1 Add mask filtering function in SNCFilterSMapH.m (flagmaskfilter)
%%      (070502) H1.0 Add trial by trial filtering (flagtrialfilter)
%%      Since here, Sunincore34H versions were made by Hisashi, based on V3.4.

%%      V3.4   The final version by Haidong
%%                  Temporally delete 'mean.txt' and 'mapinfo.txt' output, will put back later
%%                  Add 'datasource' into suninuser so can process previous saved ivf files besides block file. 
%%                  Add 'operation' into suninusmaser
%%      V3.3   Add high-pass, low-pass filtering feature for average maps	(051128-060123)		
%%                  Also add 'tDCFrames=tDCFrames+1;' around line 520 for solving the problem that some pixels are 0 in dye imaging
%%      V3.2v add vector analysis, 051025: correct one error in sp value calculating: 
%%                  previous one is scaled by number of blocks (relative shape doesn't change)
%%      V3.2 (050519-20) Changed defination of 'flagmap' to better control output maps. 1 -> 4, less maps -> more maps
%%                  Add 'blockselect' for control which block to process, delete 'startblk' and 'startstim'
%%                  Add 'blockfilenum' to separate from 'blockselectnum', the latter one is the number in 'blockselect'
%%                  Add 'DataType' to 'anapar', and modified all blk reading function to recognize different data type (12, 13, 14)
%%                  Add 'flagsaveaverageblk' for save averaged block file, 'flagquantify' to control fun.txt output
%%                  Corrected 'goodstim' usage, was shaped wrong before.
%%                  Add function that will search under '\expfolder\masks' for *.bmp or *.txt for masking (if name is not provided)
%%      V3.1  Add randomnize stim sequence option 'flagrandom', will read '_stimseq.txt' from data folder.
%%                  Correction: Modified 'operation' part: sum --> mean
%%                  Modify OIAlign2, to include method 4 (normxcorr2.m)
%%                  Modified filenames for accummap and tempmap output 
%%      V3.0  Add 'goodstim' matrix, to select only good conditions for averaging. 
%%                  add 'stimsumnum' to keep track of stim number been summed (useful in 'flagrandom' and '..' cases
%%                  unfinished: need change the way to calculate domainavg and domainstdev

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check and create directories
% FrameRate=4;    %hz
resultdriver=outputdriver;
resultfolder=[resultdriver, outputfolder, expname, runname, resultname];
runfolder=[resultdriver, outputfolder, expname, runname];
blockfolder = strcat(datadriver, datafolder, expname, runname); 
if ~isdir([resultdriver, outputfolder])
    mkdir(resultdriver, outputfolder);    
end
if ~isdir([resultdriver, outputfolder, expname])
    mkdir([resultdriver, outputfolder], expname);    
end
if ~isdir([resultdriver, outputfolder, expname, runname])
    mkdir([resultdriver, outputfolder, expname], runname);    
end
if flagmap~=0
    if ~isdir([resultdriver, outputfolder, expname, runname, resultname])
        mkdir([resultdriver, outputfolder, expname, runname], resultname);    
    end
    if (~isdir([resultfolder, 'blkavg\']) & flagmap~=0)
        mkdir(resultfolder, 'blkavg\');    
    end
    if flagmap==2 & ~isdir([resultfolder, 'blkavg\singlecond\'])
        mkdir([resultfolder, 'blkavg\'], 'singlecond\');    
    end
    if flagmap==3 | flagmap==4
        if ~isdir([resultfolder, 'singleblk\'])
            mkdir(resultfolder, 'singleblk\');    
        end
        if ~isdir([resultfolder, 'singleblk\singlecond\'])
            mkdir([resultfolder, 'singleblk'], '\singlecond\');    
        end
    end
    if tempsaveblocks~=0 && ~isdir([resultfolder, 'tempmap\'])
        mkdir(resultfolder, 'tempmap\');    
    end
    if saveaccum~=0 && ~isdir([resultfolder, 'accummap\'])
        mkdir(resultfolder, 'accummap\');    
    end
end

if flagalign~=0 && flagalign~=41
    shiftlog = strcat(resultfolder, 'shiftlog.txt');     %for output log file during processing
    fidshiftlog = fopen(shiftlog, 'w');  %for shiftlog output
end


% default value for vector analysis, 
if flagvector>0
    vectclipsd=0;	% clip image befor calculating, to remove some extreams, may not necessory, put 0 for no clipping
%     vectmask=ones(504); % mask used for Bosking normalization
    %vectmask=imread('masks\v1mask.bmp');
    if flagvectormaskfilter
        filtervectormaskname = strcat(resultdriver, outputfolder, expname, 'masks\vectormask\default.bmp');
        vectmasktemp = imread (filtervectormaskname, 'bmp');
        if size(masktemp,3)==1
            vectmask = vectmasktemp;
        else
            vectmask = vectmasktemp(:,:,1);
        end
    else
        vectmask=ones(504); % mask used for Bosking normalization
    end
%     polarclipsd=1;  % angle*map map is usually too dark, use lower clipsd to make it brighter, used 1 for 040113GarRun2
end
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
	Polarmap=OIPolarH2(singlemap, lut, vectmask, vectclipsd, lowpass, highpass, LPfiltermethod, HPfiltermethod);
%	figure; image(norm_to_uint8(Polarmap.ang));  axis equal; axis off; colormap(lut);  colorbar; %%load billcolorfix; 
%	figure; imagesc(norm_to_uint8(Polarmap.mag)); axis equal; axis off; colormap(gray(256)); colorbar; 
	filename1 = strcat('Angle', num2str(lowpass), '-', num2str(highpass), HPfiltermethod, ext, '.tif');
	filename1 = strcat(resultfolder, 'vector\', filename1);
	imwrite(norm_to_uint8(Polarmap.ang), lut, filename1, 'tiff');
%	OIWriteIVF(Polarmap.ang, [filename1(1:end-4), ext, '.ivf']);
	filename2 = strcat('Mag', num2str(lowpass), '-', num2str(highpass), HPfiltermethod, ext, '.tif');
	filename2 = strcat(resultfolder, 'vector\', filename2);
	imwrite(norm_to_uint8(OIClipH(Polarmap.mag, clipmethod, clipvalue)), filename2, 'tiff');
%	OIWriteIVF(Polarmap.mag, [filename2(1:end-4), ext, '.ivf']);
	filename3 = strcat('Polar', num2str(lowpass), '-', num2str(highpass), HPfiltermethod, ext, '.tif');
    filename3 = strcat(resultfolder, 'vector\', filename3);
    mag=norm_to_01(Polarmap.mag);
    mag=OIClipH(mag, 1, polarclipsd);   % to adjust map darkness
    switch precisiontype
        case 'single'
            ang=single(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
        case 'double'
            ang=double(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
    end
    polarmap=ang;
    polarmap(:,:,1)=ang(:,:,1).*mag;
    polarmap(:,:,2)=ang(:,:,2).*mag;
    polarmap(:,:,3)=ang(:,:,3).*mag;
    imwrite(norm_to_uint8(polarmap), lut, filename3, 'tiff');
%	OIWriteIVF(polarmap, [filename3(1:end-4), ext, '.ivf']);
	filename4 = strcat('Sum', num2str(lowpass), '-', num2str(highpass), HPfiltermethod, ext, '.tif');
    filename4 = strcat(resultfolder, 'vector\', filename4);
	imwrite(norm_to_uint8(OIClipH(Polarmap.sum, clipmethod, clipvalue)), filename4, 'tiff');
%	OIWriteIVF(Polarmap.sum, [filename4(1:end-4), ext, '.ivf']);
	fprintf('\rDone (time: %f minutes)\r', etime(clock, t1)/60);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determine data source & check availability
filtermaskname = strcat(resultdriver, outputfolder, expname, 'masks\filtermask\default.bmp');
clipmaskname = strcat(resultdriver, outputfolder, expname, 'masks\clipmask\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2

if flagmaskfilter
    filtermasktemp = (double(imread(filtermaskname, 'bmp'))/255);
    if size(filtermasktemp,3)==1
%         fprintf('size 1');
        filtermask = filtermasktemp;
    else
%         fprintf('size 3');
        filtermask = filtermasktemp(:,:,1);
    end
    clear filtermasktemp;
%     filtermask=filtermask*(-1)+1; %090204 Hisashi
%     fprintf('filtermask: max %g\r', max(max(filtermask)));
%     fprintf('filtermask: min %g\r', min(min(filtermask)));
end

if flaglpfilter==0
    lpkernelsize=0; %090204 Hisashi
end
if flaghpfilter==0
    hpkernelsize=0; %090204 Hisashi
end

if (clipmethod==2)      %use masking (to exclude blood vessels or only include certain area)
    cpmasktemp = imread (clipmaskname, 'bmp');
    %        cpmask = cpmask (cropsize+1:end-cropsize, cropsize+1:end-cropsize);
    if size(cpmasktemp,3)==1
%         fprintf('size 1');
        cpmask = cpmasktemp;
    else
%         fprintf('size 3');
        cpmask = cpmasktemp(:,:,1);
    end
    clear cpmasktemp;    
else 
    cpmask = 0;
end

% [stim]=textread(strcat(runfolder, 'conditionname.txt'), '%s') %By Hisashi 080310
% % [smap1, smap2, smap3]=textread(strcat(runfolder, 'subtractionmap.txt'), '%s,[%u %u] %q') %Hisashi 071023
% % temp=textread(strcat(runfolder, 'subtractionmap.txt'), '%s') %By Hisashi 080310
% % size(temp,1)
% fid=fopen(strcat(runfolder, 'subtractionmap.txt'));
% i=1;
% temp='aaa';
% while temp~=-1
%     temp = fgetl(fid)
%     if temp~=-1
%         smap(i,:)=temp
%         i=i+1;
%     end
% end
% fclose(fid);

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
%         for i=1:size(filename,1)
%             fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
%         end
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
        SNCFilterSMapH2;
    end
    if flagmap>0
        SNCSaveSMap;
    end
    fprintf('\rProcess finished\r');
    return;
% elseif datasource==3 %added by Hisashi
%     %check;
% %     if isempty(filename) % Search for block files if filenames are not provided
% %         if system=='v'
% %             tempfilename=struct2cell(dir([blockfolder, '*.aBLK']));
% %         elseif system=='r'
% %             fprintf('Note: Averaged block file is available only as VDAQ format\n');
% %         end
% %         filename=sort(tempfilename(1,:)');
% %         for i=1:size(filename,1)
% %             fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
% %         end
% %         fprintf('\nFound %d aBLK files(sorted, check sequence).\n', size(filename,1));
% %     end
%     
%     if system=='v'
%         tempfilename=strcat(blockfolder, ablkname, '_avg', '.aBLK');
%     elseif system=='r'
%         fprintf('Note: Averaged block file is available only as VDAQ format\n');
%     end
%     filename=tempfilename(1,:)';
%     for i=1:size(filename,1)
%         fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
%     end
%     fprintf('\nFound %d aBLK files(sorted, check sequence).\n', size(filename,1));
%     blockfilenum=size(filename, 1);      % how many blocks
end
if flagloadspvalues
    flagmap=0;
    flagvector=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read  info
anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(1), 'junk'), 'junk')), system);
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
% if isempty(superpixel)
%     flagsp=0;
% end
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

% if (isempty(blockselect) | datasource ==3)
if isempty(blockselect)
    blockselect=1:blockfilenum;
end
if isempty(stimselect)
    stimselect=1:NStim;
end
% if isempty(ivfblockselect)
%     ivfblockselect=1:blockfilenum;
% end
blockselectnum=size(blockselect,2);
stimselectnum=size(stimselect,2);
fprintf('Selected %d blk files.\r', blockselectnum);
fprintf('Selected %d stims\n', stimselectnum');
% stimselect
fprintf('\r');
%if isempty(flagsaveprcblk)
%     flagsaveprcblk=0;
% end


subtraction = 1;    % 0: division, 1: subtraction   (division to be implemented)


% === shifting ===
flagdarksub = 0;    % whether or not subtract dark frame
%cropsize = 0;   % how many pixel crop off along 4 edges, should be used only when doing image alignment (which causes 'blank' margins), (currently not used), may couse problem in clipping

% === data/map saving ===
flagsavedata=0;    % 1: save data file inaddition to the maps, note: slow
SDtype=0;          % 0: SD is across blocks, 1: SD is across different dots (only for only for 'flagspmask=2', i.e. dots from coordinates)
flagafunout=0;     % controls fun1-3 output: 0: quantify only single-condition map, 1: quantify both single-condition and analysis maps
flagsavedetailfun = 0;  % save detailed averaged (across different spots) domain values for each stim and block.  0: no detail save
dotdetail=[             % first colum is mapnum (single cond), second column is masknum
      ];     
if flagquantify~=0
    if flagspmask==0
        fprintf('Error: Need masks for quantification (flagquantify~=0 and flagspmask==0)!\n');
    end
    % open fun file here
end


% === masking/Quantifying ===
%savemaskmap=1;     % 1: save maskmap based on coordinates (only for 'flagspmask=2')
if flagspmask>0 
    maskfolder=strcat(resultdriver, outputfolder, expname, 'masks\spmask\');
    switch flagspmask
    case 1      % for map-loading, Mask bmp files are 0-1 monochrome bmp files, desired areas are 1's
        if isempty(spmaskname)        
            tempfilename=struct2cell(dir([maskfolder, '*.bmp']));
            maskfilename=sort(tempfilename(1,:)');
            masknum=size(maskfilename, 1);
            if masknum == 0
                fprintf('Error: Need spmask in spmask folder!\n');
            end
            for i=1:masknum
                spmasknametemp(i,:)=getfield(cell2struct(maskfilename(i), 'junk'), 'junk');
            end
            spmaskname=char(spmasknametemp(:,1:end-4));         % somehow this line can't be put together with the previous i loop
        else
            masknum=size(spmaskname,1);
        end
        for i=1:masknum
            strcat(maskfolder, spmaskname(i,:));
            masktemp=imread (strcat(maskfolder, spmaskname(i,:), '.bmp'), 'bmp');
            if size(masktemp,3)==1
                mask(:,:,i) = masktemp;
            else
                mask(:,:,i) = masktemp(:,:,1);
            end
        end
        mask=uint8(mask);   % need check if mask is binary.
    case 2     % masks are from coordinates, create masks in memory here
        if isempty(spmaskname)
            tempfilename=struct2cell(dir([maskfolder, '*x.txt']));  % only look for x.txt
            maskfilename=sort(tempfilename(1,:)');
            masknum=size(maskfilename, 1);
            dotsinmap=200;       % assume maximum 100 dots in each domain map, only for 'flagmask==2'
            xx=zeros(dotsinmap, masknum);       % store the dot center coordinates          
            yy=zeros(dotsinmap, masknum);
            for i=1:masknum
                spmasknametemp(i,:)=getfield(cell2struct(maskfilename(i), 'junk'), 'junk');
            end
            spmaskname=char(spmasknametemp(:,1:end-6));
        else
            masknum=size(spmaskname,1);
        end
%         if isempty(domainradius)
%             domainradius=zeros(masknum)+5;
            domainradius=zeros(masknum)+domainradiussize;
%         end                
        for i=1:masknum
            domainfidx=fopen (strcat(resultdriver, outputfolder, expname, 'masks\spmask\', spmaskname(i,:), '_x.txt'), 'r'); %Hisashi
            domainfidy=fopen (strcat(resultdriver, outputfolder, expname, 'masks\spmask\', spmaskname(i,:), '_y.txt'), 'r'); %Hisashi
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
                imwrite(uint8(maskmap), strcat(resultdriver, outputfolder, expname, 'masks\spmask_', spmaskname(i,:), ext, '.bmp'), 'bmp');
            end
        end
        clear maskmap;
    end
maskfilename
fprintf('\nFound %d masks for superpixel anslysis.\n', masknum);
end

flagpixelintensityoutpt=0;
if flagpixelintensityoutpt>0 
    pixelmaskfolder=strcat(resultdriver, outputfolder, expname, 'masks\pixelmask\');
	% for map-loading, Mask bmp files are 0-1 monochrome bmp files, desired areas are 1's
    if isempty(pixelmaskname)        
        tempfilename=struct2cell(dir([pixelfolder, '*.bmp']));
        pixelmaskfilename=sort(tempfilename(1,:)');
        pixelmasknum=size(pixelmaskfilename, 1);
        if pixelmasknum ~= 0
            for i=1:pixelmasknum
                pixelmasknametemp(i,:)=getfield(cell2struct(pixelmaskfilename(i), 'junk'), 'junk');
            end
            pixelmaskname=char(pixelmasknametemp(:,1:end-4));         % somehow this line can't be put together with the previous i loop
        else
            pixelmasknum=1;
            pixelmaskname='whole'
        end
    else
        pixelmasknum=size(pixelmaskname,1);
    end
    if pixelmasknum ~= 0
        for i=1:pixelmasknum
            strcat(pixelmaskfolder, pixelmaskname(i,:));
            pixelmasktemp=imread (strcat(pixelmaskfolder, pixelmaskname(i,:), '.bmp'), 'bmp');
            if size(pixelmasktemp,3)==1
                pixelmask(:,:,i) = pixelmasktemp;
            else
                pixelmask(:,:,i) = pixelmasktemp(:,:,1);
            end
        end
    else
        pixelmask(:,:,1) = ones(height, width)*255;
    end
    pixelmask=uint8(pixelmask);   % need check if mask is binary.
    pixelmaskfilename
    fprintf('\nFound %d masks for pixel-by-piexel intensity output.\n', masknum);
end
mapinfo=zeros(3, blockselectnum, mapnum);      %save (median, lowclip, highclip) for each maps, blknum, mapnum. /What is the purpose of this line? Hisashi 

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
nobvmask = strcat(resultdriver, outputfolder, expname, 'masks\nobvlarge.bmp');     %used when histindex(1,1)=0, which means calculate histograms on nobv rigion, if no such bmp exist, calculate the whole map and give warning by the end.
% need improve, delete nobvmask

if (flagrandom==1);
    stimseqfile=[runfolder, '_stimseq.txt'];
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
    
% if (flaggoodstim==1 & datasource~=3) % load good stim condition file
if (flaggoodstim==1) % load good stim condition file
%     if gsfilename=='';
%         gsfilename='goodstim.txt';
%     end
	fidgoodstim=fopen(strcat(runfolder, gsfilename), 'r'); %Hisashi 071023
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
%     if prod(sum(goodstim, 2))==0
    if prod(sum(goodstim, 1))==0 % By Hisashi 081111
        fprintf('Worning: too few good stims, one or more condition has no data!\n');
    end
	fprintf('\r** Good stim parameters loaded **\r');
    if flagmap>=3   % user reqires single block maps?
        fprintf('Warning: ''goodstim=1'' and ''flagmap>=3'', you may not get all single block maps, anykey to continue...');
        pause;
    end
else
    goodstim=ones(blockfilenum, NStim);
end	


% check
if (flagspmask==1 & SDtype==1)
	fprintf('ERROR: check flagspmask and SDtype\r');
end
% if flaggoodstim==1 & flagtrialfilter==1 % Commented out by Hisashi 071023
% 	fprintf('Can not perform trial-by-trial filtering if using goodstim, flagtrialfilter disabled\r')
% 	flagtrialfilter=0;
% end
finalnotice='Process Finished: ';    % Programmer: add notice into this string, it will be displayed at the end of the calculation, use <finalnotice=strcat(finalnotice, 'addhere')>
timer0=clock;       % used for calculate how much time this program runing

% Log the parameters for this process secession.
% fidmasterlog = fopen(strcat(resultfolder, 'masterlog.txt'), 'a');  
%fprintf(fidmasterlog, '%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\tstart%d-%d-%d %d:%d.%d\r',ext, blockfolder, blockselectnum, operation, flagalign, shiftrange(1), shiftrange(2), shiftrange(3), shiftrange(4), shiftstep1, clipmethod, clipvalue, sframe, fframe, bframe, flagspmask, fix(clock));

% Open & Initiate log files
%fprintf(fidshiftlog, '__________Start Date: %d %d %d    Time: %d %d %d\r', fix(clock));

%Read Header info from one block file

% Other variables used
%        DCframes[y,x,j]: Stores DC frames in one trial
%        fmShifts[j, i, k]: shifts for best alignment
%

% if do shift alignment, readin the very first frame
if flagalign==1
    frameone=OIReadFrame(strcat(blockfolder, firstframe), 'v', 1, 1);
%    shiftmask = imread(strcat(resultdriver, outputfolder, expname, 'masks\shiftmask1.bmp'), 'bmp');     % bv mask for shifting method3
end
if flagalign==11    %% Load shift parameters
    [sfname, sfblock, sfstim, sfframe, sfx, sfy, sfcoor]=textread(strcat(resultfolder, 'shiftinput.txt'), '%s %d %d %d %f %f %f');
    if size(sfname,1)~=blockselectnum*NStim*FramesPerStim
        fprintf('\r Wrong shift parameternum');
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blockselectnum
        for i=1:NStim
            for j=1:FramesPerStim
                if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim*FramesPerStim+(i-1)*FramePerStim+j).sfname
                    fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname((k-1)*NStim*FramesPerStim+(i-1)*FramePerStim+j).sfname);
                end
            end
        end
    end
end
% if (flagalign==41 & datasource~=3)  % Load shift parameters generated by flagalign==4 process (include '.1' file for all first frame shift and '.3' file for all rest frames)
if (flagalign==41)  % Load shift parameters generated by flagalign==4 process (include '.1' file for all first frame shift and '.3' file for all rest frames)
    [sfname, sfblock, sfstim, sfframe, sfx1, sfy1, sfcoor]=textread(strcat(runfolder, 'shiftinput.2.txt'), '%s %d %d %d %f %f %f'); %Hisashi 071023
    fprintf('\r** Shift input parameters loaded **\r');
    if size(sfname,1)~=blockfilenum*NStim
        fprintf('Error: flagalign===41, .2 file doesnot contain blockfilenum*NStim number of entries!\n');
        size(sfname,1)
        blockfilenum
        NStim
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blockselectnum
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
%     for k=1:blockselectnum
%         for i=1:NStim
%             for j=1:FramesPerStim
%                 if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramePerStim-1)+j).sfname
%                     fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramePerStim-1)+j).sfname);
%                 end
%             end
%         end
%     end
end

if (flaghpfilter||flaglpfilter)
    fprintf('\r** Spatial filtering **');
    if flaglpfilter
        fprintf('\rImages will be low-pass filtered with size %g \''%s\''.', lpkernelsize, lpfmethod);
    end
    if flaghpfilter
        fprintf('\rImages will be high-pass filtered with size %g \''%s\''.', hpkernelsize, hpfmethod);
    end
    if flagmaskfilter
        fprintf('\rMask will be applied for high-pass filtering.');
    end
    if flagtrialfilter
        fprintf('\rSpatial filtering will be applied trial-by-trial.');
    end
    if (flagsp && (flaghpfilter||flaglpfilter))
        fprintf('\rSpatial filtering will be applied frame-by-frame for super-pixel analysis.\r');
    end
else
    fprintf('\r** No spatial filtering **\r');
end
        

% Print out important information from the block header: width, height in terms of pixels; frame number per stim, stim numbers
fprintf('\r');
if flagmap
    fprintf('Results of imaging map will be saved in   "%s"\n', resultfolder);
end
if flagsp
    fprintf('Results of superpixel analysis will be saved in   "%s"\n', [runfolder, SPresultname]);
end
if flagvector
    fprintf('Results of vector map will be saved in   "%s"\n', [runfolder, vectorresultname]);
end
fprintf('\rFrameWidth=%d FrameHeight=%d NFramesPerStim=%d NStimuli=%d\n',FrameWidth, FrameHeight, FramesPerStim, NStim);

%*********************************************************************************************************************
%********************************        Start Process     ***********************************************************
%*********************************************************************************************************************

fprintf('\rWait...');
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

if flagsp
    if flagsaveframeave
        framenum=FramesPerStim+2;
    else
        framenum=FramesPerStim;
    end
end
if flagspcond==1
    NSPcond=size(SPcond, 1);
else
    NSPcond=0;
end
% if flagspsubcond==1
%     NSPsubcond=size(SPsubcond, 1);
% else
%     NSPsubcond=0;
% end

switch precisiontype %by Hisashi
    case 'single'
        if (flagspmask>0) 
            domainvalues=zeros(masknum, mapnum, blockselectnum, 'single');    %store average value for each mask, each stim, and each block, can expand to include analysis map values
        end
        if (flagmap~=0 && flaggoodstim ~=1) || (flagvector~=0 && flaggoodstim ~=1) || (flagmap==3) || (flagmap==4) ||(flagquantify~=0)
            maps=zeros(FrameHeight, FrameWidth, mapnum, 'single');     %for average maps
        end
        if tempsaveblocks~=0
            tempmaps=zeros(FrameHeight, FrameWidth, mapnum, 'single');     %for tempory maps
        end
        if (flaggoodstim==1)
            NewSum=zeros(FrameHeight,FrameWidth,NStim, 'single');		  % a new way to calcultate average trials (can be different number for different stim)
        end
        
%         stimsumnum=zeros(NStim,1);	%keep how many stim was summed. (in randomized & goodstim cases, it may not same for each condition).
        stimsumnum=zeros(NStim+NSPcond,1);	%keep how many stim was summed. (in randomized & goodstim cases, it may not same for each condition).
        if flagsp
            spvalues=NaN(framenum,masknum,NStim,blockselectnum, 'single');
            if (flaghpfilter||flaglpfilter)
                spvalues2=NaN(framenum,masknum,NStim,blockselectnum, 'single');
            end
        end
    case 'double'
        if (flagspmask>0) 
            domainvalues=zeros(masknum, mapnum, blockselectnum);    %store average value for each mask, each stim, and each block, can expand to include analysis map values
        end
        if (flagmap~=0 && flaggoodstim ~=1) || (flagvector~=0 && flaggoodstim ~=1)  || (flagmap==3) || (flagmap==4) ||(flagquantify~=0)
            maps=zeros(FrameHeight, FrameWidth, mapnum);     %for average maps
        end
        if tempsaveblocks~=0
            tempmaps=zeros(FrameHeight, FrameWidth, mapnum);     %for tempory maps
        end
        if (flaggoodstim==1)
            NewSum=zeros(FrameHeight,FrameWidth,NStim);		  % a new way to calcultate average trials (can be different number for different stim)
        end
 %         stimsumnum=zeros(NStim,1);	%keep how many stim was summed. (in randomized & goodstim cases, it may not same for each condition).
        stimsumnum=zeros(NStim+NSPcond,1);	%keep how many stim was summed. (in randomized & goodstim cases, it may not same for each condition).
        if flagsp
            spvalues=NaN(framenum,masknum,NStim,blockselectnum);
            if (spoperation==2 || spoperation==3)
                spvaluesdrr=NaN(framenum,masknum,NStim,blockselectnum);
%                 size(spvalues)
%                 size(spvaluesdrr)
            end
%             if (flaghpfilter||flaglpfilter)
%                 spvaluesdrr=NaN(framenum,masknum,NStim,blockselectnum);
%             end
        end
end

timertrigger=1;
ctime=clock;
ctime0=clock;
flagfirst=0;
% Start read/process block by block*************
fprintf('total block number = %d\r\n', blockselectnum);
if (flagloadspvalues~=1)
    for k=blockselect  % start block by block processing
        if blockselectnum>1
            if sum(stimsumnum)>0 && flagfirst==0 % && k==blockselect(2)    % just for estimate process time 
                blockproctime=etime(clock, ctime0);
                fprintf('\rTime for one block: %7.4f secs', blockproctime);
                temptime=blockproctime/sum(stimsumnum)*sum(sum(goodstim))*(stimselectnum/NStim)*(blockselectnum/blockfilenum);
                fprintf('\rTotal time for all (%d) blocks estimated: %5.2f minutes (%5.2f hours)\r', blockselectnum, temptime/60, temptime/3600);
                flagfirst=1;
            end
        end
        SumFrame=zeros([FrameHeight, FrameWidth, NStim]);
        fprintf('\rblock#%d(%s): ',k, getfield(cell2struct(filename(k), 'junk'), 'junk'));
        fid=fopen(strcat(blockfolder, getfield(cell2struct(filename(k), 'junk'), 'junk')),'rb','ieee-le'); % open block file
        fseek(fid, headerlength,'bof');    % skip header
        % Here start read/process each stim in block 'k'
        for i = 1:NStim
            if goodstim(k, i)==0 || sum(find(stimselect==i))==0  % skip this stim; k: blockselect; i: NStim. Since H029, skip stims unselected by 'stimselect'
            fseek(fid, stimlength,'cof');    
            else
            if flagrandom   % if stim seq in data is random, use i for reading, and use 'currentstimid' for averageing
                currentstimid=stimseq(k, i);
            else
                currentstimid=i;
            end
            fprintf(' stim%d ', currentstimid);
            switch precisiontype %by Hisashi
            case 'single'
                switch system
                case 'v'    % for VDAQ data 
                    switch DataType
                    case 12
                        tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>float32');   % read a long vector of integers from block file, ('t' is for temp, may need shifting); By Hisashi on 080130
                    case 13
                        tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint32=>float32');   % read a long vector of integers from block file, ('t' is for temp, may need shifting); By Hisashi on 080130
                    case 14
                        tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'float32=>float32');   % read a long vector of integers from block file, ('t' is for temp, may need shifting); By Hisashi on 080130
                    otherwise
                        fprintf('Error: data type must be 12, 13 or 14!\n');
                        return;
                    end
                    tDCFrames=reshape(tDCFrames, FrameWidth, FrameHeight, FramesPerStim);     % reshape the long vector into 3-D matrix
                    tDCFrames=permute(tDCFrames, [2,1,3]); % Since matlab read column first. Same purpose as the old method (flip&rotate). Note, the width and height changed
                case 'r'    % for RedShirt data
                    tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>float32');  % Redshirt use 'uint16' data type; By Hisashi on 080130
                    tBNCFrames=fread(fid,[8*FramesPerStim], 'uint16');   %NBNC=8 Read 8 BNCs
                    tDarkFrames=fread(fid,[FrameWidth*FrameHeight],'uint16');    % read dark frames
                    tBNCDark=fread(fid, [8], 'uint16');  % read BNC dark 
                    tDCFrames=reshape(tDCFrames, FramesPerStim, FrameWidth, FrameHeight); % reshape into 3-D matrix
                    tDCFrames=permute(tDCFrames, [3,2,1]);    % so the dimension order is (height, width, frame)
                case 'd'
                     tDCFrames=fread(fid,[FrameWidth*FrameHeight*FramesPerStim],'uint16=>float32');   % read a long vector of integers from block file, ('t' is for temp, may need shifting); By Hisashi on 080130
                    tDCFrames=reshape(tDCFrames, FrameWidth, FrameHeight, FramesPerStim);     % reshape the long vector into 3-D matrix
                    tDCFrames=permute(tDCFrames, [2,1,3]); % Since matlab read column first. Same purpose as the old method (flip&rotate). Note, the width and height changed
                end
            case 'double'
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
%             if datasource~=3
            switch flagalign
            case 0 %if (flagalign==0)   % no shifting
                DCFrames = tDCFrames;   % how much time will it cost for this operation?
            case 1 %elseif (flagalign==1)        % shifting 1: all aligned with 'frameone'
                pixoff=[0 0]; 
                maxcor=0;
                % note: fmshifted frame is the same size as before shifting, may contain zeros due to the shifting
                for j = 1:FramesPerStim            
                    if timertrigger==1|timertrigger==2  % estimate time
                        [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blockselectnum, NStim*FramesPerStim);
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
                            [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blockselectnum, NStim*FramesPerStim);
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
%             else
%                 DCFrames = tDCFrames;   % By Hisashi
%             end
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
            case 2      % 2: dR/R0: percentage change, R0 is the first frame of each stim condition; default
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
%             if flagtrialfilter && (flaghpfilter||flaglpfilter) && datasource~=3 %added by Hisashi on 080229
            if flagtrialfilter && (flaghpfilter||flaglpfilter) %added by Hisashi on 080229
                if flagmaskfilter == 0
                    SumFrame(:,:,currentstimid)=OIEasyFilterH2(SumFrame(:,:,currentstimid), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize); %modified by Hisashi	
                else
                    SumFrame(:,:,currentstimid)=OIEasyFilterH2wMask(SumFrame(:,:,currentstimid), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize, filtermask); %modified by Hisashi	
                end
            end  
%             if flagalign==41 && datasource~=3 % Modified by Hisashi on 080229
            if flagalign==41 % Modified by Hisashi on 080229
                SumFrame(:,:,currentstimid)=OIShift(SumFrame(:,:,currentstimid), sfx1((k-1)*NStim+i), sfy1((k-1)*NStim+i));
            end

            if (flaggoodstim==1)
                NewSum(:,:,currentstimid)=NewSum(:,:,currentstimid)+SumFrame(:,:,currentstimid);    % Add this sum because of use of 'goodstim'
            end
            stimsumnum(currentstimid)=stimsumnum(currentstimid)+1;     % record how many repeats for each stim

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Pixel time course recording
            if (flagsp>0)
%                 if (flaghpfilter||flaglpfilter) && datasource~=3 %added by Hisashi on 080307
                if (flaghpfilter||flaglpfilter)
                    % Note :: For superpixel analysis, filtering is applied trial-by-trial (080307 Hisashi)
                    for f=1:FramesPerStim
                        if  (spoperation==1 || spoperation==3)
                            if flagmaskfilter == 0
                                DCFrames3(:,:,f)=OIEasyFilterH2(DCFrames(:,:,f), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize); %modified by Hisashi
                            else
                                DCFrames3(:,:,f)=OIEasyFilterH2wMask(DCFrames(:,:,f), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize, filtermask); %modified by Hisashi
                            end
                        end
                        if  (spoperation==2 || spoperation==3)
                            DCFramesfframe=mean(DCFrames(:,:,fframe),3);
                            DCFrames3drr(:,:,f)=(DCFrames(:,:,f)-DCFramesfframe)./DCFramesfframe; %dRR
                            if flagmaskfilter == 0
                                DCFrames3drr(:,:,f)=OIEasyFilterH2(DCFrames3drr(:,:,f), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize); %modified by Hisashi
                            else
                                DCFrames3drr(:,:,f)=OIEasyFilterH2wMask(DCFrames3drr(:,:,f), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize, filtermask); %modified by Hisashi
                            end
                        end
                    end
                else
                    DCFrames3=DCFrames;
                end
%                 if flagalign==41 && datasource~=3 %added by Hisashi on 080307
                if flagalign==41 %added by Hisashi on 080307
                    for j = 1:FramesPerStim
                        if  ~(flaghpfilter||flaglpfilter) || ((flaghpfilter||flaglpfilter) && (spoperation==1 || spoperation==3))
                            DCFrames3(:,:,j) = OIShift(DCFrames3(:,:,j), sfx1((k-1)*NStim+i), sfy1((k-1)*NStim+i));
                        end
                        if (flaghpfilter||flaglpfilter) && (spoperation==2 || spoperation==3)
                            DCFrames3drr(:,:,j) = OIShift(DCFrames3drr(:,:,j), sfx1((k-1)*NStim+i), sfy1((k-1)*NStim+i));
                        end
                    end
                end
                for j=1:masknum
                    for f=1:framenum
                        if  ~(flaghpfilter||flaglpfilter) || ((flaghpfilter||flaglpfilter) && (spoperation==1 || spoperation==3))
                            if f<=FramesPerStim
                                switch precisiontype
                                    case 'single'
                                        maskedmap= single(DCFrames3(:,:,f)).*single(mask(:, :, j));
                                        spvalues(f,j,currentstimid,stimsumnum(currentstimid))=sum(sum(maskedmap))/sum(sum(single(mask(:,:,j))));
                                    case 'double'
                                        maskedmap= double(DCFrames3(:,:,f)).*double(mask(:, :, j));
                                        spvalues(f,j,currentstimid,stimsumnum(currentstimid))=sum(sum(maskedmap))/sum(sum(double(mask(:,:,j))));
                                end
                            elseif f==FramesPerStim+1
                                spvalues(f,j,currentstimid,stimsumnum(currentstimid))=mean(squeeze(spvalues(fframe,j,currentstimid,stimsumnum(currentstimid))));
                            elseif f==FramesPerStim+2
                                spvalues(f,j,currentstimid,stimsumnum(currentstimid))=mean(squeeze(spvalues(sumrange,j,currentstimid,stimsumnum(currentstimid))));
                            end
                        end
                        if (flaghpfilter||flaglpfilter) && (spoperation==2 || spoperation==3)
                            if f<=FramesPerStim
                                switch precisiontype
                                    case 'single'
                                        maskedmap2= single(DCFrames3drr(:,:,f)).*single(mask(:, :, j));
                                        spvaluesdrr(f,j,currentstimid,stimsumnum(currentstimid))=sum(sum(maskedmap2))/sum(sum(single(mask(:,:,j))));
                                    case 'double'
                                        maskedmap2= double(DCFrames3drr(:,:,f)).*double(mask(:, :, j));
                                        spvaluesdrr(f,j,currentstimid,stimsumnum(currentstimid))=sum(sum(maskedmap2))/sum(sum(double(mask(:,:,j))));
                                end
                            elseif f==FramesPerStim+1
                                spvaluesdrr(f,j,currentstimid,stimsumnum(currentstimid))=mean(squeeze(spvaluesdrr(fframe,j,currentstimid,stimsumnum(currentstimid))));
                            elseif f==FramesPerStim+2
                                spvaluesdrr(f,j,currentstimid,stimsumnum(currentstimid))=mean(squeeze(spvaluesdrr(sumrange,j,currentstimid,stimsumnum(currentstimid))));
                            end
                        end
                    end              
                end
                if  ~(flaghpfilter||flaglpfilter) || ((flaghpfilter||flaglpfilter) && (spoperation==1 || spoperation==3))
                    clear DCFrames3;
                end
                if (flaghpfilter||flaglpfilter) && (spoperation==2 || spoperation==3)
                    clear DCFrames3drr;
                end
            end    %% (flagsp>0)
          end     %% if goodstim(k, i)==0
        end    %%i = 1:NStim    
        status = fclose(fid);
%         if (flagalign~=0 & datasource~=3)
        if (flagalign~=0 && flagalign~=41)
            fclose(fidshiftlog);
            fidshiftlog = fopen(shiftlog, 'a'); %close and reopen to save current results
        end


        if (flagmap~=0 && flaggoodstim ~=1) || (flagvector~=0 && flaggoodstim ~=1)  || (flagmap==3) || (flagmap==4) ||(flagquantify~=0)
            singlemap(:,:,1:NStim)=SumFrame(:,:,1:NStim);   % use singlemap to save both single condition and subtraction map from this block
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Generate subtraction maps    
        if (flagmap~=0)
            for j=1:size(Smap, 1)						% note: not differentiate goodstim=1 or 0, get single block map anyway
                Sum1=zeros(FrameHeight, FrameWidth);
                Sum2=zeros(FrameHeight, FrameWidth);
                condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
%                 j
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
                if (flagmap~=0 && flaggoodstim ~=1) || (flagmap==3) || (flagmap==4) ||(flagquantify~=0)
                    singlemap(:,:,NStim+j)=Sum1-Sum2;  
                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sum over blocks
    %         if (flaggoodstim==1)    % do not sum subtraction map (will re-calculate later), only sum single condition map
        % 	    for i = 1:NStim
        %             maps(:,:,i) =  maps(:,:,i) + singlemap(:,:,i).*goodstim(k, i);  % only good stim are sumed (need work on avg map.
        %         end
    %         else
            if (flagmap~=0 && flaggoodstim ~=1) || (flagvector~=0 && flaggoodstim ~=1)  || (flagmap==3) || (flagmap==4) ||(flagquantify~=0)
                for i = 1:mapnum   
                    maps(:,:,i) =  maps(:,:,i) + singlemap(:,:,i); 
                    if tempsaveblocks~=0
                        tempmaps(:,:,i) =  tempmaps(:,:,i) + singlemap(:,:,i); 
                    end
                end
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  save single block maps for each blocks
            if flagmap==3
                for i=NStim+1:mapnum
                    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
                    [singlemap(:,:,i), mapinfo(1, k, i), mapinfo(2, k, i), mapinfo(3, k, i)]=OIClipH(singlemap(:,:,i), clipmethod, clipvalue, cpmask);       % mapinfo will be write into 'mapinfo.txt'
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
                    [singlemap(:,:,i), mapinfo(1, k, i), mapinfo(2, k, i), mapinfo(3, k, i)]=OIClipH(singlemap(:,:,i), clipmethod, clipvalue, cpmask);       % mapinfo will be write into 'mapinfo.txt'
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
                        switch precisiontype
                            case 'single'
                                maskedmap= singlemap(:,:,i).*single(mask(:, :, j));
                                domainvalues(j,i,k)=sum(sum(maskedmap))/sum(sum(single(mask(:,:,j)))); % avg value for block 'k', stim 'i', and mask 'j'.
                            case 'double'
                                maskedmap= singlemap(:,:,i).*double(mask(:, :, j));
                                domainvalues(j,i,k)=sum(sum(maskedmap))/sum(sum(double(mask(:,:,j)))); % avg value for block 'k', stim 'i', and mask 'j'.
                        end
                    end
                end
                clear maskedmap;
            end
        end    %for if (flagmap~=0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save temp sum maps & data (for every x blocks, then clear matrix), useful for animal status changing)
        if (mod(k, tempsaveblocks)==0)  % set tempsaveblocks=0 or 999 to avoid save tempmaps
            for i=1:mapnum    
               mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
               if (flagsavedata) 
                   tempname = strcat(resultfolder, 'tempdata\', num2str(i), '_k', num2str(k), '.txt');
                   dlmwrite(tempname, tempmaps(:,:,i), '\t'); 
               end
               [tempmaps(:,:,i), framemedian, lowClip, highClip]=OIClipH(tempmaps(:,:,i), clipmethod, clipvalue, cpmask);
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
                [tempaccu, framemedian, lowClip, highClip]=OIClipH(maps(:,:,i), clipmethod, clipvalue, cpmask);
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
    end % for k1:blockselectnum
end % if 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Done with file reading and initial processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tempmaps;
clear DCFrames;
% if (flagalign~=0 && datasource~=3)
if (flagalign~=0  && flagalign~=41)
    fclose(fidshiftlog);
end

% map
if ~((flagmap~=0 && flaggoodstim ~=1) || (flagvector~=0 && flaggoodstim ~=1) || (flagmap==3) || (flagmap==4) ||(flagquantify~=0))
    switch precisiontype %by Hisashi
        case 'single'
            maps=zeros(FrameHeight, FrameWidth, mapnum, 'single');
        case 'double'
            maps=zeros(FrameHeight, FrameWidth, mapnum);     %for average maps
    end
end
if (flagmap~=0) || (flagvector~=0)     
    if (flaggoodstim==1)        % if use goodstim, avg map is from averaged single condition map,  if goodstim=0, use averaged subtraction map, should be the same. 
    	for i=1:NStim
    		maps(:,:,i)=NewSum(:,:,i)/stimsumnum(i); %maps(:,:,i)/sum(goodstim(k,:));
    	end
    	% for 'flaggoodstim==1', need calculate subtraction map based on selected stims
        if flagmap~=0
            SNCCalculateSMap;
        end
    else 	%         
        for i=1:mapnum
            maps(:,:,i)=maps(:,:,i)/blockselectnum;
        end
    end
end
clear NewSum;
clear SumFrame;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output mapinfo
fprintf('\r\r********************\r');
fprintf('**** Start analyses ****\r');
fprintf('********************\r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output histograms (pixel gray value distributions) for each map
if (histout==2) % need implement histout==1 (or delete)
    fprintf('\r** Output histograms **\r');
    SNCOutputHist
end			
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Pixel time course output, only for 'flagsp>0' and 'flagspmask>0'
if (flagsp>0)
    fprintf('\r** Super Pixel time course output **\r');
    SNCSPTimeCourseH
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Super Pixel time course all trials output, only for 'flagsavespdataall>0' and 'flagspmask>0'
if (flagsp>0)
    if (flagsavespdataall>0  || flagnoiseanalysis>0) && (flagsp || flagloadspvalues)
        fprintf('\r** Super Pixel time course all trials output **\r');
        SNCSPAllTrialsOutputH
        if flagnoiseanalysis>0 % for Noise analysis
            fprintf('\r** Super Pixel time course all trials - Noise Analysis **\r');
            SNCSPNoiseAnalysis
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T-test for Super Pixel time course, only for 'flagsp>0' and 'flagspttest>0'
if (flagsp>0)
    if flagspttest==1 && flagsp>0
        fprintf('\r** T-test for Super Pixel time course **\r');
        SNCSPTtestH
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save average values across domains, only for 'flagspmask>0' 
if (flagquantify>0)
    fprintf('\r** Save average values across domains **\r');
    SNCQuantify
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vector analysis
if flagvector==1
    fprintf('\r** Vector analysis **\r');
    SNCVectorMap
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering & Save sum maps
if flagmap>0
    
    if flaglpfilter | flaghpfilter>0
        if flagtrialfilter == 0
%             if flagmaskfilter == 0
%                 SNCFilterSMapH2; % modified by Hisashi
%             else
%                 SNCFilterSMapH2wMask; % modified by Hisashi
%             end
            SNCFilterSMapH2; % modified by Hisashi
        end
    end
    SNCSaveSMapH; % modified by Hisashi
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vector analysis
flagpixelintensityoutput=0;
if flagpixelintensityoutput==1
    fprintf('\r** Pixel Intensity Output **\r');
    SNCPixelIntensityOutputH
end


% finishing up
% fprintf(fidmasterlog, '\tend%d-%d-%d %d:%d.%d\n', fix(clock));
% fclose(fidmasterlog);
totaltime=etime(clock, timer0);
fprintf('\r%s', finalnotice);
fprintf('\n\n\nTotal time used: %5.2f minutes (%5.2f hours).\r\r', totaltime/60, totaltime/3600);

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