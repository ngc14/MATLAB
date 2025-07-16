%% Alignuser.m  Optical Imaging image alignment     HDL 060417 modified from 'suninalign.m'
%% This file is not a function. 
%% This file will be called by 'sunincore', and use parameters defined in 'suninuser' 
%% This file also can be used alone (set manual==1 and provide parameters)
%% This file calls 'OIAlignCore.m' which aligns two frames.
%% input: 
%%      Optical Imaging Data location
%%      Align type (see below)
%%      Obj function (see below)
%%      Base frame  (for certain alignment need a common base frame)
%%      Preprocess (see below)
%% output: 
%%      Shiftlog (dx, dy to shift frame back)
%%      goodstim.txt: evaluation of  alignment results (0 for good conditions), note this depends on 'Aligntype'.

%% aligntype:
%%      1:  General alignment, align all frames to one selected frame (xxx), output 'shiftlog.1.txt'
%%      2:  First frame alignment, align all first frames (can be a range) in each stim to a selected frame 
%%          (usually a good frame in the middle of the imaging session,e.g. first frame in block 100), result: 'shiftlog.2.txt', 'goodstim-ff.txt'
%%      3:  Within condition alignment, align all subsquent frame to the first frame selected in each stim presentation. result: 'shiftlog.3.txt', 'goodstim-within.txt'
%%      4:  combination of 2 and 3
%% objfun: Objective function
%%      1: fast correlation, use 'normxcorr2.m', precision=1 pixel
%%      2: slow correlation, a combination of fast correlation (1) and sub-pixel correlation, default precision=0.1 pixel (change 'precision) to modify this default.  
%%      3: use stdev of difference image (default precision=1 pixel)
%%      4: min-difference score (default precision=1 pixel)
%%      5: use mutual info to align (note: only uint8 image can be use here, precision=1 pixel)
%% Pre-Process:
%%      *  illumination bkground subtraction
%%      *  low-pass filtering
%%      *  cropping
%%      1: masking, need implement
%%      2: contrast enhancement, need implement
%%      3: dark sub, need implement

function Align0(datadrive)
manual=1;   % set manual=0 when it is called by 'suninuser.m'
if manual==1;
    clear all;              % do not use this if this function is not used alone
    system='v';             % 'v' for VDAQ, 'r' for RedShirt
    datadrive = 'j:\';     % Data disk name
    datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in 'expresult'
    expname = '070605LapV_RGLum\'; % Exp folder name (in both data folder and result folder)
    runname = 'run00\';      % Run foler name (in both data folder and result folder)
    blkfilename={
    };
    aligntype=4;    % 1: general, 2: fframe, 3:within stim, 4: combine of 2&3
    objfun = 2;     % 1: fast correlation, 2: slow correlation, 3: stdev of difference, 4: min-difference, 5: mutual info
    shiftrange = [-20 20 -20 20];  % how many pixels to shift to search for best fit, towards left, right, up, down
    baseframefile={'Lap_E01B000.BLK', 1, 1};  % filename, stim number, frame number for 'aligntype=2 or 4'
    ffrange=[1];      % the range of first frame to be aligned in 'aligntyp=2 or 4'

    LPKernel=2;       % low-pass kernel, part of pre-process befor alignment, put 0 if no need 
    crop1=[153 222 200 306];  % [x1, y1, x2, y2] of clip region, part of pre-process, put 0 if no need
    flagillubksub=0;  % illum bkground sub, calculated once per block with 50kernel, part of pre-porcess, 0 for none
%     percentage2=0.9;	% percent of trials to be included after shift2 process
%     percentage3=0.9;	% percent of trials to be included after shift3 process

%     corrthreshold1=0.95; %Threshold for the correlation coefficient in shift2 process
    corrthreshold1=0.8; %Threshold for the correlation coefficient in shift2 process
    shiftthreshold1=10; %Threshold for shift in pixel: if mean pixel shift for framerange exceeds this value, the trial is marked as '0'
    
    corrthreshold2=0.8; %Threshold for the correlation coefficient in shift3 process
    shiftthreshold2=3; %Threshold for shift in pixel: if mean pixel shift for framerange exceeds this value, the trial is marked as '0'
end
% end of input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if aligntype~=2 
    if size(ffrange, 2)~=1
        fprintf('Error: Check ffrange, only aligntype==2 can have more than one first frame selected!\n');
    else
        fframe=ffrange(1,1);    % used in aligntype=3
    end    
end
    
blockfolder = strcat(datadrive, datafolder, expname, runname); 
ctime=fix(clock);
% A shift log with date and time as file name will always saved in data folder
if aligntype==1 | aligntype==2 | aligntype==3
%     shiftlog = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.', num2str(aligntype), '.txt');
    shiftlog = strcat(blockfolder, '_shiftlog','.', num2str(aligntype), '.txt');
    fidshiftlog = fopen(shiftlog, 'w');  %for shiftlog output
end
if aligntype==4
%     shiftlog1 = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.2.txt');
    shiftlog1 = strcat('_shiftlog.2.txt');
    fidshiftlog1 = fopen(shiftlog1, 'w');  %for shiftlog output
%     shiftlog2 = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.3.txt');
    shiftlog2 = strcat('_shiftlog.3.txt');
    fidshiftlog2 = fopen(shiftlog2, 'w');  %for shiftlog output
end
fidshiftinput=fopen(strcat('shiftinput.2.txt'), 'w');

if isempty(blkfilename)
    if system=='v'
        tempfilename=struct2cell(dir([blockfolder, '*.blk']));
    elseif system=='r'
        fprintf('Note: you may need delete non-block "*.da" files in data folder\n');
        tempfilename=struct2cell(dir([blockfolder, '*.da']));
    end
    blkfilename=sort(tempfilename(1,:)');
    for i=1:size(blkfilename,1)
        fprintf('''%s''\n', getfield(cell2struct(blkfilename(i), 'junk'), 'junk'));
    end
    fprintf('\nfound %d blk files(sorted, check sequence).\n', size(blkfilename,1));
end
blocknum=size(blkfilename, 1);      % how many blocks

% Read head info
anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(blkfilename(1), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;

baseframe=OIReadFrame(strcat(blockfolder, getfield(cell2struct(baseframefile(1), 'junk'), 'junk')), system, getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'));

% Preprocess
if flagillubksub
    baseillubk=OIMeanFilt(baseframe, 50);
    imbase=OIAPreProcess(baseframe, baseillubk, LPKernel, [0 0], crop1);
else
    imbase=OIAPreProcess(baseframe, 0, LPKernel, [0 0], crop1);    
end

if system=='v'
    headerlength = 1716;    % byte
elseif system=='r'
    headerlength = 5120;
end
timertrigger=1;
ctime=clock;
for k=1:blocknum  % Start read/process block by block
    fprintf('\rblock=%d  ',k);
    filename=getfield(cell2struct(blkfilename(k), 'junk'), 'junk');
    if flagillubksub    % calculate illumination bk
        illubk=OIReadFrame(strcat(blockfolder,filename), system, 3,1);  % note: 3 is just for selecting a middle frame to calculate illumination background
        illubk=OIMeanFilt(illubk, 50);              
    else
        illubk=0;
    end    
    for i=1:NStim
        switch aligntype
        case 1  % General alignment
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            for j = 1:FramesPerStim            
                if timertrigger==1|timertrigger==2  % estimate time
                    [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*FramesPerStim);
                end
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(baseframe, imshift, shiftrange, objfun); % find the best alignment
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
            end
        case 2        % First frame alignment
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            if timertrigger==1|timertrigger==2  % estimate time
                [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*size(ffrange, 2));
            end
            for j=ffrange
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1);
                [imshifted, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, shiftrange, objfun);
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
%                imwrite(nc(imshifted-imbase), strcat('tempd\', num2str(k), '-', num2str(i), '.bmp')); 
            end
        case 3        % Within condition alignment
            if timertrigger==1|timertrigger==2  % estimate time
                 [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim);
            end            
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            imfframe=OIAPreProcess(tDCFrames(:,:,fframe), illubk, LPKernel, [0 0], crop1);        % note: a bug here, fframe maybe more than one
            OIAWriteShiftLog(fidshiftlog, filename, k, i, 1, 0, 0, 1, filename, i, 1, 1);          % just for consistency
            for j = 2:FramesPerStim
                if timertrigger==1|timertrigger==2  % estimate time
                     [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*(FramesPerStim-1));
                end            
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(imfframe, imshift, shiftrange, objfun); 
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, filename, i, 1, currentobj);
            end
        case 4        % combine 2 and 3
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            if timertrigger==1|timertrigger==2  % estimate time
                 [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim);
            end            
            for j=ffrange
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1);
                [imshifted, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, shiftrange, objfun);
                OIAWriteShiftLog(fidshiftlog1, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
                fprintf(fidshiftinput, '%s \t%d \t%d \t%d \t%3.1f \t%3.1f \t%6.6f\r\n', filename, k, i, j, pixoff(1), pixoff(2), bestobj);
            end
            imfframe=OIAPreProcess(tDCFrames(:,:,fframe), illubk, LPKernel, [0 0], crop1);        % note: a bug here, fframe maybe more than one
            OIAWriteShiftLog(fidshiftlog2, filename, k, i, 1, 0, 0, 1, filename, i, 1, 1);          % just for consistency
            for j = 2:FramesPerStim
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(imfframe, imshift, shiftrange, 1);        % Note: here forced 'objfun' to 1
                OIAWriteShiftLog(fidshiftlog2, filename, k, i, j, pixoff(1), pixoff(2), bestobj, filename, i, 1, currentobj);
            end
        otherwise
            fprintf('error: wrong alignment method\n');
        end        
        % close & reopen files
		switch aligntype
		case 1
            fclose(fidshiftlog);
		case 2
            fclose(fidshiftlog);
		case 3
            fclose(fidshiftlog);
		case 4
            fclose(fidshiftlog1);
            fclose(fidshiftlog2);
            fclose(fidshiftinput);
		end
		switch aligntype
		case 1
            fidshiftlog = fopen(shiftlog, 'a');
		case 2
            fidshiftlog = fopen(shiftlog, 'a'); 
		case 3
            fidshiftlog = fopen(shiftlog, 'a'); 
		case 4
            fidshiftlog1 = fopen(shiftlog1, 'a'); 
            fidshiftlog2 = fopen(shiftlog2, 'a'); 
            fidshiftinput= fopen('shiftinput.2.txt', 'a');
		end
    end % for i=1:NStim
end % for k1:blocknum

switch aligntype
case 1
    fclose(fidshiftlog);
%    OIAshift2goodstim1(shiftlog);  % not implemented yet
case 2
    fclose(fidshiftlog);
%     OIAshift2goodstim2(shiftlog1, 0.0001, NStim, percentage2);
    OIAshift2goodstim2H(shiftlog1, corrthreshold1, shiftthreshold1, NStim);
case 3
    fclose(fidshiftlog);
%     OIAshift2goodstim3(shiftlog2, percentage3, 1, NStim, FramesPerStim, [5:16]);
    OIAshift2goodstim3H(shiftlog2, corrthreshold2, shiftthreshold2, NStim, FramesPerStim, [5:16]);
case 4
    fclose(fidshiftlog1);
% 	OIAshift2goodstim2(shiftfile2, threshold, NStim, percentage)
    OIAshift2goodstim2H(shiftlog1, corrthreshold1, shiftthreshold1, NStim);
    fclose(fidshiftlog2);
%   OIAshift2goodstim3(shiftfile3, percentage, shiftthreshold, NStim, NFrame, framerange)    
%     OIAshift2goodstim3(shiftlog2, percentage3, 1, NStim, FramesPerStim, [5:16]); %commented out by Hisashi
    OIAshift2goodstim3H(shiftlog2, corrthreshold2, shiftthreshold2, NStim, FramesPerStim, [5:16]);
    fclose(fidshiftinput);
otherwise
end


return