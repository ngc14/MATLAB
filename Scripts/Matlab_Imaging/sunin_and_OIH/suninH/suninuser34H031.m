%% Sunin parameter file
%%
%%  (080301) Hisashi: This parameter file works with sunincore34H026.m or
%%  newer version.
%%
%% Last modification 5/2/06, add trial by trial filtering
%% Last modification 3/24/06, add analysis based on ivf files. 
clear all;
system='v';             % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'G:\';     % Data disk name
datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in outputfolder
outputdriver = 'D:\My Documents\Magic Briefcase\';     % Data disk name
outputfolder = 'expt1\'; % Output folder name on data disk, results will be saved here
expname = '100208JerL2_Atn_S\'; % Exp folder name (in both data folder and result folder)
runname = 'run99\';      % Run foler name (in both data folder and result folder)
resultname = 'sunin\sunin_noSF_0C00020_F0306\';
SPresultname = 'suninSP\suninSP00_Md100_bmp\';
vectorresultname = 'suninVector\suninVector_Mn15Md100_mC00002\';

filename={  % leave blank for atuo search
};

% Condition names
stim(1)={'Atn_Blank'};  
stim(2)={'Atn_M2K_025'};  
stim(3)={'Atn_M2K_025'}; 
stim(4)={'Atn_M2K_025'};  
stim(5)={'Atn_M2K_025'}; 
stim(6)={'Atn_M2K_025'};  
stim(7)={'Atn_M2K_025'}; 
stim(8)={'Atn_M2K_025'};  
stim(9)={'Atn_M2K_025'}; 
stim(10)={'Psv_Blank'};  
stim(11)={'Psv_M2K_025'};  
stim(12)={'Psv_M2K_025'}; 
stim(13)={'Psv_M2K_025'};  
stim(14)={'Psv_M2K_025'}; 
stim(15)={'Psv_M2K_025'};  % no data
stim(16)={'Psv_M2K_025'}; 
stim(17)={'Psv_M2K_025'};  
stim(18)={'Psv_M2K_025'}; 


% Subtraction maps
Smap(1, :)={'All_M2K_025-Blank', [2:9 11:14 16:18], [1 10]};
Smap(2, :)={'Atn_M2K_025-Blank', [2:9], [1]};
Smap(3, :)={'Psv_M2K_025-Blank', [11:14 16:18], [10]};
Smap(4, :)={'All_M2K_025-0', [2:9 11:14 16:18], []};
Smap(5, :)={'Atn_M2K_025-0', [2:9], []};
Smap(6, :)={'Psv_M2K_025-0', [11:14 16:18], []};
Smap(7, :)={'Atn-Psv_M2K_025', [2:9], [11:14 16:18]};
Smap(8, :)={'Psv-Atn_M2K_025', [11:14 16:18], [2:9]};

% Basic parameters
flagmap=1;          % 0: do not generate any maps, 1: generate only block-averaged maps, 2: generate both block-averaged maps and single-block maps
blockselect =[];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
stimselect=[]; %select stims for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]) -- Since H029
FrameRate=2;        % Hz, usually 4Hz for VDAQ and 7Hz for RedShirt, not critical for analysis, just make output sp file nicer
fframe = [1];     % first frame range e.g. [], [1], [1, 2] or [1:3]; [1:3] for RedShirt
sumrange=[3:6];    % sum frame range e.g. [5 6 7 8 9 10 11 12] or [5:12]; [9:28] for RedShirt
operation=2;        % specify type of pixel value: 0: raw DC, 1: dR, 2:dR/R, 3:dR/Rb, 4:dR/R0end
ext = '_00';        % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)

% Filtering
flagmaskfilter=0; % 0: no mask filtering; 1: mask filtering is applied to reduce vessel artifacts. A mask bmp file have to be saved as '/masks/filtermask/default.bmp'. By Hisashi on 080118
%filtermaskname='filtermask01.bmp'; %This mask bmp file have to be saved in '/masks/filtermask/' for mask filtering.
flagtrialfilter=1;	% 0: no trial by trial filtering (default); 1: filtering subtraction map for each trial before summation (take much time).
flaglpfilter=0;		% for average maps only 0: no low-pass filter, 1: lowpass filter based on data
lpfmethod='fastmedianH'; % filter kernel type, see 'filtermethod' below in vector section, default 'gaussian'
lpkernelsize=15;		% diameter in pixels. 5-10 depends on high frequency noise size. 
flaghpfilter=0;		% for average maps only (takes time) 0: no high-pass filter, 1: highpass filter based on data
hpfmethod='fastmedianH';	% filter kernel type, see 'filtermethod' below in vector section, default 'fastmean'
hpkernelsize=100;	% diameter in pixels. 50-200 depends on baseline-noise size. slow for large kernel. For 'ribot' filter, it's order of polynomial (usually 2 or 3)
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, apply median filter, then expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', filtering using FFT (Fast Fourier Transform).

% Clipping
clipmethod = 6;
	% 0: no clipping; 
	% 1: clipping at median+-SD (value);  
	% 2: clipping at median+-SD (value) with a mask (masks\clipmask\default.bmp);  
	% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
	% 4: clipping at median+-intensity change (value), Ex. 0.0005 or 0.0008;  
	% 5: clipping at 0+-SD (value);
	% 6: clipping at 0+-intensity change (value);
clipvalue  = 0.0020; %0.001 %1.5 %1.5;     % How many SD to be used for clipping, usually =1 (range is plus and minus 1SD on both sides of median)	% clipvalue = [10, 10; 20, 20; 0.8, 0];   %for method 3 only (window cliping), represents [x1, y1; x2, y2; sd, nouse]

%---------------- More Advanced Options (default: 0 or []) -------------%
datasource=1;        %1: from original OI block files, 2: from previously saved single-condition 'ivf' files (save time). 
                            %3: from averaged block file, which have been processed for 'goodstim', 'shiftinput', 'operation'. 'flaggoodstim', 'flagalign', and 'flagsaveblk' are ignored.
precisiontype='double';
flagrandom = 0;     % 1: if is awake & randomnized data (requires '_stimseq.txt' at block file folder. Usually 0
%flagperform= 0;     % 1: if is awake data (requires '_performance.txt' at block file folder), usually 0
flaggoodstim = 1; % 1: will looking for 'goodstim.txt' which contains '1's for good stim conditions, and used these for average map & quantification
gsfilename = 'goodstim99.txt'; % Specify goodstim file. The default is 'goodstim.txt'.
saveaccum=0;       % 1: save the temporary accumulate maps, usually =0
accummapnum=[];  % Smap number to be saved accumulatly, no single condition map will be saved this way to save space.
tempsaveblocks = 0;     % how often save the temp analysis maps, eg. 5 means save results every 5 blocks. % set tempsaveblocks=999 to avoid save tempmaps
flagsaveivf=0;      % 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
% ivfblockselect=[]; % select blocks for making ivf file, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
flagsaveaverageblk=0;     % 1: save averaged data into a block file (under testing). 0: not save
% ablkname='Lap_E01_test3'; % This parameter is not used currently.
flagquantify=0;    % 1: quantify map value based on masks (fun1-fun3.txt output), 2: detailed quantify (blk by blk, slow), 0: no quantification

% Alignment
flagalign= 0; % default 41;       % 0: no shift/align;  1: shift/align all frames with asigned frame; 2: shift/align within stim frames, 11: load shift parameters (shiftinput.txt in run folder);
% shiftmethod = 41;   % 1: coorleation method, 2: min difference method, 3:masked, 4:normxcorr2 (fast)
% shiftrange = [-9 9 9 9];  % how many pixels to shift to search for best fit, towards left, right, up, down
% shiftstep1 = 3;    % Shift by 'step1' first, then use finer shift 'step2' to find better alignment around the result from first shift.
% shiftstep2 = 0;    % ==0 if no fine shift is needed
% firstframe = 'Lap_E01B063.BLK';    % read the very first frame for shifting (some case the first block is not in file array), or you want to use a frame in the middle. only for flag 1,2

%---------------- Superpixel and Time course analysis --------------------------------%
flagsp=0;           % '0': no superpixel time course output, '1': for map loading, '2' for coordinates loading (for '1' and '2' need change 'flagspmask' to same value
flagsavespvalues=0; % 1: Values related to superpixel analysis are saved as MAT files in datafolder/expname/runname/. This value is ignored when flagloadspvalues==1
flagloadspvalues=0; % 1: Values related to superpixel analysis are loaded from MAT files in datafolder/expname/runname/. When flagloadspvalues==1, flagmap and flagvector become 0 automatically
flagspmask=flagsp;      % 0: no masking in average value calculation, '1' for map loading from 'masks', '2' for coordinates loading
% flagsptrialfilter=1;	% 0: no trial by trial filtering for SP; 1: filtering map for each trial before summation for SP (take much time).
spmaskname=[          % Names of the masks (need to be same length), these names will be added '.bmp' for loading mask map or add '_x.txt'/'_y.txt' for loading coordinates.
];
domainradiussize=5; %default 5 pixels (80 micro)
savemaskmap=0;
flagsavespdataall=0; % 1: save trial-by-trial raw data in text file; default 0; (since H027)
flagsaveframeave=1; % 1: calculate the averages of frames in fframe and sumrange; default 0; (since H027)
spoperation=2;        % specify type of superpixel value: 1: raw DC; 2:dR/R; 3: raw DC & dR/R; default=3, but that take more time if flagsptrialfilter==1. (since H028)
flagnoiseanalysis=0; %1: ; default=0 (since H028)
% SPresultname = 'suninSP\suninSP00_Md100_bmp\';
% domainBar=[      % Can specify different size for different domains (e.g.
% V2 orientation domains are larger than V1's), default 5 pixels (leave empty)
% ];
% superpixel=[0];

% % Calculate average or differece of superpixel response of conditions, if flagspcond==1.
% % If the 2nd [] is usually empty, calculate average of 1st []. If the 2nd [] is not usually empty, calculate differece of average between 1st [] and 2nd [].
flagspcond=1;
sumtype=1; % Default=1; 1: averaging or subtraction are done when at least one of conditions is success trial. 2: averaging or subtraction are done only when all of conditions are success trial.
SPcond(1, :)={'RG_A_Patch', [3], []};
SPcond(2, :)={'RG_O_Patch', [5], []};
SPcond(3, :)={'Lum_A_Patch', [7], []};
SPcond(4, :)={'Lum_O_Patch', [9], []};
SPcond(5, :)={'RG_A_Patch-B', [3], [1]};
SPcond(6, :)={'RG_O_Patch-B', [5], [1]};
SPcond(7, :)={'Lum_A_Patch-B', [7], [1]};
SPcond(8, :)={'Lum_O_Patch-B', [9], [1]};
SPcond(9, :)={'RG_A&O_Patch', [3 5], []};
SPcond(10, :)={'Lum_A&O_Patch', [7 9], []};
SPcond(11, :)={'A_RG&Lum_Patch', [3 7], []};
SPcond(12, :)={'O_RG&Lum_Patch', [5 9], []};
SPcond(13, :)={'RG_A&O_Patch-B', [3 5], [1]};
SPcond(14, :)={'Lum_A&O_Patch-B', [7 9], [1]};
SPcond(15, :)={'A_RG&Lum_Patch-B', [3 7], [1]};
SPcond(16, :)={'O_RG&Lum_Patch-B', [5 9], [1]};


% Calculate superpixel response across several conditions, if flagspttest==1
flagspttest=1;
tail='both';            % tail = 'both' specifies the alternative A~=B (default). tail = 'right' specifies the alternative A>B. tail = 'left' specifies the alternative A<B., Note since OI signal is negative, A>B meens A has smaller response
SPttestcond(1, :)={'RGvsLum_A_Patch', [3],  [7]};
SPttestcond(2, :)={'RGvsLum_O_Patch', [5],  [9]};
SPttestcond(3, :)={'RGvsLum_AO_Patch', [3 5],  [7 9]};
SPttestcond(4, :)={'AvsO_RG_Patch', [3], [5]};
SPttestcond(5, :)={'AvsO_Lum_Patch', [7], [9]};
SPttestcond(6, :)={'AvsO_RGLum_Patch', [3 7], [5 9]};
%---------------- Vector analysis --------------------------------%
flagvector=0;   % 0: no vector analysis, 1: vector map from raw data, 2: vector map from previous analyzed single condition maps (i.e. only if you have run vector analysis before, and you have single condition map in 'vector\' folder).
% specify single condition maps used for vector analysis:
vect(1,:)={'000', [3 7 13 17], []};     % e.g. 'Horizontal orientation' has [1 5 9]: left, right, both eye horizontal. the 2nd [] is usually empty.
vect(2,:)={'045', [4 8 14 18], []};
vect(3,:)={'090', [5 9 15 19], []};
vect(4,:)={'135', [6 10 16 20], []};

flagabsoluteforvector=0; %default=0; if 1, calculate absolute value of subtractions shown above. By Hisashi on 100709
% Since H031, the parameters for filtering of vector maps follows those for difference maps.
% lowpass=15;     % low-pass filter kernel size (diameter in pixels, typical: ~7s-20 for size 300-500)
% highpass=100;	% high-pass filter kenel size; (typical: ~80 for 504x504)
% LPfiltermethod='gaussian';
% HPfiltermethod='fastmedianH';
polarclipmethod = 1; %default 1 for polar map
polarclipvalue  = 1.0; % angle*map map is usually too dark, use lower polarclipvalue for magnitude map to make it brighter, used 1 for 040113GarRun2
% polarclipsd=1;  % angle*map map is usually too dark, use lower clipsd to make it brighter, used 1 for 040113GarRun2
flagvectormaskfilter=flagmaskfilter; % 0: no mask filtering; 1: mask filtering is applied to reduce vessel artifacts. A mask bmp file have to be saved as '/masks/filtermask/default.bmp'. By Hisashi on 080118
flagsaveivfvector=0; % 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
% lutfile='rygcbm.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi                              
lutfile='bwfix.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi  
% lutfile='rygb.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi  
% vectorresultname = 'suninVector\suninVector_Mn15Md100_mC00002\';
%---------------- Do not change following lines ------------------%
sunincore34H031;
% sound(sin([1:100]), 100);   % give a sound notification
clear all
return;
