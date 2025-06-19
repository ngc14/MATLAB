%% Sunin parameter file
%% Last modification 5/2/06, add trial by trial filtering
%% Last modification 3/24/06, add analysis based on ivf files. 
clear all;
system='v';             % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'c:\';     % Data disk name
datafolder = 'expt\';   % Data folder name on data disk, results will be saved in 'expresult'
expname = '060323Lap\'; % Exp folder name (in both data folder and result folder)
runname = 'run08\';      % Run foler name (in both data folder and result folder)

filename={
};

% Condition name
stim(1)={'LRGH1'};  

% Subtraction maps
Smap(1, :)={'LRall', [1:16], [17:32]};


datasource=2;       % 1: from original OI block files, 2: from previously saved single-condition 'ivf' files.
fframe = [1:2];     % first frame range e.g. [], [1], [1, 2] or [1:3]
sumrange=[5:16];    % sum frame range e.g. [5 6 7 8 9 10 11 12] or [5:12]
operation=2;        % specify type of pixel value: 0: raw DC, 1: dR, 2:dR/R, 3:dR/Rb, 4:dR/R0end
clipmethod = 1;     % 0: no clipping; 1: clipping to +-SD (value);  2: clipping to +-SD (value) with a mask (default.bmp);  3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
clipvalue  = 1.5;     % How many SD to be used for clipping, usually =1 (range is plus and minus 1SD on both sides of median)	% clipvalue = [10, 10; 20, 20; 0.8, 0];   %for method 3 only (window cliping), represents [x1, y1; x2, y2; sd, nouse]
flagmasking=1;      % 0: no masking in average value calculation, '1' for map loading, '2' for coordinates loading
maskname=[          % Names of the masks (need to be same length), these names will be added '.bmp' for loading mask map or add '_x.txt'/'_y.txt' for loading coordinates.
];
domainradius=[      % Can specify different size for different domains (e.g. V2 orientation domains are larger than V1's), default 5 pixels (leave empty)
];
flagmap=1;          % 0: do not generate any maps, 1: generate only block-averaged maps, 2: generate both block-averaged maps and single-block maps
flagsp=0;           % '0': no superpixel time course output, '1': for map loading, '2' for coordinates loading (for '1' and '2' need change 'glagmasking' to same value
superpixel=[0];
FrameRate=4;        % Hz, usually 4Hz for VDAQ and 7Hz for RedShirt, not critical for analysis, just make output sp file nicer
ext = '.t0';        % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)

%---------------- More Advanced Options (default: 0 or []) -------------%
flagrandom = 0;     % 1: if is awake & randomnized data (requires '_stimseq.txt' at block file folder. Usually 0
flagperform= 0;     % 1: if is awake data (requires '_performance.txt' at block file folder), usually 0
flaggoodstim = 0; % 1: will looking for 'goodstim.txt' which contains '1's for good stim conditions, and used these for average map & quantification
saveaccum=0;       % 1: save the temporary accumulate maps, usually =0
accummapnum=[];  % Smap number to be saved accumulatly, no single condition map will be saved this way to save space.
tempsaveblocks = 0;     % how often save the temp analysis maps, eg. 5 means save results every 5 blocks. % set tempsaveblocks=999 to avoid save tempmaps
blockselect =[];  % select blocks for processing, eg. [ 1 3 4 5], leave empty for all block being selected ([]).
flagsaveivf=1;      % 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
flagsaveblk=0;     % 1: save averaged data into a block file (under testing). 0: not save
flagquantify=0;    % 1: quantify map value based on masks (fun1-fun3.txt output), 2: detailed quantify (blk by blk, slow), 0: no quantification
flagtrialfilter=0;	% 0: no trial by trial filtering; 1: filtering subtraction map for each trial before summation.
flaghpfilter=1;		% for average maps only (takes time) 0: no high-pass filter, 1: highpass filter based on data
hpfmethod='ribot'	% filter kernel type, see 'filtermethod' below in vector section, default 'fastmean'
hpkernelsize=2;	% diameter in pixels. 50-200 depends on baseline-noise size. slow for large kernel. For 'ribot' filter, it's order of polynomial (usually 2 or 3)
flaglpfilter=1;		% for average maps only 0: no low-pass filter, 1: lowpass filter based on data
lpfmethod='gaussian' % filter kernel type, see 'filtermethod' below in vector section, default 'gaussian'
lpkernelsize=5;		% diameter in pixels. 5-10 depends on high frequency noise size. 

% Alignment
flagalign=0;       % 0: no shift/align;  1: shift/align all frames with asigned frame; 2: shift/align within stim frames, 11: load shift parameters (shiftinput.txt in run folder);
shiftmethod = 41;   % 1: coorleation method, 2: min difference method, 3:masked, 4:normxcorr2 (fast)
shiftrange = [-9 9 9 9];  % how many pixels to shift to search for best fit, towards left, right, up, down
shiftstep1 = 3;    % Shift by 'step1' first, then use finer shift 'step2' to find better alignment around the result from first shift.
shiftstep2 = 0;    % ==0 if no fine shift is needed
firstframe = 'els_E00B050.blk';    % read the very first frame for shifting (some case the first block is not in file array), or you want to use a frame in the middle. only for flag 1,2

%---------------- Vector analysis --------------------------------%
flagvector=0;   % 0: no vector analysis, 1: vector map from raw data, 2: vector map from previous analyzed single condition maps (i.e. only if you have run vector analysis before, and you have single condition map in 'vector\' folder).
% specify single condition maps used for vector analysis:
vect(1,:)={'H', [1 5 9], []};     % e.g. 'Horizontal orientation' has [1 5 9]: left, right, both eye horizontal. the 2nd [] is usually empty.
vect(2,:)={'A', [2 6 10], []};
vect(3,:)={'V', [3 7 11], []};
vect(4,:)={'O', [4 8 12], []};
lowpass=15;     % low-pass filter kernel size (diameter in pixels, typical: ~7s-20 for size 300-500)
highpass=80;	% high-pass filter kenel size; (typical: ~80 for 504x504)
filtermethod='gaussian';    % choose from 'fastmean', 'slowmean', 'gaussian'
                            % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
                            % 'slowmean' 'disk-like' but better at edge
                            % 'gaussian' Gaussian filter with half sd
                            % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
                              
%---------------- Do not change following lines ------------------%
sunincore34;
sound(sin([1:100]), 100);   % give a sound notification
clear all
return;
