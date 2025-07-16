%% Sunin parameter file
%%
%%  This parameter file works with sunincore33H033.m
%%
%%  Ver 3.4.2 (150224)	sumtype=3 was introduced. just for check. 
%%						Frame-by-frame processes were improved. "flagmaskfilter" was omitted. Hisashi
%%  Ver 3.4 (150221) flagalign = 1 now perfroms image shift frame by frame. need .4 file. dalsalineremove was introduced. Hisashi
%%  Ver 3.3.1 (150220) When flagblanksubtraction==2, blank subtraction will be performed after averaging all blank condition. Hisashi
%%
clear all;
system='v';             % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'K:\';     % Data disk name
datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in outputfolder
outputdriver = 'D:\Vandy\';     % Data disk name
outputfolder = 'expt1\'; % Output folder name on data disk, results will be saved here
expname = '101231JerL2G_Atn\'; % Exp folder name (in both data folder and result folder)
runname = 'run00\';      % Run foler name (in both data folder and result folder)
resultname = 'sunin\sunin_H34\'; % Used only if flagmap=1.
SPresultname = 'suninSP\suninSP_H34\'; % Used only if flagsp=1.
vectorresultname = 'suninVector\suninVector_H34\'; % Used only if flagvector=1.

filename={  % leave blank for atuo search
};

% Condition names
stim(1)={'AtnIn_Blank'};  
stim(2)={'AtnIn_R2K_070_SacUp'};  
stim(3)={'AtnIn_G2K_070_SacUp'}; 
stim(4)={'AtnIn_R2K_160_SacUp'};  
stim(5)={'AtnIn_G2K_160_SacUp'}; 
stim(6)={'AtnIn_R2K_070_SacDw'};  
stim(7)={'AtnIn_G2K_070_SacDw'}; 
stim(8)={'AtnIn_R2K_160_SacDw'};  
stim(9)={'AtnIn_G2K_160_SacDw'}; 
stim(10)={'AtnOut_Blank'};  
stim(11)={'AtnOut_R2K_070_SacUp'};  
stim(12)={'AtnOut_G2K_070_SacUp'}; 
stim(13)={'AtnOut_R2K_160_SacUp'};  
stim(14)={'AtnOut_G2K_160_SacUp'}; 
stim(15)={'AtnOut_R2K_070_SacDw'};  
stim(16)={'AtnOut_G2K_070_SacDw'}; 
stim(17)={'AtnOut_R2K_160_SacDw'};  
stim(18)={'AtnOut_G2K_160_SacDw'}; 


% Basic parameters
blockselect =[12 16:41];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
% blockselect =[1:4];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
stimselect=[]; %select stims for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]) -- Since H029
FrameRate=4;        % Hz, usually 4Hz for VDAQ and 7Hz for RedShirt, not critical for analysis, just make output sp file nicer
fframe = [1:2];     % first frame range e.g. [], [1], [1, 2] or [1:3]; [1:3] for RedShirt
sumrange=[5:16];    % sum frame range e.g. [5 6 7 8 9 10 11 12] or [5:12]; [9:28] for RedShirt

% Goodstim
flaggoodstim = 1; % 1: will looking for 'goodstim.txt' which contains '1's for good stim conditions, and used these for average map & quantification
gsfilename = 'goodstim.txt'; % Specify goodstim file. The default is 'goodstim.txt'.

% Alignment
flagalign= 1; % default=1 or 2: no shift/align; 1: shift/align using the first frame on trial-by-trial. need 'shiftinput.2.txt'; 2: shift/align using all frames on trial-by-trial. need 'shiftinput.4.txt'

% Filtering
flagtrialfilter=1;	% 0: no trial by trial filtering (default); 1: filtering subtraction map for each trial before summation (take much time).
lpfmethod='fastmean'; % filter kernel type, see 'filtermethod' below in vector section, default 'gaussian'
lpkernelsize=0;		% diameter in pixels. 5-10 depends on high frequency noise size. 
hpfmethod='fastmean';	% filter kernel type, see 'filtermethod' below in vector section, default 'fastmean'
hpkernelsize=0;	% diameter in pixels. 50-200 depends on baseline-noise size. slow for large kernel. For 'ribot' filter, it's order of polynomial (usually 2 or 3)
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, apply median filter, then expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', filtering using FFT (Fast Fourier Transform).

dalsalineremove=0;	% It's better to set 1 when flagtrialfilter=2: 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels. 2: apply LP filter (size 3) to remove line artifact.

%---------------- More Advanced Options (default: 0 or []) -------------%
datasource=1;        % default=1; 1: from original OI block files, 2: from previously saved single-condition 'ivf' files (save time). 
                            %3: from averaged block file, which have been processed for 'goodstim', 'shiftinput', 'operation'. 'flaggoodstim', 'flagalign', and 'flagsaveblk' are ignored.
precisiontype='double'; %'double' or 'single'
flagrandom = 0;     % default=0; 1: if is awake & randomnized data (requires '_stimseq.txt' at block file folder. Usually 0
%flagperform= 0;     % 1: if is awake data (requires '_performance.txt' at block file folder), usually 0
saveaccum=0;       % default=0;  1: save the temporary accumulate maps
accummapnum=[];  % default=[];  Smap number to be saved accumulatly, no single condition map will be saved this way to save space.
tempsaveblocks = 0;     % default=0; how often save the temp analysis maps, eg. 5 means save results every 5 blocks. % set tempsaveblocks=999 to avoid save tempmaps
flagsaveivf=0;      % default=0; 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
% ivfblockselect=[]; % select blocks for making ivf file, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
flagsaveaverageblk=0;     % default=0; 1: save averaged data into a block file (under testing). 0: not save
% ablkname='Lap_E01_test3'; % This parameter is not used currently.
flagquantify=0;    % default=0; 1: quantify map value based on masks (fun1-fun3.txt output), 2: detailed quantify (blk by blk, slow), 0: no quantification

%---------------- Make single-condition and difference maps --------------------------------%
flagmap=1; %    % 0: do not generate any maps, 1: generate only block-averaged maps, 2: generate both block-averaged maps and single-block maps
outputimageformat  = 'png'; % The format for subtraction maps. % Choose from these format: 'bmp', 'gif', 'jpg', 'pcx', 'tif', etc. Since H33
ext = '';        % this name will be add to each output bmp maps (e.g. 10_OD_k01.a.bmp)
operation=2;        % specify type of pixel value: 0: raw DC, 1: dR, 2:dR/R, 3:dR/Rb, 4:dR/R0end
% Clipping
clipmethod = 1;
	% 0: no clipping; 
	% 1: clipping at median+-SD (value);  
	% 2: clipping at median+-SD (value) with a mask (masks\clipmask\default.bmp);  
	% 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
	% 4: clipping at median+-intensity change (value), Ex. 0.0005 or 0.0008;  
	% 5: clipping at 0+-SD (value);
	% 6: clipping at 0+-intensity change (value);
clipvalue  = 1.5; %0.001 %1.5 %1.5;     % How many SD to be used for clipping, usually =1 (range is plus and minus 1SD on both sides of median)	% clipvalue = [10, 10; 20, 20; 0.8, 0];   %for method 3 only (window cliping), represents [x1, y1; x2, y2; sd, nouse]

% Color-corded map, since H33
flagcolorcordedmap = 0; % If 1, make color-coded map
	ccoutputfolder='cc_converted\';
	cclutfile = 'jetH2r.lut';
	ccprefix = 'lut_';
	ccLPMethod='median';    % Low-pass filtering method, usually 'gaussian'
	ccLPKernel=0;%10; %15;             % set as 0 if no LP filter
	ccHPMethod='fastmedianH';       % High-pass filtering method, 'ribot' filter is fastest
	ccHPKernel=0; %100;             % For 'ribot' filter, it's the order of polynomial fit, usually set as 2, set as 0 if no HP filter
	ccflagoutputcolortable = 0;

% Difference maps
Smap(1, :)={'AtnIn-AtnOut_Cocktail', [2:9], [11:18]};
Smap(2, :)={'AtnOut-AtnIn_Cocktail', [11:18], [2:9]};
Smap(3, :)={'AtnIn-PsvSac_Cocktail', [2:9], [20:27]};
Smap(4, :)={'PsvSac-AtnIn_Cocktail', [20:27], [2:9]};
Smap(5, :)={'AtnOut-PsvSac_Cocktail', [11:18], [20:27]};
Smap(6, :)={'PsvSac-AtnOut_Cocktail', [20:27], [11:18]};



flagblanksubtraction=0; % default=0; %Since H32
if flagblanksubtraction
%     SmapBSub(1, :)={'AtnIn-AtnOut_Cocktail', [1], [10]};
%     SmapBSub(2, :)={'AtnOut-AtnIn_Cocktail', [10], [1]};
%     SmapBSub(3, :)={'AtnIn-PsvSac_Cocktail', [1], [19]};
%     SmapBSub(4, :)={'PsvSac-AtnIn_Cocktail', [19], [1]};
%     SmapBSub(5, :)={'AtnOut-PsvSac_Cocktail', [10], [19]};
%     SmapBSub(6, :)={'PsvSac-AtnOut_Cocktail', [19], [10]};
end

%---------------- Superpixel and Time course analysis --------------------------------%
flagsp=0;           % '0': no superpixel time course output, '1': for map loading, '2' for coordinates loading (for '1' and '2' need change 'flagspmask' to same value
flagsavespvalues=0; % default=0; 1: Values related to superpixel analysis are saved as MAT files in datafolder/expname/runname/. This value is ignored when flagloadspvalues==1
flagloadspvalues=0; % default=0; 1: Values related to superpixel analysis are loaded from MAT files in datafolder/expname/runname/. When flagloadspvalues==1, flagmap and flagvector become 0 automatically
flagspmask=flagsp;      % default=0; 0: no masking in average value calculation, '1' for map loading from 'masks', '2' for coordinates loading
% flagsptrialfilter=1;	% default=0; 0: no trial by trial filtering for SP; 1: filtering map for each trial before summation for SP (take much time).
spmaskname=[          % Names of the masks (need to be same length), these names will be added '.bmp' for loading mask map or add '_x.txt'/'_y.txt' for loading coordinates.
];
domainradiussize=5; %default 5 pixels (80 micro)
savemaskmap=1;
flagsavespdataall=0; % default=0; 1: save trial-by-trial raw data in text file; default 0; (since H027)
flagsaveframeave=0; % default=0; 1: calculate the averages of frames in fframe and sumrange; default 0; (since H027)
spoperation=3;        % default=3; specify type of superpixel value: 1: raw DC; 2:dR/R; 3: raw DC & dR/R; default=3, but that take more time if flagsptrialfilter==1. (since H028)
flagnoiseanalysis=1; % default=0 (since H028)
% SPresultname = 'suninSP\suninSP00_Md100_bmp\';
% domainBar=[      % Can specify different size for different domains (e.g.
% V2 orientation domains are larger than V1's), default 5 pixels (leave empty)
% ];
% superpixel=[0];

flagspcond=1;
% % Calculate average or differece of superpixel response of conditions, if flagspcond==1.
% % If the 2nd [] is usually empty, calculate average of 1st []. If the 2nd [] is not usually empty, calculate differece of average between 1st [] and 2nd [].
sumtype=3; % Default=3; 
% 1: Averaging or subtraction are done when at least one of conditions is success trial in a block. 
% 2: Averaging or subtraction are done only when all of conditions are success trial in a block.
% 3: Subtraction will be done when at least one of trials is success in all blocks. Not for statistics.
% SPcond(1, :)={'AtnIn-AtnOut_RGAO', [2:9], [11:18]};
% 
SPcond(1, :)={'AtnIn_RGAO', [2:9], []};
SPcond(2, :)={'AtnOut_RGAO', [11:18], []};
SPcond(3, :)={'PsvSac_RGAO', [20:27], []};
SPcond(4, :)={'AtnIn-AtnOut_RGAO', [2:9], [11:18]};
SPcond(5, :)={'AtnIn-PsvSac_RGAO', [2:9], [20:27]};
SPcond(6, :)={'AtnOut-PsvSac_RGAO', [11:18], [20:27]};


flagSPblanksubtraction=0; % default = 0 (Since H32)
% 0 ... no blank subtraction; 1 ... subtraction within each block file; 2 ... subtraction by the average of blank response across all blocks
if flagSPblanksubtraction
% SPcondBSub(1, :)={'AtnIn_RGAO', [1], []};
% SPcondBSub(2, :)={'AtnOut_RGAO', [10], []};
% SPcondBSub(3, :)={'PsvSac_RGAO', [19], []};
% SPcondBSub(4, :)={'AtnIn-AtnOut_RGAO', [1], [10]};
% SPcondBSub(5, :)={'AtnIn-PsvSac_RGAO', [1], [19]};
% SPcondBSub(6, :)={'AtnOut-PsvSac_RGAO', [10], [19]};
end

flagspttest=1;
% Calculate superpixel response across several conditions, if flagspttest==1
tail='both';            % tail = 'both' specifies the alternative A~=B (default). tail = 'right' specifies the alternative A>B. tail = 'left' specifies the alternative A<B., Note since OI signal is negative, A>B meens A has smaller response
SPttestcond(1, :)={'AtnIn-AtnOut_RGAO', [2:9], [11:18]};
SPttestcond(2, :)={'AtnIn-PsvSac_RGAO', [2:9], [20:27]};
SPttestcond(3, :)={'AtnOut-PsvSac_RGAO', [11:18], [20:27]};


flagSPTtestblanksubtraction=0; % default = 0 (Since H32)
% 0 ... no blank subtraction; 1 ... subtraction within each block file; 2 ... subtraction by the average of blank response across all blocks
if flagSPTtestblanksubtraction
% SPttestcondBSub(1, :)={'AtnIn-AtnOut_RGAO', [1], [10]};
% SPttestcondBSub(2, :)={'AtnIn-PsvSac_RGAO', [1], [19]};
% SPttestcondBSub(3, :)={'AtnOut-PsvSac_RGAO', [10], [19]};

end

%---------------- Vector analysis --------------------------------%
flagvector=0;   % 0: no vector analysis, 1: vector map from raw data, 2: vector map from previous analyzed single condition maps (i.e. only if you have run vector analysis before, and you have single condition map in 'vector\' folder).
% specify single condition maps used for vector analysis:
vect(1,:)={'R', [2 3], []};     % e.g. 'Horizontal orientation' has [1 5 9]: left, right, both eye horizontal. the 2nd [] is usually empty.
vect(2,:)={'Y', [4 5], []};
vect(3,:)={'G', [6 7], []};
vect(4,:)={'C', [8 9], []};
vect(5,:)={'B', [10 11], []};
vect(6,:)={'M', [12 13], []};
% Since H031, the parameters for filtering of vector maps follows those for difference maps.
% lowpass=15;     % low-pass filter kernel size (diameter in pixels, typical: ~7s-20 for size 300-500)
% highpass=100;	% high-pass filter kenel size; (typical: ~80 for 504x504)
% LPfiltermethod='gaussian';
% HPfiltermethod='fastmedianH';
polarclipmethod = 4;
polarclipvalue  = 0.001; % angle*map map is usually too dark, use lower clipsd to make it brighter, used 1 for 040113GarRun2
% polarclipsd=1;  % angle*map map is usually too dark, use lower clipsd to make it brighter, used 1 for 040113GarRun2
flagvectormask=0; % 0: no vector mask; 1: vector mask is applied for Bosking normalization. A mask bmp file have to be saved as '/masks/vectormask/default.bmp'.
flagsaveivfvector=0; % 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
% lutfile='rygcbm.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi                              
lutfile='bwfix.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi  
% lutfile='rygb.lut'; %The file for color table, which should be save in the same folder as this file. By Hisashi  
%---------------- Do not change following lines ------------------%
sunincore34H034;
% sound(sin([1:100]), 100);   % give a sound notification
clear all
return;
