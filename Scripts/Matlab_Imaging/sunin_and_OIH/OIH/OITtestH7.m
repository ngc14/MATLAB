function OITtestH7()
% OITtestH
% T-Test for Optical Imaging Image analysis, written by Hisashi
% Ver 7.1.5: (150224) Frame-by-frame processes were improved.
% Ver 7.0: (150221) dalsalineremove=1 were introduced. flagalign=2 were introduced. Hisashi
% Ver 6.1: (150219) If flagblanksubtraction==2, blank subtraction will be performed after averaging all blank condition.
% Ver 6.0: (100611) Non-parametric statistics were added. 
% Ver 5.0: Blank subtraction before the analysis was introduced. 
% Ver 4.1.1: a small fix 
% Ver 4.1: Band-pass filtering using FFT was added. A bug in Mask filtering was fixed. 
% Ver 4.0: Data from Each trial is treated as one data for statistics. In
%   the previous versions, averaged data across several trials within a block
%   was treated as one data for statistics. This revision requires more
%   process time, but improves statistical power.
% Ver 3.1: 'outputdriver' and 'gsfilename' were introduced
% Ver 3.0: 
% Ver 2.2: use Single-Precision to reduce memory usage. 080124 Hisashi
% Ver 2.1: Analysis time and memory usage was improved. Please fill 'stimeselect' to specify
%               stim ID to read. -080123
% Ver 2.0: Comparison can be done between multiple pairs. 071025 Hisashi
%             Goodstim can be analyzed.
% Ver 1.0: Shift align can be analyzed.
% modified from blockview -HDL 041001-050729

clear all;
% __________________ Start User Input ____________________
system='v';     % 'v' or 'r'
datadriver = 'C:\Users\OMAR\Documents\MATLAB\';     % Data disk name    /Users/omargharbawie/Documents/Exp/110725BooneMotor/run05
datafolder = 'Exp\';   % Data folder name on data disk, results will be saved in outputfolder


outputdriver = 'C:\Users\OMAR\Documents\MATLAB\ExpResults\';     % Data disk name   
outputfolder = 'Results\'; % Output folder name on data disk, results will be saved here
expname = 'MsHowell\Right_Chamber\Awake_Grasping\2014_09_12_MsHowell\'; % Exp folder name (in both data folder and result folder)
runname = 'run0\';      % Run foler name (in both data folder and result folder)
resultname = 'OITtestH7\NewAlign_BLK33\Frames20-30\LP2_HP300\FFrame1\Masked\';

lutfile = 'statcolorBR.lut'




filename = {            % blk file names, leave empty for automatic searching
};
TestType=2;             % 1: Within condition comparison (fn vs. f0). In thmap, the second stim group will be ignored. 2: Compare 2 conditions (stim group 1 vs. stim group 2 in thmap) after first frame subtraction
tail='both';            % tail = 'both' specifies the alternative A~=B (default). tail = 'right' specifies the alternative A>B. tail = 'left' specifies the alternative A<B., Note since OI signal is negative, A>B meens A has smaller response
parametric=1;           % 1: parametric 2: non-parametric
prefix1='';


Smap(1, :)={'Button-Blank_', [2], [1]};
Smap(2, :)={'LargeSphere-Blank_', [3], [1]};
Smap(3, :)={'SmallSphere-Blank_', [4], [1]};
Smap(4, :)={'Large_Sphere-Button', [3], [2]};
Smap(5, :)={'Small_Sphere-Button_', [4], [2]};
Smap(6, :)={'Small_Sphere-Large_Sphere_', [4], [3]};


flagblanksubtraction=0; % default=0; 0 ... no blank subtraction; 1 ... subtraction within each block file; 2 ... subtraction by the average of blank response across all blocks
if flagblanksubtraction
% Bmap(1, :)={'R2K-G2K_All', [], []};
% Bmap(2, :)={'070-160_All', [], []};
% Bmap(3, :)={'AtnIn-AtnOut_All', [1], [10]};

end


fframe=[1];
framerange=[20:30];
blockselect =[];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).

threshold=[0.0001 0.001 0.01 0.05];	% threshold p values, just for generating thresholded p maps, each value will have a thresholded map (i.e. 0.05 will generate a p<0.05 map)
pmax=0.0001;
pmin=0.0000001;

% clipsd=1;				% for map output, usually 1
LPMethod='gaussian';    % Low-pass filtering method, usually 'gaussian'
LPKernel=2;             % set as 0 if no LP filter
HPMethod='fastmean';       % High-pass filtering method, 'ribot' filter is fastest
HPKernel=300;             % Set as 0 if no HP filter. For 'ribot' filter, it's the order of polynomial fit, usually set as 2. For 'fft' filter, it's the period of cycle for cutoff (in pixel). Setting as the half of the image size means that the cutoff is 2.
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, apply median filter, then expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', band-pass filtering using FFT (Fast Fourier Transform).
flagfiltermask=1; % default=0; 0: no mask filtering; 1: mask filtering is applied to reduce vessel artifacts. A mask bmp file have to be saved as '/masks/filtermask/default.bmp'. By Hisashi on 080118
Binning=1;             % note: it is the x times of original size (i.e. 0.5 is equal to 2x2 binning)
flaggoodstim=0;     % similar to the 'goodstim' in sunin, if 'goodstim=1', will search for 'goodstim.txt' in block folder
gsfilename='goodstim.txt';
flagalign=2;        % 0: no shift/align; 1: shift/align using the first frame on trial-by-trial. need 'shiftinput.2.txt'; 2: shift/align using all frames on trial-by-trial. need 'shiftinput.4.txt'
flagrandom=0;       % default=0; similar to the 'flagrandom' in sunin, if 'flagrandom=1', will search for 'stimseq.txt' in block folder
flagsaveivf=0;      % default=0; 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
dalsalineremove=0;	% It's better to set 1 when flagtrialfilter=2: 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels

% ___________________ end of user input _______________
resultdriver=outputdriver; % updated at Ver 3.1
if isempty(resultname)
    resultname=strcat('HP',num2str(HPKernel), '_LP', num2str(LPKernel), '_Bin', num2str(Binning), '\');
elseif 	resultname(end)~='\';
	resultname=[resultname, '\'];	
end
resultfolder=[resultdriver, outputfolder, expname, runname, resultname];
runfolder=[resultdriver, outputfolder, expname, runname];
blockfolder=[datadriver, datafolder, expname, runname]; % updated at Ver 3.1
% filtermask = strcat(resultdriver, 'expt1\', expname, 'masks\filtermask\default.bmp');

if 	blockfolder(end)~='\';
	blockfolder=[blockfolder, '\'];	
end

if ~isdir(resultfolder)
    mkdir(resultfolder);    
end

if isempty(filename)
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

blockfilenum=size(filename, 1);
if isempty(blockselect)
    blockselect=1:blockfilenum;
end
blockselectnum=size(blockselect,2);
% stimselectnum=size(stimselect,2);

fprintf('total block number = %d\r\n', blockselectnum);

if flagfiltermask==1
	filtermaskname = strcat(resultdriver, outputfolder, expname, runname,'masks\filtermask\default.bmp');
    filtermasktemp = imread(filtermaskname, 'bmp');
    if size(filtermasktemp,3)==1
        filtermask = double(filtermasktemp)/255;
    else
        filtermask = double(filtermasktemp(:,:,1))/255;
    end
    clear filtermasktemp;
    if ~((strcmp(LPMethod, 'fastmean') || strcmp(LPMethod, 'gaussian')) && (strcmp(HPMethod, 'fastmean') || strcmp(HPMethod, 'gaussian')))
        fprintf('Mask filtering is available only for "fastmean" and "gaussian" now. The other methods work without mask filtering.\r');
    end
else
    filtermask = ones(504,504,'double');
end


anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(1), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;
DataType=anapar.DataType;

if blockselectnum==0;
    blockselectnum=NStim;
end
% if isempty(stimselect)
%     stimselect=1:NStim;
% else
% 	if max(stimselect)>NStim
%         fprintf('\rstim ID you specified exceed actual stim ID\r');
%     end
% end

if flaggoodstim
    goodstim=textread(strcat(runfolder, gsfilename), '%d');
    if isempty(goodstim)|size(goodstim,1)~=blockfilenum*NStim
        fprintf('Error, "goodstim.txt" does not contain right number of conditions\r');
    end
	if sum(sum(goodstim))~=sum(sum(goodstim.*goodstim))	% only 0*0=0 and 1*1=1
		fprintf('Error: ''goodstim.txt'' should contain only "1" and "0"'\n);
	end
    goodstim=reshape(goodstim, [NStim, blockfilenum]);
    goodstim=goodstim';
    if prod(sum(goodstim, 2))==0
        fprintf('Worning: too few good stims, one or more condition has no data!\n');
    end
    fprintf('\rGood stim parameters loaded\r');
else
	goodstim=ones(blockfilenum, NStim);	% if no goodstim.txt provided, use all conditions.    
end

if flagalign == 1
    [sfname, sfblock, sfstim, sfframe, sfx1, sfy1, sfcoor, junk1, junk2, junk3, junk4]=textread(strcat(runfolder, 'shiftinput.2.txt'), '%s %d %d %d %f %f %f %s %d %d %f');
    if size(sfname,1)~=blockfilenum*NStim
        fprintf('Error: flagalign=1, .2 file doesnot contain blockfilenum*NStim number of entries!\r');
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blockfilenum
        for i=1:NStim
            if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim+i).sfname
                fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstruct(k*NStim*NFrames).sfname);
            end
        end
    end
elseif flagalign == 2
    [sfname, sfblock, sfstim, sfframe, sfx3, sfy3, sfcoor]=textread(strcat(runfolder, 'shiftinput.4.txt'), '%s %d %d %d %f %f %f'); %Hisashi 071023
    fprintf('\r** Shift input parameters loaded **\r');
    if size(sfname,1)~=blockfilenum*NStim*FramesPerStim
        fprintf('Error: flagalign===2, .4 file doesnot contain blockfilenum*NStim number of entries!\n');
        size(sfname,1)
        blockfilenum
        NStim
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blockselectnum
        for i=1:NStim
        	for fr=1:FramesPerStim
                if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramesPerStim-1)+fr).sfname
                    fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname((k-1)*NStim*(FramesPerStim-1)+(i-1)*(FramesPerStim-1)+fr).sfname);
                end
            end
        end
    end
end

if (flagrandom==1);
    file=[runfolder, '_stimseq.txt'];
    stimseq=textread([blockfolder, '_stimseq.txt'], '%d');
	if size(stimseq,1)~=blockfilenum*NStim
        fprintf('Error: number of stim in ''_stimseq.txt'' is %d\r', size(stimseq));
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
else 
    stimseq=repmat([1:NStim], blockfilenum, 1);    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading frames %%%%%%%%%%%%%%%%%%%%%%%
newsize=round(FrameWidth*Binning);  % assume width=height here
numSmap=size(Smap,1);   
% Sample1=zeros(newsize, newsize, blockselectnum, numSmap, 'single');  
% Sample2=zeros(newsize, newsize, blockselectnum, numSmap, 'single');
Range1=zeros(newsize, newsize, NStim, 'single');  
Range2=zeros(newsize, newsize, NStim, 'single');  
Rangetemp1=zeros(newsize, newsize, 'single');
Rangetemp1=zeros(newsize, newsize, 'single');

if Binning ==1
    if flagfiltermask~=0
        filtermask2=filtermask;
    end
else
    if LPKernel
        LPKernel=round(LPKernel*Binning);
    end
    if HPMethod
        switch HPMethod
            case 'ribot'
                LPKernel=LPKernel;
            otherwise
                HPKernel=round(HPKernel*Binning);
        end
    end
    if flagfiltermask~=0
        filtermask2=imresize(filtermask, Binning, 'bilinear');
    end
end


if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
    firfilter=OIMakFIRfilterH(newsize,LPMethod, LPKernel, HPMethod, HPKernel);
end    

if flagblanksubtraction == 2
	Bcounter1=zeros(size(Bmap,1));
	Bcounter2=zeros(size(Bmap,1));
	blankgs1=zeros(size(Bmap,1));
	blankgs2=zeros(size(Bmap,1));
% 	Bsample1=zeros(newsize, newsize, 'single');  
% 	Bsample2=zeros(newsize, newsize, 'single');
	Bsample1=zeros(newsize, newsize, size(Bmap,1), 'single');  
	Bsample2=zeros(newsize, newsize, size(Bmap,1), 'single');
	for k=blockselect
        blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
        fprintf('%s:', blkfilename);

		for j=1:size(Bmap,1)
			Bcondlist1=getfield(cell2struct(Bmap(j, 2), 'junk'), 'junk');
			for m=Bcondlist1; %1:NStim
		%             counter0=counter0+1;
				if goodstim(k, m)
					Bcounter1(j)=Bcounter1(j)+1;
					if flagrandom
						stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
					else
						stimloc=m;
					end
		%             stimloc2(m)=counter0;
					fprintf(' stim%d ', stimloc);
					AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);

					if dalsalineremove == 1
						for fr = [fframe framerange]  % 150221 HT
							 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
						end
					end
					if flagalign == 2
						for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
							 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
						end
					end      
					
					if Binning ==1
						Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
						Rangetemp2=mean(AllFrames(:,:,framerange), 3);
					else
						Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
						Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
					end

					Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

					if LPKernel || HPKernel
						if flagfiltermask==1
							Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
						else
							Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end
					end

					if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
						Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
					end

					if flagalign == 1
						if Binning ==1
							Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
						else
							Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
						end
					end

					Bsample1(:,:,j)=Bsample1(:,:,j)+Sampletemp;
				end
			end % for m

			Bcondlist2=getfield(cell2struct(Bmap(j, 3), 'junk'), 'junk');        	
			for m=Bcondlist2; %1:NStim
		%             counter0=counter0+1;
				if goodstim(k, m)
					Bcounter2(j)=Bcounter2(j)+1;
					if flagrandom
						stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
					else
						stimloc=m;
					end
		%             stimloc2(m)=counter0;
					fprintf(' stim%d ', stimloc);
					AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);

					if dalsalineremove == 1
						for fr = [fframe framerange]  % 150221 HT
							 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
						end
					end
					if flagalign == 2
						for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
							 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
						end
					end      

					if Binning ==1
						Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
						Rangetemp2=mean(AllFrames(:,:,framerange), 3);
					else
						Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
						Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
					end

					Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

					if LPKernel || HPKernel
						if flagfiltermask==1
							Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
						else
							Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end
					end

					if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
						Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
					end

					if flagalign == 1
						if Binning ==1
							Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
						else
							Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
						end
					end

					Bsample2(:,:,j)=Bsample2(:,:,j)+Sampletemp;
				end
			end % for m
		end % for j=1:size(Bmap,1)
		fprintf('\r');
	end % for k=blockselect
	
	for j=1:size(Bmap,1)
		if Bcounter1(j)
			Bsample1(:,:,j)=Bsample1(:,:,j)/Bcounter1(j);
			blankgs1(j)=1;
		else
			Bsample1(:,:,j)=0;
			blankgs1(j)=0;
		end
		if Bcounter2(j)
			Bsample2(:,:,j)=Bsample2(:,:,j)/Bcounter2(j);
			blankgs2(j)=1;
		else
			Bsample2(:,:,j)=0;
			blankgs2(j)=0;
		end
	end
end % fi flagblanksubtraction == 2

ctime0=clock;
flagfirst=0;
% counter1=0;
% counter2=0;
counter5=0;
for j=1:size(Smap,1)
    fprintf('Statistical map: %s\r',getfield(cell2struct(Smap(j, 1),'junk'),'junk'));
    counter1=0;
    counter2=0;
    NStimToRead=size(Smap,2)+size(Smap,3);
    Sample1=zeros(newsize, newsize, blockselectnum*NStimToRead, 'single');  
    Sample2=zeros(newsize, newsize, blockselectnum*NStimToRead, 'single');
    for k=blockselect
        blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
        if blockselectnum>1
            if counter5>0 && flagfirst==0 % && k==blockselect(2)    % just for estimate process time 
                blockproctime=etime(clock, ctime0);
                fprintf('\rTime for one block: %7.4f secs', blockproctime);
                temptime=blockproctime*blockselectnum*size(Smap,1);
                fprintf('\rTime for all (%d) blocks: %5.2f minutes (%5.2f hours)\r', blockselectnum, temptime/60, temptime/3600);
                flagfirst=1;
            end
        end
        fprintf('%s:', blkfilename);
%         counter1=counter1+1;
        counter5=counter5+1;
        counter0=0;
    %     fid=fopen(strcat(blockfolder, blkfilename),'rb','ieee-le'); % open block file
        if flagblanksubtraction == 1
            Bcounter1=0;
            Bsample1=zeros(newsize, newsize, 'single');  
            Bcondlist1=getfield(cell2struct(Bmap(j, 2), 'junk'), 'junk');
            for m=Bcondlist1; %1:NStim
    %             counter0=counter0+1;
                if goodstim(k, m)
                    Bcounter1=Bcounter1+1;
                    if flagrandom
                        stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
                    else
                        stimloc=m;
                    end
        %             stimloc2(m)=counter0;
                    fprintf(' stim%d ', stimloc);
                    AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                    
					if dalsalineremove == 1
						for fr = [fframe framerange]  % 150221 HT
							 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
						end
					end
					if flagalign == 2
						for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
							 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
						end
					end      
                    
                    if Binning ==1
                        Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
                        Rangetemp2=mean(AllFrames(:,:,framerange), 3);
                    else
                        Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
                        Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
                    end

                    Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                    if LPKernel || HPKernel
                    	if flagfiltermask==1
                    		Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
                    	else
							Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end
					end

                    if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                        Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
                    end

                    if flagalign == 1
                        if Binning ==1
                            Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        else
                            Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                        end
                    end

                    Bsample1=Bsample1+Sampletemp;
                end
            end % for m
            if Bcounter1
                Bsample1=Bsample1/Bcounter1;
                blankgs1=1;
            else
                Bsample1=0;
                blankgs1=0;
            end
            
            Bcounter2=0;
            Bsample2=zeros(newsize, newsize, 'single');
            Bcondlist2=getfield(cell2struct(Bmap(j, 3), 'junk'), 'junk');
            for m=Bcondlist2; %1:NStim
    %             counter0=counter0+1;
                if goodstim(k, m)
                    Bcounter2=Bcounter2+1;
                    if flagrandom
                        stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
                    else
                        stimloc=m;
                    end
        %             stimloc2(m)=counter0;
                    fprintf(' stim%d ', stimloc);
                    AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                    
					if dalsalineremove == 1
						for fr = [fframe framerange]  % 150221 HT
							 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
						end
					end
					if flagalign == 2
						for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
							 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
						end
					end      
                    
                    if Binning ==1
                        Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
                        Rangetemp2=mean(AllFrames(:,:,framerange), 3);
                    else
                        Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
                        Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
                    end

                    Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                    if LPKernel || HPKernel
                    	if flagfiltermask==1
                    		Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
                    	else
							Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end
					end

                    if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                        Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
                    end

                    if flagalign == 1
                        if Binning ==1
                            Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        else
                            Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                        end
                    end

                    Bsample2=Bsample2+Sampletemp;
                end
            end % for m
            if Bcounter2
                Bsample2=Bsample2/Bcounter2;
                blankgs2=1;
            else
                Bsample2=0;
                blankgs2=0;
            end
        end % if flagblanksubtraction == 1

		if flagblanksubtraction  
			Bcondlist1num=size(getfield(cell2struct(Bmap(j, 2), 'junk'), 'junk'),2);
			Bcondlist2num=size(getfield(cell2struct(Bmap(j, 3), 'junk'), 'junk'),2);
		end


        condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
        for m=condlist1; %1:NStim
%             counter0=counter0+1;
            if goodstim(k, m) && (flagblanksubtraction == 0 || (flagblanksubtraction == 1 && blankgs1) || (flagblanksubtraction == 2 && blankgs1(j)) || (flagblanksubtraction == 1 && Bcondlist1num==0) || (flagblanksubtraction == 2 && Bcondlist1num==0))
                counter1=counter1+1;
                if flagrandom
                    stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
                else
                    stimloc=m;
                end
    %             stimloc2(m)=counter0;
                fprintf(' stim%d ', stimloc);
                AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                
				if dalsalineremove == 1
					for fr = [fframe framerange]  % 150221 HT
						 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
					end
				end
				if flagalign == 2
					for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
						 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
					end
				end      
                
                if Binning ==1
                    Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
                    Rangetemp2=mean(AllFrames(:,:,framerange), 3);
                else
                    Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
                    Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
                end

                Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                if LPKernel || HPKernel
                	if flagfiltermask==1
                		Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
                	else
						Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
					end
				end
                
                if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                    Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
                end
                
                if flagalign == 1
                    if Binning ==1
                        Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    else
                        Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                    end
                end

                if flagblanksubtraction == 1
                    Sample1(:,:,counter1)=Sampletemp-Bsample1;
                elseif flagblanksubtraction == 2 && Bcondlist1num~=0
                    Sample1(:,:,counter1)=Sampletemp-Bsample1(:,:,j);
                else
                    Sample1(:,:,counter1)=Sampletemp;
                end
            end
        end % for m
        
        condlist2=getfield(cell2struct(Smap(j, 3), 'junk'), 'junk');
        for m=condlist2; %1:NStim
%             counter0=counter0+1;
            if goodstim(k, m) && (flagblanksubtraction == 0 || (flagblanksubtraction == 1 && blankgs2) || (flagblanksubtraction == 2 && blankgs2(j)) || (flagblanksubtraction == 1 && Bcondlist2num==0) || (flagblanksubtraction == 2 && Bcondlist2num==0))
                counter2=counter2+1;
                if flagrandom
                    stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
                else
                    stimloc=m;
                end
    %             stimloc2(m)=counter0;
                fprintf(' stim%d ', stimloc);
                AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                
				if dalsalineremove == 1
					for fr = [fframe framerange]  % 150221 HT
						 AllFrames(:,:,fr) = OIDalsalineremover(AllFrames(:,:,fr));
					end
				end
				if flagalign == 2
					for fr = [fframe framerange]  % 150221 HT, for shift image frame by frame
						 AllFrames(:,:,fr) = OIShift(AllFrames(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
					end
				end      
                
                if Binning ==1
                    Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
                    Rangetemp2=mean(AllFrames(:,:,framerange), 3);
                else
                    Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
                    Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
                end

                Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                if LPKernel || HPKernel
                	if flagfiltermask==1
                		Sampletemp=OIEasyFilterH2wMask(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel,filtermask2);
                	else
						Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
					end
				end
                
                if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                    Sampletemp=ifft2(fft2(Sampletemp).*firfilter);
                end

                if flagalign == 1
                    if Binning ==1
                        Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    else
                        Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                    end
                end
                
                if flagblanksubtraction == 1
                    Sample2(:,:,counter2)=Sampletemp-Bsample2;
                elseif flagblanksubtraction == 2 && Bcondlist2num~=0
                    Sample2(:,:,counter2)=Sampletemp-Bsample2(:,:,j);
                else
                    Sample2(:,:,counter2)=Sampletemp;
                end
            end
        end % for m
        fprintf('\r');
    end	% end of "k"
    fprintf('\r');


% 	tempSample1=squeeze(Sample1(:,:,:,j));
% 	tempSample2=squeeze(Sample2(:,:,:,j));
    if parametric==1
        tempSample1=squeeze(Sample1(:,:,1:counter1));
        tempSample2=squeeze(Sample2(:,:,1:counter2));
        if TestType==1
            [H,p,CI,STATS] = ttest(tempSample1, tempSample2, 0.05, tail, 3);    % [H,p,CI,STATS]=TTEST(X,Y,ALPHA,TAIL,DIM) % for paired
        else
            [H,p,CI,STATS] = ttest2(tempSample1, tempSample2, 0.05, tail, 'equal', 3);    %[H,P,CI,STATS] = TTEST2(X,Y,ALPHA,TAIL,VARTYPE,DIM); % for unpaired

        end
    elseif parametric==2    
        Bcounter1=0;
        Bcounter2=0;
        p=ones(newsize,newsize);
        for x=1:newsize
            for y=1:newsize
                Bcounter1=Bcounter1+1;
                tempSample1=squeeze(Sample1(x,y,1:counter1));
                tempSample2=squeeze(Sample2(x,y,1:counter2));
                if TestType==1
                    p(x,y) = signrank(tempSample1,tempSample2,'alpha',0.05);
                    if Bcounter1>=newsize*newsize*0.1*(Bcounter2+1)
                        fprintf('\rWilcoxon signed rank test process: %g%% completed', 10*(Bcounter2+1));
                        Bcounter2=Bcounter2+1;
                    end
                elseif TestType==2
                    p(x,y) = ranksum(tempSample1,tempSample2,'alpha',0.05);
                    if Bcounter1>=newsize*newsize*0.1*(Bcounter2+1)
                        fprintf('\rMann-Whitney U-test process: %g%% completed', 10*(Bcounter2+1));
                        Bcounter2=Bcounter2+1;
                    end
                end
            end
        end
        clear tempSample1;
        clear tempSample2;
        tempSample1=squeeze(Sample1(:,:,1:counter1));
        tempSample2=squeeze(Sample2(:,:,1:counter2));        
    end
    
    map1=squeeze(mean(tempSample1,3)); 
    map2=squeeze(mean(tempSample2,3));


%     fprintf('map1 min: %g   max: %g   median: %g\r', min(min(map1)), max(max(map1)), median(median(map1)));
% 	fprintf('map2 min: %g   max: %g   median: %g\r', min(min(map2)), max(max(map2)), median(median(map2)));
    submap=map1-map2;
%     fprintf('submap min: %g   max: %g   median: %g\r', min(min(submap)), max(max(submap)), median(median(submap)));
	% signmap=double(submap>=0)*2-1;
	signmap=double(submap>=0)*(-2)+1;
%     fprintf('signmap min: %g   max: %g   median: %g\r', min(min(signmap)), max(max(signmap)), median(median(signmap)));
	% max(p)
	% min(p)
	% median(p)
	
	%p_ns=(p>=0.05);
	
	Templ= p > pmin; % logical operation 0 or 1; locations of lower value
	Tempu= p < pmax;% logical operation 0 or 1; locations of higher value
	Tempul= (Templ.*Tempu).*p;% bewteen low and high clips
	Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
	Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
	p2=Tempul2;
	% max(p2)
	% min(p2)
	% median(p2)
	p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));
	% max(p)
	% min(p2_log)
	
	pmap=uint8(round(p2_log.*signmap*127+128));
%     fprintf('pmap min: %g   max: %g   median: %g\r', min(min(pmap)), max(max(pmap)), median(median(pmap)));
	% max(pmap)
	% min(pmap)
	% median(pmap)
	% min(round(p2_log.*signmap*127+128))
    % Smap(j,1);
    mapname=getfield(cell2struct(Smap(j,1), 'junk'), 'junk');
	
    if (j<=9)
        prefix=strcat('0', num2str(j), '_');
    else
        prefix=strcat(num2str(j), '_');
    end
    
	imwrite(pmap, strcat(resultfolder, prefix1, prefix, mapname, '_p-', num2str(pmax,'%1.5f'), '_pmap.bmp'));
	lut=textread(lutfile);  % this color table should be in sunin folder
	imwrite(pmap,lut, strcat(resultfolder, prefix1, prefix, mapname, '_p-', num2str(pmax,'%1.5f'), '_pmapcolor.bmp'));
	
	[maptemp, framemedian, lowClip, highClip] = OIClipH(p, 1, 1);
    imwrite(norm_to_uint8(maptemp), strcat(resultfolder, prefix1, prefix, mapname, '_p-all.bmp'));
	for i=threshold
	    b=double(p<i);
	    imwrite(uint8(b*255), strcat(resultfolder, prefix1, prefix, mapname, '_p-', num2str(i,'%1.5f'), '.bmp'));
	end
	if flagsaveivf
        OIWriteIVF(p, strcat(resultfolder, prefix1, prefix, mapname, '_p.ivf'));             % p values
        OIWriteIVF(STATS.tstat, strcat(resultfolder, prefix1, prefix, mapname, '_t.ivf'));   % t statstic
        OIWriteIVF(STATS.df, strcat(resultfolder, prefix1, prefix, mapname, '_df.ivf'));     % degree of freedom
        OIWriteIVF(STATS.sd, strcat(resultfolder, prefix1, prefix, mapname, '_sd.ivf'));     % stdandard divation of A-B
       	% also output subtraction map for comparison
		OIWriteIVF(submap, strcat(resultfolder, prefix1, prefix, mapname, '_submap.ivf'));
    end
    [maptemp, framemedian, lowClip, highClip] = OIClipH(submap, 1, 1);
	imwrite(norm_to_uint8(maptemp), strcat(resultfolder, prefix1, prefix, mapname, '_submap.bmp'));
    
%     imwrite(norm_to_uint8(OIClip(map1,1,1)), strcat(resultfolder, mapname, '_map1.bmp'));
%     imwrite(norm_to_uint8(OIClip(map2,1,1)), strcat(resultfolder, mapname, '_map2.bmp'));
%     imwrite(norm_to_uint8(map1), strcat(resultfolder, mapname, '_map1.bmp'));
%     imwrite(norm_to_uint8(map2), strcat(resultfolder, mapname, '_map2.bmp'));

%     [maptemp, framemedian, lowClip, highClip] = OIClipH(map1, 1, 1);
% 	imwrite(norm_to_uint8(maptemp), strcat(resultfolder, mapname, '_map1.bmp'));
%     [maptemp, framemedian, lowClip, highClip] = OIClipH(map2, 1, 1);
% 	imwrite(norm_to_uint8(maptemp), strcat(resultfolder, mapname, '_map2.bmp'));
    
%     fprintf('norm_to_uint8(maptemp) min: %g   max: %g\r', min(min(norm_to_uint8(maptemp))), max(max(norm_to_uint8(maptemp))));

    clear Sample1;
    clear Sample2;
    clear tempSample1;
    clear tempSample2;

end %for n=1:size(Smap,1)  

% output a color table
colorbarmap=zeros(60, 256, 3);
for i=1:30
    colorbarmap(i, :, :)=256.*lut(:,:);
end
imwrite(uint8(colorbarmap), strcat(resultfolder, 'colortable.tif'), 'tiff'); 

totaltime=etime(clock, ctime0);
fprintf('\rProcess Finished: ');
fprintf('\n\n\nTotal time used: %8.4f secs (%8.4f hours).\r\r', totaltime, totaltime/3600);

return;

