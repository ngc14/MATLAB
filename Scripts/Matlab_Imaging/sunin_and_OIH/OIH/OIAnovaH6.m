function OIAnovaH6()
% OIAnovaH: perform one-way ANOVA pixel-by-pixel
% Ver 6.1.2: (150224) Frame-by-frame processes were improved. "flagmaskfilter" was omitted. Hisashi
% Ver 6.0: (150221) dalsalineremove=1 were introduced. flagalign=2 were introduced. Hisashi
% Ver 5.1.1: (150219) When flagblanksubtraction==2, blank subtraction will be performed after averaging all blank condition. Hisashi
% Ver 5.0: (111004) Blank subtraction before the analysis was introduced.  Hisashi
% Ver 4.2.2: (100525) 'gsfilename' was introduced. Hisashi
% Ver 4.2.1: (100525) The default filtermask location was chnaged Hisashi
% Ver 4.1: (091224) three-way ANOVA was added. Hisashi
% Ver 4.0: (091224) two-way ANOVA was added. A bug was fixed in one-way ANOVA Hisashi
% Ver 3.0: 
% Ver 2.0: use Single-Precision to reduce memory usage. 080124 Hisashi
%             can use Goodstim.
% Ver 1.0: can use Shift align.
% modified from blockview -HDL 041001-050729

clear all;
% __________________ Start User Input ____________________
system='v';     % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'K:\';     % Data disk name
datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in outputfolder
outputdriver = 'D:\Vandy\';     % Data disk name
outputfolder = 'expt1\'; % Output folder name on data disk, results will be saved here
expname = '101231JerL2G_Atn\'; % Exp folder name (in both data folder and result folder)
runname = 'run98\';      % Run foler name (in both data folder and result folder)
resultfolder='OIAnova\OIAnovaH6\';        % specify a founder name results to be saved. if not specified, program will generate one like "H2L2B0.5', 
lutfile = 'statcolorBR.lut';

filename = {            % blk file names, leave empty for automatic searching
};

anovatype=3; % 1: one-way ANOVA; 2: two-way ANOVA

if anovatype==1
    % Stim condition list for one-way ANOVA
    stim(1,:)=[2 6 11 15];
    stim(2,:)=[4 8 13 17]; 
    stim(3,:)=[3 7 12 16];
    stim(4,:)=[5 9 14 18];  
elseif anovatype==2
    % Stim condition list for two-way ANOVA
    stim(1,1,:)=[2 6 11 15];
    stim(1,2,:)=[4 8 13 17];
    stim(2,1,:)=[3 7 12 16];
    stim(2,2,:)=[5 9 14 18];
elseif anovatype==3
    % Stim condition list for three-way ANOVA
    stim(1,1,1,:)=[2 6]; %(color,orientation,attended location)
    stim(1,2,1,:)=[4 8];
    stim(2,1,1,:)=[3 7];
    stim(2,2,1,:)=[5 9];
    stim(1,1,2,:)=[11 15];
    stim(1,2,2,:)=[13 17];
    stim(2,1,2,:)=[12 16];
    stim(2,2,2,:)=[14 18];
end

flagblanksubtraction=0; % default = 0
% 0 ... no blank subtraction; 1 ... subtraction within each block file; 2 ... subtraction by the average of blank response across all blocks
if flagblanksubtraction
%     Bstim(1,1,1,:)=[1]; %(color,orientation,attended location)
%     Bstim(1,2,1,:)=[1];
%     Bstim(2,1,1,:)=[1];
%     Bstim(2,2,1,:)=[1];
%     Bstim(1,1,2,:)=[10];
%     Bstim(1,2,2,:)=[10];
%     Bstim(2,1,2,:)=[10];
%     Bstim(2,2,2,:)=[10];

end


fframe=[1:2];
framerange=[5:16];
blockselect =[];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).

threshold=[0.0001 0.001 0.01 0.05 0.1];	% threshold p values, just for generating thresholded p maps, each value will have a thresholded map (i.e. 0.05 will generate a p<0.05 map)
pmax=0.05;
pmin=0.00001;

% clipsd=1;				% for map output, usually 1
LPMethod='fastmean';    % Low-pass filtering method, usually 'gaussian'
LPKernel=0;             % set as 0 if no LP filter
HPMethod='fastmean';       % High-pass filtering method, 'ribot' filter is fastest
HPKernel=0;             % For 'ribot' filter, it's the order of polynomial fit, usually set as 2, set as 0 if no HP filter
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, apply median filter, then expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', band-pass filtering using FFT (Fast Fourier Transform).
Binning=0.5;             % note: it is the x times of original size (i.e. 0.5 is equal to 2x2 binning)

flaggoodstim=1;     % similar to the 'goodstim' in sunin, if 'goodstim=1', will search for 'goodstim.txt' in block folder
gsfilename='goodstim.txt';
flagalign=0;        % 0: no shift/align; 1: shift/align using the first frame on trial-by-trial. need 'shiftinput.2.txt'; 2: shift/align using all frames on trial-by-trial. need 'shiftinput.4.txt'
flagrandom=0;       % default=0; similar to the 'flagrandom' in sunin, if 'flagrandom=1', will search for 'stimseq.txt' in block folder
flagsaveivf=0;      % default=0; 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
dalsalineremove=0;	% It's better to set 1 when flagtrialfilter=2: 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels

%___________________ end of user input _______________
resultdriver=outputdriver;
if isempty(resultfolder)
    resultfolder=strcat('HP',num2str(HPKernel), '_LP', num2str(LPKernel), '_Bin', num2str(Binning), '\');
elseif 	resultfolder(end)~='\';
	resultfolder=[resultfolder, '\'];	
end

resultfolder=[resultdriver, outputfolder, expname, runname, resultfolder];
runfolder=[resultdriver, outputfolder, expname, runname];
blockfolder=[datadriver, datafolder, expname, runname];
% filtermaskname = strcat(resultdriver, outputfolder, expname, 'masks\filtermask\default.bmp');
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
%     for i=1:size(filename,1)
%         fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
%     end
    fprintf('\nFound %d blk files(sorted, check sequence).\n', size(filename,1));
end
blockfilenum=size(filename, 1);      % how many blocks
if isempty(blockselect)
    blockselect=1:blockfilenum;
end
blockselectnum=size(blockselect,2);
stimselectnum=size(stim,2)*size(stim,2);

fprintf('total block number = %d\r\n', blockselectnum);

% if flagfiltermask==1
% 	filtermaskname = strcat(runfolder, 'masks\filtermask\default.bmp');
%     filtermasktemp = imread(filtermaskname, 'bmp');
%     if size(filtermasktemp,3)==1
%         filtermask = double(filtermasktemp)/255;
%     else
%         filtermask = double(filtermasktemp(:,:,1))/255;
%     end
%     clear filtermasktemp;
% %     if ~(strcmp(LPMethod, 'fastmean') && strcmp(LPMethod, 'gaussian') && strcmp(HPMethod, 'fastmean') && strcmp(HPMethod, 'gaussian'))
% %         fprintf('Mask filtering is available only for "fastmean" and "gaussian" currently. The other methods work without mask filtering.\r\r');
% %     end
% else
%     filtermask = ones(504,504,'double');
% end

anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(5), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;


% if flaggoodstim
%     goodstim=textread(strcat(runfolder, 'goodstim.txt'), '%d');
%     if isempty(goodstim)|size(goodstim,1)~=blockfilenum*NStim
%         fprintf('Error, "goodstim.txt" does not contain right number of conditions\r');
%     end
%     goodstim=reshape(goodstim, [NStim, blockfilenum]);
%     goodstim=goodstim';
% else
% 	goodstim=ones(blockfilenum, NStim);	% if no goodstim.txt provided, use all conditions.    
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
    file=[blockfolder, '_stimseq.txt'];
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading frames %%%%%%%%%%%%%%%%%%%%%%%
if anovatype==1
    ngroup=size(stim,1);
    nstimpercondition=size(stim,2);
    if flagblanksubtraction
    	nblankpercondition=size(Bstim,2);
    end
    newWidth=round(FrameWidth*Binning);
    newHeight=round(FrameHeight*Binning);
    Obsrv=NaN(newWidth, newHeight, blockselectnum*ngroup*nstimpercondition, 'single');
    Group=NaN(newWidth, newHeight, blockselectnum*ngroup*nstimpercondition, 'single');
    p=zeros(newWidth, newHeight);
    fprintf('*** One-way ANOVA ***\r');
elseif anovatype==2
    ngroup1=size(stim,1);
    ngroup2=size(stim,2);
    nstimpercondition=size(stim,3);
    if flagblanksubtraction
    	nblankpercondition=size(Bstim,3);
    end
    newWidth=round(FrameWidth*Binning);
    newHeight=round(FrameHeight*Binning);
    Obsrv=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*nstimpercondition, 'single');
    Group1=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*nstimpercondition, 'single');
    Group2=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*nstimpercondition, 'single');
    p=zeros(newWidth, newHeight, 3);
    fprintf('*** Two-way ANOVA ***\r');
elseif anovatype==3
    ngroup1=size(stim,1);
    ngroup2=size(stim,2);
    ngroup3=size(stim,3);
    nstimpercondition=size(stim,4);
    if flagblanksubtraction
    	nblankpercondition=size(Bstim,4);
    end    
    newWidth=round(FrameWidth*Binning);
    newHeight=round(FrameHeight*Binning);
    Obsrv=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*ngroup3*nstimpercondition, 'single');
    Group1=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*ngroup3*nstimpercondition, 'single');
    Group2=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*ngroup3*nstimpercondition, 'single');
    Group3=NaN(newWidth, newHeight, blockselectnum*ngroup1*ngroup2*ngroup3*nstimpercondition, 'single');
    p=zeros(newWidth, newHeight, 7);
    fprintf('*** Three-way ANOVA ***\r');
end    
    
if Binning ==1
%     filtermask2=filtermask;
else
    LPKernel=round(LPKernel*Binning);
    switch HPMethod
        case 'ribot'
            LPKernel=LPKernel;
        otherwise
            HPKernel=round(HPKernel*Binning);
    end
%     filtermask2=imresize(filtermask, Binning, 'bilinear');
end

if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
    fftfilter=OIMakFFTfilterH(newsize,LPMethod, LPKernel, HPMethod, HPKernel);
end    

if flagblanksubtraction == 2
	if anovatype==1
		counter3=zeros(ngroup);
		blanksubtemp=zeros(newWidth, newHeight, ngroup);
		for k=blockselect
			blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
			fprintf('%s:', blkfilename);		
			for o=1:ngroup			
				for j=1:nblankpercondition 
					if flagrandom
						stimloc=find(stimseq(k, :)==Bstim(o,j));  % location of stim(o,j) in this randomnized trial.
					else
						stimloc=Bstim(o,j);
					end
		
					if goodstim(k, stimloc)
						fprintf(' Bstim%d ', stimloc);
						counter3(o)=counter3(o)+1;
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

						Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

						if LPKernel || HPKernel
							Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end

						if flagalign == 1
							if Binning ==1
								Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
							else
								Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
							end
						end
						blanksubtemp(:,:,o)=blanksubtemp(:,:,o)+Obsrvtemp;
					end % goodstim(k, stimloc)
				end % j=1:nblankpercondition
			end % o=1:ngroup
            fprintf('\r');
		end % k=blockselect
		for o=1:ngroup
			blanksubtemp(:,:,o)=blanksubtemp(:,:,o)/counter3(o);
		end
    elseif anovatype==2
		counter3=zeros(ngroup1,ngroup2);
		blanksubtemp=zeros(newWidth,newHeight,ngroup1,ngroup2);
		for k=blockselect
			blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
			fprintf('%s:', blkfilename);		
			for o=1:ngroup1
				for q=1:ngroup2
					for j=1:nblankpercondition 
						if flagrandom
							stimloc=find(stimseq(k, :)==Bstim(o,q,j));  % location of stim(o,j) in this randomnized trial.
						else
							stimloc=Bstim(o,q,j);
						end
		
						if goodstim(k, stimloc)
							fprintf(' Bstim%d ', stimloc);
							counter3(o,q)=counter3(o,q)+1;
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

							Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

							if LPKernel || HPKernel
								Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
							end

							if flagalign == 1
								if Binning ==1
									Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
								else
									Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
								end
							end
							blanksubtemp(:,:,o,q)=blanksubtemp(:,:,o,q)+Obsrvtemp;
						end % goodstim(k, stimloc)
					end % j=1:nblankpercondition
				end % q=1:ngroup2
			end % o=1:ngroup1
            fprintf('\r');
		end % k=blockselect
		for o=1:ngroup1
			for q=1:ngroup2
				blanksubtemp(:,:,o,q)=blanksubtemp(:,:,o,q)/counter3(o,q);
			end
		end
	elseif anovatype==3
		counter3=zeros(ngroup1,ngroup2,ngroup3);
		blanksubtemp=zeros(newWidth,newHeight,ngroup1,ngroup2,ngroup3);
		for k=blockselect
			blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
			fprintf('%s:', blkfilename);		
			for o=1:ngroup1
				for q=1:ngroup2
					for r=1:ngroup3
						for j=1:nblankpercondition 
							if flagrandom
								stimloc=find(stimseq(k, :)==Bstim(o,q,r,j));  % location of stim(o,j) in this randomnized trial.
							else
								stimloc=Bstim(o,q,r,j);
							end
		
							if goodstim(k, stimloc)
								fprintf(' Bstim%d ', stimloc);
								counter3(o,q,r)=counter3(o,q,r)+1;
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

								Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

								if LPKernel || HPKernel
									Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
								end

								if flagalign == 1
									if Binning ==1
										Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
									else
										Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
									end
								end
								blanksubtemp(:,:,o,q,r)=blanksubtemp(:,:,o,q,r)+Obsrvtemp;
							end % goodstim(k, stimloc)
						end % j=1:nblankpercondition
					end % r=1:ngroup3
				end % q=1:ngroup2
			end % o=1:ngroup1
            fprintf('\r');
		end % k=blockselect
		for o=1:ngroup1
			for q=1:ngroup2
				for r=1:ngroup3
					blanksubtemp(:,:,o,q,r)=blanksubtemp(:,:,o,q,r)/counter3(o,q,r);
				end
			end
		end
	end %if anovatype
end % flagblanksubtraction == 2

counter1=0;
counter2=0;
ctime=clock;
ctime0=clock;
flagfirst=0;

for k=blockselect
    if blockselectnum>1
        if counter2>0 && flagfirst==0 % && k==blockselect(2)    % just for estimate process time 
            blockproctime=etime(clock, ctime0);
            fprintf('\rTime for one block: %7.4f secs', blockproctime);
            temptime=blockproctime/counter2*sum(sum(goodstim))*blockselectnum/blockfilenum*stimselectnum/NStim+170;
            fprintf('\rTime for all (%d) blocks: %5.2f min (%5.2f hours)\r', blockselectnum, temptime/60, temptime/3600);
%             fprintf('\rTime for all (%d) blocks: %5.2f hours\r', blockselectnum, blockproctime/counter2*sum(sum(goodstim))*blockselectnum/blockfilenum*stimselectnum/NStim/3600);
            flagfirst=1;
        end
    end
    blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
    fprintf('%s:', blkfilename);
    counter1=counter1+1;
    if anovatype==1
        for o=1:ngroup
        	if flagblanksubtraction == 1
        		counter3=0;
                blanksubtemp=zeros(newWidth, newHeight);
				for j=1:nblankpercondition 
					if flagrandom
						stimloc=find(stimseq(k, :)==Bstim(o,j));  % location of stim(o,j) in this randomnized trial.
					else
						stimloc=Bstim(o,j);
					end
					
					if goodstim(k, stimloc)
						fprintf(' Bstim%d ', stimloc);
                        counter3=counter3+1;
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
	
						Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;
	
						if LPKernel || HPKernel
							Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end
	
						if flagalign == 1
							if Binning ==1
								Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
							else
								Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
							end
						end
	
						blanksubtemp=blanksubtemp+Obsrvtemp;
					end % goodstim(k, stimloc)
				end % j=1:nblankpercondition
				blanksubtemp=blanksubtemp/counter3;
        	end % flagblanksubtraction
            for j=1:nstimpercondition 
                if flagrandom
                    stimloc=find(stimseq(k, :)==stim(o,j));  % location of stim(o,j) in this randomnized trial.
                else
                    stimloc=stim(o,j);
                end
                
                if goodstim(k, stimloc) && (flagblanksubtraction == 0  || (flagblanksubtraction == 1 && counter3) || flagblanksubtraction == 2)
                    fprintf(' stim%d ', stimloc);
                    counter2=counter2+1;
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

                    Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                    if LPKernel || HPKernel
						Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
					end

                    if flagalign == 1
                        if Binning ==1
                            Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        else
                            Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                        end
                    end

%                     Obsrv(:,:,(o+ngroup*(counter1-1)))=Obsrv(:,:,(o+ngroup*(counter1-1)))+Obsrvtemp;
					if flagblanksubtraction == 1
						Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp;
                        clear blanksubtemp;
                        counter3=0;
					elseif flagblanksubtraction == 2
						Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp(:,:,o);
					else
						Obsrv(:,:,(counter2))=Obsrvtemp;
					end                    
                    Group(:,:,(counter2))=o;
                end
            end % j=1:nstimpercondition
        end	% o=1:ngroup
    elseif anovatype==2
        for o=1:ngroup1
            for q=1:ngroup2
				if flagblanksubtraction == 1
					counter3=0;
                    blanksubtemp=zeros(newWidth, newHeight);
					for j=1:nblankpercondition
						if flagrandom
							stimloc=find(stimseq(k, :)==Bstim(o,q,j));  % location of stim(o,j) in this randomnized trial.
						else
							stimloc=Bstim(o,q,j);
						end
						
						if goodstim(k, stimloc)
                            fprintf(' Bstim%d ', stimloc);
							counter3=counter3+1;
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
	
							Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;
	
							if LPKernel || HPKernel
								Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
							end
	
							if flagalign == 1
								if Binning ==1
									Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
								else
									Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
								end
							end
	
							blanksubtemp=blanksubtemp+Obsrvtemp;
						end
					end
					blanksubtemp=blanksubtemp/counter3;
            	end %if flagblanksubtraction == 1
                for j=1:nstimpercondition
                    if flagrandom
                        stimloc=find(stimseq(k, :)==stim(o,q,j));  % location of stim(o,j) in this randomnized trial.
                    else
                        stimloc=stim(o,q,j);
                    end
                    
                    if goodstim(k, stimloc) && (flagblanksubtraction == 0  || (flagblanksubtraction == 1 && counter3) || flagblanksubtraction == 2)
                        fprintf(' stim%d ', stimloc);
                        counter2=counter2+1;
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

                        Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                        if LPKernel || HPKernel
							Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
						end

                        if flagalign == 1
                            if Binning ==1
                                Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                            else
                                Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                            end
                        end

%                         Obsrv(:,:,o,(q+ngroup*(counter1-1)))=Obsrv(:,:,o,(q+ngroup2*(counter1-1)))+Obsrvtemp;
						if flagblanksubtraction == 1
							Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp;
                            clear blanksubtemp;
                            counter3=0;
						elseif flagblanksubtraction == 2
							Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp(:,:,o,q);
						else
							Obsrv(:,:,(counter2))=Obsrvtemp;
						end         						
                        Group1(:,:,(counter2))=o;
                        Group2(:,:,(counter2))=q;
                    end
                end
            end
        end
    elseif anovatype==3
        for o=1:ngroup1
            for q=1:ngroup2
                for r=1:ngroup3
					if flagblanksubtraction == 1
						counter3=0;
                        blanksubtemp=zeros(newWidth, newHeight);
						for j=1:nblankpercondition
							if flagrandom
								stimloc=find(stimseq(k, :)==Bstim(o,q,r,j));  % location of stim(o,j) in this randomnized trial.
							else
								stimloc=Bstim(o,q,r,j);
							end
 							
							if goodstim(k, stimloc)
								fprintf(' Bstim%d ', stimloc);
                                counter3=counter3+1;
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
	
								Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;
	
								if LPKernel || HPKernel
									Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
								end
	
								if flagalign == 1
									if Binning ==1
										Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
									else
										Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
									end
								end
								blanksubtemp=blanksubtemp+Obsrvtemp;
%                                 counter3
							end %if goodstim
						end %for j
						blanksubtemp=blanksubtemp/counter3;
					end %if flagblanksubtraction == 1
                    for j=1:nstimpercondition
                        if flagrandom
                            stimloc=find(stimseq(k, :)==stim(o,q,r,j));  % location of stim(o,j) in this randomnized trial.
                        else
                            stimloc=stim(o,q,r,j);
                        end
                        
                        if goodstim(k, stimloc) && (flagblanksubtraction == 0  || (flagblanksubtraction == 1 && counter3) || flagblanksubtraction == 2)
                            fprintf(' stim%d ', stimloc);
                            counter2=counter2+1;
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

                            Obsrvtemp=(Rangetemp2-Rangetemp1)./Rangetemp1;

                            if LPKernel || HPKernel
								Obsrvtemp=OIEasyFilterH2(Obsrvtemp, LPMethod, LPKernel, HPMethod, HPKernel);
							end

                            if flagalign == 1
                                if Binning ==1
                                    Obsrvtemp=OIShift(Obsrvtemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                                else
                                    Obsrvtemp=OIShift(Obsrvtemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                                end
                            end

							if flagblanksubtraction == 1
								Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp;
                                clear blanksubtemp;
                                counter3=0;
							elseif flagblanksubtraction == 2
								Obsrv(:,:,(counter2))=Obsrvtemp-blanksubtemp(:,:,o,q,r);
							else
								Obsrv(:,:,(counter2))=Obsrvtemp;
							end         						
                            Group1(:,:,(counter2))=o;
                            Group2(:,:,(counter2))=q;
                            Group3(:,:,(counter2))=r;
                        end
                    end
                end    
            end
        end
    end
    fprintf('\r');
end	% end of "k"
% fprintf('\r');

%temp=zeros(size(stim,1), blockselectnum);
%[H,p,CI,STATS] = ttest(Obsrv, Group, 0.05, tail, 3);    % [H,p,CI,STATS]=TTEST(X,Y,ALPHA,TAIL,DIM)  % for unpaired: %[H,P,CI,STATS] = TTEST2(X,Y,ALPHA,TAIL,VARTYPE); 
counter3=0;
counter4=0;
for x=1:newWidth
    for y=1:newHeight
        if anovatype==1
            counter3=counter3+1;
            temp1=squeeze(Obsrv(x,y,1:counter2));
            temp2=squeeze(Group(x,y,1:counter2));
            ptemp = anova1(temp1,temp2,'off');
            p(x,y) = ptemp;
            if counter3>=newWidth*newHeight*0.1*(counter4+1)
                fprintf('\rANOVA process: %g%% completed', 10*(counter4+1));
                counter4=counter4+1;
            end
        elseif anovatype==2
            counter3=counter3+1;
            temp1=squeeze(Obsrv(x,y,1:counter2));
            temp2=squeeze(Group1(x,y,1:counter2));
            temp3=squeeze(Group2(x,y,1:counter2));
            ptemp = anovan(temp1,{temp2 temp3},'model','interaction','display','off');
            p(x,y,:) = ptemp;
%             p(x,y,2) = ptemp(2);
%             p(x,y,3) = ptemp(3);
            if counter3>=newWidth*newHeight*0.1*(counter4+1)
                fprintf('\rANOVA process: %g%% completed', 10*(counter4+1));
                counter4=counter4+1;
            end
        elseif anovatype==3
            counter3=counter3+1;
            temp1=squeeze(Obsrv(x,y,1:counter2));
            temp2=squeeze(Group1(x,y,1:counter2));
            temp3=squeeze(Group2(x,y,1:counter2));
            temp4=squeeze(Group3(x,y,1:counter2));
            ptemp = anovan(temp1,{temp2 temp3 temp4},'model','full','display','off');
            p(x,y,:) = ptemp;
%             p(x,y,2) = ptemp(2);
%             p(x,y,3) = ptemp(3);
%             p(x,y,4) = ptemp(4);
%             p(x,y,5) = ptemp(5);
%             p(x,y,6) = ptemp(6);
%             p(x,y,7) = ptemp(7);
            if counter3>=newWidth*newHeight*0.1*(counter4+1)
                fprintf('\rANOVA process: %g%% completed', 10*(counter4+1));
                counter4=counter4+1;
            end
        end
    end
end
fprintf('\rANOVA process: 100%% completed\r');

if anovatype==1
    Templ= p > pmin; %logical operation 0 or 1; locations of lower value
    Tempu= p < pmax; %logical operation 0 or 1; locations of higher value
    Tempul= (Templ.*Tempu).*p;% bewteen low and high clips
    Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
    Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
    p2=Tempul2;
    p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));

    pmap=uint8(round(p2_log*127+128));

    imwrite(pmap, [resultfolder, 'anovapmap', '_p-', num2str(pmax,'%1.5f'),'.bmp']);
    lut=textread(lutfile);  % this color table should be in sunin folder
    imwrite(pmap,lut, [resultfolder, 'anovapmapcolor', '_p-', num2str(pmax,'%1.5f'), '.bmp']);

    imwrite(norm_to_uint8(OIClip(p,1,1)), [resultfolder, 'p-all.bmp']);
    for i=threshold
        b=double(p<i);
        imwrite(uint8(b*255), [resultfolder, 'p-', num2str(i,'%f'), '.bmp']);
    end

    if flagsaveivf
        OIWriteIVF(p, [resultfolder, 'p.ivf']);             % p values
    end
elseif anovatype==2
    for j=1:3
        switch j
            case 1
                suffix='Factor1';
            case 2
                suffix='Factor2';
            case 3
                suffix='Intrct';
        end
        Templ= p(:,:,j) > pmin; %logical operation 0 or 1; locations of lower value
        Tempu= p(:,:,j) < pmax; %logical operation 0 or 1; locations of higher value
        Tempul= (Templ.*Tempu).*p(:,:,j);% bewteen low and high clips
        Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
        Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
        p2=Tempul2;
        p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));

        pmap=uint8(round(p2_log*127+128));

        imwrite(pmap, [resultfolder, 'Anovapmap_', suffix, '_p-', num2str(pmax,'%1.5f'),'.bmp']);
        lut=textread(lutfile);  % this color table should be in sunin folder
        imwrite(pmap,lut, [resultfolder, 'Anovapmapcolor_', suffix, '_p-', num2str(pmax,'%1.5f'), '.bmp']);

        imwrite(norm_to_uint8(OIClip(p(:,:,j),1,1)), [resultfolder, suffix, '_p-all.bmp']);
        for i=threshold
            b=double(p(:,:,j)<i);
            imwrite(uint8(b*255), [resultfolder, suffix, '_p-', num2str(i,'%f'), '.bmp']);
        end

        if flagsaveivf
            OIWriteIVF(p(:,:,j), [resultfolder, suffix, '_p.ivf']);             % p values
        end        
    end
elseif anovatype==3
    for j=1:7
        switch j
            case 1
                suffix='Factor1';
            case 2
                suffix='Factor2';
            case 3
                suffix='Factor3';
            case 4
                suffix='Intrct1-2';
            case 5
                suffix='Intrct1-3';
            case 6
                suffix='Intrct2-3';
            case 7
                suffix='Intrct1-2-3';
        end
        Templ= p(:,:,j) > pmin; %logical operation 0 or 1; locations of lower value
        Tempu= p(:,:,j) < pmax; %logical operation 0 or 1; locations of higher value
        Tempul= (Templ.*Tempu).*p(:,:,j);% bewteen low and high clips
        Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
        Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
        p2=Tempul2;
        p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));

        pmap=uint8(round(p2_log*127+128));

        imwrite(pmap, [resultfolder, 'Anovapmap_', suffix, '_p-', num2str(pmax,'%1.5f'),'.bmp']);
        lut=textread(lutfile);  % this color table should be in sunin folder
        imwrite(pmap,lut, [resultfolder, 'Anovapmapcolor_', suffix, '_p-', num2str(pmax,'%1.5f'), '.bmp']);

        imwrite(norm_to_uint8(OIClip(p(:,:,j),1,1)), [resultfolder, suffix, '_p-all.bmp']);
        for i=threshold
            b=double(p(:,:,j)<i);
            imwrite(uint8(b*255), [resultfolder, suffix, '_p-', num2str(i,'%f'), '.bmp']);
        end

        if flagsaveivf
            OIWriteIVF(p(:,:,j), [resultfolder, suffix, '_p.ivf']);             % p values
        end        
    end
end

totaltime=etime(clock, ctime0);
fprintf('\rProcess Finished: ');
fprintf('\n\n\nTotal time used: %5.2f min (%5.2f hours).\r\r', totaltime/60, totaltime/3600);


return;
% if run into memory problem, try bin more or try modify this program to do pixel-wise test

