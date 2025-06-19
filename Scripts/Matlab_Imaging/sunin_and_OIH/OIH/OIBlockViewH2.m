function OIBlockViewH2()
%
% OIBlockViewH2
% (150221) Ver 3.1.1    dalsalineremove=1 were introduced. flagalign=2 were introduced. Hisashi
% (121118) Ver 2.4    'outputimageformat' was introduced. Hisashi
% (090207) Ver 2.3    Bug fix for the case when fframe=[] 
% (090207) Ver 2.2    Comaprison between groups of blocks is available when TestType=3. 
% (090207) Ver 2.1    Critical bugs for spatial filtering were fixed
% (081106) Ver 2    Make a output file frame by frame
% 
%
% blockview: view all frames in one stim or a subtraction.     -HDL 041001-050729

clear all;
% __________________ Start User Input ____________________
system='v';             % 'v' for VDAQ, 'r' for RedShirt
datadriver = 'C:\Users\OMAR\Documents\MATLAB\';     % Data disk name    /Users/omargharbawie/Documents/Exp/110725BooneMotor/run05
datafolder = 'Exp\';   % Data folder name on data disk, results will be saved in outputfolder

outputdriver = 'C:\Users\OMAR\Documents\MATLAB\ExpResults\';     % Data disk name   
outputfolder = 'Results\'; % Output folder name on data disk, results will be saved here
expname = 'MsHowell\Right_Chamber\Awake_Grasping\2014_09_18_MsHowell\'; % Exp folder name (in both data folder and result folder)
runname = 'run0\';      % Run foler name (in both data folder and result folder)

resultfolder='OIBlockViewH2\';        % specify a founder name results to be saved. if not specified, program will generate one like "H2L2B0.5', 
resultname='Frames20-25\LP2_HP400\FFrame1\NoAlign\Masked\';



outputimageformat  = 'png'; % The format for subtraction maps. % Choose from these format: 'bmp', 'gif', 'jpg', 'pcx', 'tif', etc. Since H33

filename = {            % blk file names, leave empty for automatic searching
};
TestType=1;             
% 1: Within condition comparison (fn vs. f0), 
% 2: Compare 2 conditions (stim1 vs. stim2), all first-frame subtracted first
% 3: Compare 2 groups of blocks in a condition .

stim1=[1];    % 10 11 12 13 14 15 16          % stim set 1, this need to be defined for both TestType=1, TestType=2, and , TestType=3
stim2=[];              % stim set 2, this is for TestType=2 only, no use for TestType=1 and TestType=3

fframe=[]; % If fframe==[], no first frame subtraction
framerange=[1]; % This will be ignored when flagframebyframe==1
flagframebyframe=0; %default 0; 1: show results brame by frame

blockselect =[1];  % select blocks for processing, eg. [ 1 3 4 5], leave empty for all block being selected ([]).
blockselect2 =[]; % select blocks for comparsion between groups of blocks, TestType=3.

numblockbinded=100; % number of blocks to be binded. If numblockbinded > the number of actual blocks, all block are binded.

% Clipping
    % 0: no clipping; 
    % 1: clipping at median+-SD (value);  
    % 2: clipping at median+-SD (value) with a mask (default.bmp);  
    % 3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
    % 4: clipping at median+-intensity change (value), Ex. 0.0005 or 0.0008;  
    % 5: clipping at 0+-SD (value);
    % 6: clipping at 0+-intensity change (value);
clipmethod = 1;   
clipvalue  = 1.5; %1.5; %0.001 %1.5 %1.5;     % How many SD to be used for clipping, usually =1 (range is plus and minus 1SD on both sides of median)	% clipvalue = [10, 10; 20, 20; 0.8, 0];   %for method 3 only (window cliping), represents [x1, y1; x2, y2; sd, nouse]
flagshowframeStd=0;

% Filtering
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, apply median filter, then expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', band-pass filtering using FFT (Fast Fourier Transform).
LPMethod='fastmean';    % Low-pass filtering method, usually 'gaussian'
LPKernel=0; %15;             % set as 0 if no LP filter
HPMethod='fastmean';       % High-pass filtering method, 'ribot' filter is fastest
HPKernel=0; %100;             % For 'ribot' filter, it's the order of polynomial fit, usually set as 2, set as 0 if no HP filter
flagmaskfilter=0; % 0: no mask filtering; 1: mask filtering is applied to reduce vessel artifacts. A mask bmp file have to be saved as '/masks/filtermask/default.', outputimageformat. By Hisashi on 080118

Binning=1;             % note: it is the x times of original size (i.e. 0.5 is equal to 2x2 binning)

flaggoodstim=0;     % similar to the 'goodstim' in sunin, if 'goodstim=1', will search for 'goodstim.txt' in block folder
gsfilename = 'goodstim.txt'; % Specify goodstim file. The default is 'goodstim.txt'.
flagalign=0;        % 0: no shift/align; 1: shift/align using the first frame on trial-by-trial. need 'shiftinput.2.txt'; 2: shift/align using all frames on trial-by-trial. need 'shiftinput.4.txt'
flagrandom=0;       % default=0; similar to the 'flagrandom' in sunin, if 'flagrandom=1', will search for 'stimseq.txt' in block folder
dalsalineremove=0;	% It's better to set 1 when flagtrialfilter=2: 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels

%___________________ end of user input _______________
% resultdriver=datadriver;
if isempty(resultfolder)
    resultfolder=strcat('HP',num2str(HPKernel), '_LP', num2str(LPKernel), '_Bin', num2str(Binning), '\');
elseif 	resultfolder(end)~='\';
	resultfolder=[resultfolder, '\'];	
end
resultfolder=[outputdriver, outputfolder, expname, runname, resultfolder];
runfolder=[outputdriver, outputfolder, expname, runname];
blockfolder=[datadriver, datafolder, expname, runname];
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
if isempty(blockselect)
    blockselectnum=size(filename, 1); 
    blockselect=1:blockselectnum;
else
    blockselectnum=size(blockselect,2);
end
if TestType==3
    if isempty(blockselect2)
        fprintf('Testtype==3 and blockselect2 is empty. Set blockselect2 for comparsion.\n');
        blockselectnum2=size(filename, 1); 
        blockselect2=1:blockselectnum;
    else
        blockselectnum2=size(blockselect2,2);
    end
end
blockfilenum=size(filename, 1);      % how many blocks

clipmask = strcat(outputdriver, outputfolder, expname, 'masks\clipmask\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2
% filtermaskname = strcat(outputdriver, outputfolder, expname, 'masks\filtermask\default.bmp');

if flagmaskfilter==1
	filtermaskname = strcat(resultdriver, outputfolder, expname, 'masks\filtermask\default.bmp');
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
fprintf('\nNumber of frames = %d\n', FramesPerStim);
fprintf('Number of Stims = %d\n', NStim);

if flaggoodstim
    strcat(runfolder, gsfilename)
    goodstim=textread(strcat(runfolder, gsfilename), '%d');
    if isempty(goodstim)|size(goodstim,1)~=blockfilenum*NStim
        fprintf('Error, "goodstim.txt" does not contain right number of conditions\r');
    end
	if sum(sum(goodstim))~=sum(sum(goodstim.*goodstim))	% only 0*0=0 and 1*1=1
		fprintf('Error: ''goodstim.txt'' should contain only "1" and "0"'\n);
	end
    goodstim=reshape(goodstim, [NStim, blockfilenum]);
    goodstim=goodstim'; %'
%     if prod(sum(goodstim, 2))==0
%         fprintf('Worning: too few good stims, one or more condition has no data!\n');
%     end
    fprintf('Good stim parameters loaded\r\r');
else
	goodstim=ones(blockfilenum, NStim);	% if no goodstim.txt provided, use all conditions.    
end

if flagalign ==1
    [sfname, sfblock, sfstim, sfframe, sfx1, sfy1, sfcoor]=textread(strcat(runfolder, 'shiftinput.2.txt'), '%s %d %d %d %f %f %f'); %Hisashi 071023
    if size(sfname,1)~=blockfilenum*NStim
        fprintf('Error: flagalign=1, .2 file doesnot contain blockfilenum*NStim number of entries!\r');
    end
    shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
    for k=1:blockfilenum
        for i=1:NStim
            if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NStim+i).sfname
                fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstructname(k*NStim*NFrames).sfname);
            end
        end
    end
elseif flagalign ==2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading frames %%%%%%%%%%%%%%%%%%%%%%%
newsize=round(FrameWidth*Binning);  % assume width=height here
if flagframebyframe==0
	sample1=zeros(newsize, newsize);  
	sample2=zeros(newsize, newsize);
else
	if TestType==1
		sample1=zeros(newsize, newsize, FramesPerStim);  
		sample2=zeros(newsize, newsize);
	elseif TestType==2
		sample1=zeros(newsize, newsize, FramesPerStim);  
		sample2=zeros(newsize, newsize, FramesPerStim);
	elseif TestType==3
		sample1=zeros(newsize, newsize, FramesPerStim);  
		sample2=zeros(newsize, newsize, FramesPerStim);
	end
end

if Binning ==1
    if flagmaskfilter~=0
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
    if flagmaskfilter~=0
        filtermask2=imresize(filtermask, Binning, 'bilinear');
    end
end

if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
    firfilter=OIMakFIRfilterH(newsize,LPMethod, LPKernel, HPMethod, HPKernel);
end

numfframe=size(fframe,2);
blkbindedcount=0;
gscount1=0;
gscount2=0;
gstemp=0;
ctime0=clock;
flagfirst=0;
counter=0;
for k=blockselect
    blkbindedcount=blkbindedcount+1;
    if blkbindedcount==1
    	firstblk=k;
    end
    blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
    if blockselectnum>1
        if counter>0 && flagfirst==0 % && k==blockselect(2)    % just for estimate process time 
            blockproctime=etime(clock, ctime0);
            fprintf('\rTime for one block: %7.4f secs', blockproctime);
            for t=stim1
                gstemp=gstemp+sum(goodstim(:,t));
            end
            if TestType==2
                for t=stim2
                    gstemp=gstemp+sum(goodstim(:,t));
                end
            end
            temptime=blockproctime/counter*gstemp;
            fprintf('\rTime for all blocks: %5.2f minutes (%5.2f hours)\r', temptime/60, temptime/3600);
            flagfirst=1;
        end
    end
    fprintf('%s\r', blkfilename);
    
    for j=1:size(stim1,2)
        if flagrandom
            stimloc=find(stimseq(k, :)==stim1(j));  % location of stim1(j) in this randomnized trial.
        else
            stimloc=stim1(j);
        end
        if goodstim(k, stimloc)
            gscount1=gscount1+1;
            counter=counter+1;
            Stimtemp=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
            
            if dalsalineremove == 1
				for fr = 1:FramesPerStim  % 150221 HT
					 Stimtemp(:,:,fr) = OIDalsalineremover(Stimtemp(:,:,fr));
				end
            end
			if flagalign == 2
				for fr = 1:FramesPerStim  % 150221 HT, for shift image frame by frame
					 Stimtemp(:,:,fr) = OIShift(Stimtemp(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
				end
			end                 

            if flagframebyframe==0
                if numfframe
                    Range1=mean(Stimtemp(:,:,fframe), 3);
                end
                Range2=mean(Stimtemp(:,:,framerange), 3);
                if flagalign == 1
                    if numfframe
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    Range2=OIShift(Range2, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                end
                if Binning ~=1
                    if numfframe
                        Range1=imresize(Range1, Binning, 'bilinear'); % assume 'bilinear' is the best
                    end
                    Range2=imresize(Range2, Binning, 'bilinear');
                end
                if TestType==1
                    sample2(:,:)=sample2(:,:)+0;
                end
                if flagmaskfilter==0
                    if numfframe
                        sampletemp=OIEasyFilterH2((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                    else
                        sampletemp=OIEasyFilterH2(Range2, LPMethod, LPKernel, HPMethod, HPKernel);
                    end
                else
                    if numfframe
                        sampletemp=OIEasyFilterH2wMask((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                    else
                        sampletemp=OIEasyFilterH2wMask(Range2, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);                    
                    end
                end
                
                if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                    sampletemp=ifft2(fft2(sampletemp).*firfilter);
                end
                
                sample1(:,:)=sample1(:,:)+sampletemp;
                clear sampletemp;
            else % if flagframebyframe==0
                if numfframe
                    Range1=mean(Stimtemp(:,:,fframe), 3);
                    if flagalign == 1
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        Range1=imresize(Range1, Binning, 'bilinear'); % assume 'bilinear' is the best
                    end
                end
                if TestType==1
                    sample2(:,:)=sample2(:,:)+0;
                end
                for f=1:FramesPerStim
                    if flagalign == 1
                        Stimtemp(:,:,f)=OIShift(Stimtemp(:,:,f), sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        sampletemp=imresize(Stimtemp(:,:,f), Binning, 'bilinear');
                    else
                        sampletemp=Stimtemp(:,:,f);
                    end
                    if numfframe
                        if flagmaskfilter==0
                            sampletemp=OIEasyFilterH2((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                        else
                            sampletemp=OIEasyFilterH2wMask((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                        end
                    else
                        if flagmaskfilter==0
                            sampletemp=OIEasyFilterH2(sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
                        else
                            sampletemp=OIEasyFilterH2wMask(sampletemp, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                        end                        
                    end
                    
                    if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                        sampletemp=ifft2(fft2(sampletemp).*firfilter);
                    end
                    
                    sample1(:,:,f)=sample1(:,:,f)+sampletemp;
                    clear sampletemp;
                end
            end % if flagframebyframe==0
        end %end for 'goodstim(k, stimloc)'
    end % end for j
    if TestType==2
        for j=1:size(stim2,2)
            if flagrandom
                stimloc=find(stimseq(k, :)==stim2(j));  % location of stim1(j) in this randomnized trial.
            else
                stimloc=stim2(j);
            end
            if goodstim(k, stimloc)
                gscount2=gscount2+1;
                Stimtemp=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                
				if dalsalineremove == 1
					for fr = 1:FramesPerStim  % 150221 HT
						 Stimtemp(:,:,fr) = OIDalsalineremover(Stimtemp(:,:,fr));
					end
				end                
				if flagalign == 2
					for fr = 1:FramesPerStim  % 150221 HT, for shift image frame by frame
						 Stimtemp(:,:,fr) = OIShift(Stimtemp(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
					end
				end
				                
                if flagframebyframe==0
                    Range1=mean(Stimtemp(:,:,fframe), 3);  
                    Range2=mean(Stimtemp(:,:,framerange), 3);
                    if flagalign == 1
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        Range2=OIShift(Range2, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        Range1=imresize(Range1, Binning, 'bilinear');  % assume 'bilinear' is the best
                        Range2=imresize(Range2, Binning, 'bilinear');
                    end
                    if flagmaskfilter==0
                        sampletemp=OIEasyFilterH2((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                    else
                        sampletemp=OIEasyFilterH2wMask((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                    end
                    
                    if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                        sampletemp=ifft2(fft2(sampletemp).*firfilter);
                    end
                    
                    sample2(:,:)=sample2(:,:)+sampletemp;
                    clear sampletemp;
                else % if flagframebyframe==0
                    Range1=mean(Stimtemp(:,:,fframe), 3);  
                    if flagalign == 1
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        Range1=imresize(Range1, Binning, 'bilinear');  % assume 'bilinear' is the best
                    end
                    for fr=1:FramesPerStim
                        if flagalign == 1
                            Stimtemp(:,:,f)=OIShift(Stimtemp(:,:,f), sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        end
                        if Binning ~=1
                            sampletemp=imresize(Stimtemp(:,:,fr), Binning, 'bilinear');
                        else
                            sampletemp=Stimtemp(:,:,fr);
                        end
                        if flagmaskfilter==0
                            sampletemp=OIEasyFilterH2((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                        else
                            sampletemp=OIEasyFilterH2wMask((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                        end
                        
                        if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                            sampletemp=ifft2(fft2(sampletemp).*firfilter);
                        end
                        
                        sample2(:,:,fr)=sample2(:,:,fr)+sampletemp;
                        clear sampletemp;
                    end
                end %if flagframebyframe==0
            end %end for 'if goodstim(k, stimloc)'
        end % end for j
    end %end for Testtype==2
    if ((blkbindedcount == numblockbinded) | (blkbindedcount == blockselectnum)) && TestType ~= 3
        if TestType==1
            gscount2=gscount1;
        end
    	if flagframebyframe==0
            if numblockbinded > 0
                if blockfilenum<100
                    outputname=strcat(resultfolder, resultname, '_Blk', twodigitsH(firstblk), '_', twodigitsH(k), '.', outputimageformat);
                else
                    outputname=strcat(resultfolder, resultname, '_Blk', threedigitsH(firstblk), '_', threedigitsH(k), '.', outputimageformat);
                end
            else
                fprintf('Error: numblockbinded should be more than 0\r');
            end
            if (gscount1 ~= 0)  && (gscount2 ~= 0)
                map1=sample1/gscount1;
                map2=sample2/gscount2;
                submap=map1-map2;
            else
                submap=zeros(newsize, newsize);
            end
            [maptemp, framemedian, lowClip, highClip] = OIClipH(submap,clipmethod,clipvalue);
			imwrite(norm_to_uint8(maptemp), outputname);
            if (clipmethod==1 || clipmethod==2 ||   clipmethod==5) && flagshowframeStd==1
                if flagshowframeStd==1
                    fprintf('\nframeStd of map %g: %8.8f.', k, (highClip-framemedian)/clipvalue);
                end
            end
			sample1=sample1*0;
			sample2=sample2*0;
		else
			for fr=1:FramesPerStim
                if numblockbinded > 0
                    if blockfilenum<100
                        outputname=strcat(resultfolder, resultname, '_Blk', twodigitsH(firstblk), '_', twodigitsH(k), '_F', twodigitsH(fr), '.', outputimageformat);
                    else
                        outputname=strcat(resultfolder, resultname, '_Blk', threedigitsH(firstblk), '_', threedigitsH(k), '_F', twodigitsH(fr), '.', outputimageformat);
                    end
                else
                    fprintf('Error: numblockbinded should be more than 0\r');
                end
				if TestType==1
					if (gscount1 ~= 0)  & (gscount2 ~= 0)
						map1=sample1(:,:,fr) /gscount1; 
						map2=sample2/gscount2;
						submap=map1-map2;
					else
						submap=zeros(newsize, newsize);
					end       
				elseif TestType==2
					if (gscount1 ~= 0)  & (gscount2 ~= 0)
						map1=sample1(:,:,fr)/gscount1; 
						map2=sample2(:,:,fr)/gscount2;
						submap=map1-map2;
					else
						submap=zeros(newsize, newsize);
                    end
                end
                [maptemp, framemedian, lowClip, highClip] = OIClipH(submap,clipmethod,clipvalue);
                imwrite(norm_to_uint8(maptemp), outputname);
                if (clipmethod==1 || clipmethod==2 || clipmethod==5) && flagshowframeStd==1
                    if flagshowframeStd==1
                        fprintf('\nframeStd of map %g: %8.8f.', k, (highClip-framemedian)/clipvalue);
                    end
                end
			end
			sample1=sample1*0;
			sample2=sample2*0;
		end
		blkbindedcount=0;
		gscount1=0;
		gscount2=0;        
    end
end	% end of "k"
%fprintf('\r');
if TestType==3
    blkbindedcount2=0;
%     gscount1b=0;
    gscount2=0;
    for k=blockselect2
        blkbindedcount2=blkbindedcount2+1;
        if blkbindedcount2==1
            firstblk=k;
        end
        blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
        fprintf('%s\r', blkfilename);

        for j=1:size(stim1,2)
            if flagrandom
                stimloc=find(stimseq(k, :)==stim1(j));  % location of stim1(j) in this randomnized trial.
            else
                stimloc=stim1(j);
            end
            if goodstim(k, stimloc)
                gscount2=gscount2+1;
                Stimtemp=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
                
				if dalsalineremove == 1
					for fr = 1:FramesPerStim  % 150221 HT
						 Stimtemp(:,:,fr) = OIDalsalineremover(Stimtemp(:,:,fr));
					end
				end
				if flagalign == 2
					for fr = 1:FramesPerStim  % 150221 HT, for shift image frame by frame
						 Stimtemp(:,:,fr) = OIShift(Stimtemp(:,:,fr), sfx3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr), sfy3((k-1)*NStim*FramesPerStim+(stimloc-1)*FramesPerStim+fr));
					end
				end    
                if flagframebyframe==0
                    Range1=mean(Stimtemp(:,:,fframe), 3);  
                    Range2=mean(Stimtemp(:,:,framerange), 3);
                    if flagalign == 1
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        Range2=OIShift(Range2, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        Range1=imresize(Range1, Binning, 'bilinear'); % assume 'bilinear' is the best
                        Range2=imresize(Range2, Binning, 'bilinear');
                    end
                    if flagmaskfilter==0
                        sampletemp=OIEasyFilterH2((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                    else
                        sampletemp=OIEasyFilterH2wMask((Range2-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                    end
                    
                    if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                        sampletemp=ifft2(fft2(sampletemp).*firfilter);
                    end
                    
                    sample2(:,:)=sample2(:,:)+sampletemp;
                    
                    clear sampletemp;
                else % if flagframebyframe==0
                    Range1=mean(Stimtemp(:,:,fframe), 3);
                    if flagalign == 1
                        Range1=OIShift(Range1, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                    end
                    if Binning ~=1
                        Range1=imresize(Range1, Binning, 'bilinear'); % assume 'bilinear' is the best
                    end
                    for fr=1:FramesPerStim
                        if flagalign == 1
                            Stimtemp(:,:,fr)=OIShift(Stimtemp(:,:,fr), sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                        end
                        if Binning ~=1
                            sampletemp=imresize(Stimtemp(:,:,fr), Binning, 'bilinear');
                        else
                            sampletemp=Stimtemp(:,:,fr);
                        end
                        if flagmaskfilter==0
                            sampletemp=OIEasyFilterH2((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel);
                        else
                            sampletemp=OIEasyFilterH2wMask((sampletemp-Range1)./Range1, LPMethod, LPKernel, HPMethod, HPKernel, filtermask);
                        end
                        
                        if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
                            sampletemp=ifft2(fft2(sampletemp).*firfilter);
                        end
                        
                        sample2(:,:,fr)=sample2(:,:,fr)+sampletemp;
                        clear sampletemp;
                    end
                end % if flagframebyframe==0
            end %end for 'goodstim(k, stimloc)'
        end % end for j
        if (blkbindedcount2 == blockselectnum2)
            if flagframebyframe==0
                if numblockbinded > 0
                    if blockfilenum<100
                        outputname=strcat(resultfolder, resultname, '_Blk', twodigitsH(firstblk), '_', twodigitsH(k), '.', outputimageformat);
                    else
                        outputname=strcat(resultfolder, resultname, '_Blk', threedigitsH(firstblk), '_', threedigitsH(k), '.', outputimageformat);
                    end
                else
                    fprintf('Error: numblockbinded should be more than 0\r');
                end
                if (gscount1 ~= 0)  & (gscount2 ~= 0)
                    map1=sample1/gscount1;
                    map2=sample2/gscount2;
                    submap=map1-map2;
                else
                    submap=zeros(newsize, newsize);
                end
                [maptemp, framemedian, lowClip, highClip] = OIClipH(submap,clipmethod,clipvalue);
                imwrite(norm_to_uint8(maptemp), outputname);
                if (clipmethod==1 || clipmethod==2 ||   clipmethod==5) && flagshowframeStd==1
                    if flagshowframeStd==1
                        fprintf('\nframeStd of map %g: %8.8f.', k, (highClip-framemedian)/clipvalue);
                    end
                end
                sample1=sample1*0;
                sample2=sample2*0;
            else
                for fr=1:FramesPerStim
                    if numblockbinded > 0
                        if blockfilenum<100
                            outputname=strcat(resultfolder, resultname, '_Blk', twodigitsH(firstblk), '_', twodigitsH(k), '_F', twodigitsH(fr), '.', outputimageformat);
                        else
                            outputname=strcat(resultfolder, resultname, '_Blk', threedigitsH(firstblk), '_', threedigitsH(k), '_F', twodigitsH(fr), '.', outputimageformat);
                        end
                    else
                        fprintf('Error: numblockbinded should be more than 0\r');
                    end
                    if TestType==1
                        if (gscount1 ~= 0)  & (gscount2 ~= 0)
                            map1=sample1(:,:,fr) /gscount1; 
                            map2=sample2/gscount2;
                            submap=map1-map2;
                        else
                            submap=zeros(newsize, newsize);
                        end       
                    elseif TestType==2 || TestType==3
                        if (gscount1 ~= 0)  & (gscount2 ~= 0)
                            map1=sample1(:,:,fr)/gscount1; 
                            map2=sample2(:,:,fr)/gscount2;
                            submap=map1-map2;
                        else
                            submap=zeros(newsize, newsize);
                        end
                    end
                    [maptemp, framemedian, lowClip, highClip] = OIClipH(submap,clipmethod,clipvalue);
                    imwrite(norm_to_uint8(maptemp), outputname);
                    if (clipmethod==1 || clipmethod==2 || clipmethod==5) && flagshowframeStd==1
                        if flagshowframeStd==1
                            fprintf('\nframeStd of map %g: %8.8f.', k, (highClip-framemedian)/clipvalue);
                        end
                    end
                end
                sample1=sample1*0;
                sample2=sample2*0;
            end
            blkbindedcount2=0;
    %         gscount1
    %         gscount2
            gscount1=0;
            gscount2=0;        
        end
    end	% end of "k"
end

return;
% if run into memory problem, try bin more or try modify this program to do pixel-wise test
