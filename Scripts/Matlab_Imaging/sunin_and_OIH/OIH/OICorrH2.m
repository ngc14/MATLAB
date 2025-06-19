function OICorrH2()
% OICorrH
% Correlation analysis for Optical Imaging Image analysis, written by Hisashi

% Ver 2.1.1: (150224) Frame-by-frame processes were improved. "flagmaskfilter" was omitted. Hisashi
% Ver 2.0: (150221) dalsalineremove=1 were introduced. flagalign=2 were introduced. Hisashi
% Ver 1.1: (101011) Sample number is outputted.
% Ver 1.0: (100411)

clear all;
% __________________ Start User Input ____________________
system='v';     % 'v' or 'r'
datadriver = 'K:\';     % Data disk name
datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in outputfolder
outputdriver = 'D:\Vandy\';     % Data disk name
outputfolder = 'expt1\'; % Output folder name on data disk, results will be saved here
expname = '101231JerL2G_Atn\'; % Exp folder name (in both data folder and result folder)
runname = 'run98\';      % Run foler name (in both data folder and result folder)
resultname='OICorr\OICorrH1\'        % specify a founder name results to be saved. if not specified, program will generate one like "H2L2B0.5', 
statlutfile = 'statcolorBR.lut';
rlutfile = 'jetH2.lut';
prefix1='';

filename = {            % blk file names, leave empty for automatic searching
};

% ___________________


Smap(1, :)={'All', [2:9 11:18]};
Smap(2, :)={'AtnIn', [2:9]};
Smap(3, :)={'AtnOut', [11:18]};



fframe=[1:2]; % leave empty if you want to use raw DC value, not dR/R value.
framerange=[5:16];
blockselect =[];  % select blocks for processing, eg. [1 3 4 5] or [1:5], leave empty for all block being selected ([]).
% precisiontype='single';

% ___________________
TestType=1;             % 1: Pearson's linear correlation coefficient; 2: 'Spearman' computes Spearman's rho
tail='both';            % 'both': Correlation is not zero (the default); 'right': Correlation is greater than zero; 'left': Correlation is less than zero
threshold=[0.0001 0.001 0.01 0.05 0.1];	% threshold p values, just for generating thresholded p maps, each value will have a thresholded map (i.e. 0.05 will generate a p<0.05 map)
pmax=0.05;
pmin=0.00001;

rrange=1;

% ___________________
corrmasktype=2; % '1': for map loading, '2' for coordinates loading (for '1' and '2' need change 'flagcorrmask' to same value
corrmaskname=[          % Names of the masks (need to be same length), these names will be added '.bmp' for loading mask map or add '_x.txt'/'_y.txt' for loading coordinates.
];
domainradiussize=5; %default 5 pixels (80 micro)
savemaskmap=1;      % default=1; 

% ___________________
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
flagsaveivf=1;      % default=0; 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)

dalsalineremove=0;	% It's better to set 1 when flagtrialfilter=2: 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels

% 
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

fprintf('total block number = %d\r\n', blockselectnum);

% if flagfiltermask==1
% 	filtermaskname = strcat(resultdriver, outputfolder, expname, runname, 'masks\filtermask\default.bmp');
%     filtermasktemp = imread(filtermaskname, 'bmp');
%     if size(filtermasktemp,3)==1
%         filtermask = double(filtermasktemp)/255;
%     else
%         filtermask = double(filtermasktemp(:,:,1))/255;
%     end
%     clear filtermasktemp;
%     if ~(strcmp(LPMethod, 'fastmean') && strcmp(LPMethod, 'gaussian') && strcmp(HPMethod, 'fastmean') && strcmp(HPMethod, 'gaussian'))
%         fprintf('Mask filtering is available only for "fastmean" and "gaussian" currently. The other methods work without mask filtering.\r\r');
%     end
% else
%     filtermask = ones(504,504,'double');
% end

anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(1), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;
DataType=anapar.DataType;


if corrmasktype>0 
    corrmaskfolder=strcat(resultdriver, outputfolder, expname, runname, 'masks\corrmask\') %100409 HT
    switch corrmasktype
    case 1      % for map-loading, Mask bmp files are 0-1 monochrome bmp files, desired areas are 1's
        if isempty(corrmaskname)        
            tempfilename=struct2cell(dir([corrmaskfolder, '*.bmp']));
            corrmaskfilename=sort(tempfilename(1,:)');
            corrmasknum=size(corrmaskfilename, 1);
            if corrmasknum == 0
                fprintf('Error: Need corrmask in corrmask folder!\n');
            end
            for i=1:corrmasknum
                corrmasknametemp(i,:)=getfield(cell2struct(corrmaskfilename(i), 'junk'), 'junk');
            end
            corrmaskname=char(corrmasknametemp(:,1:end-4)); 
        else
            corrmasknum=size(corrmaskname,1);
        end
        for i=1:corrmasknum
            strcat(corrmaskfolder, corrmaskname(i,:));
            masktemp=imread (strcat(corrmaskfolder, corrmaskname(i,:), '.bmp'), 'bmp');
            if size(masktemp,3)==1
                mask(:,:,i) = masktemp;
            else
                mask(:,:,i) = masktemp(:,:,1);
            end
        end
        mask=uint8(mask);   % need check if mask is binary.
    case 2     % masks are from coordinates, create masks in memory here
        if isempty(corrmaskname)
            tempfilename=struct2cell(dir([corrmaskfolder, '*x.txt']));  % only look for x.txt
            corrmaskfilename=sort(tempfilename(1,:)');
            corrmasknum=size(corrmaskfilename, 1);
            dotsinmap=200;       % assume maximum 100 dots in each domain map, only for 'flagmask==2'
            xx=zeros(dotsinmap, corrmasknum);       % store the dot center coordinates          
            yy=zeros(dotsinmap, corrmasknum);
            for i=1:corrmasknum
                corrmasknametemp(i,:)=getfield(cell2struct(corrmaskfilename(i), 'junk'), 'junk');
            end
            corrmaskname=char(corrmasknametemp(:,1:end-6));
        else
            corrmasknum=size(corrmaskname,1);
        end
%         if isempty(domainradius)
%             domainradius=zeros(corrmasknum)+5;
            domainradius=zeros(corrmasknum)+domainradiussize;
%         end                
        for i=1:corrmasknum
            domainfidx=fopen (strcat(resultdriver, outputfolder, expname, runname, 'masks\corrmask\', corrmaskname(i,:), '_x.txt'), 'r'); %Hisashi
            domainfidy=fopen (strcat(resultdriver, outputfolder, expname, runname, 'masks\corrmask\', corrmaskname(i,:), '_y.txt'), 'r'); %Hisashi
            xtemp=fscanf(domainfidx, '%f');
            ytemp=fscanf(domainfidy, '%f');
            domaindotnum(i)=size(xtemp, 1);
            for j=1:domaindotnum(i)
                xx(j, i)=xtemp(j);      % j is dots, i is mask type
                yy(j, i)=ytemp(j);
            end
            fclose(domainfidx);
            fclose(domainfidy);        
            corrmask(:,:,i)=OIMaskMake(0, xx(:,i), yy(:,i), domainradius(i), domaindotnum(i), FrameWidth, FrameHeight); % maskmake(manual, xs, ys, pixradius, dotnum, xdim, ydim)
            if (savemaskmap==1)     % output memory masks for evaluating if selected 
                corrmaskmap=corrmask(:,:,i)*255;
                imwrite(uint8(corrmaskmap), strcat(resultdriver, outputfolder, expname, runname, 'masks\corrmask_', corrmaskname(i,:), '.bmp'), 'bmp');
            end
        end
        clear maskmap;
    end
	corrmaskfilename
	fprintf('\nFound %d masks for superpixel anslysis.\n', corrmasknum);
else
	fprintf('\nError: specify corrmasktype correctly.\n');
end

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
        	for fr=[fframe framerange]
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
% newsize=round(FrameWidth*Binning);  % assume width=height here
newWidth=round(FrameWidth*Binning);
newHeight=round(FrameHeight*Binning);
numSmap=size(Smap,1);   
Range1=zeros(newWidth, newHeight, NStim);  
% Range2=zeros(newsize, newsize, NStim, 'single');  
Rangetemp1=zeros(newWidth, newHeight);
Rangetemp2=zeros(newWidth, newHeight);

if Binning ~=1
    LPKernel=round(LPKernel*Binning);
    switch HPMethod
        case 'ribot'
            LPKernel=LPKernel;
        otherwise
            HPKernel=round(HPKernel*Binning);
    end
end

if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
    fftfilter=OIMakFFTfilterH(newWidth,LPMethod, LPKernel, HPMethod, HPKernel);
end    


ctime0=clock;
flagfirst=0;
% counter1=0;
% counter2=0;
counter5=0;
for j=1:size(Smap,1)
    fprintf('\r\rMap for correlation analysis: %s\r',getfield(cell2struct(Smap(j, 1),'junk'),'junk'));
    counter1=0;
    counter2=0;
    NStimToRead=size(Smap,2); %+size(Smap,3);
    Sample1=zeros(newWidth, newHeight, blockselectnum*NStimToRead);  
%     Sample2=zeros(newsize, newsize, blockselectnum*NStimToRead, 'single');
	for k=blockselect
% 		k
		blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
		if blockselectnum>1
			if counter5>0 && flagfirst==0 % && k==blockselect(2)    % just for estimate process time 
				blockproctime=etime(clock, ctime0);
				fprintf('Time for one block: %7.4f secs', blockproctime);
				temptime=blockproctime*blockselectnum*size(Smap,1);
				fprintf('\rTime for all (%d) blocks: %5.2f minutes (%5.2f hours)\r\r', blockselectnum, temptime/60, temptime/3600);
				flagfirst=1;
			end
		end
		fprintf('%s:', blkfilename);
%         counter1=counter1+1;
		counter5=counter5+1;
		counter0=0;
	%     fid=fopen(strcat(blockfolder, blkfilename),'rb','ieee-le'); % open block file
		condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
		for m=condlist1; %1:NStim
% 			m
%             counter0=counter0+1;
			if goodstim(k, m)
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

				if (LPKernel || HPKernel)
					Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
				end
				
				if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
					Sampletemp=ifft2(fft2(Sampletemp).*fftfilter);
				end
				
				if flagalign == 1
					if Binning ==1
						Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
					else
						Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
					end
				end

				Sample1(:,:,counter1)=Sampletemp;
				
                for i=1:corrmasknum
					if Binning ==1
						corrmask2=corrmask(:,:,i);
					else
						corrmask2=imresize(corrmask(:,:,i), Binning, 'bilinear');
					end
					
					maskedmap2= Sample1(:,:,counter1).*corrmask2;
					maskedvalues(i,counter1)=sum(sum(maskedmap2))/sum(sum(corrmask2));
				end								
			end
		end % for m
% 		counter1
		
%         condlist2=getfield(cell2struct(Smap(j, 3), 'junk'), 'junk');
%         for m=condlist2; %1:NStim
% %             counter0=counter0+1;
%             if goodstim(k, m)
%                 counter2=counter2+1;
%                 if flagrandom
%                     stimloc=find(stimseq(k, :)==m);  % location of stim1(j) in this randomnized trial.
%                 else
%                     stimloc=m;
%                 end
%     %             stimloc2(m)=counter0;
%                 fprintf(' stim%d ', stimloc);
%                 AllFrames=OIReadStim(strcat(blockfolder, blkfilename), stimloc, system);
%                 if Binning ==1
%                     Rangetemp1=mean(AllFrames(:,:,fframe), 3);  
%                     Rangetemp2=mean(AllFrames(:,:,framerange), 3);
%                 else
%                     Rangetemp1=imresize(mean(AllFrames(:,:,fframe), 3), Binning, 'bilinear');  
%                     Rangetemp2=imresize(mean(AllFrames(:,:,framerange), 3), Binning, 'bilinear');
%                 end
% 
%                 Sampletemp=(Rangetemp2-Rangetemp1)./Rangetemp1;
% 
%                 if (LPKernel || HPKernel)
%                     Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
%                 end
%                 
%                 if strcmp(LPMethod, 'fft') || strcmp(HPMethod, 'fft')
%                     Sampletemp=ifft2(fft2(Sampletemp).*fftfilter);
%                 end
% 
%                 if flagalign
%                     if Binning ==1
%                         Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
%                     else
%                         Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
%                     end
%                 end
% 
%                 Sample2(:,:,counter2)=Sampletemp;
%             end
%         end % for m
        fprintf('\r');
    end	% end of "k"
    fprintf('\r');

	for i=1:corrmasknum
        
        
        fprintf('\rSite to analyze: %s', corrmaskname(i,:));
        counter3=0;
        counter4=0;
        
        for x=1:newWidth
            for y=1:newHeight


                counter3=counter3+1;
                temp1=squeeze(Sample1(x,y,1:counter1));
                temp2=squeeze(maskedvalues(i,1:counter1))';
                if TestType==1
                    [RHO,PVAL] = corr(temp1, temp2, 'type', 'Pearson');
                else
                    [RHO,PVAL] = corr(temp1, temp2, 'type', 'Spearman');
                end
                p(x,y) = PVAL;
                r(x,y) = RHO;
                if counter3>=newWidth*newHeight*0.1*(counter4+1)
                    fprintf('\rCorrelation analysis process: %g%% completed', 10*(counter4+1));
                    counter4=counter4+1;
                end

            end
        end
        fprintf('\rCorrelation analysis process: 100%% completed\r');
        
        signmap=double(r>=0)*2-1;

        Templ= p > pmin; % logical operation 0 or 1; locations of lower value
        Tempu= p < pmax;% logical operation 0 or 1; locations of higher value
        Tempul= (Templ.*Tempu).*p;% bewteen low and high clips
        Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
        Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
        p2=Tempul2;
        p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));

        pmap=uint8(round(p2_log.*signmap*127+128));
        mapname=getfield(cell2struct(Smap(j,1), 'junk'), 'junk');

%         if (j<=9)
%             prefix=strcat('0', num2str(j), '_');
%         else
%             prefix=strcat(num2str(j), '_');
%         end

        prefix='';

        sitefolder=[resultdriver, outputfolder, expname, runname, resultname, mapname, '\'];
        if ~isdir(sitefolder)
            mkdir(sitefolder);    
        end
        
        imwrite(pmap, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p-', num2str(pmax,'%1.5f'), '_pmap.bmp'));
        lutstat=textread(statlutfile);  % this color table should be in sunin folder
        imwrite(pmap,lutstat, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p-', num2str(pmax,'%1.5f'), '_pmapcolor.bmp'));

%         [maptemp, framemedian, lowClip, highClip] = OIClipH(p, 1, 1);
%         imwrite(norm_to_uint8(maptemp), strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p-all.bmp'));
        for q=threshold
            b=double(p<q);
            imwrite(uint8(b*255), strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p-', num2str(q,'%1.5f'), '.bmp'));
        end
        
        if rrange==999
            if min(min(r))*-1 < max(max(r))
                rrange=max(max(r));
            else
                rrange=min(min(r))*-1;
            end
        end

        cliptextfid= fopen(strcat(sitefolder, strcat('r_map_clipping_range', '.txt')), 'a');
        fprintf(cliptextfid, 'r_map %s was clipped at 0 +- %g.\r', strcat(corrmaskname(i,:)), rrange);                
        fclose(cliptextfid); 
        
        rhomaptemp = norm_to_uint8b(r,rrange*-1,rrange);
        
        imwrite(rhomaptemp, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_rmap.bmp'));
        lutrho=textread(rlutfile);  % this color table should be in sunin folder
        imwrite(rhomaptemp, lutrho, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)),  '_rmapcolor.bmp'));

        if flagsaveivf
%             OIWriteIVF(p, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p.ivf'));             % p values
            OIWriteIVF(r, strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_rmap.ivf'));             % r values
        end
       
			
	end %for i=1:corrmasknum
    if flagsaveivf
        samplenumfid= fopen(strcat(sitefolder, prefix1, prefix, 'sample_number.txt'), 'a');
        fprintf(samplenumfid, '%g\r', counter1);
        fclose(samplenumfid); 
    end
end %for n=1:size(Smap,1)  

% output a color table
colorbarmap=zeros(60, 256, 3);
for i=1:30
    colorbarmap(i, :, :)=256.*lutstat(:,:);
end
imwrite(uint8(colorbarmap), lutstat, strcat(resultfolder, 'statcolortable.tif'), 'tiff'); 

for i=1:30
    colorbarmap(i, :, :)=256.*lutrho(:,:);
end
imwrite(uint8(colorbarmap), lutrho, strcat(resultfolder, 'rcolortable.tif'), 'tiff'); 

totaltime=etime(clock, ctime0);
fprintf('\rProcess Finished: ');
fprintf('\n\n\nTotal time used: %8.4f secs (%8.4f hours).\r\r', totaltime, totaltime/3600);

return;

