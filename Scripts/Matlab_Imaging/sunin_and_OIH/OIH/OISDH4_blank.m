function OISDH4()
% OISDH: Visualize the SD of signals pixel by pixel. Hisashi

%% Ver 4.1.1 (150224) Frame-by-frame processes were improved. Hisashi
%% Ver 4.0 (150221) dalsalineremove=1 were introduced. flagalign=2 were introduced. Hisashi
%% Ver 3.3 (121118) 'outputimageformat' was introduced. Hisashi
%% Ver 3.2 (091227) 'flagspfilterforSD' was introduced. Hisashi
%% Ver 3.1 (090817) Spatial filtering is applied before calculating SD. Hisashi

clear all;
% __________________ Start User Input ____________________
system='v';     % 'v' or 'r'
datadriver = 'K:\';     % Data disk name
datafolder = 'expt0\';   % Data folder name on data disk, results will be saved in outputfolder
outputdriver = 'D:\Vandy\';     % Data disk name
outputfolder = 'expt1\'; % Output folder name on data disk, results will be saved here
expname = '110328VitRG_Atn\'; % Exp folder name (in both data folder and result folder)
runname = 'run99\';      % Run foler name (in both data folder and result folder)
resultfolder='OISD\OISDH4_blank\';        % specify a founder name results to be saved. if not specified, program will generate one like "H2L2B0.5', 
mapname='blank'
outputimageformat  = 'png'; % The format for subtraction maps. % Choose from these format: 'bmp', 'gif', 'jpg', 'pcx', 'tif', etc. Since H33

filename = {            % blk file names, leave empty for automatic searching
};
%List of stim groups for calucalation
stim(1,:)=[1 ]; %Stim ID for calculation of SD


fframe=[1:2];
framerange=[5:16];
blockselect =[];  % select blocks for processing, eg. [ 1 3 4 5], leave empty for all block being selected ([]).

%threshold=[0.00001 0.0001 0.001 0.01 0.05];	% threshold p values, just for generating thresholded p maps, each value will have a thresholded map (i.e. 0.05 will generate a p<0.05 map)
% SDthresholdpecentile=0.95; % determine percentile which is used for thresholding, usually 0.95
% SDthreshold=0.000171;
% precisiontype='double';

%clipsd=1;				% for map output, usually 1
Binning=1; %0.5;             % note: it is the x times of original size (i.e. 0.5 is equal to 2x2 binning)
LPMethod='fastmean';    % Low-pass filtering method, usually 'gaussian'
LPKernel=0; %10;             % set as 0 if no LP filter
HPMethod='fastmean';       % High-pass filtering method, 'ribot' filter is fastest
HPKernel=0;             % For 'ribot' filter, it's the order of polynomial fit, usually set as 2, set as 0 if no HP filter
        % choose from 'fastmean', 'slowmean', 'gaussian', 'ribot'
        % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
        % 'slowmean' 'disk-like' but better at edge
        % 'gaussian' Gaussian filter with half sd
        % 'fastmedian' Median filter 
        % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
                           
flaggoodstim = 1; % 1: will looking for the file named by gsfilename parameter in block folder. The file contains '1's for good stim conditions, and used these for average map & quantification
gsfilename='goodstim.txt';
flagalign=0;        % 0: no shift/align; 1: shift/align using the first frame on trial-by-trial. need 'shiftinput.2.txt'; 2: shift/align using all frames on trial-by-trial. need 'shiftinput.4.txt'
flagrandom=0;       % default=0; similar to the 'flagrandom' in sunin, if 'flagrandom=1', will search for 'stimseq.txt' in block folder
% flaggroupaveraging=1; %averaging data in each stim group. Default= 1
flagspfilterforSD=1;% default=1; 1: apply spatial filters to SD map
flagsaveivf=0;      % default=0; 1: save ivf format inaddition to bmps (only for average maps), ivf is folat data type, can only be viewed by WinMix, also will be useful for subsquent data analysis (no need to read from source block files)
flagSDmask=0;       % default=0; 0: Only center portion (250X250) will be used for calculation % the portion for calculation will be determined by a mask file ('expresult\masks\SDmask\default.bmp')
flaghist=0;         % default=0; 1: make a histogram
dalsalineremove=0;	% 1: remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels

%
%
% if run into memory problem, try bin more or reduce the number of blocks
% for calculation.
%


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
if isempty(blockselect)
    blockselectnum=size(filename, 1); 
    blockselect=1:blockselectnum;
else
    blockselectnum=size(blockselect,2);
end
blockfilenum=size(filename, 1);


anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(filename(5), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;


if flaggoodstim
    goodstim=textread(strcat(runfolder, gsfilename), '%d');
    if isempty(goodstim)|size(goodstim,1)~=blockfilenum*NStim
        fprintf('Error, "goodstim.txt" does not contain right number of conditions\r');
    end
	if sum(sum(goodstim))~=sum(sum(goodstim.*goodstim))	% only 0*0=0 and 1*1=1
		fprintf('Error: ''goodstim.txt'' should contain only "1" and "0"\n');
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

if (flagSDmask==1)      %use masking (to exclude blood vessels or only include certain area)
    maskfilename = strcat(resultdriver, outputfolder, expname, 'masks\SDmask\default.bmp');     % to determine which pixel is used for SD calculation.
    SDmask = imread (maskfilename, 'bmp') / 255;
%     maskpixelnum=sum(sum(SDmask));
else 
    SDmask = ones(FrameWidth, FrameHeight);
    SDmask(1:50,:)=0; %Only center portion will be used for calculation
    SDmask(454:504,:)=0;
    SDmask(:,1:50)=0;
    SDmask(:,454:504)=0;
%     maskpixelnum=sum(sum(SDmask));
end

if Binning ==1
    SDmask2=SDmask;
else
    LPKernel=round(LPKernel*Binning);
    switch HPMethod
        case 'ribot'
        otherwise
            HPKernel=round(HPKernel*Binning);
    end
    SDmask2=imresize(SDmask, Binning, 'nearest');
end

maskpixelnum=sum(sum(SDmask2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading frames %%%%%%%%%%%%%%%%%%%%%%%
%ngroup=size(stim,1);
nstimselect=size(stim,2)
newWidth=round(FrameWidth*Binning);
newHeight=round(FrameHeight*Binning);
% sample1=zeros(newWidth, newHeight, blockselectnum*nstimselect);
% sample2=zeros(newWidth, newHeight, blockselectnum*nstimselect);
sample=zeros(newWidth, newHeight, blockselectnum);
% samplegs=ones(nstimselect, blockselectnum);

blkcount=0;
ctime0=clock;
for k=blockselect
%     blkcount=blkcount+1;
    blkfilename=getfield(cell2struct(filename(k), 'junk'), 'junk');
    fprintf('%s\r', blkfilename);
    if k==blockselect(2)    % just for estimate process time
        blockproctime=etime(clock, ctime0);
        fprintf('\rTime for one block: %7.4f secs', blockproctime);
        fprintf('\rTime for all (%d) blocks: %5.2f hours\r', blockselectnum, (blockproctime*blockselectnum+190)/3600);
    end
    samplegs=0;
    flagfirst=1;
    for j=1:nstimselect
        if flagrandom
            stimloc=find(stimseq(k, :)==stim(1,j));  % location of stim(o,j) in this randomnized trial.
        else
            stimloc=stim(1,j);
        end
%         goodstim(k, stimloc)
        if goodstim(k, stimloc)
            if flagfirst==1
                blkcount=blkcount+1;
                flagfirst=0;
            end
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
            
%             Range1=mean(Stimtemp(:,:,fframe), 3);  
%             Range2=mean(Stimtemp(:,:,framerange), 3);
%             sample(:,:,(nstimselect*(blkcount-1)+j))=(Range2-Range1)./Range1;
            Sampletemp=((Rangetemp2-Rangetemp1)./Rangetemp1);

            if flagspfilterforSD==0
                Sampletemp=OIEasyFilterH2(Sampletemp, LPMethod, LPKernel, HPMethod, HPKernel);
            end

            if flagalign == 1
                if Binning ==1
                    Sampletemp=OIShift(Sampletemp, sfx1((k-1)*NStim+stimloc), sfy1((k-1)*NStim+stimloc));
                else
                    Sampletemp=OIShift(Sampletemp, (sfx1((k-1)*NStim+stimloc)*Binning), (sfy1((k-1)*NStim+stimloc)*Binning));
                end
            end
                
            sample(:,:, blkcount)=sample(:,:, blkcount)+Sampletemp;

%             fprintf('%g, %g, Act\r', j, k);
%             median(median(sample(:,:, j, blkcount)))
            samplegs=samplegs+1;
%         else
% %             sample(:,:,(nstimselect*(blkcount-1)+j))=NaN;
%             sample(:,:, j, blkcount)=NaN;
%             samplegs(j, blkcount)=0;
% %             fprintf('%g, %g, NaN\r', j, k);
        end
    end % end of "j"
    if flagfirst==0
        sample(:,:, blkcount)=sample(:,:, blkcount)/samplegs;
    end
end	% end of "k"
fprintf('\r');


SDmaptemp=zeros(newWidth, newHeight);

counter3=0;
counter4=0;
% for j=1:nstimselect
%     temp2=squeeze(samplegs(j,:))';
    for x=1:newWidth
        for y=1:newHeight
            counter3=counter3+1;
%             temp1=squeeze(sample(x,y,:));
% %             size(temp2);
% %             size(temp1);
%             combinetemp=[temp1, temp2];
%             combinetemp=flipud(sortrows(combinetemp,2));
% %             sum(combinetemp(:,2));
%             combinetemp2=squeeze(combinetemp(1:sum(combinetemp(:,2)),1));
%             sdtemp = std(combinetemp2);
            sdtemp = std(sample(x,y,:));
%             sdtemp = var(temp1);
            SDmaptemp(x,y) = sdtemp;
            if counter3>=newWidth*newHeight*0.1*(counter4+1)
                fprintf('\rOISDH analysis: %g%% completed', 10*(counter4+1));
                counter4=counter4+1;
            end
        end
    end
% end

% SDmap=mean(SDmaptemp,3);
SDmap=SDmaptemp;

if flagspfilterforSD==1
    SDmap=OIEasyFilterH2(SDmap, LPMethod, LPKernel, HPMethod, HPKernel);
end

clear sample;

%% Output
% A=reshape(SDmap, prod(size(SDmap)), 1);
% B=reshape(double(SDmask), prod(size(SDmap)), 1);
combined=[reshape(SDmap, prod(size(SDmap)), 1), reshape(double(SDmask2), prod(size(SDmap)), 1)];    %following 3 lines for pick up useful pixels
combined=flipud(sortrows(combined,2));
usefulpixels=squeeze(combined(1:maskpixelnum, 1));
usefulpixels_sort=sortrows(usefulpixels,1);

SDmax=max(usefulpixels);
SDmin=min(usefulpixels);
SDmedian=median(usefulpixels);
SDmean=mean(usefulpixels);
SDSD=std(usefulpixels);
SD5thPercentile=usefulpixels_sort(round(maskpixelnum*0.05));
SD95thPercentile=usefulpixels_sort(round(maskpixelnum*0.95));
SDCI95lower=sqrt((blockfilenum-1)*SDmean*SDmean/chi2inv(0.975, blockfilenum-1));
SDCI95upper=sqrt((blockfilenum-1)*SDmean*SDmean/chi2inv(0.025, blockfilenum-1));
SDCI99lower=sqrt((blockfilenum-1)*SDmean*SDmean/chi2inv(0.995, blockfilenum-1));
SDCI99upper=sqrt((blockfilenum-1)*SDmean*SDmean/chi2inv(0.005, blockfilenum-1));

fprintf('\r\rmax of SD: %1.6f\r', SDmax);
fprintf('min of SD: %1.6f\r', SDmin);
fprintf('mean of SD: %1.6f\r', SDmean);
fprintf('SD of SD: %1.6f\r', SDSD);
fprintf('median of SD: %1.6f\r', SDmedian);
fprintf('the 5-th percentile of SD: %1.6f\r', SD5thPercentile);
fprintf('the 95-th percentile of SD: %1.6f\r', SD95thPercentile);
fprintf('95%%-confidence interval lower limit: %1.6f\r', SDCI95lower);
fprintf('95%%-confidence interval upper limit: %1.6f\r', SDCI95upper);
fprintf('99%%-confidence interval lower limit: %1.6f\r', SDCI99lower);
fprintf('99%%-confidence interval upper limit: %1.6f\r', SDCI99upper);


SDfid= fopen(strcat(resultfolder, mapname, strcat('_SD_values', '.txt')), 'w');

fprintf(SDfid, 'max of SD: %1.6f\r', SDmax);
fprintf(SDfid, 'min of SD: %1.6f\r', SDmin);
fprintf(SDfid, 'mean of SD: %1.6f\r', SDmean);
fprintf(SDfid, 'SD of SD: %1.6f\r', SDSD);
fprintf(SDfid, 'median of SD: %1.6f\r', SDmedian);
fprintf(SDfid, 'the 5-th percentile of SD: %1.6f\r', SD5thPercentile);
fprintf(SDfid, 'the 95-th percentile of SD: %1.6f\r', SD95thPercentile);
fprintf(SDfid, '95%%-confidence interval lower limit: %1.6f\r', SDCI95lower);
fprintf(SDfid, '95%%-confidence interval upper limit: %1.6f\r', SDCI95upper);
fprintf(SDfid, '99%%-confidence interval lower limit: %1.6f\r', SDCI99lower);
fprintf(SDfid, '99%%-confidence interval upper limit: %1.6f\r', SDCI99upper);
fclose(SDfid);

SDfid2= fopen(strcat(resultfolder, mapname, strcat('_SD_values_raw', '.txt')), 'w');
for i=1:maskpixelnum
    fprintf(SDfid2, '%d\r', usefulpixels(i));
end
fclose(SDfid2);

SDthreshold=SDCI95upper;


%% flaghist
if flaghist
    % hist(usefulpixels,0:0.0000005:0.00016)
    % [n,xout] = hist(usefulpixels,0:0.0000005:0.00016);
    histdiv=0.00001;
%     histmax=0.003;
    histmax=floor((SDmedian+SDSD*2)/histdiv)*histdiv
    % hist(usefulpixels,0:histdiv:histmax)
    [n,xout] = hist(usefulpixels,0:histdiv:histmax);
    
    n=n(1:int16(histmax/histdiv));
    xout=xout(1:int16(histmax/histdiv));
%     bar(xout,n)

    modevalue=xout(find(n == max(n)));
    
%     temp_threshold = modevalue+(modevalue-SD5thPercentile)-histdiv/2;
    temp_threshold = SDthreshold;
    
    fprintf('Mode of SD= %g \rThreshold= %g\r', modevalue, temp_threshold);
    
    SDfid= fopen(strcat(resultfolder, mapname, strcat('_SD_values', '.txt')), 'a');
    fprintf(SDfid, 'Mode of SD= %g \rThreshold= %g\r', modevalue, temp_threshold);
    fclose(SDfid);
    
%     SDthreshold=temp_threshold;
    
    
%     plot(xout,n,'o')
    bar(xout,n)
    hold on
    hold off
    
end

%% flagfithist
flagfithist=0;
if flagfithist
    % hist(usefulpixels,0:0.0000005:0.00016)
    % [n,xout] = hist(usefulpixels,0:0.0000005:0.00016);
    histdiv=0.00002;
    histmax=0.003;
    % hist(usefulpixels,0:histdiv:histmax)
    [n,xout] = hist(usefulpixels,0:histdiv:histmax);
    n=n(1:histmax/histdiv);
    xout=xout(1:histmax/histdiv);
    bar(xout,n)

    maxvalue=xout(find(n == max(n)));
    
    temp_threshold = xout(find(n == max(n))+4);
    fprintf('maxvalue= %g \rtemp_threshold= %g\r', maxvalue, temp_threshold);
    outliers = excludedata(xout,n,'domain',[temp_threshold Inf]);

    opt = fitoptions('Method', 'NonlinearLeastSquares','Lower', [-Inf -Inf 0], 'Exclude', ~outliers, 'Normalize', 'on');
    %opt = fitoptions('Method', 'NonlinearLeastSquares','Lower', [-Inf -Inf 0], 'Normalize', 'on');
    [cfun,gof]=fit(xout(:), n(:), 'gauss1', opt);

    [cprd,ypred] = predint(cfun,0:histdiv:histmax,0.95, 'obs', 'off');
    % [cprd2,ypred2] = predint(cfun,0:histdiv:histmax,0.95, 'obs', 'off')

    gof

    plot(xout,n,'o')
    hold on
    plot(cfun,'m')
    hold on
    plot(0:histdiv:histmax, cprd, 'k-.')
    % hold on
    % plot(0:histdiv:histmax, cprd2, 'r:.')
    hold off

    check=find(n == max(n))+2;
    while n(check) < cprd(check, 2)
        fprintf('%g %g %g\r', xout(check), n(check), cprd(check, 2));
        check=check+1
    end

%     SDthreshold=xout(check-1);
    fprintf('SDthreshold= %g \r', SDthreshold);


    save(strcat(resultfolder, mapname, strcat('_SD_hist', '.mat')), 'n', 'xout')

    goffid= fopen(strcat(resultfolder, mapname, strcat('_SD_hist_gof', '.txt')), 'w');
    fprintf(goffid, 'sse: %6.4g', gof.sse);
    fprintf(goffid, 'r-square: %6.4g', gof.rsquare);
    fprintf(goffid, 'dfe: %6.4g', gof.dfe);
    fprintf(goffid, 'adjrsquare: %6.4g', gof.adjrsquare);
    fprintf(goffid, 'rmse: %6.4g', gof.rmse);
    fclose(goffid);

    %SDthreshold=usefulpixels_sort(round(maskpixelnum*SDthresholdpecentile));
end

%%
TempU= SDmap < SDthreshold; %logical operation 0 or 1; locations of higher value
SDmapUTemp = TempU.*SDmap;% bewteen low and high clips
SDmapU=SDmapUTemp + (SDthreshold*(~TempU));%clip all avules higher than highClip

min(min(SDmap))
min(min(SDmapU))
SDmapBMP=(SDmap-min(min(SDmap)))/(max(max(SDmap))-min(min(SDmap)));
SDmapUBMP=(SDmapU-min(min(SDmapU)))/(max(max(SDmapU))-min(min(SDmapU)));

imwrite(uint8(SDmapBMP*255), strcat(resultfolder, mapname, '_SD_map.', outputimageformat));
imwrite(uint8(SDmapUBMP*255), strcat(resultfolder, mapname, '_SD_map_clipped.', outputimageformat));
imwrite(uint8(TempU*255), strcat(resultfolder, mapname, '_SD_map_threshold.', outputimageformat));

% SDLog=OIEasyFilterH2(SDmap, 'log', LPKernel, HPMethod, 0);
% imwrite(norm_to_uint8(SDLog), [resultfolder, mapname, '_SD_map_log.bmp']);
% 
% imwrite(uint8(medfilt2(TempU, [5 5])*255), [resultfolder, mapname, '_SD_map_threshold_median05.bmp']);
% imwrite(uint8(medfilt2(TempU, [10 10])*255), [resultfolder, mapname, '_SD_map_threshold_median10.bmp']);
% imwrite(uint8(medfilt2(TempU, [15 15])*255), [resultfolder, mapname, '_SD_map_threshold_median15.bmp']);
% imwrite(uint8(medfilt2(TempU, [20 20])*255), [resultfolder, mapname, '_SD_map_threshold_median20.bmp']);
% imwrite(uint8(medfilt2(TempU, [25 25])*255), [resultfolder, mapname, '_SD_map_threshold_median25.bmp']);
% imwrite(uint8(medfilt2(TempU, [30 30])*255), [resultfolder, mapname, '_SD_map_threshold_median30.bmp']);

if flagsaveivf
	OIWriteIVF(SDmap, [resultfolder, mapname, '_SD_map.ivf']);             % p values
end


return;


