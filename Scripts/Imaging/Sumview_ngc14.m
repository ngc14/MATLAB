function SumviewRMF()

% A Matlab program does simple analysis on optical imaging data
% Basic function: Generate a time-sequence frames averaged over trials.
% m files called;  OIHeadRead.m, OIReadStim.m, OIClipH2.m, norm_to_uint8.m
% See readme for details.			070420

clear all; close all;
%__________________ User Input Area ________________________

system='v';  % 'r' for redshirt, 'v' for VDAQ


datafolder = '';   % Data folder name on data disk, results will be saved in outputfolder


% blockfolder='f:\expt\060818Jem\run0\';  % where '.blk' or '.da' file are
% located. variable specified later in code.

animal = 'Gilligan';
date = '12_13_2018';
run = 'run00';


datadriver = ['S:\Lab\',animal,'\All Data\'];     % Data disk name
outputdriver = ['S:\Lab\',animal,'\All Data\'];     % Data disk name

%outputfolder = 'Results\'; % Output folder name on data disk, results will be saved here
expname = [animal,'_',date,'\Imaging\']; % Exp folder name (in both data folder and result folder)
runname = [run,'\'];      % Run folder name (in both data folder and result folder)
outputfolder = [expname,runname,'Results\'];
%resultfolder='SumView\LP2_HP300\FFrame1\';        % specifies a FILE where results to be saved. 


filename={}; % blk file names, leave empty for automatic searching


outfile='Rest_.bmp';                  % part of output file name
% % blockfolder='C:\expt\Exp061207\run1\';   % where '.blk' or '.da' files
% located, program will auto search based on these default file extension.

style=2;    % 1: simple DC, 2:ffsub (need specify fframe), 3 subtraction (stim1 -stim2), no ffsub, 4: subtraction (stim1-stim2) with ffsub
fframe=[1];   % specify first frame (or reference frames), used in 'style=2' or 'style=4'
stim1=[4];  % specify stimulus list 1, used in all styles, if more than one stim, will take average
stim2=[0];  % specify stimulus list 2, used in style=3 and style=4, if more than one stim, will take average
flagsingle=1;   % whether or not to save single-block maps (1=save, 0=no save)
flagsummap=1;   % whether or not to save the average of the block maps (1=save, 0=no save)

stims = containers.Map([0 2 4 5],{'Rest','ExtraSmallSphere','LargeSphere', 'Photocell'});
filename={ %if empty will find all data files in RUN folder
%'file1.da'
%'file2.da'
};

% Filtering
    % choose filtering method from the followings
    % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
    % 'slowmean' 'disk-like' but better at edge
    % 'gaussian' Gaussian filter with half sd
    % 'fastmedian' median filter
    % 'fastmedianH' faster median filter developed by Hisashi; shrink the image to the half size first, 
    % apply median filter, then 
    % expand the image to the original size.
    % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM) kernel is order of polynomial (usually 2 or 3)
    % 'fft', band-pass filtering using FFT (Fast Fourier Transform).
    
LPMethod='gaussian';    % Low-pass filtering method, usually 'gaussian'
LPKernel=5; %15;             % set as 0 if no LP filter
HPMethod='fastmedian';       % High-pass filtering method, 'ribot' filter is fastest
HPKernel=550; %100;             % For 'ribot' filter, it's the order of polynomial fit, usually set as 2, 
%set as 0 if no HP filter
flagtrialfilter=0;	% 0: no trial by trial filtering; 1: filtering every frame in every block (slow)
flagsumfilter=0;    % whether or not to filter average block maps (1=save, 0=no save). 

clipMethod = 6;
sd=.3;       % number of standard deviations in map clipping

% following 3 lines for layout of the output bmp map
colnum=5;     % number of frames in a row, number of rows will be calculated
spacepixx=3; % space between two frames in one row in pixel
spacepixy=3;  % space between two frames in one column in pixel

BNC=0;     % BNC=1 will save voltage value at 8 BNC port (in BNCs.txt), BNC=0: no save, note for redshirt data only
% not tested

if style==3 || style==4
    stimName = ['[',stims(stim1),' - ',stims(stim2),']'];
else
    stimName = ['[',stims(stim1),']'];
end

if flagtrialfilter || flagsumfilter
    resultfolder= ['SumView\',stimName,'\LP',num2str(LPKernel),'HP',num2str(HPKernel),'C',num2str(sd)]; 
else
    resultfolder= ['SumView\',stimName,'\LP0HP0C',num2str(sd)]; 
end

if style==2 || style==4
    resultfolder = [resultfolder,'\FFrame',num2str(fframe)]; 
else
    resultfolder = [resultfolder,'\NoFF']; 
end

%___________________ End of User Input (do not edit below)_______________

hpfmethod=HPMethod;	% haidong's = hisashi's 
hpkernelsize=HPKernel;     % haidong's = hisashi's 
lpfmethod=LPMethod; % haidong's = hisashi's 
lpkernelsize=LPKernel;		% haidong's = hisashi's 

%_____Folder Specification/Generation______

if isempty(resultfolder)
    resultfolder=strcat('HP',num2str(HPKernel), '_LP', num2str(LPKernel), '_Bin', num2str(Binning), '/');
    elseif 	resultfolder(end)~='/'
        resultfolder=[resultfolder, '/'];	
end

blockfolder=strcat(datadriver, datafolder, expname, runname);
runfolder=strcat(outputdriver, outputfolder, expname, runname);
resultfolder=strcat(outputdriver, outputfolder, resultfolder);
if ~isdir(resultfolder)
    mkdir(resultfolder);    
end
cd(blockfolder)
if isempty(dir('*.BLK'))
    blockfolder = ['D:\Data\',animal,'_SqM\',hemi,'Hemisphere\',date,'\',run,'\'];
end


%___________ Locating Block Data Files ___________

if isempty(filename)
	if system=='v'
	    tempfilename=struct2cell(dir([blockfolder, '*.BLK']));
	    BNC=0;
	elseif system=='r'
	    fprintf('Note: you may need delete non-block "*.da" files in data folder\n');
	    tempfilename=struct2cell(dir([blockfolder, '*.da']));
	end
	filename=sort(tempfilename(1,:)');
	for i=1:size(filename,1)
	    fprintf('''%s''\n', getfield(cell2struct(filename(i), 'junk'), 'junk'));
	end
	fprintf('\nFound %d blk files(sorted, check the sequence).\n', size(filename,1));
end

blockfilenum=size(filename, 1);      % number of blocks (files)

%___________Acquiring Data File Information______
blocks = dir('*.BLK');
[~,indx] = natsort({blocks.name});
blocks = blocks(indx);

anapar = OIHeadRead((blocks(1).name), 'v');
Width = anapar.FrameWidth;
Height = anapar.FrameHeight;             % anapar finds the information about the block files
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
NBlocks = length(blocks);
conds = [];
for n =1:size(blocks,1)
    anapar = OIHeadRead((blocks(n).name), 'v');
    conds = [conds, anapar.Cond];
end
blocks = blocks(conds==stim1(1));
numConds = unique(conds);
%blocks = blocks(find(conds==stim1));

blocks2 = dir('*.BLK');
[~,indx] = natsort({blocks2.name});
blocks2 = blocks2(indx);
%blocks2 = blocks2(find(conds==stim2));

rownum=ceil(NFrames/colnum);

if clipMethod~=2
    clipMask = logical(ones(Height,Width));
else
    clipMask = imread(['\\univ.pitt.edu\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_',date, '\Imaging\',run,'\clipMask',run(end-1:end),'.bmp']);
    clipMask = clipMask(:,:,1)>0;
end

%___________Initializing Map Files______


singlemaps=zeros(Height, Width, NFrames); % hold single block maps
maps=zeros(Height, Width, NFrames); % hold final single frames
mapsAll = zeros(Height, Width, NFrames, length(blocks));

singlemapsf=zeros(Height, Width, NFrames); % hold single block maps
mapsf=zeros(Height, Width, NFrames); % hold final single frames
mapsfs=zeros(Height, Width, NFrames); % hold final single frames

bigframe=uint8(zeros((Height+spacepixy)*rownum, (Width+spacepixx)*colnum)); %final output
bigframef=uint8(zeros((Height+spacepixy)*rownum, (Width+spacepixx)*colnum)); %final output
bigframefs=uint8(zeros((Height+spacepixy)*rownum, (Width+spacepixx)*colnum)); %final output

%___________start reading blocks/frames loop____________

for k=1:length(blocks)
    blkfilename=blocks(k).name;
    outfile = blkfilename;
    blkfilename2 = blocks2(k).name;
    
    fprintf('%s\r', blkfilename);
    Frames=zeros(Height, Width, NFrames, size(stim1,2));    % for stim set 1
    Frames2=zeros(Height, Width, NFrames, size(stim2,2));   % for stim set 2
    for j=1:size(stim1,2)   % read file
        Frames(:,:,:,j)=OIReadStim([blocks(k).folder,'\', blkfilename], 1,system);
    end
	if style==3||style==4

        for j=1:size(stim2,2)
            Frames2(:,:,:,j)=OIReadStim([blocks2(k).folder,'\', blkfilename2], 1,system);
        end
    end
    
    switch style
    case 1      % simple DC image
        singlemaps=mean(Frames, 4);
        maps=maps+singlemaps;
        mapsAll(:,:,:,k) = singlemaps;
        
    case 2      % ffsub image
        for j=1:size(stim1,2)
            fframeavg=mean(Frames(:,:,fframe,j), 3);
            for i=1:NFrames
                Frames(:,:,i,j)=(Frames(:,:,i,j)-fframeavg)./fframeavg;
            end
        end
        singlemaps=mean(Frames, 4);
        maps=maps+singlemaps;
        mapsAll(:,:,:,k) = singlemaps;
        
    case 3      % difference image, no ffsub
        singlemaps=mean(Frames,4)-mean(Frames2,4);
        maps=maps+singlemaps;
        mapsAll(:,:,:,k) = singlemaps;
        
    case 4		% difference map, with ffsub
        for j=1:size(stim1,2)
            fframeavg=mean(Frames(:,:,fframe,j), 3);
            for i=1:NFrames
                Frames(:,:,i,j)=(Frames(:,:,i,j)-fframeavg)./fframeavg;
            end
        end
        for j=1:size(stim2,2)
            fframeavg=mean(Frames2(:,:,fframe,j), 3);
            for i=1:NFrames
                Frames2(:,:,i,j)=(Frames2(:,:,i,j)-fframeavg)./fframeavg;
            end
        end
        singlemaps=mean(Frames,4)-mean(Frames2,4);
        maps=maps+singlemaps;
        mapsAll(:,:,:,k) = singlemaps;
        
        % RMF Add filter option flagtrialfilter;	
        if flagtrialfilter==1                
            fprintf('\rFiltering: \rframe');
            fprintf('%d ', i);
            singlemapsf=OIEasyFilter(singlemaps, LPMethod, LPKernel, HPMethod, HPKernel);
            mapsf=mapsf+singlemapsf; 
        end
        fprintf('\r');
    end    
    
    if flagsingle==1	% save map for each file
        
        %RMF move BNC before flagsingle
        if BNC==1
            % output BNC value
            fid=fopen(strcat(blockfolder, blkfilename),'rb','ieee-le'); % open block file
            for stim=1:NConds    
                fseek(fid,5120,-1); % skip 5120 byte header for RedShirt system
                fseek(fid, 2*(Width*Height*NFrames+8*NFrames+Width*Height+8)*(stim-1), 'cof');
                fseek(fid, 2*Width*Height*NFrames,'cof'); 
                tBNCFrames=fread(fid,[8*NFrames], 'uint16');   %NBNC=8 Read 8 BNCs
                tBNCFrames=reshape(tBNCFrames, NFrames, 8); 
                BNCs(:,:,stim,k)=permute(tBNCFrames, [2,1]);    % so the dimension order is (height, width, frame)
            end
            fidBNC=fopen(strcat(resultfolder, 'BNCs-', num2str(k), '.txt'), 'w');
            fprintf(fidBNC, 'This file contains voltage timecourse of 8 BNCs\r\n\t');
            fprintf(fidBNC, '%f\t', 1:NFrames);
            fprintf(fidBNC, '\r\n');
            for stim=1:NConds
                fprintf(fidBNC, 'Condition#%d \r\n', stim);
                for i=1:8
                    fprintf(fidBNC, 'BNC#%d\t', i);
                    fprintf(fidBNC, '%f \t', BNCs(i,:,stim,k)./10000);
                    fprintf(fidBNC, '\r\n');
                end
            end
            fclose(fidBNC);
        end
        
        % to plot frames on single-block big map
		framenum=0;
		for i=1:rownum %col
            for j=1:colnum  %row
                x=(Width+spacepixx)*(j-1)+1;    %x/y location from top left cornor
                y=(Height+spacepixy)*(i-1)+1;
                framenum=framenum+1;
                if framenum<=NFrames
                    
                    bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemaps(:,:,framenum),clipMethod,sd,clipMask));
                    
                    %{
                    switch style
                    case 1
                        bigframe(y:y+Height-1, x:x+Width-1)=singlemaps(:,:,framenum);
                    case 2
                        if size(fframe,2)==1 && framenum==fframe
                        else
                            bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemaps(:,:,framenum),1,sd));
                        end
                    case 3
                        bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemaps(:,:,framenum),1,sd));
                    case 4
                        bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemaps(:,:,framenum),1,sd));
                    end
                    %}
                    
                   
                    if flagtrialfilter==1
                        switch style
                        case 1
                            bigframef(y:y+Height-1, x:x+Width-1)=singlemapsf(:,:,framenum);
                        case 2
                            if size(fframe,2)==1 && framenum==fframe
                            else
                             bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemapsf(:,:,framenum),1,sd));
                            end
                        case 3
                            bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemapsf(:,:,framenum),1,sd));
                        case 4
                            bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(singlemapsf(:,:,framenum),1,sd));
                        end
                    end
                    
                    
                end
            end
		end
        blkoutfile=strcat(resultfolder,  num2str(k),' - ',outfile(1:end-4), '.bmp');
		imwrite(imresize(bigframe,0.25), blkoutfile , 'bmp');
        
        if flagtrialfilter==1
            blkoutfile=strcat(resultfolder, outfile(1:end-4), 'f-', num2str(k), '.bmp');
            imwrite(bigframef, blkoutfile , 'bmp');
        end
        
    end
end

% to plot average frames on the big map

if flagsummap==1
    
bigframe=uint8(zeros((Height+spacepixy)*rownum, (Width+spacepixx)*colnum)); %final output
bigframeSTD=zeros((Height+spacepixy)*rownum, (Width+spacepixx)*colnum); %final output

maps = maps./length(blocks);
mapsSTD = std(mapsAll,0,4);
if flagtrialfilter==1
    mapsf=mapsf./length(blocks);   
end 

mapsZ = zeros(Height,Width,NFrames,length(blocks));
for f = 1:NFrames
    for n = 1:length(blocks)
        
        % abs of z score of each block for each frame
        tempMapsZ = abs((mapsAll(:,:,f,n) - maps(:,:,f)) ./ mapsSTD(:,:,f));
        tempMapsZ(~clipMask) = NaN;
        mapsZ(:,:,f,n) = tempMapsZ;
        
        % avg abs z score for each image
        mapsZmean(f,n) = nanmean(reshape(mapsZ(:,:,f,n),[Height*Width 1]));
    end
end
blockAvgZ = nanmean(mapsZmean,1);

fig1 = figure;
plot(blockAvgZ,'o-'); grid on;
xlabel('Block #');
ylabel('Avg Abs[Z-Score] for Each Block');

framenum=0;
for i=1:rownum %col
    for j=1:colnum  %row
        x=(Width+spacepixx)*(j-1)+1;    %x/y location from top left cornor
        y=(Height+spacepixy)*(i-1)+1;
        framenum=framenum+1;
        if framenum<=NFrames
            
            bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH(maps(:,:,framenum),clipMethod,sd,clipMask));
            
            mapsSTD_temp = mapsSTD(:,:,framenum);
            mapsSTD_temp(~clipMask) = 0;
            bigframeSTD(y:y+Height-1, x:x+Width-1) = mapsSTD_temp;
            
            %{
            switch style
            case 1
                bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(maps(:,:,framenum),1,sd));
            case 2
                if size(fframe,2)==1 && framenum==fframe
                else
                    bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(maps(:,:,framenum),1,sd));
                end
            case 3
                bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(maps(:,:,framenum),1,sd));
            case 4
                bigframe(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(maps(:,:,framenum),1,sd));
                bigframeSTD(y:y+Height-1, x:x+Width-1)=norm_to_uint8(mapsSTD(:,:,framenum));
            end
            %}
            
            if flagtrialfilter==1
                switch style
                case 1
                    bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsf(:,:,framenum),1,sd));
                case 2
                    if size(fframe,2)==1 && framenum==fframe
                    else
                        bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsf(:,:,framenum),1,sd));
                    end
                case 3
                    bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsf(:,:,framenum),1,sd));
                case 4
                    bigframef(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsf(:,:,framenum),1,sd));
                end
            end
            
            
            if flagsumfilter==1
                
                % filtering
                fprintf('\rFiltering: \rframe');
                fprintf('%d ', i);
                if style==2 | style==4 
                    if size(fframe,2)==1 & i==fframe % for style 2 and 4, if there is a black frame, no filtering.
                    else
                    singlemapsf=OIEasyFilter(singlemaps(framenum), LPMethod, LPKernel, HPMethod, HPKernel);
                    mapsfs=OIEasyFilter(maps(framenum), LPMethod, LPKernel, HPMethod, HPKernel);
                    end
                end
                fprintf('\r');

                singlemapsfmean=mean(singlemapsf(:));
                singlemapsfstd=std(singlemapsf(:));
                mapsfsmean=mean(mapsfs(:));
                mapsfsstd=std(mapsfs(:));              

                
                
                switch style
                case 1
                    bigframefs(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsfs(framenum),1,sd));
                case 2
                    if size(fframe,2)==1 && framenum==fframe
                    else
                        bigframefs(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsfs(framenum),1,sd));
                    end
                case 3
                    bigframefs(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsfs(framenum),1,sd));
                case 4
                    bigframefs(y:y+Height-1, x:x+Width-1)=norm_to_uint8(OIClipH2(mapsfs(framenum),1,sd));
                end
            end
            
            
        end
    end
end

% plot a average map at bottom right corner (note: it will overwrite a map if there is one frame
%if style==1
%    bigframe(end-Height+1-spacepixy:end-spacepixy, end-Width+1-spacepixx:end-spacepixx)=mean(maps, 3);
%else
%    bigframe(end-Height+1-spacepixy:end-spacepixy, end-Width+1-spacepixx:end-spacepixx)=norm_to_uint8(OIClipH2(mean(maps, 3),1,sd));
%end

imwrite(imresize(bigframe,0.25), strcat(resultfolder,  outfile(1:end-4), '-avg.bmp'), 'bmp');
imwrite(imresize(norm_to_uint8(bigframeSTD),0.25), jet, strcat(resultfolder,  outfile(1:end-4), '-std.bmp'), 'bmp');
saveas(fig1,[resultfolder,'Z scores.fig']);

if flagtrialfilter==1
    imwrite(bigframef, strcat(resultfolder,  outfile(1:end-4), 'f-avg.bmp'), 'bmp');
end

if flagsumfilter==1
    imwrite(bigframefs, strcat(resultfolder,  outfile(1:end-4), 'f-avg.bmp'), 'bmp');
end

end

if BNC==1
    BNCmean=mean(BNCs, 4);
    fidBNC=fopen(strcat(resultfolder, 'BNCs-avg', '.txt'), 'w');
    fprintf(fidBNC, 'This file contains voltage timecourse of 8 BNCs\r\n\t');
    fprintf(fidBNC, '%f\t', 1:NFrames);
    fprintf(fidBNC, '\r\n');
    for stim=1:NConds
        fprintf(fidBNC, 'Condition#%d \r\n', stim);
        for i=1:8
            fprintf(fidBNC, 'BNC#%d\t', i);
            fprintf(fidBNC, '%f \t', BNCmean(i,:,stim)./10000);
            fprintf(fidBNC, '\r\n');
        end
    end
    fclose(fidBNC);
end
return
