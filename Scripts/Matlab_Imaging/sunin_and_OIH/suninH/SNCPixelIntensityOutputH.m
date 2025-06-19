% This file is a section of Sunin code for OI data process
% Purpose: Output intensity at each pixel in maps for making histograms
% 
% Originally written by Hisashi

if flagmap==1 % by Hisashi on 080229
    startmap=NStim+1;
else
    startmap=1;
end

for k=startmap:mapnum
    piofid1= fopen(strcat(runfolder, strcat('OIPixelIntensityOutputs\', 'superpixel_raw1', '.txt')), 'w');
    for i=1:pixelmasknum
        
        maps(:,:,k)


    
    
    
    
end





if k==-1            % use non-blood-vessel region
    histomask=imread(nobvmask, 'bmp');
elseif k==0         % use the whole map
    histomask=ones(FrameHeight,FrameWidth);
else                % use specific domains
    histomask=mask(:,:,k);    
end
for i=1:mapnum     % get the values from the masked maps
    switch precisiontype
        case 'single'
            pixnum(i)=sum(sum(single(histomask)));              % following 6 lines put masked pixels into 'maskary', any better way?
            maskmap= maps(:,:,i) .* single(histomask);
            maskmap= maskmap + (single(histomask)-1)*(-999);    %% NOTE: -999 need change if the elements of the map can be that value or higher
        case 'double'
            pixnum(i)=sum(sum(double(histomask)));              % following 6 lines put masked pixels into 'maskary', any better way?
            maskmap= maps(:,:,i) .* double(histomask);
            maskmap= maskmap + (double(histomask)-1)*(-999);    %% NOTE: -999 need change if the elements of the map can be that value or higher
    end
    maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
    maskmap= sort(maskmap, 1);
    maskary= maskmap(1:pixnum(i), 1);
    maphist(:,i)=(hist(maskary, bins))';
end
% output
histfid1= fopen(strcat(resultfolder, strcat('hist1', '.txt')), 'w');  % histo type 1: across all maps but using one mask
fprintf(histfid1, '#\tbin\t');  
for i=1:mapnum
   fprintf(histfid1, '%s\t', getfield(cell2struct(mapnames(i), 'junk'), 'junk'));
end
%    fprintf(histfid1, '\r\t\t');
%    for i=1:mapnum
%        fprintf(histfid1, '%d\t', pixnum(i));
%    end
fprintf(histfid1, '\r');
for j=1:binnum
    fprintf(histfid1, '%d\t%f\t', j, bins(j));
    for i=1:mapnum
        fprintf(histfid1, '%f\t', maphist(j, i)/pixnum(i));
    end
    fprintf(histfid1, '\r');
end
clear maphist;
clear pixnum;
fclose(histfid1);

% 2. hist type 2: output pixel value distrubution from one map but for each masks
i=histindex(1,2);
for j=1:masknum
    switch precisiontype
        case 'single'
            pixnum(j)=sum(sum(single(mask(:,:,j))));         % following 6 lines put masked pixels into 'maskary', better way?
            maskmap= maps(:,:,i) .* single(mask(:, :, j));
            maskmap= maskmap + (single(mask(:, :, j))-1)*(-999);    
        case 'double'
            pixnum(j)=sum(sum(double(mask(:,:,j))));         % following 6 lines put masked pixels into 'maskary', better way?
            maskmap= maps(:,:,i) .* double(mask(:, :, j));
            maskmap= maskmap + (double(mask(:, :, j))-1)*(-999);    
    end
    maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
    maskmap= sort(maskmap, 1);
    maskary= maskmap(1:pixnum(j), 1);
    maphist(:,j)=(hist(maskary, bins))';
end
% output
histfid2= fopen(strcat(resultfolder, strcat('hist2', '.txt')), 'w');  
fprintf(histfid2, '#\tbin\t');
for i=1:masknum
    fprintf(histfid2, '%s\t',spmaskname(i,:));
end
%    fprintf(histfid2, '\r\t\t');
%    for i=1:masknum
%        fprintf(histfid2, '%d\t', pixnum(i));
%    end
fprintf(histfid2, '\r');
for j=1:binnum 
    fprintf(histfid2, '%d\t%f\t', j, bins(j));
    for i=1:masknum
        fprintf(histfid2, '%f\t', maphist(j, i)/pixnum(i));
    end
    fprintf(histfid2, '\r');
end
clear maphist;
clear pixnum;
fclose(histfid2);

% 3. hist type 3: output gray value histograms from provided map-mask pairs
pairs=size(hist3index, 1);  % number of pairs
for i=1:pairs
    imap =hist3index(i, 1); %index of map
    imask=hist3index(i, 2);
    switch precisiontype
        case 'single'
            pixnum(i)=sum(sum(single(mask(:,:,imask))));         % following 6 lines put masked pixels into 'maskary', better way?
            maskmap= maps(:,:,imap) .* single(mask(:, :, imask));
            maskmap= maskmap + (single(mask(:, :, imask))-1)*(-999);    
        case 'double'
            pixnum(i)=sum(sum(double(mask(:,:,imask))));         % following 6 lines put masked pixels into 'maskary', better way?
            maskmap= maps(:,:,imap) .* double(mask(:, :, imask));
            maskmap= maskmap + (double(mask(:, :, imask))-1)*(-999);    
    end
    maskmap= reshape(maskmap, FrameHeight*FrameWidth,1);
    maskmap= sort(maskmap, 1);
    maskary= maskmap(1:pixnum(i), 1);
    maphist(:,i)=(hist(maskary, bins))';
end
% output
histfid3= fopen(strcat(resultfolder, strcat(getfield(cell2struct(filename(1), 'junk'), 'junk'), 'hist3', '.txt')), 'w');  %here also used file names
fprintf(histfid3, '\t\t');
for i=1:pairs
    mapname=getfield(cell2struct(mapnames(hist3index(i,1)), 'junk'), 'junk');
    fprintf(histfid3, '%s\t',mapname);
end
fprintf(histfid3, '\r');
fprintf(histfid3, '#\tbin\t');
for i=1:pairs
    fprintf(histfid3, '%s\t',spmaskname(hist3index(i,2),:));
end
fprintf(histfid3, '\r');
for j=1:binnum 
    fprintf(histfid3, '%d\t%f\t', j, bins(j));
    for i=1:pairs
        fprintf(histfid3, '%f\t', maphist(j, i)/pixnum(i));
    end
    fprintf(histfid3, '\r');
end
clear maphist;
clear pixnum;
fclose(histfid3);



