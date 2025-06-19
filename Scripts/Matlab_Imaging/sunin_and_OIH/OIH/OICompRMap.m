function OIDiffRMap()
%
% Ver 1.0: (101011 Hisashi)     Compare between R-value maps
%                               This program works with OICorrH
%
%
clear all;
% __________________ Start User Input ____________________

OIcorrfolder = 'D:\expt1\100511JerL2_AtnS\run99\OICorr\OICorrH2\';
imagefolder1 = 'Atn_M2K_025\'; %  
imagefolder2 = 'Psv_M2K_025\'; % 
outputfolder = 'Diff_Atn-Psv_M2K_025\'; %Output files will be saved

% filename = {            % blk file names, leave empty for automatic searching
% };
statlutfile = 'statcolorBR.lut';
lutfile = 'jetH2.lut';
prefix = 'Comp_Atn-Psv_';
flagoutputcolortable = 0;
inputimageformat = 'ivf'; %'bmp' or 'ivf'
outputimageformat = 'bmp';

inputrrange=1; %1 means r value range -1~+1; This is used only when inputimageformat='bmp'
outputrrange=0.5; %1 means r value range -1~+1

% ___________________
TestType=1;             % 1: Pearson's linear correlation coefficient; 2: 'Spearman' computes Spearman's rho
tail=2;             %1: one sided; 2: both sizded
threshold=[0.0001 0.001 0.01 0.05 0.1];	% threshold p values, just for generating thresholded p maps, each value will have a thresholded map (i.e. 0.05 will generate a p<0.05 map)
pmax=0.05;
pmin=0.00001;

%___________________ end of user input _______________
if 	imagefolder1(end)~='\';
	imagefolder1=[imagefolder1, '\'];	
end
imagefolder1path=[OIcorrfolder, imagefolder1];

if 	imagefolder2(end)~='\';
	imagefolder2=[imagefolder2, '\'];	
end
imagefolder2path=[OIcorrfolder, imagefolder2];

if 	outputfolder(end)~='\';
	outputfolder=[outputfolder, '\'];	
end
outputfolderpath=[OIcorrfolder, outputfolder];

if ~isdir(outputfolderpath)
    mkdir(outputfolderpath);    
end

% if isempty(filename)
tempfilename1=struct2cell(dir([imagefolder1path, '*', 'rmap.', inputimageformat]));
filename1=sort(tempfilename1(1,:)');

tempfilename2=struct2cell(dir([imagefolder2path, '*', 'rmap.', inputimageformat]));
filename2=sort(tempfilename2(1,:)');

for i=1:size(filename1,1)
    fprintf('''%s\\%s - %s\\%s''\n', imagefolder1, getfield(cell2struct(filename1(i), 'junk'), 'junk'), imagefolder2, getfield(cell2struct(filename2(i), 'junk'), 'junk'));
end

if size(filename1,1) == size(filename2,1)
    fprintf('\nFound %d r-map pairs (sorted, check sequence).\n', size(filename1,1));
else
    fprintf('\nThe number of r-maps should be the same between both folder.\n');
end
% end

lut=textread(lutfile);  % this color table should be in sunin folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading frames %%%%%%%%%%%%%%%%%%%%%%%


for k=1:size(filename1,1)
    filename1b=getfield(cell2struct(filename1(k), 'junk'), 'junk');
    filename2b=getfield(cell2struct(filename2(k), 'junk'), 'junk');
%     fprintf('%s\r', bmpfilename1);
    imagename1=strcat(imagefolder1path, filename1b);
    imagename2=strcat(imagefolder2path, filename2b);
    if strcmp(inputimageformat,'bmp')
        image1=double(imread(imagename1)-127)/128*inputrrange;
        image2=double(imread(imagename2)-127)/128*inputrrange;
    elseif strcmp(inputimageformat,'ivf')
        image1=OIReadIVF(imagename1);
        image2=OIReadIVF(imagename2);
        
        samplenumfid=fopen (strcat(imagefolder1path, 'sample_number.txt'), 'r'); %Hisashi
        samplenum1=fscanf(samplenumfid, '%f');
        fclose(samplenumfid);
        
        samplenumfid=fopen (strcat(imagefolder2path, 'sample_number.txt'), 'r'); %Hisashi
        samplenum2=fscanf(samplenumfid, '%f');
        fclose(samplenumfid);
    else
        fprintf('\nThe parameter inputimageformat should be bmp or ivf file.\n');
    end
    diffimage=image1-image2;
%     max(max(diffimage))
    diffimage255=norm_to_uint8b(diffimage, outputrrange*-1, outputrrange);
    outputname=strcat(outputfolderpath, prefix, filename1b);
    imwrite(diffimage255, lut, strcat(outputname(1:end-3), outputimageformat));
    
    z1=log((1+image1)./(1-image1))/2; %Fisher transformation
    z2=log((1+image2)./(1-image2))/2; %Fisher transformation
    
    if TestType==1
        z=(z1-z2)./sqrt((1./(samplenum1-3)+1./(samplenum2-3))); % z-score for 
    elseif TestType==2
        z=(z1-z2)./sqrt((1.06./(samplenum1-3)+1.06./(samplenum2-3)));
    else
        fprintf('\nThe parameter TestType should be 1 or 2.\n');
    end    
    
    p=tail*(1-normcdf(abs(z)));
    
    
    % Statistics
    signmap=double(z>=0)*2-1;
    
    Templ= p > pmin; % logical operation 0 or 1; locations of lower value
    Tempu= p < pmax;% logical operation 0 or 1; locations of higher value
    Tempul= (Templ.*Tempu).*p;% bewteen low and high clips
    Tempul2=Tempul + (pmax*(~Tempu));%clip all avules higher than highClip
    Tempul2=Tempul2+ (pmin*(~Templ));%clip all avules lower than lowClip
    p2=Tempul2;
    p2_log=(log(p2)-log(pmax))/(log(pmin)-log(pmax));

    pmap=uint8(round(p2_log.*signmap*127+128));

    imwrite(pmap, strcat(outputname(1:end-3), '_p-', num2str(pmax,'%1.5f'), '_pmap.bmp'));
    lutstat=textread(statlutfile);  % this color table should be in sunin folder
    imwrite(pmap,lutstat, strcat(outputname(1:end-3), '_p-', num2str(pmax,'%1.5f'), '_pmapcolor.bmp'));

%         [maptemp, framemedian, lowClip, highClip] = OIClipH(p, 1, 1);
%         imwrite(norm_to_uint8(maptemp), strcat(sitefolder, prefix1, prefix, strcat(corrmaskname(i,:)), '_p-all.bmp'));
    for q=threshold
        b=double(p<q);
        imwrite(uint8(b*255), strcat(outputname(1:end-3), '_p-', num2str(q,'%1.5f'), '.bmp'));
    end            
end

fprintf('\rConverted image files were saved in %s\n\r', outputfolder);

% output a color table
if flagoutputcolortable
    colorbarmap=zeros(60, 256, 3);
    for i=1:30
        colorbarmap(i, :, :)=256.*lutstat(:,:);
    end
    imwrite(uint8(colorbarmap), lutstat, strcat(resultfolder, 'statcolortable.tif'), 'tiff'); 

    for i=1:30
        colorbarmap(i, :, :)=256.*lutrho(:,:);
    end
    imwrite(uint8(colorbarmap), lutrho, strcat(resultfolder, 'rcolortable.tif'), 'tiff'); 
end


%%%%% Reference %%%%%
% Spearman's rank correlation coefficient
% http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
% http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html
% http://faculty.vassar.edu/lowry/rdiff.html
%
% Pearson product-moment correlation coefficient
% http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
% http://support.sas.com/kb/24/995.html
%
% Fieller, E.C. et al (1957) Tests for rank correlation coefficients :I. Biometrika 44, 470–481


return;
