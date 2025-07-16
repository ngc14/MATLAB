% This file is a section of Sunin code for OI data process
% Purpose: Vector analysis
% 
% (091230): fixed in 'vectormask' function by Hisashi
% (090807): norm_to_uint8 was changed to norm_to_uint8b by Hisashi
% (081125): Separated from Core program by Hisashi
% Originally written by HDL

fprintf('\r\rCreating polar image, vector calculating...\r');
t1=clock;
Nvect=size(vect, 1);
% calculate vector map
for j=1:Nvect
    Sum1=zeros(FrameHeight, FrameWidth);
    Sum2=zeros(FrameHeight, FrameWidth);
    condlist1=getfield(cell2struct(vect(j, 2), 'junk'), 'junk');
    for n=condlist1
        Sum1=Sum1+maps(:,:,n);
    end
    Sum1=Sum1/size(condlist1,2);
    
    condlist2=getfield(cell2struct(vect(j, 3), 'junk'), 'junk');
    for n=condlist2
        Sum2=Sum2+maps(:,:,n);
    end
    if size(condlist2,2)~=0
        Sum2=Sum2/size(condlist2,2);
    end
    vmaps(:,:,j)=Sum1-Sum2;
end
if ~isdir([runfolder, vectorresultname])
    mkdir(runfolder, vectorresultname);    
end

% Save single vector map

for i=1:Nvect
%     [maptemp, framemedian, lowClip, highClip] = OIClipH(vmaps(:,:,i), clipmethod, clipvalue, immask);
%     [maptemp, framemedian, lowClip, highClip] = OIClipH(vmaps(:,:,i), 1, 1, 1.5); %by Hisashi on 090118

    if flagtrialfilter == 0
        if flagvectormaskfilter==0
            maptemp=OIEasyFilterH2(vmaps(:,:,i), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize);

            if strcmp(lpfmethod, 'fft') || strcmp(hpfmethod, 'fft')
                if FrameHeight == FrameWidth
                    fftfilter=OIMakFIRfilterH(FrameHeight,LPfiltermethod, lowpass, HPfiltermethod, highpass);
                else
                    fprintf('Error: For FFT filtering, maps should be square.\n');
                end     
                maptemp=ifft2(fft2(maptemp).*fftfilter);
            end        
        else
            maptemp=OIEasyFilterH2wMask(vmaps(:,:,i), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize, vectmask);

            if strcmp(lpfmethod, 'fft') || strcmp(hpfmethod, 'fft')
                fprintf('Error: For FFT filtering, maskfiltering can not be used.\n');
            end
        end
    else
        maptemp=vmaps(:,:,i);
    end   
    
    
    
    [maptemp2, framemedian, lowClip, highClip] = OIClipH(maptemp, clipmethod, clipvalue, cpmask); %by Hisashi on 090118
    %     maptemp = norm_to_uint8(maptemp);
    maptemp2 = norm_to_uint8b(maptemp2, lowClip, highClip); %fixed by Hisashi on 8/7/2009
    if i<10
        savefilename=strcat('vect0', num2str(i), ext, '.bmp');  
    else
        savefilename=strcat('vect', num2str(i), ext, '.bmp');  
    end
    
    savefilename=strcat(runfolder, vectorresultname, savefilename);
    imwrite(maptemp2, savefilename, 'bmp');
    if flagsaveivfvector
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(maptemp, savefilename);	
    end
end

clear maptemp;
clear maptemp2;

% Generate vector map

singlemap = vmaps; % by Hisashi
%     if flagsaveivfvector
%         for i = 1:Nvect
%             if k<10
%                 filename=strcat('vect0', num2str(i), ext, '.ivf');  
%             else
%                 filename=strcat('vect', num2str(i), ext, '.ivf');  
%             end
%             filename=strcat(runfolder, vectorresultname, filename);
%             b=OIReadIVF(filename);
%     %         if (i==1)
%     %             [height, width]=size(b);
%     %             singlemap = zeros(height, width, Nvect); 
%     %         end
%     %         singlemap(:, :, i) = b; %norm_to_01(b);
%         end
%     end

% 	lut=textread('bwfix.lut');
lut=textread(lutfile); % by Hisashi
% Polarmap=OIPolarH2(singlemap, lut, vectmask, vectclipsd, lowpass, highpass, LPfiltermethod, HPfiltermethod); %modified by Hisashi 071022


if flagmaskfilter==1
    if ~isdir(strcat(outputdriver, outputfolder, expname, runname, 'masks\vectormask'))
        vectmaskname = strcat(outputdriver, outputfolder, expname, runname, 'masks\filtermask\default.bmp');
    else
        vectmaskname = strcat(outputdriver, outputfolder, expname, runname, 'masks\vectormask\default.bmp');
    end
    vectmasktemp = (double(imread(vectmaskname, 'bmp'))/255);
    if size(vectmasktemp,3)==1
        vectmask = vectmasktemp;
    else
        vectmask = vectmasktemp(:,:,1);
    end
    clear vectmasktemp;
end

if flagtrialfilter == 0
    Polarmap=OIPolarH2(singlemap, lut, clipmethod, clipvalue, lpfmethod, lpkernelsize, hpfmethod, hpkernelsize, flagmaskfilter, vectmask); %modified by Hisashi 090118
else
    Polarmap=OIPolarH2(singlemap, lut, clipmethod, clipvalue, lpfmethod, 0, hpfmethod, 0, flagmaskfilter, vectmask); %modified by Hisashi 090118
end    
    %	figure; image(norm_to_uint8(Polarmap.ang));  axis equal; axis off; colormap(lut);  colorbar; %%load billcolorfix; 
%	figure; imagesc(norm_to_uint8(Polarmap.mag)); axis equal; axis off; colormap(gray(256)); colorbar; 

%%
filename1 = strcat('Angle_', lpfmethod, '_', num2str(lpkernelsize), '_', hpfmethod, '_', num2str(hpkernelsize), ext, '.tif');
filename1 = strcat(runfolder, vectorresultname, filename1);
imwrite(norm_to_uint8(Polarmap.ang), lut, filename1, 'tiff');

% filename1b = strcat('Angle_grayscale', LPfiltermethod, '_', num2str(lowpass), '_', HPfiltermethod, '_', num2str(highpass), ext, '.tif');
% filename1b = strcat(runfolder, vectorresultname, filename1b);
% imwrite(norm_to_uint8(Polarmap.ang), filename1b, 'tiff');
% 
% filename1c = strcat('Angle_contour', LPfiltermethod, '_', num2str(lowpass), '_', HPfiltermethod, '_', num2str(highpass), ext, '.tif');
% filename1c = strcat(runfolder, vectorresultname, filename1c);
% 
% I = norm_to_uint8(Polarmap.ang);
% figure('Position', [101 101 604 604]);
% [C,h] = imcontour(I,Nvect*2);
% set(h, 'LineWidth',2);
% clabel(C,h);
% print('-dtiff ', 'aaa');



% imwrite(norm_to_uint8(Polarmap.ang), filename1c, 'tiff');

%	OIWriteIVF(Polarmap.ang, [filename1(1:end-4), ext, '.ivf']);
%%
filename2 = strcat('Mag_', lpfmethod, '_', num2str(lpkernelsize), '_', hpfmethod, '_', num2str(hpkernelsize), ext, '.tif');
filename2 = strcat(runfolder, vectorresultname, filename2);
% imwrite(norm_to_uint8(OIClipH(Polarmap.mag, clipmethod, clipvalue)), filename2, 'tiff');
imwrite(norm_to_uint8(OIClipH_lowMin(Polarmap.mag, polarclipmethod, polarclipvalue)), filename2, 'tiff'); %by Hisashi on 090118
if flagsaveivfvector
    OIWriteIVF(Polarmap.mag, [filename2(1:end-4), ext, '.ivf']);
end

%%
filename3 = strcat('Polar_', lpfmethod, '_', num2str(lpkernelsize), '_', hpfmethod, '_', num2str(hpkernelsize), ext, '.tif');
filename3 = strcat(runfolder, vectorresultname, filename3);
mag=norm_to_01(Polarmap.mag);
% mag=OIClipH(mag, 1, polarclipsd);   % to adjust map darkness  
mag=OIClipH_lowMin(mag, polarclipmethod, polarclipvalue);   % to adjust map darkness   %by Hisashi on 090118
switch precisiontype
    case 'single'
        ang=single(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
    case 'double'
        ang=double(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
end
polarmap=ang;
polarmap(:,:,1)=ang(:,:,1).*mag;
polarmap(:,:,2)=ang(:,:,2).*mag;
polarmap(:,:,3)=ang(:,:,3).*mag;
imwrite(norm_to_uint8(polarmap), lut, filename3, 'tiff');
%	OIWriteIVF(polarmap, [filename3(1:end-4), ext, '.ivf']);

%%
filename4 = strcat('Sum_', lpfmethod, '_', num2str(lpkernelsize), '_', hpfmethod, '_', num2str(hpkernelsize), ext, '.tif');
filename4 = strcat(runfolder, vectorresultname, filename4);
% imwrite(norm_to_uint8(OIClipH(Polarmap.sum, clipmethod, clipvalue)), filename4, 'tiff');
% imwrite(norm_to_uint8(OIClipH(Polarmap.sum, 1, 1.5)), filename4, 'tiff');
imwrite(norm_to_uint8(OIClipH(Polarmap.sum, clipmethod, clipvalue)), filename4, 'tiff');
%	OIWriteIVF(Polarmap.sum, [filename4(1:end-4), ext, '.ivf']);

% output a color table
colorbarmap=zeros(60, 256, 3);
for i=1:30
    colorbarmap(i, :, :)=256.*lut(:,:);
end
for i=1:Nvect
    x=(i-1)*floor(256/Nvect)+1;
    colorbarmap(36:60, x:x+9, 1)=256*ones(25, 10)*lut(x,1);
    colorbarmap(36:60, x:x+9, 2)=256*ones(25, 10)*lut(x,2);
    colorbarmap(36:60, x:x+9, 3)=256*ones(25, 10)*lut(x,3);
end
imwrite(uint8(colorbarmap), strcat(runfolder, vectorresultname, 'colortable.tif'), 'tiff'); 

fprintf('Done (time for vector analysis: %f minutes)\r', etime(clock, t1)/60);


