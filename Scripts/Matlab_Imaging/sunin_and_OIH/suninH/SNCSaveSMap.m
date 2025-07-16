% This file is a section of Sunin code for OI data process
% Purpose: Save final block-averaged maps (subtraction maps (and single condition maps)). 
% 
% Note
% HDL 060324

if flagmap==1
    startmap=NStim+1;
else
    startmap=1;
end
for i=startmap:mapnum
    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
    if (flagsavedata)
        dlmwrite(strcat(resultfolder, 'data\', num2str(i),'.txt'), tempmaps(:,:,i), '\t'); 
    end
    [maptemp, framemedian, lowClip, highClip] = OIClip(maps(:,:,i), clipmethod, clipvalue, immask);
    maptemp = norm_to_uint8(maptemp);         
    if (i<=9)
        savefilename=strcat('0', num2str(i), '_', mapname, ext, '.bmp');
    else    
        savefilename=strcat(num2str(i), '_', mapname, ext, '.bmp');
    end
    if i<=NStim
        savefilename=strcat(resultfolder, 'blkavg\singlecond\', savefilename);
    else
        savefilename=strcat(resultfolder, 'blkavg\', savefilename);
    end
    imwrite(maptemp, savefilename, 'bmp');
    %% IVF image output for constructiong Polar
    if flagsaveivf
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(maps(:,:,i), savefilename); %note: width and height are switched;
    end
end
