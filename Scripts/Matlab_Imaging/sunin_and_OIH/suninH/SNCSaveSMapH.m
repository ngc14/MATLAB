% This file is a section of Sunin code for OI data process
% Purpose: Save final block-averaged maps (subtraction maps (and single condition maps)). 
% 
% Note
% HIsashi 090807: norm_to_uint8 was changed to norm_to_uint8b
% HDL 060324

if flagmap==1
    startmap=NStim+1;
else
    startmap=1;
end
maxframeStd=0;
fprintf('\n');
for i=startmap:mapnum
    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
    if (flagsavedata)
        dlmwrite(strcat(resultfolder, 'data\', num2str(i),'.txt'), tempmaps(:,:,i), '\t'); 
    end
    [maptemp, framemedian, lowClip, highClip] = OIClipH(maps(:,:,i), clipmethod, clipvalue, cpmask);
    if clipmethod
        cliptextfid= fopen(strcat(resultfolder, strcat('map_clipping_range', '.txt')), 'a');
        if (clipmethod==1 || clipmethod==2 || clipmethod==3)
    %         if flagshowframeStd==1
    %             fprintf('\nmap %s was clipped at median %8.8f ± %1.3f X SD %8.8f.', mapname, framemedian, clipvalue, (highClip-framemedian)/clipvalue);
                fprintf(cliptextfid, 'Map %s was clipped at median %g +- %g X SD %g.\n', mapname, framemedian, clipvalue, (highClip-framemedian)/clipvalue);
    %         end
        elseif (clipmethod==4)
                fprintf(cliptextfid, 'Map %s was clipped at median %g +- %g.\n', mapname, framemedian, (highClip-framemedian));
        elseif (clipmethod==5)
                fprintf(cliptextfid, 'Map %s was clipped at %g +- %g X SD %g.\n', mapname, framemedian, clipvalue, (highClip-framemedian)/clipvalue);
        elseif (clipmethod==6)
                fprintf(cliptextfid, 'Map %s was clipped at %g +- %g.\n', mapname, framemedian, (highClip-framemedian));                
        end
        fclose(cliptextfid); 
    end
%     maptemp = norm_to_uint8(maptemp);
    maptemp = norm_to_uint8b(maptemp, lowClip, highClip); %fixed by Hisashi on 8/7/2009
    if i<=NStim
        if (i<=9)
            savefilename=strcat('0', num2str(i), '_', mapname, ext, '.bmp');
        else    
            savefilename=strcat(num2str(i), '_', mapname, ext, '.bmp');
        end
    else
        j=i-NStim;
        if (j<=9)
            savefilename=strcat('0', num2str(j), '_', mapname, ext, '.bmp');
        else    
            savefilename=strcat(num2str(j), '_', mapname, ext, '.bmp');
        end
    end
    if i<=NStim
        savefilename=strcat(resultfolder, 'blkavg\singlecond\', savefilename);
    else
        savefilename=strcat(resultfolder, 'blkavg\', savefilename);
    end
    imwrite(maptemp, savefilename, 'bmp');
%     imwrite(maptemp, strcat(savefilename(1:end-4), outputimageformat)); %since sunincore34H028
    %% IVF image output for constructiong Polar
    if flagsaveivf
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(maps(:,:,i), savefilename); %note: width and height are switched;
    end
end
% fprintf('\nMaximum frameStd: %8.8f.\r', maxframeStd);