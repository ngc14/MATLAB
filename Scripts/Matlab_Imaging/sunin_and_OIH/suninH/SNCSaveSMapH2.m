% This file is a section of Sunin code for OI data process
% Purpose: Save final block-averaged maps (subtraction maps (and single condition maps)). 
% 
% Note
% Hisashi 121118: The outout format can be changed. The format should be defined by 'outputimageformat' in suninuser. Since H33.
% Hisashi 090807: norm_to_uint8 was changed to norm_to_uint8b
% HDL 060324

if exist('outputimageformat','var') == 0
	outputimageformat='bmp';
end
if exist('flagcolorcordedmap','var') == 0
	flagcolorcordedmap=0;
end
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
            savefilename=strcat('0', num2str(i), '_', mapname, ext);
        else    
            savefilename=strcat(num2str(i), '_', mapname, ext);
        end
    else
        j=i-NStim;
        if (j<=9)
            savefilename=strcat('0', num2str(j), '_', mapname, ext);
        else    
            savefilename=strcat(num2str(j), '_', mapname, ext);
        end
    end
    if i<=NStim
%         savefilename2=strcat(resultfolder, 'blkavg\singlecond\', savefilename, '.', outputimageformat);
        savefilename2=strcat(resultfolder, 'singlecond\', savefilename, '.', outputimageformat);
    else
%         savefilename2=strcat(resultfolder, 'blkavg\', savefilename, '.', outputimageformat);
        savefilename2=strcat(resultfolder, savefilename, '.', outputimageformat);
    end
    imwrite(maptemp, savefilename2, outputimageformat);
%     imwrite(maptemp, strcat(savefilename(1:end-4), outputimageformat)); %since sunincore34H028
    %% IVF image output for constructiong Polar
    if flagsaveivf
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(maps(:,:,i), savefilename); %note: width and height are switched;
    end
    if flagcolorcordedmap && i>NStim % Since H33
		if ccLPKernel>0 || ccHPKernel>0
			maptemp2=double(maptemp);
			maptemp2=OIEasyFilterH2(maptemp2, ccLPMethod, ccLPKernel, ccHPMethod, ccHPKernel);
			maptemp2=norm_to_uint8b(maptemp2, 0, 255);
		else
			maptemp2=maptemp;
		end
		if ~isdir(strcat(resultfolder, ccoutputfolder))
			mkdir(strcat(resultfolder, ccoutputfolder));    
		end
	    savefilename3=strcat(resultfolder, ccoutputfolder,  ccprefix, savefilename, '.', outputimageformat);
	    cclut=textread(cclutfile);  % this color table should be in sunin folder
		imwrite(maptemp2, cclut, savefilename3);
		if ccflagoutputcolortable
			colorbarmap=zeros(60, 256, 3);
			for i=1:30
				colorbarmap(i, :, :)=256.*cclut(:,:);
			end
			imwrite(uint8(colorbarmap), cclut, strcat(resultfolder, ccoutputfolder, 'colortable', '.', outputimageformat)); 
		end   	
    end
end
% fprintf('\nMaximum frameStd: %8.8f.\r', maxframeStd);