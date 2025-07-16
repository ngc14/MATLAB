fidmapinfo=fopen(strcat(resultfolder, strcat('mapinfo', ext, '.txt')), 'w');
fprintf(fidmapinfo, 'this file contains mapinfo\r');
%if (flagmap<2)
%    fprintf(fidmapinfo, '!Note: flagmap < 2, may not have all the outputs. \r');
%end    
for i=1:mapnum
    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
    fprintf(fidmapinfo, '%d%s\t', i,mapname);
    for j=1:blocknum
        fprintf(fidmapinfo, '%f\t', mapinfo(1, j, i));
    end
    fprintf(fidmapinfo, '\r');
end
fprintf(fidmapinfo, '\r\rHighClip-Median\r');
for i=1:mapnum
    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
    fprintf(fidmapinfo, '%d%s\t', i,mapname);
    for j=1:blocknum
        fprintf(fidmapinfo, '%f\t', mapinfo(2, j, i)-mapinfo(1, j, i));
    end
    fprintf(fidmapinfo, '\r');
end    



meanfid = fopen(strcat(resultfolder, 'mean.txt'), 'a');   %save mean value


fprintf(fidmapinfo, 'averaged map\r');
if flagmap==1
    startmap=NStim+1;
else
    startmap=1;
end
fprintf(meanfid, 'min\tmax\tmedian\tmean\tstd\n');
for i=startmap:mapnum
    mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
    % save min, max, mean, median, stdev for each map
    temp1=reshape(maps(:,:,i), prod(size(maps(:,:,i))), 1);
    fprintf(meanfid, 'map%d\t%10.7f\t%10.7f\t%10.7f\t%10.7f\t%10.7f\n',i, min(temp1), max(temp1), median(temp1), mean(temp1), std(temp1));
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
    if ((flagsaveivf))
        savefilename=strcat(savefilename(1:end-4), '.ivf');
        OIWriteIVF(maps(:,:,i), savefilename); %note: width and height are switched;
    end
    fprintf(fidmapinfo, '%d(%s)\t%f \t%f\r', i, mapname, lowClip, highClip);
end


fprintf(meanfid, '\r\r');
fclose(meanfid);
fclose(fidmapinfo);
