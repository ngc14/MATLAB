% This file is a section of Sunin code for OI data process
% Purpose: Super Pixel time course all trials output 
% 
% 081125: Separated from Core program by Hisashi
% Originally written by Hisashi

if flagloadspvalues==1
    load(strcat(resultfolder, 'spvalues.mat'));
    load(strcat(resultfolder, 'stimsumnum.mat'));
end

if ~isdir([runfolder, SPresultname, 'superpixel_data_all'])
    mkdir([runfolder, SPresultname], 'superpixel_data_all');    
end
xlabel=(1:FramesPerStim)/FrameRate;     % in unit of second

%% raw
if spoperation==1 || spoperation==3
    for (j=1:masknum)
        spfidalla= fopen(strcat(runfolder, SPresultname, strcat('superpixel_data_all\', 'superpixel_raw_data_all_', spmaskname(j, :), '.txt')), 'w'); 
    %         spfidalla= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
        fprintf(spfidalla, 'raw value\t');
        fprintf(spfidalla, '%f\t', xlabel);
        fprintf(spfidalla, '\r');
        for (i=1:NStim)
            fprintf(spfidalla, 'stim%d\t', i);
            fprintf(spfidalla, 'f%d\t', 1:FramesPerStim);
            if flagsaveframeave
                fprintf(spfidalla, 'fframe\t');
                fprintf(spfidalla, 'sumrange\t');
            end
    %             fprintf(spfidalla, 'gs%d\t');
            fprintf(spfidalla, '\r');
            counter=0;
            for k=blockselect  % start block by block processing
                fprintf(spfidalla, 'block%d\t', k);
                if goodstim(k, i)~=0
                    counter=counter+1;
                    for (f=1:framenum)
                        fprintf(spfidalla, '%6.5f\t', spvalues(f, j, i, counter));			% original value
                    end
    %                     fprintf(spfidalla, '%d', 1);
                else
                    for (f=1:framenum)
                        fprintf(spfidalla, '\t');
                    end
    %                     fprintf(spfidalla, '%d', 0);
                end
                fprintf(spfidalla, '\r');
            end  
            fprintf(spfidalla, '\r');
        end
        fclose(spfidalla);
    end
end
    
%% dRR 	
if spoperation==2 || spoperation==3
    for (j=1:masknum)
        spfidallb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_data_all\', 'superpixel_dRR_data_all_', spmaskname(j, :), '.txt')), 'w');
    %         spfidallb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
        fprintf(spfidallb, 'dRR value\t');
        fprintf(spfidallb, '%f\t', xlabel);
        fprintf(spfidallb, '\r');
        for (i=1:NStim)
            fprintf(spfidallb, 'stim%d\t', i);
            fprintf(spfidallb, 'f%d\t', 1:FramesPerStim);
            if flagsaveframeave
                fprintf(spfidallb, 'fframe\t');
                fprintf(spfidallb, 'sumrange\t');
            end
    %             fprintf(spfidallb, 'gs%d\t');
            fprintf(spfidallb, '\r');
            counter=0;
            for k=blockselect  % start block by block processing
                fprintf(spfidallb, 'block%d\t', k);
                if goodstim(k, i)~=0
                    counter=counter+1;
                    for (f=1:framenum)
                        fprintf(spfidallb, '%3.8f\t', spvaluesdrr(f, j, i, counter));			% original value
                    end
    %                     fprintf(spfidallb, '%d', 1);
                else
                    for (f=1:framenum)
                        fprintf(spfidallb, '\t');
                    end
    %                     fprintf(spfidallb, '%d', 0);
                end
                fprintf(spfidallb, '\r');
            end  
            fprintf(spfidallb, '\r');
        end
        fclose(spfidallb);
    end
end



