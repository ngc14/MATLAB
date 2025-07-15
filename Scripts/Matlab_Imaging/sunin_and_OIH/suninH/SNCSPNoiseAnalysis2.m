% This file is a section of Sunin code for OI data process
% Purpose: Super Pixel time course all trials output - Noise Analysis
% 
% 150225: It's found that the calculation of these values for SPcond is not correct. Need fix. Hisashi
% 100213: Variance, Fano (Variance-to-Mean Ratio), CV (Coefficient of Variation) are calculated by Hisashi
% 081205: Originally written by Hisashi

if ~isdir([runfolder, SPresultname, 'superpixel_noise_analysis'])
    mkdir([runfolder, SPresultname], 'superpixel_noise_analysis');    
end
xlabel=(1:FramesPerStim)/FrameRate;     % in unit of second

if spoperation==1 || spoperation==3
    spvaluesad=zeros(framenum, masknum, NStim+NSPcond, blockselectnum);
    spvaluesnad=zeros(framenum, masknum, NStim+NSPcond, blockselectnum);
    spvaluesds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared
    spvaluesnds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared, devided by mean
    spvaluesnnds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared, devided by mean^2
    
    %% All trials
    %% Deviation for raw data
	for (j=1:masknum)
		spfidallda= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw_data_all_Deviation_', spmaskname(j, :), '.txt')), 'w'); 
	%         spfidallda= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidallda, 'Dev_raw value\t');
		fprintf(spfidallda, '%f\t', xlabel);
		fprintf(spfidallda, '\r');
		for (i=1:NStim)
			fprintf(spfidallda, 'stim%d\t', i);
			fprintf(spfidallda, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidallda, 'fframe\t');
				fprintf(spfidallda, 'sumrange\t');
			end
	%             fprintf(spfidallda, 'gs%d\t');
			fprintf(spfidallda, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidallda, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
						fprintf(spfidallda, '%6.5f\t', spvalues(f, j, i, counter)-mean(spvalues(f, j, i, 1:stimsumnum(i))));			% original value
					end
	%                     fprintf(spfidallda, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidallda, '\t');
					end
	%                     fprintf(spfidallda, '%d', 0);
				end
				fprintf(spfidallda, '\r');
			end  
			fprintf(spfidallda, '\r');
		end
		fclose(spfidallda);
	end
	
	%% Absolute Deviation for raw data
	for (j=1:masknum)
		spfidallada= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw_data_all_AbsoluteDeviation_', spmaskname(j, :), '.txt')), 'w'); 
	%         spfidallada= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidallada, 'AD_raw value\t');
		fprintf(spfidallada, '%f\t', xlabel);
		fprintf(spfidallada, '\r');
		for (i=1:NStim)
			fprintf(spfidallada, 'stim%d\t', i);
			fprintf(spfidallada, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidallada, 'fframe\t');
				fprintf(spfidallada, 'sumrange\t');
			end
	%             fprintf(spfidallada, 'gs%d\t');
			fprintf(spfidallada, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidallada, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
						fprintf(spfidallada, '%6.5f\t', abs(spvalues(f, j, i, counter)-mean(spvalues(f, j, i, 1:stimsumnum(i)))));			% original value
                        spvaluesad(f, j, i, counter)=abs(spvalues(f, j, i, counter)-mean(spvalues(f, j, i, 1:stimsumnum(i))));
                        spvaluesds(f, j, i, counter)=(spvalues(f, j, i, counter)-mean(spvalues(f, j, i, 1:stimsumnum(i))))^2; % deviation squared
                    end
	%                     fprintf(spfidallada, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidallada, '\t');
					end
	%                     fprintf(spfidallada, '%d', 0);
				end
				fprintf(spfidallada, '\r');
			end  
			fprintf(spfidallada, '\r');
		end
		fclose(spfidallada);
	end
	
	%% Normalized Absolute Deviation for raw data
	for (j=1:masknum)
		spfidallnada= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw_data_all_NormalizedAbsoluteDeviation_', spmaskname(j, :), '.txt')), 'w'); 
	%         spfidallnada= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidallnada, 'NAD_raw value\t');
		fprintf(spfidallnada, '%f\t', xlabel);
		fprintf(spfidallnada, '\r');
		for (i=1:NStim)
			fprintf(spfidallnada, 'stim%d\t', i);
			fprintf(spfidallnada, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidallnada, 'fframe\t');
				fprintf(spfidallnada, 'sumrange\t');
			end
	%             fprintf(spfidallnada, 'gs%d\t');
			fprintf(spfidallnada, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidallnada, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
                        tempmean1=mean(spvalues(f, j, i, 1:stimsumnum(i)));
                        if tempmean1~=0
                            fprintf(spfidallnada, '%6.5f\t', abs(spvalues(f, j, i, counter)-tempmean1)/tempmean1);			% original value
                            spvaluesnad(f, j, i, counter)=abs(spvalues(f, j, i, counter)-tempmean1)/tempmean1;
                            spvaluesnds(f, j, i, counter)=((spvalues(f, j, i, counter)-tempmean1)^2)/tempmean1; % deviation squared, devided by mean
                            spvaluesnnds(f, j, i, counter)=((spvalues(f, j, i, counter)-tempmean1)^2)/tempmean1/tempmean1; % deviation squared, devided by mean ^2
                        else
                            fprintf(spfidallnada, '%6.5f\t', 0);
                            spvaluesnad(f, j, i, counter)=0;
                            spvaluesnds(f, j, i, counter)=0;
                            spvaluesnnds(f, j, i, counter)=0;
                        end
					end
	%                     fprintf(spfidallnada, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidallnada, '\t');
					end
	%                     fprintf(spfidallnada, '%d', 0);
				end
				fprintf(spfidallnada, '\r');
			end  
			fprintf(spfidallnada, '\r');
		end
		fclose(spfidallnada);
    end
    
    %% Average, SD, Variance, Fano factor, CV

    spnoiseadavg=zeros(framenum, masknum, NStim+NSPcond);
    spnoisenadavg=zeros(framenum, masknum, NStim+NSPcond);
    
    spnoiseadstd=zeros(framenum, masknum, NStim+NSPcond);
    spnoisenadstd=zeros(framenum, masknum, NStim+NSPcond);

    spnoisevar=zeros(framenum, masknum, NStim+NSPcond); % variance
    spnoisefano=zeros(framenum, masknum, NStim+NSPcond); % fano factor?
    spnoisecv=zeros(framenum, masknum, NStim+NSPcond); % coefficient of variation
   
    for i=1:NStim
        if stimsumnum(i)>0
            spnoiseadavg(:,:,i)=sum(spvaluesad(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
            spnoisenadavg(:,:,i)=sum(spvaluesnad(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different

            spnoiseadstd(:,:,i)=std(spvaluesad(:,:,i,1:stimsumnum(i)),0,4);
            spnoisenadstd(:,:,i)=std(spvaluesnad(:,:,i,1:stimsumnum(i)),0,4);
            
            spnoisevar(:,:,i)=sum(spvaluesds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1); % unbiased variance
            spnoisefano(:,:,i)=sum(spvaluesnds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1); % unbiased variance devided by mean, fano factor?
            spnoisecv(:,:,i)=sqrt(sum(spvaluesnnds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1)); % standard deviation devided by mean, coefficient of variation (CV)
       end
    end

    %% Absolute Deviation for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesad(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoiseadavg(:,:,j)=SPcondsum/SPcondsumnum;
                spnoiseadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    
    %% Normalized Absolute Deviation for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnad(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisenadavg(:,:,j)=SPcondsum/SPcondsumnum;
                spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    
    %% Variance for SPcond by Hisashi on 081003 and 090129
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisevar(:,:,j)=SPcondsum/(SPcondsumnum-1); % unbiased variance
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
     %% Fano factor for SPcond by Hisashi on 081003 and 090129
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisefano(:,:,j)=SPcondsum/(SPcondsumnum-1);
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
     %% Coefficient of Variation (CV) for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnnds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisecv(:,:,j)=sqrt(SPcondsum/(SPcondsumnum-1));
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    
    %% Output of Absolute Deviation for raw data
    spfiada1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw1_AbsoluteDeviation', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfiada2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw2_AbsoluteDeviation', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfiada1, 'AD_raw value\t');
    fprintf(spfiada1, '%f\t', xlabel);
    fprintf(spfiada1, '\r');
    for (j=1:masknum)
        fprintf(spfiada1, '%s\t', spmaskname(j, :));
        fprintf(spfiada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada1, 'fframe\t');
            fprintf(spfiada1, 'sumrange\t');
        end
        fprintf(spfiada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfiada1, '%6.5f\t', spnoiseadavg(f, j, i));			% original value
            end
            fprintf(spfiada1, '\r');
        end
        fprintf(spfiada1, '\r');
    end
    fprintf(spfiada1, '\r\r_____SD_____\r\t');
    fprintf(spfiada1, '%f\t', xlabel);
    fprintf(spfiada1, '\r');
    for (j=1:masknum)
        fprintf(spfiada1, '%s\t', spmaskname(j, :));
        fprintf(spfiada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada1, 'fframe\t');
            fprintf(spfiada1, 'sumrange\t');
        end
        fprintf(spfiada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfiada1, '%6.5f\t', spnoiseadstd(f, j, i));
            end
            fprintf(spfiada1, '\r');
        end
        fprintf(spfiada1, '\r');
    end
    fprintf(spfiada1, '\r\r_____Num_____\r\t');
    fprintf(spfiada1, '%f\t', xlabel);
    fprintf(spfiada1, '\r');
    for j=1:masknum
        fprintf(spfiada1, '%s\t', spmaskname(j, :));
        fprintf(spfiada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada1, 'fframe\t');
            fprintf(spfiada1, 'sumrange\t');
        end
        fprintf(spfiada1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfiada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfiada1, '%d\t', stimsumnum(i));
            end
            fprintf(spfiada1, '\r');
        end
        fprintf(spfiada1, '\r');
    end
    fprintf(spfiada1, '\r\r_____SEM_____\r\t');
    fprintf(spfiada1, '%f\t', xlabel);
    fprintf(spfiada1, '\r');
    for (j=1:masknum)
        fprintf(spfiada1, '%s\t', spmaskname(j, :));
        fprintf(spfiada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada1, 'fframe\t');
            fprintf(spfiada1, 'sumrange\t');
        end
        fprintf(spfiada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if stimsumnum(i)>0
                    fprintf(spfiada1, '%6.5f\t', spnoiseadstd(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfiada1, '\r');
        end
        fprintf(spfiada1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfiada2, 'AD_raw value\t');
    fprintf(spfiada2, '%f\t', xlabel);
    fprintf(spfiada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada2, 'fframe\t');
            fprintf(spfiada2, 'sumrange\t');
        end
        fprintf(spfiada2, '\r');
        for (j=1:masknum)
            fprintf(spfiada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfiada2, '%6.5f\t', spnoiseadavg(f, j, i));			% original value
            end
            fprintf(spfiada2, '\r');
        end
        fprintf(spfiada2, '\r');
    end
    fprintf(spfiada2, '\r\r_____SD_____\r\t');
    fprintf(spfiada2, '%f\t', xlabel);
    fprintf(spfiada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada2, 'fframe\t');
            fprintf(spfiada2, 'sumrange\t');
        end
        fprintf(spfiada2, '\r');
        for (j=1:masknum)
            fprintf(spfiada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfiada2, '%6.5f\t', spnoiseadstd(f, j, i));
                end
            end
            fprintf(spfiada2, '\r');
        end
        fprintf(spfiada2, '\r');
    end
    fprintf(spfiada2, '\r\r_____Num_____\r\t');
    fprintf(spfiada2, '%f\t', xlabel);
    fprintf(spfiada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada2, 'fframe\t');
            fprintf(spfiada2, 'sumrange\t');
        end
        fprintf(spfiada2, '\r');
        for (j=1:masknum)
            fprintf(spfiada2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfiada2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfiada2, '\r');
        end
        fprintf(spfiada2, '\r');
    end
    fprintf(spfiada2, '\r\r_____SEM_____\r\t');
    fprintf(spfiada2, '%f\t', xlabel);
    fprintf(spfiada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiada2, 'fframe\t');
            fprintf(spfiada2, 'sumrange\t');
        end
        fprintf(spfiada2, '\r');
        for (j=1:masknum)
            fprintf(spfiada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if stimsumnum(i)>0
                        fprintf(spfiada2, '%6.5f\t', spnoiseadstd(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfiada2, '\r');
        end
        fprintf(spfiada2, '\r');
    end
    fclose(spfiada1);
    fclose(spfiada2);
    
    %% Output of Normalized Absolute Deviation for raw data
    spfinada1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw1_NormalizedAbsoluteDeviation', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfinada2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw2_NormalizedAbsoluteDeviation', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfinada1, 'NAD_raw value\t');
    fprintf(spfinada1, '%f\t', xlabel);
    fprintf(spfinada1, '\r');
    for (j=1:masknum)
        fprintf(spfinada1, '%s\t', spmaskname(j, :));
        fprintf(spfinada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada1, 'fframe\t');
            fprintf(spfinada1, 'sumrange\t');
        end
        fprintf(spfinada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfinada1, '%6.5f\t', spnoisenadavg(f, j, i));			% original value
            end
            fprintf(spfinada1, '\r');
        end
        fprintf(spfinada1, '\r');
    end
    fprintf(spfinada1, '\r\r_____SD_____\r\t');
    fprintf(spfinada1, '%f\t', xlabel);
    fprintf(spfinada1, '\r');
    for (j=1:masknum)
        fprintf(spfinada1, '%s\t', spmaskname(j, :));
        fprintf(spfinada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada1, 'fframe\t');
            fprintf(spfinada1, 'sumrange\t');
        end
        fprintf(spfinada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfinada1, '%6.5f\t', spnoisenadstd(f, j, i));
            end
            fprintf(spfinada1, '\r');
        end
        fprintf(spfinada1, '\r');
    end
    fprintf(spfinada1, '\r\r_____Num_____\r\t');
    fprintf(spfinada1, '%f\t', xlabel);
    fprintf(spfinada1, '\r');
    for j=1:masknum
        fprintf(spfinada1, '%s\t', spmaskname(j, :));
        fprintf(spfinada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada1, 'fframe\t');
            fprintf(spfinada1, 'sumrange\t');
        end
        fprintf(spfinada1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfinada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfinada1, '%d\t', stimsumnum(i));
            end
            fprintf(spfinada1, '\r');
        end
        fprintf(spfinada1, '\r');
    end
    fprintf(spfinada1, '\r\r_____SEM_____\r\t');
    fprintf(spfinada1, '%f\t', xlabel);
    fprintf(spfinada1, '\r');
    for (j=1:masknum)
        fprintf(spfinada1, '%s\t', spmaskname(j, :));
        fprintf(spfinada1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada1, 'fframe\t');
            fprintf(spfinada1, 'sumrange\t');
        end
        fprintf(spfinada1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinada1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinada1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if stimsumnum(i)>0
                    fprintf(spfinada1, '%6.5f\t', spnoisenadstd(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfinada1, '\r');
        end
        fprintf(spfinada1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfinada2, 'NAD_raw value\t');
    fprintf(spfinada2, '%f\t', xlabel);
    fprintf(spfinada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada2, 'fframe\t');
            fprintf(spfinada2, 'sumrange\t');
        end
        fprintf(spfinada2, '\r');
        for (j=1:masknum)
            fprintf(spfinada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfinada2, '%6.5f\t', spnoisenadavg(f, j, i));			% original value
            end
            fprintf(spfinada2, '\r');
        end
        fprintf(spfinada2, '\r');
    end
    fprintf(spfinada2, '\r\r_____SD_____\r\t');
    fprintf(spfinada2, '%f\t', xlabel);
    fprintf(spfinada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada2, 'fframe\t');
            fprintf(spfinada2, 'sumrange\t');
        end
        fprintf(spfinada2, '\r');
        for (j=1:masknum)
            fprintf(spfinada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfinada2, '%6.5f\t', spnoisenadstd(f, j, i));
                end
            end
            fprintf(spfinada2, '\r');
        end
        fprintf(spfinada2, '\r');
    end
    fprintf(spfinada2, '\r\r_____Num_____\r\t');
    fprintf(spfinada2, '%f\t', xlabel);
    fprintf(spfinada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada2, 'fframe\t');
            fprintf(spfinada2, 'sumrange\t');
        end
        fprintf(spfinada2, '\r');
        for (j=1:masknum)
            fprintf(spfinada2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfinada2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfinada2, '\r');
        end
        fprintf(spfinada2, '\r');
    end
    fprintf(spfinada2, '\r\r_____SEM_____\r\t');
    fprintf(spfinada2, '%f\t', xlabel);
    fprintf(spfinada2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinada2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinada2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinada2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinada2, 'fframe\t');
            fprintf(spfinada2, 'sumrange\t');
        end
        fprintf(spfinada2, '\r');
        for (j=1:masknum)
            fprintf(spfinada2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if stimsumnum(i)>0
                        fprintf(spfinada2, '%6.5f\t', spnoisenadstd(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfinada2, '\r');
        end
        fprintf(spfinada2, '\r');
    end
    fclose(spfinada1);
    fclose(spfinada2);    
    %% Output of Variance for raw data
    spfivara1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw1_Variance', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfivara2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw2_Variance', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfivara1, 'Var_raw value\t');
    fprintf(spfivara1, '%f\t', xlabel);
    fprintf(spfivara1, '\r');
    for (j=1:masknum)
        fprintf(spfivara1, '%s\t', spmaskname(j, :));
        fprintf(spfivara1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivara1, 'fframe\t');
            fprintf(spfivara1, 'sumrange\t');
        end
        fprintf(spfivara1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfivara1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfivara1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfivara1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfivara1, '%6.5f\t', spnoisevar(f, j, i));			% original value
            end
            fprintf(spfivara1, '\r');
        end
        fprintf(spfivara1, '\r');
    end

    fprintf(spfivara1, '\r\r_____Num_____\r\t');
    fprintf(spfivara1, '%f\t', xlabel);
    fprintf(spfivara1, '\r');
    for j=1:masknum
        fprintf(spfivara1, '%s\t', spmaskname(j, :));
        fprintf(spfivara1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivara1, 'fframe\t');
            fprintf(spfivara1, 'sumrange\t');
        end
        fprintf(spfivara1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfivara1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfivara1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfivara1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfivara1, '%d\t', stimsumnum(i));
            end
            fprintf(spfivara1, '\r');
        end
        fprintf(spfivara1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfivara2, 'Var_dRR value\t');
    fprintf(spfivara2, '%f\t', xlabel);
    fprintf(spfivara2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfivara2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfivara2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfivara2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfivara2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivara2, 'fframe\t');
            fprintf(spfivara2, 'sumrange\t');
        end
        fprintf(spfivara2, '\r');
        for (j=1:masknum)
            fprintf(spfivara2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfivara2, '%6.5f\t', spnoisevar(f, j, i));			% original value
            end
            fprintf(spfivara2, '\r');
        end
        fprintf(spfivara2, '\r');
    end
    fprintf(spfivara2, '\r\r_____Num_____\r\t');
    fprintf(spfivara2, '%f\t', xlabel);
    fprintf(spfivara2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfivara2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfivara2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfivara2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfivara2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivara2, 'fframe\t');
            fprintf(spfivara2, 'sumrange\t');
        end
        fprintf(spfivara2, '\r');
        for (j=1:masknum)
            fprintf(spfivara2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfivara2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfivara2, '\r');
        end
        fprintf(spfivara2, '\r');
    end
    fclose(spfivara1);
    fclose(spfivara2);
    
    %% Output of Fano factor for raw data
    spfifanoa1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw1_Fano', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfifanoa2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw2_Fano', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfifanoa1, 'Fano_raw value\t');
    fprintf(spfifanoa1, '%f\t', xlabel);
    fprintf(spfifanoa1, '\r');
    for (j=1:masknum)
        fprintf(spfifanoa1, '%s\t', spmaskname(j, :));
        fprintf(spfifanoa1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanoa1, 'fframe\t');
            fprintf(spfifanoa1, 'sumrange\t');
        end
        fprintf(spfifanoa1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfifanoa1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfifanoa1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfifanoa1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfifanoa1, '%6.5f\t', spnoisefano(f, j, i));			% original value
            end
            fprintf(spfifanoa1, '\r');
        end
        fprintf(spfifanoa1, '\r');
    end

    fprintf(spfifanoa1, '\r\r_____Num_____\r\t');
    fprintf(spfifanoa1, '%f\t', xlabel);
    fprintf(spfifanoa1, '\r');
    for j=1:masknum
        fprintf(spfifanoa1, '%s\t', spmaskname(j, :));
        fprintf(spfifanoa1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanoa1, 'fframe\t');
            fprintf(spfifanoa1, 'sumrange\t');
        end
        fprintf(spfifanoa1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfifanoa1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfifanoa1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfifanoa1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfifanoa1, '%d\t', stimsumnum(i));
            end
            fprintf(spfifanoa1, '\r');
        end
        fprintf(spfifanoa1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfifanoa2, 'Fano_raw value\t');
    fprintf(spfifanoa2, '%f\t', xlabel);
    fprintf(spfifanoa2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfifanoa2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfifanoa2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfifanoa2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfifanoa2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanoa2, 'fframe\t');
            fprintf(spfifanoa2, 'sumrange\t');
        end
        fprintf(spfifanoa2, '\r');
        for (j=1:masknum)
            fprintf(spfifanoa2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfifanoa2, '%6.5f\t', spnoisefano(f, j, i));			% original value
            end
            fprintf(spfifanoa2, '\r');
        end
        fprintf(spfifanoa2, '\r');
    end
    fprintf(spfifanoa2, '\r\r_____Num_____\r\t');
    fprintf(spfifanoa2, '%f\t', xlabel);
    fprintf(spfifanoa2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfifanoa2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfifanoa2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfifanoa2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfifanoa2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanoa2, 'fframe\t');
            fprintf(spfifanoa2, 'sumrange\t');
        end
        fprintf(spfifanoa2, '\r');
        for (j=1:masknum)
            fprintf(spfifanoa2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfifanoa2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfifanoa2, '\r');
        end
        fprintf(spfifanoa2, '\r');
    end
    fclose(spfifanoa1);
    fclose(spfifanoa2);
    
    %% Output of Coefficient of Variation (CV) for raw data
    spficva1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw1_CV', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spficva2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_raw2_CV', '.txt')), 'w');    

    % superpix_original1
    fprintf(spficva1, 'CV_raw value\t');
    fprintf(spficva1, '%f\t', xlabel);
    fprintf(spficva1, '\r');
    for (j=1:masknum)
        fprintf(spficva1, '%s\t', spmaskname(j, :));
        fprintf(spficva1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficva1, 'fframe\t');
            fprintf(spficva1, 'sumrange\t');
        end
        fprintf(spficva1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spficva1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spficva1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spficva1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spficva1, '%6.5f\t', spnoisecv(f, j, i));			% original value
            end
            fprintf(spficva1, '\r');
        end
        fprintf(spficva1, '\r');
    end

    fprintf(spficva1, '\r\r_____Num_____\r\t');
    fprintf(spficva1, '%f\t', xlabel);
    fprintf(spficva1, '\r');
    for j=1:masknum
        fprintf(spficva1, '%s\t', spmaskname(j, :));
        fprintf(spficva1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficva1, 'fframe\t');
            fprintf(spficva1, 'sumrange\t');
        end
        fprintf(spficva1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spficva1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spficva1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spficva1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spficva1, '%d\t', stimsumnum(i));
            end
            fprintf(spficva1, '\r');
        end
        fprintf(spficva1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spficva2, 'CV_raw value\t');
    fprintf(spficva2, '%f\t', xlabel);
    fprintf(spficva2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spficva2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spficva2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spficva2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spficva2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficva2, 'fframe\t');
            fprintf(spficva2, 'sumrange\t');
        end
        fprintf(spficva2, '\r');
        for (j=1:masknum)
            fprintf(spficva2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spficva2, '%6.5f\t', spnoisecv(f, j, i));			% original value
            end
            fprintf(spficva2, '\r');
        end
        fprintf(spficva2, '\r');
    end
    fprintf(spficva2, '\r\r_____Num_____\r\t');
    fprintf(spficva2, '%f\t', xlabel);
    fprintf(spficva2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spficva2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spficva2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spficva2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spficva2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficva2, 'fframe\t');
            fprintf(spficva2, 'sumrange\t');
        end
        fprintf(spficva2, '\r');
        for (j=1:masknum)
            fprintf(spficva2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spficva2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spficva2, '\r');
        end
        fprintf(spficva2, '\r');
    end
    fclose(spficva1);
    fclose(spficva2);
    
end

if spoperation==2 || spoperation==3
    spvaluesad=zeros(framenum, masknum, NStim+NSPcond, blockselectnum);
    spvaluesnad=zeros(framenum, masknum, NStim+NSPcond, blockselectnum);
    spvaluesds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared
    spvaluesnds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared, devided by mean
    spvaluesnnds=zeros(framenum, masknum, NStim+NSPcond, blockselectnum); % deviation squared, devided by mean^2

    %% All trials
    %% Deviation for dR/R
	for (j=1:masknum)
		spfidalldb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR_data_all_Deviation_', spmaskname(j, :), '.txt')), 'w');
	%         spfidalldb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidalldb, 'Dev_dRR value\t');
		fprintf(spfidalldb, '%f\t', xlabel);
		fprintf(spfidalldb, '\r');
		for (i=1:NStim)
			fprintf(spfidalldb, 'stim%d\t', i);
			fprintf(spfidalldb, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidalldb, 'fframe\t');
				fprintf(spfidalldb, 'sumrange\t');
			end
	%             fprintf(spfidalldb, 'gs%d\t');
			fprintf(spfidalldb, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidalldb, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
						fprintf(spfidalldb, '%3.8f\t', spvaluesdrr(f, j, i, counter)-mean(spvaluesdrr(f, j, i, 1:stimsumnum(i))));			% dR/R value
                    end
	%                     fprintf(spfidalldb, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidalldb, '\t');
					end
	%                     fprintf(spfidalldb, '%d', 0);
				end
				fprintf(spfidalldb, '\r');
			end  
			fprintf(spfidalldb, '\r');
		end
		fclose(spfidalldb);
    end
    
    %% Absolute Deviation for dR/R
	for (j=1:masknum)
		spfidalladb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR_data_all_AbsoluteDeviation_', spmaskname(j, :), '.txt')), 'w');
	%         spfidalladb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidalladb, 'AD_dRR value\t');
		fprintf(spfidalladb, '%f\t', xlabel);
		fprintf(spfidalladb, '\r');
		for (i=1:NStim)
			fprintf(spfidalladb, 'stim%d\t', i);
			fprintf(spfidalladb, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidalladb, 'fframe\t');
				fprintf(spfidalladb, 'sumrange\t');
			end
	%             fprintf(spfidalladb, 'gs%d\t');
			fprintf(spfidalladb, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidalladb, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
						fprintf(spfidalladb, '%3.8f\t', abs(spvaluesdrr(f, j, i, counter)-mean(spvaluesdrr(f, j, i, 1:stimsumnum(i)))));			% dR/R value
                        spvaluesad(f, j, i, counter)=abs(spvaluesdrr(f, j, i, counter)-mean(spvaluesdrr(f, j, i, 1:stimsumnum(i))));
                        spvaluesds(f, j, i, counter)=(spvaluesdrr(f, j, i, counter)-mean(spvaluesdrr(f, j, i, 1:stimsumnum(i))))^2; % deviation squared
                    end
	%                     fprintf(spfidalladb, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidalladb, '\t');
					end
	%                     fprintf(spfidalladb, '%d', 0);
				end
				fprintf(spfidalladb, '\r');
			end  
			fprintf(spfidalladb, '\r');
		end
		fclose(spfidalladb);
	end
	
	%% Normalized Absolute Deviation for dR/R
	for (j=1:masknum)
		spfidallnadb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR_data_all_NormalizedAbsoluteDeviation_', spmaskname(j, :), '.txt')), 'w');
	%         spfidallnadb= fopen(strcat(runfolder, SPresultname, strcat('superpixel_raw_all.txt')), 'w'); 
		fprintf(spfidallnadb, 'NAD_dRR value\t');
		fprintf(spfidallnadb, '%f\t', xlabel);
		fprintf(spfidallnadb, '\r');
		for (i=1:NStim)
			fprintf(spfidallnadb, 'stim%d\t', i);
			fprintf(spfidallnadb, 'f%d\t', 1:FramesPerStim);
			if flagsaveframeave
				fprintf(spfidallnadb, 'fframe\t');
				fprintf(spfidallnadb, 'sumrange\t');
			end
	%             fprintf(spfidallnadb, 'gs%d\t');
			fprintf(spfidallnadb, '\r');
			counter=0;
			for k=blockselect  % start block by block processing
				fprintf(spfidallnadb, 'block%d\t', k);
				if goodstim(k, i)~=0
					counter=counter+1;
					for (f=1:framenum)
                        tempmean1=mean(spvaluesdrr(f, j, i, 1:stimsumnum(i)));
                        if tempmean1~=0
                            fprintf(spfidallnadb, '%3.8f\t', abs(spvaluesdrr(f, j, i, counter)-tempmean1)/tempmean1);			% dR/R value
                            spvaluesnad(f, j, i, counter)=abs(spvaluesdrr(f, j, i, counter)-tempmean1)/tempmean1;
                            spvaluesnds(f, j, i, counter)=((spvaluesdrr(f, j, i, counter)-tempmean1)^2)/tempmean1; % deviation squared, devided by mean
                            spvaluesnnds(f, j, i, counter)=((spvaluesdrr(f, j, i, counter)-tempmean1)^2)/tempmean1/tempmean1; % deviation squared, devided by mean ^2
                        else
                            fprintf(spfidallnadb, '%3.8f\t', 0);
                            spvaluesnad(f, j, i, counter)=0;
                            spvaluesnds(f, j, i, counter)=0;
                            spvaluesnnds(f, j, i, counter)=0;
                        end
					end
	%                     fprintf(spfidallnadb, '%d', 1);
				else
					for (f=1:framenum)
						fprintf(spfidallnadb, '\t');
					end
	%                     fprintf(spfidallnadb, '%d', 0);
				end
				fprintf(spfidallnadb, '\r');
			end  
			fprintf(spfidallnadb, '\r');
		end
		fclose(spfidallnadb);
    end
   
    %% Average, SD, Variance, Fano factor, CV

    spnoiseadavg=zeros(framenum, masknum, NStim+NSPcond);
    spnoisenadavg=zeros(framenum, masknum, NStim+NSPcond);
    
    spnoiseadstd=zeros(framenum, masknum, NStim+NSPcond);
    spnoisenadstd=zeros(framenum, masknum, NStim+NSPcond);

    spnoisevar=zeros(framenum, masknum, NStim+NSPcond); % variance
    spnoisefano=zeros(framenum, masknum, NStim+NSPcond); % fano factor?
    spnoisecv=zeros(framenum, masknum, NStim+NSPcond); % coefficient of variation
    
    for i=1:NStim
        if stimsumnum(i)>0
            spnoiseadavg(:,:,i)=sum(spvaluesad(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
            spnoisenadavg(:,:,i)=sum(spvaluesnad(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different

            spnoiseadstd(:,:,i)=std(spvaluesad(:,:,i,1:stimsumnum(i)),0,4);
            spnoisenadstd(:,:,i)=std(spvaluesnad(:,:,i,1:stimsumnum(i)),0,4);
            
            spnoisevar(:,:,i)=sum(spvaluesds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1); % unbiased variance
            spnoisefano(:,:,i)=sum(spvaluesnds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1); % unbiased variance devided by mean, fano factor?
            spnoisecv(:,:,i)=sqrt(sum(spvaluesnnds(:,:,i,1:stimsumnum(i)),4)/(stimsumnum(i)-1)); % standard deviation devided by mean, coefficient of variation (CV)
       end
    end

    %% Absolute Deviation for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesad(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoiseadavg(:,:,j)=SPcondsum/SPcondsumnum;
                spnoiseadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    
    %% Normalized Absolute Deviation for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnad(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisenadavg(:,:,j)=SPcondsum/SPcondsumnum;
                spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    %% Variance for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisevar(:,:,j)=SPcondsum/(SPcondsumnum-1); % unbiased variance
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    %% Fano factor for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisefano(:,:,j)=SPcondsum/(SPcondsumnum-1);
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    %% Coefficient of Variation (CV) for SPcond
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
%             SPsubblockcounter=0;
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum);
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
%                         squeeze(spvalues(:,:,i,1:stimsumnum(i)));
%                         SPcondtemp1=SPcondtemp1+squeeze(spvalues(:,:,m,cumgs));
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesnnds(:,:,m,cumgs));
                    end
                end
                if sumtype==1 || sumtype==3 % 150225 HT
                    if SPcondcounter1>0
                        if SPcondcounter2>0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                elseif sumtype==2
                    if SPcondcounter1==size(condlist1,2)
                        if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
                            end                            
                                SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                                SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        elseif size(condlist2,2)==0
                            for n=1:SPcondcounter1
                                SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                            end                          
                            SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
                            SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end   
    %                     SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
    %                     SPcondsumnum=SPcondsumnum+SPcondcounter1;
                    end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
                clear SPcondtemp1;
                clear SPcondtemp2;
    %             end        		
            end % for k=blockselect
%             stimsumnum(j)=SPsubblockcounter;
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
%                 spavg(:,:,j)=SPcondsum/SPsubblockcounter;
                spnoisecv(:,:,j)=sqrt(SPcondsum/(SPcondsumnum-1));
%                 spnoisenadstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
            end
            clear SPcondtemp;
        end
    end
    
    %% Output of Absolute Deviation for raw data
    spfiadb1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR1_AbsoluteDeviation', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfiadb2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR2_AbsoluteDeviation', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfiadb1, 'AD_dRR value\t');
    fprintf(spfiadb1, '%f\t', xlabel);
    fprintf(spfiadb1, '\r');
    for (j=1:masknum)
        fprintf(spfiadb1, '%s\t', spmaskname(j, :));
        fprintf(spfiadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb1, 'fframe\t');
            fprintf(spfiadb1, 'sumrange\t');
        end
        fprintf(spfiadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfiadb1, '%6.5f\t', spnoiseadavg(f, j, i));			% original value
            end
            fprintf(spfiadb1, '\r');
        end
        fprintf(spfiadb1, '\r');
    end
    fprintf(spfiadb1, '\r\r_____SD_____\r\t');
    fprintf(spfiadb1, '%f\t', xlabel);
    fprintf(spfiadb1, '\r');
    for (j=1:masknum)
        fprintf(spfiadb1, '%s\t', spmaskname(j, :));
        fprintf(spfiadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb1, 'fframe\t');
            fprintf(spfiadb1, 'sumrange\t');
        end
        fprintf(spfiadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfiadb1, '%6.5f\t', spnoiseadstd(f, j, i));
            end
            fprintf(spfiadb1, '\r');
        end
        fprintf(spfiadb1, '\r');
    end
    fprintf(spfiadb1, '\r\r_____Num_____\r\t');
    fprintf(spfiadb1, '%f\t', xlabel);
    fprintf(spfiadb1, '\r');
    for j=1:masknum
        fprintf(spfiadb1, '%s\t', spmaskname(j, :));
        fprintf(spfiadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb1, 'fframe\t');
            fprintf(spfiadb1, 'sumrange\t');
        end
        fprintf(spfiadb1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfiadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfiadb1, '%d\t', stimsumnum(i));
            end
            fprintf(spfiadb1, '\r');
        end
        fprintf(spfiadb1, '\r');
    end
    fprintf(spfiadb1, '\r\r_____SEM_____\r\t');
    fprintf(spfiadb1, '%f\t', xlabel);
    fprintf(spfiadb1, '\r');
    for (j=1:masknum)
        fprintf(spfiadb1, '%s\t', spmaskname(j, :));
        fprintf(spfiadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb1, 'fframe\t');
            fprintf(spfiadb1, 'sumrange\t');
        end
        fprintf(spfiadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfiadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfiadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if stimsumnum(i)>0
                    fprintf(spfiadb1, '%6.5f\t', spnoiseadstd(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfiadb1, '\r');
        end
        fprintf(spfiadb1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfiadb2, 'AD_dRR value\t');
    fprintf(spfiadb2, '%f\t', xlabel);
    fprintf(spfiadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb2, 'fframe\t');
            fprintf(spfiadb2, 'sumrange\t');
        end
        fprintf(spfiadb2, '\r');
        for (j=1:masknum)
            fprintf(spfiadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfiadb2, '%6.5f\t', spnoiseadavg(f, j, i));			% original value
            end
            fprintf(spfiadb2, '\r');
        end
        fprintf(spfiadb2, '\r');
    end
    fprintf(spfiadb2, '\r\r_____SD_____\r\t');
    fprintf(spfiadb2, '%f\t', xlabel);
    fprintf(spfiadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb2, 'fframe\t');
            fprintf(spfiadb2, 'sumrange\t');
        end
        fprintf(spfiadb2, '\r');
        for (j=1:masknum)
            fprintf(spfiadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfiadb2, '%6.5f\t', spnoiseadstd(f, j, i));
                end
            end
            fprintf(spfiadb2, '\r');
        end
        fprintf(spfiadb2, '\r');
    end
    fprintf(spfiadb2, '\r\r_____Num_____\r\t');
    fprintf(spfiadb2, '%f\t', xlabel);
    fprintf(spfiadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb2, 'fframe\t');
            fprintf(spfiadb2, 'sumrange\t');
        end
        fprintf(spfiadb2, '\r');
        for (j=1:masknum)
            fprintf(spfiadb2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfiadb2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfiadb2, '\r');
        end
        fprintf(spfiadb2, '\r');
    end
    fprintf(spfiadb2, '\r\r_____SEM_____\r\t');
    fprintf(spfiadb2, '%f\t', xlabel);
    fprintf(spfiadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfiadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfiadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfiadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfiadb2, 'fframe\t');
            fprintf(spfiadb2, 'sumrange\t');
        end
        fprintf(spfiadb2, '\r');
        for (j=1:masknum)
            fprintf(spfiadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if stimsumnum(i)>0
                        fprintf(spfiadb2, '%6.5f\t', spnoiseadstd(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfiadb2, '\r');
        end
        fprintf(spfiadb2, '\r');
    end
    fclose(spfiadb1);
    fclose(spfiadb2);

    %% Output of Normalized Absolute Deviation for raw data
    spfinadb1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR1_NormalizedAbsoluteDeviation', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfinadb2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR2_NormalizedAbsoluteDeviation', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfinadb1, 'NAD_dRR value\t');
    fprintf(spfinadb1, '%f\t', xlabel);
    fprintf(spfinadb1, '\r');
    for (j=1:masknum)
        fprintf(spfinadb1, '%s\t', spmaskname(j, :));
        fprintf(spfinadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb1, 'fframe\t');
            fprintf(spfinadb1, 'sumrange\t');
        end
        fprintf(spfinadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfinadb1, '%6.5f\t', spnoisenadavg(f, j, i));			% original value
            end
            fprintf(spfinadb1, '\r');
        end
        fprintf(spfinadb1, '\r');
    end
    fprintf(spfinadb1, '\r\r_____SD_____\r\t');
    fprintf(spfinadb1, '%f\t', xlabel);
    fprintf(spfinadb1, '\r');
    for (j=1:masknum)
        fprintf(spfinadb1, '%s\t', spmaskname(j, :));
        fprintf(spfinadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb1, 'fframe\t');
            fprintf(spfinadb1, 'sumrange\t');
        end
        fprintf(spfinadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfinadb1, '%6.5f\t', spnoisenadstd(f, j, i));
            end
            fprintf(spfinadb1, '\r');
        end
        fprintf(spfinadb1, '\r');
    end
    fprintf(spfinadb1, '\r\r_____Num_____\r\t');
    fprintf(spfinadb1, '%f\t', xlabel);
    fprintf(spfinadb1, '\r');
    for j=1:masknum
        fprintf(spfinadb1, '%s\t', spmaskname(j, :));
        fprintf(spfinadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb1, 'fframe\t');
            fprintf(spfinadb1, 'sumrange\t');
        end
        fprintf(spfinadb1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfinadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfinadb1, '%d\t', stimsumnum(i));
            end
            fprintf(spfinadb1, '\r');
        end
        fprintf(spfinadb1, '\r');
    end
    fprintf(spfinadb1, '\r\r_____SEM_____\r\t');
    fprintf(spfinadb1, '%f\t', xlabel);
    fprintf(spfinadb1, '\r');
    for (j=1:masknum)
        fprintf(spfinadb1, '%s\t', spmaskname(j, :));
        fprintf(spfinadb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb1, 'fframe\t');
            fprintf(spfinadb1, 'sumrange\t');
        end
        fprintf(spfinadb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfinadb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfinadb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if stimsumnum(i)>0
                    fprintf(spfinadb1, '%6.5f\t', spnoisenadstd(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfinadb1, '\r');
        end
        fprintf(spfinadb1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfinadb2, 'NAD_dRR value\t');
    fprintf(spfinadb2, '%f\t', xlabel);
    fprintf(spfinadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb2, 'fframe\t');
            fprintf(spfinadb2, 'sumrange\t');
        end
        fprintf(spfinadb2, '\r');
        for (j=1:masknum)
            fprintf(spfinadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfinadb2, '%6.5f\t', spnoisenadavg(f, j, i));			% original value
            end
            fprintf(spfinadb2, '\r');
        end
        fprintf(spfinadb2, '\r');
    end
    fprintf(spfinadb2, '\r\r_____SD_____\r\t');
    fprintf(spfinadb2, '%f\t', xlabel);
    fprintf(spfinadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb2, 'fframe\t');
            fprintf(spfinadb2, 'sumrange\t');
        end
        fprintf(spfinadb2, '\r');
        for (j=1:masknum)
            fprintf(spfinadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfinadb2, '%6.5f\t', spnoisenadstd(f, j, i));
                end
            end
            fprintf(spfinadb2, '\r');
        end
        fprintf(spfinadb2, '\r');
    end
    fprintf(spfinadb2, '\r\r_____Num_____\r\t');
    fprintf(spfinadb2, '%f\t', xlabel);
    fprintf(spfinadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb2, 'fframe\t');
            fprintf(spfinadb2, 'sumrange\t');
        end
        fprintf(spfinadb2, '\r');
        for (j=1:masknum)
            fprintf(spfinadb2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfinadb2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfinadb2, '\r');
        end
        fprintf(spfinadb2, '\r');
    end
    fprintf(spfinadb2, '\r\r_____SEM_____\r\t');
    fprintf(spfinadb2, '%f\t', xlabel);
    fprintf(spfinadb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfinadb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfinadb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfinadb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfinadb2, 'fframe\t');
            fprintf(spfinadb2, 'sumrange\t');
        end
        fprintf(spfinadb2, '\r');
        for (j=1:masknum)
            fprintf(spfinadb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if stimsumnum(i)>0
                        fprintf(spfinadb2, '%6.5f\t', spnoisenadstd(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfinadb2, '\r');
        end
        fprintf(spfinadb2, '\r');
    end
    fclose(spfinadb1);
    fclose(spfinadb2);
    
    %% Output of Variance for raw data
    spfivarb1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR1_Variance', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfivarb2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR2_Variance', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfivarb1, 'Var_dRR value\t');
    fprintf(spfivarb1, '%f\t', xlabel);
    fprintf(spfivarb1, '\r');
    for (j=1:masknum)
        fprintf(spfivarb1, '%s\t', spmaskname(j, :));
        fprintf(spfivarb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivarb1, 'fframe\t');
            fprintf(spfivarb1, 'sumrange\t');
        end
        fprintf(spfivarb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfivarb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfivarb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfivarb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfivarb1, '%6.5f\t', spnoisevar(f, j, i));			% original value
            end
            fprintf(spfivarb1, '\r');
        end
        fprintf(spfivarb1, '\r');
    end

    fprintf(spfivarb1, '\r\r_____Num_____\r\t');
    fprintf(spfivarb1, '%f\t', xlabel);
    fprintf(spfivarb1, '\r');
    for j=1:masknum
        fprintf(spfivarb1, '%s\t', spmaskname(j, :));
        fprintf(spfivarb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivarb1, 'fframe\t');
            fprintf(spfivarb1, 'sumrange\t');
        end
        fprintf(spfivarb1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfivarb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfivarb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfivarb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfivarb1, '%d\t', stimsumnum(i));
            end
            fprintf(spfivarb1, '\r');
        end
        fprintf(spfivarb1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfivarb2, 'Var_dRR value\t');
    fprintf(spfivarb2, '%f\t', xlabel);
    fprintf(spfivarb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfivarb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfivarb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfivarb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfivarb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivarb2, 'fframe\t');
            fprintf(spfivarb2, 'sumrange\t');
        end
        fprintf(spfivarb2, '\r');
        for (j=1:masknum)
            fprintf(spfivarb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfivarb2, '%6.5f\t', spnoisevar(f, j, i));			% original value
            end
            fprintf(spfivarb2, '\r');
        end
        fprintf(spfivarb2, '\r');
    end
    fprintf(spfivarb2, '\r\r_____Num_____\r\t');
    fprintf(spfivarb2, '%f\t', xlabel);
    fprintf(spfivarb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfivarb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfivarb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfivarb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfivarb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfivarb2, 'fframe\t');
            fprintf(spfivarb2, 'sumrange\t');
        end
        fprintf(spfivarb2, '\r');
        for (j=1:masknum)
            fprintf(spfivarb2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfivarb2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfivarb2, '\r');
        end
        fprintf(spfivarb2, '\r');
    end
    fclose(spfivarb1);
    fclose(spfivarb2);
    
    %% Fano factor for raw data
    spfifanob1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR1_Fano', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfifanob2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR2_Fano', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfifanob1, 'Fano_dRR value\t');
    fprintf(spfifanob1, '%f\t', xlabel);
    fprintf(spfifanob1, '\r');
    for (j=1:masknum)
        fprintf(spfifanob1, '%s\t', spmaskname(j, :));
        fprintf(spfifanob1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanob1, 'fframe\t');
            fprintf(spfifanob1, 'sumrange\t');
        end
        fprintf(spfifanob1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfifanob1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfifanob1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfifanob1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfifanob1, '%6.5f\t', spnoisefano(f, j, i));			% original value
            end
            fprintf(spfifanob1, '\r');
        end
        fprintf(spfifanob1, '\r');
    end

    fprintf(spfifanob1, '\r\r_____Num_____\r\t');
    fprintf(spfifanob1, '%f\t', xlabel);
    fprintf(spfifanob1, '\r');
    for j=1:masknum
        fprintf(spfifanob1, '%s\t', spmaskname(j, :));
        fprintf(spfifanob1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanob1, 'fframe\t');
            fprintf(spfifanob1, 'sumrange\t');
        end
        fprintf(spfifanob1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfifanob1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfifanob1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfifanob1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfifanob1, '%d\t', stimsumnum(i));
            end
            fprintf(spfifanob1, '\r');
        end
        fprintf(spfifanob1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfifanob2, 'Fano_dRR value\t');
    fprintf(spfifanob2, '%f\t', xlabel);
    fprintf(spfifanob2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfifanob2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfifanob2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfifanob2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfifanob2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanob2, 'fframe\t');
            fprintf(spfifanob2, 'sumrange\t');
        end
        fprintf(spfifanob2, '\r');
        for (j=1:masknum)
            fprintf(spfifanob2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfifanob2, '%6.5f\t', spnoisefano(f, j, i));			% original value
            end
            fprintf(spfifanob2, '\r');
        end
        fprintf(spfifanob2, '\r');
    end
    fprintf(spfifanob2, '\r\r_____Num_____\r\t');
    fprintf(spfifanob2, '%f\t', xlabel);
    fprintf(spfifanob2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfifanob2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfifanob2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfifanob2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfifanob2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfifanob2, 'fframe\t');
            fprintf(spfifanob2, 'sumrange\t');
        end
        fprintf(spfifanob2, '\r');
        for (j=1:masknum)
            fprintf(spfifanob2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfifanob2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfifanob2, '\r');
        end
        fprintf(spfifanob2, '\r');
    end
    fclose(spfifanob1);
    fclose(spfifanob2);
    
    %% Coefficient of Variation (CV) for raw data
    spficvb1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR1_CV', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spficvb2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_noise_analysis\', 'superpixel_dRR2_CV', '.txt')), 'w');    

    % superpix_original1
    fprintf(spficvb1, 'CV_dRR value\t');
    fprintf(spficvb1, '%f\t', xlabel);
    fprintf(spficvb1, '\r');
    for (j=1:masknum)
        fprintf(spficvb1, '%s\t', spmaskname(j, :));
        fprintf(spficvb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficvb1, 'fframe\t');
            fprintf(spficvb1, 'sumrange\t');
        end
        fprintf(spficvb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spficvb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spficvb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spficvb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spficvb1, '%6.5f\t', spnoisecv(f, j, i));			% original value
            end
            fprintf(spficvb1, '\r');
        end
        fprintf(spficvb1, '\r');
    end

    fprintf(spficvb1, '\r\r_____Num_____\r\t');
    fprintf(spficvb1, '%f\t', xlabel);
    fprintf(spficvb1, '\r');
    for j=1:masknum
        fprintf(spficvb1, '%s\t', spmaskname(j, :));
        fprintf(spficvb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficvb1, 'fframe\t');
            fprintf(spficvb1, 'sumrange\t');
        end
        fprintf(spficvb1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spficvb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spficvb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spficvb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spficvb1, '%d\t', stimsumnum(i));
            end
            fprintf(spficvb1, '\r');
        end
        fprintf(spficvb1, '\r');
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spficvb2, 'CV_dRR value\t');
    fprintf(spficvb2, '%f\t', xlabel);
    fprintf(spficvb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spficvb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spficvb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spficvb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spficvb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficvb2, 'fframe\t');
            fprintf(spficvb2, 'sumrange\t');
        end
        fprintf(spficvb2, '\r');
        for (j=1:masknum)
            fprintf(spficvb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spficvb2, '%6.5f\t', spnoisecv(f, j, i));			% original value
            end
            fprintf(spficvb2, '\r');
        end
        fprintf(spficvb2, '\r');
    end
    fprintf(spficvb2, '\r\r_____Num_____\r\t');
    fprintf(spficvb2, '%f\t', xlabel);
    fprintf(spficvb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spficvb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spficvb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spficvb2, '%s\t',
%             getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spficvb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spficvb2, 'fframe\t');
            fprintf(spficvb2, 'sumrange\t');
        end
        fprintf(spficvb2, '\r');
        for (j=1:masknum)
            fprintf(spficvb2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spficvb2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spficvb2, '\r');
        end
        fprintf(spficvb2, '\r');
    end
    fclose(spficvb1);
    fclose(spficvb2);
    
end