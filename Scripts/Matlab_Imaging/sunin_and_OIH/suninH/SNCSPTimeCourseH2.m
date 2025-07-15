% This file is a section of Sunin code for OI data process
% Purpose: Super Pixel time course output
% 
% Ver 3.0 (150224) sumtype = 3 was introduced. Subtraction will be done when at least one of trials is success in all blocks. Just for check
% Ver 1.0 (081125) Separated from Core program by Hisashi
%
% Originally written by Hisashi

% spvalues=reshape(spvalues, FramesPerStim+2,masknum,NStim,max(stimsumnum));     % don't know why this is necessory

flagspmedian=1;

if flagsavespvalues==1 && flagloadspvalues~=1
    save(strcat(resultfolder, 'spvalues.mat'), 'spvalues');
    save(strcat(resultfolder, 'stimsumnum.mat'), 'stimsumnum');
    save(strcat(resultfolder, 'goodstim.mat'), 'goodstim');
end
if flagloadspvalues==1
    load(strcat(resultfolder, 'spvalues.mat'));
    load(strcat(resultfolder, 'stimsumnum.mat'));
    load(strcat(resultfolder, 'goodstim.mat'));
end

if ~isdir([runfolder, SPresultname])
    mkdir(runfolder, SPresultname);    
end
if ~isdir([runfolder, SPresultname, 'superpixel_timecourse'])
    mkdir([runfolder, SPresultname], 'superpixel_timecourse');    
end

xlabel=(1:FramesPerStim)/FrameRate;     % in unit of second
% if flagspcond==1
%     NSPcond=size(SPcond, 1);
% else
%     NSPcond=0;
% end
% if flagspcond==1
%     NSPcond=size(SPcond, 1);
% else
%     NSPcond=0;
% end

if spoperation==1 || spoperation==3
    spavg=zeros(framenum, masknum, NStim+NSPcond);
    spstd=zeros([framenum, masknum, NStim+NSPcond]);	
	if flagspmedian
		spmed=zeros(framenum, masknum, NStim+NSPcond);
	end
    
    for i=1:NStim
        if stimsumnum(i)>0
            spavg(:,:,i)=sum(spvalues(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
        end
    end
    for i=1:NStim
        if stimsumnum(i)>0
            spstd(:,:,i)=std(spvalues(:,:,i,1:stimsumnum(i)),0,4);
        end
    end
    if flagspmedian
		for i=1:NStim
			if stimsumnum(i)>0
				spmed(:,:,i)=median(spvalues(:,:,i,1:stimsumnum(i)),4);     % average over blocks, note: if use goodstim, then number of sum maybe different
			end
		end
    end



%% SP for SPcond by Hisashi on 081003 and 090129
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            if exist('flagSPblanksubtraction','var') == 0
                flagSPblanksubtraction=0;
            end
            if flagSPblanksubtraction
                condlistBSub1=getfield(cell2struct(SPcondBSub(j-NStim, 2), 'junk'), 'junk');
                condlistBSub2=getfield(cell2struct(SPcondBSub(j-NStim, 3), 'junk'), 'junk');
            end
            
            if flagSPblanksubtraction == 2 % 150220 HT
				SPcondcounterBSub1=0;
				SPcondcounterBSub2=0;
				SPcondtempBSub1=zeros(framenum,masknum);
				SPcondtempBSub2=zeros(framenum,masknum);
				for k=blockselect  % start block by block processing                	
                    for m=condlistBSub1
                        if  goodstim(k, m)==1
                            SPcondcounterBSub1=SPcondcounterBSub1+1;
                            cumgs=0;
                            for l=blockselect
                                if l<=k && goodstim(l,m)==1
                                    cumgs=cumgs+1;
                                end
                            end
                            SPcondtempBSub1=SPcondtempBSub1+squeeze(spvalues(:,:,m,cumgs));
                        end
                    end                    
                    for m=condlistBSub2
                        if  goodstim(k, m)==1
                            SPcondcounterBSub2=SPcondcounterBSub2+1;
                            cumgs=0;
                            for l=blockselect
                                if l<=k && goodstim(l,m)==1
                                    cumgs=cumgs+1;
                                end
                            end
                            SPcondtempBSub2=SPcondtempBSub2+squeeze(spvalues(:,:,m,cumgs));
                        end
                    end
				end % for k=blockselect
				spavgBsub1=SPcondtempBSub1/SPcondcounterBSub1;
				spavgBsub2=SPcondtempBSub2/SPcondcounterBSub2;					
			end % if flagSPblanksubtraction == 2
			
            for k=blockselect  % start block by block processing
                SPcondcounter1=0;
                SPcondcounter2=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                SPcondtemp2=zeros(framenum,masknum); %100419 HT
                if flagSPblanksubtraction == 1 %110301 HT
                	SPcondcounterBSub1=0;
                	SPcondcounterBSub2=0;
                	SPcondtempBSub1=zeros(framenum,masknum);
                	SPcondtempBSub2=zeros(framenum,masknum);
                	
                    for m=condlistBSub1
                        if  goodstim(k, m)==1
                            SPcondcounterBSub1=SPcondcounterBSub1+1;
                            cumgs=0;
                            for l=blockselect
                                if l<=k && goodstim(l,m)==1
                                    cumgs=cumgs+1;
                                end
                            end
                            SPcondtempBSub1=SPcondtempBSub1+squeeze(spvalues(:,:,m,cumgs));
                        end
                    end
                    spavgBsub1=SPcondtempBSub1/SPcondcounterBSub1;
                    for m=condlistBSub2
                        if  goodstim(k, m)==1
                            SPcondcounterBSub2=SPcondcounterBSub2+1;
                            cumgs=0;
                            for l=blockselect
                                if l<=k && goodstim(l,m)==1
                                    cumgs=cumgs+1;
                                end
                            end
                            SPcondtempBSub2=SPcondtempBSub2+squeeze(spvalues(:,:,m,cumgs));
                        end
                    end
                    spavgBsub2=SPcondtempBSub2/SPcondcounterBSub2;
                end % if flagSPblanksubtraction == 1
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvalues(:,:,m,cumgs));
                        if flagSPblanksubtraction %110301 HT
                        	SPcondtemp1(:,:,SPcondcounter1)=SPcondtemp1(:,:,SPcondcounter1)-spavgBsub1;
                        end                        
                    end
                    SPcondcounter1
                end
                for m=condlist2 %100419 HT                    
                    if  goodstim(k, m)==1
                        SPcondcounter2=SPcondcounter2+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
                        SPcondtemp2=SPcondtemp2+squeeze(spvalues(:,:,m,cumgs));
                        if flagSPblanksubtraction %110301 HT
                        	SPcondtemp2=SPcondtemp2-spavgBsub2;
                        end       
                    end
                     SPcondcounter2
               end
				if sumtype==1
					if flagSPblanksubtraction == 0 %110301 HT
						if SPcondcounter1>0
							if SPcondcounter2>0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
								end                            
							elseif size(condlist2,2)==0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
								end                          
							end   
							SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
							SPcondsumnum=SPcondsumnum+SPcondcounter1;
						end
					else
						if SPcondcounter1>0 && SPcondcounterBSub1>0
							if SPcondcounter2>0 && SPcondcounterBSub2>0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2;
								end                            
							elseif size(condlist2,2)==0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
								end                          
							end   
							SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
							SPcondsumnum=SPcondsumnum+SPcondcounter1;
						end
					end
				elseif sumtype==2
                    if flagSPblanksubtraction == 0 %110301 HT
                        if SPcondcounter1==size(condlist1,2)
                            if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                                for n=1:SPcondcounter1
                                    SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2; %need to be improved 100419 HT
                                end                            
                            elseif size(condlist2,2)==0
                                for n=1:SPcondcounter1
                                    SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
                                end                          
                            end   
							SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
							SPcondsumnum=SPcondsumnum+SPcondcounter1;
                        end
                    else
						if SPcondcounter1==size(condlist1,2) && SPcondcounterBSub1 > 0
							if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0 && SPcondcounterBSub2 >0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2; %need to be improved 100419 HT
								end                            
							elseif size(condlist2,2)==0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n);
								end                          
							end   
							SPcondsum=SPcondsum+sum(SPcondtemp(:,:,SPcondsumnum+1:SPcondsumnum+SPcondcounter1),3);
							SPcondsumnum=SPcondsumnum+SPcondcounter1;
						end
					end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)
            end % for k=blockselect
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
                spavg(:,:,j)=SPcondsum/SPcondsumnum;
                spstd(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
                if flagspmedian
                	spmed(:,:,j)=median(SPcondtemp(:,:,1:SPcondsumnum), 3); %100419 HT
				end					
            end
            clear SPcondtemp;
        end
    end



    %% Output

    spavg=reshape(spavg, [framenum, masknum, NStim+NSPcond]);
    spstd=reshape(spstd, [framenum, masknum, NStim+NSPcond]);



    spfida1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_timecourse\', 'superpixel_raw1', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfida2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_timecourse\', 'superpixel_raw2', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfida1, 'raw value\t');
    fprintf(spfida1, '%f\t', xlabel);
    fprintf(spfida1, '\r');
    for (j=1:masknum)
        fprintf(spfida1, '%s\t', spmaskname(j, :));
        fprintf(spfida1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida1, 'fframe\t');
            fprintf(spfida1, 'sumrange\t');
        end
        fprintf(spfida1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfida1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfida1, '%6.5f\t', spavg(f, j, i));			% original value
            end
            fprintf(spfida1, '\r');
        end
        fprintf(spfida1, '\r');
    end
    fprintf(spfida1, '\r\r_____Median_____\r\t');
    fprintf(spfida1, '%f\t', xlabel);
    fprintf(spfida1, '\r');
    for (j=1:masknum)
        fprintf(spfida1, '%s\t', spmaskname(j, :));
        fprintf(spfida1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida1, 'fframe\t');
            fprintf(spfida1, 'sumrange\t');
        end
        fprintf(spfida1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfida1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfida1, '%6.5f\t', spmed(f, j, i));
            end
            fprintf(spfida1, '\r');
        end
        fprintf(spfida1, '\r');
    end
    fprintf(spfida1, '\r\r_____SD_____\r\t');
    fprintf(spfida1, '%f\t', xlabel);
    fprintf(spfida1, '\r');
    for (j=1:masknum)
        fprintf(spfida1, '%s\t', spmaskname(j, :));
        fprintf(spfida1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida1, 'fframe\t');
            fprintf(spfida1, 'sumrange\t');
        end
        fprintf(spfida1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfida1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfida1, '%6.5f\t', spstd(f, j, i));
            end
            fprintf(spfida1, '\r');
        end
        fprintf(spfida1, '\r');
    end
    fprintf(spfida1, '\r\r_____Num_____\r\t');
    fprintf(spfida1, '%f\t', xlabel);
    fprintf(spfida1, '\r');
    for j=1:masknum
        fprintf(spfida1, '%s\t', spmaskname(j, :));
        fprintf(spfida1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida1, 'fframe\t');
            fprintf(spfida1, 'sumrange\t');
        end
        fprintf(spfida1, '\r');
        for i=1:NStim+NSPcond
            if i<=NStim
                fprintf(spfida1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for f=1:framenum
                fprintf(spfida1, '%d\t', stimsumnum(i));
            end
            fprintf(spfida1, '\r');
        end
        fprintf(spfida1, '\r');
    end
    fprintf(spfida1, '\r\r_____SEM_____\r\t');
    fprintf(spfida1, '%f\t', xlabel);
    fprintf(spfida1, '\r');
    for (j=1:masknum)
        fprintf(spfida1, '%s\t', spmaskname(j, :));
        fprintf(spfida1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida1, 'fframe\t');
            fprintf(spfida1, 'sumrange\t');
        end
        fprintf(spfida1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfida1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfida1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if stimsumnum(i)>0
                    fprintf(spfida1, '%6.5f\t', spstd(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfida1, '\r');
        end
        fprintf(spfida1, '\r');
    end



    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfida2, 'raw value\t');
    fprintf(spfida2, '%f\t', xlabel);
    fprintf(spfida2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfida2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfida2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida2, 'fframe\t');
            fprintf(spfida2, 'sumrange\t');
        end
        fprintf(spfida2, '\r');
        for (j=1:masknum)
            fprintf(spfida2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfida2, '%6.5f\t', spavg(f, j, i));			% original value
            end
            fprintf(spfida2, '\r');
        end
        fprintf(spfida2, '\r');
    end
    fprintf(spfida2, '\r\r_____Median_____\r\t');
    fprintf(spfida2, '%f\t', xlabel);
    fprintf(spfida2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfida2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfida2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida2, 'fframe\t');
            fprintf(spfida2, 'sumrange\t');
        end
        fprintf(spfida2, '\r');
        for (j=1:masknum)
            fprintf(spfida2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfida2, '%6.5f\t', spmed(f, j, i));
                end
            end
            fprintf(spfida2, '\r');
        end
        fprintf(spfida2, '\r');
    end
    fprintf(spfida2, '\r\r_____SD_____\r\t');
    fprintf(spfida2, '%f\t', xlabel);
    fprintf(spfida2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfida2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfida2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida2, 'fframe\t');
            fprintf(spfida2, 'sumrange\t');
        end
        fprintf(spfida2, '\r');
        for (j=1:masknum)
            fprintf(spfida2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfida2, '%6.5f\t', spstd(f, j, i));
                end
            end
            fprintf(spfida2, '\r');
        end
        fprintf(spfida2, '\r');
    end
    fprintf(spfida2, '\r\r_____Num_____\r\t');
    fprintf(spfida2, '%f\t', xlabel);
    fprintf(spfida2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfida2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfida2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida2, 'fframe\t');
            fprintf(spfida2, 'sumrange\t');
        end
        fprintf(spfida2, '\r');
        for (j=1:masknum)
            fprintf(spfida2, '%d', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfida2, '%6.5f\t', stimsumnum(i));
                end
            end
            fprintf(spfida2, '\r');
        end
        fprintf(spfida2, '\r');
    end
    fprintf(spfida2, '\r\r_____SEM_____\r\t');
    fprintf(spfida2, '%f\t', xlabel);
    fprintf(spfida2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfida2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfida2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfida2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfida2, 'fframe\t');
            fprintf(spfida2, 'sumrange\t');
        end
        fprintf(spfida2, '\r');
        for (j=1:masknum)
            fprintf(spfida2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if stimsumnum(i)>0
                        fprintf(spfida2, '%6.5f\t', spstd(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfida2, '\r');
        end
        fprintf(spfida2, '\r');
    end
    fclose(spfida1);
    fclose(spfida2);
end

%%% dRR 
if spoperation==2 || spoperation==3
    %     spfid3a= fopen(strcat(resultfolder, strcat('superpixeldRR3', '.txt')), 'w');  
    %     spvaluesdrr=spvalues; %commented out by Hisashi on 081112

%     if ~(flaghpfilter||flaglpfilter) ?? commented out on 150223 HT
% %         clear spvaluesdrr;
%         spvaluesdrr=spvalues;
% %         fprintf('dR/R\r');
% %         size(spvaluesdrr)
%         for i=1:NStim
%             for k=1:stimsumnum(i)
%                 for j=1:masknum
%                     spvaluesdrr(:, j, i, k)=(spvalues(:,j,i,k)-mean(spvalues(fframe,j,i,k),1))/mean(spvalues(fframe,j,i,k),1); %dR/R
%                 end
%             end
%         end
%     end

% % 	clear spvaluesdrr;
% 	spvaluesdrr=spvalues;
% %         fprintf('dR/R\r');
% %         size(spvaluesdrr)
% 	for i=1:NStim
% 		for k=1:stimsumnum(i)
% 			for j=1:masknum
% 				spvaluesdrr(:, j, i, k)=(spvalues(:,j,i,k)-mean(spvalues(fframe,j,i,k),1))/mean(spvalues(fframe,j,i,k),1); %dR/R
% 			end
% 		end
% 	end
    
	spavgdrr=zeros(framenum, masknum, NStim+NSPcond);
    spstddrr=zeros(framenum, masknum, NStim+NSPcond);
    if flagspmedian
    	spmeddrr=zeros(framenum, masknum, NStim+NSPcond);
	end

    
    for i=1:NStim
        if stimsumnum(i)>0
            spavgdrr(:,:,i)=sum(spvaluesdrr(:,:,i,1:stimsumnum(i)),4)/stimsumnum(i);     % average over blocks, note: if use goodstim, then number of sum maybe different
        end
    end
    for i=1:NStim
        if stimsumnum(i)>0
            spstddrr(:,:,i)=std(spvaluesdrr(:,:,i,1:stimsumnum(i)),0,4);
        end
    end
    if flagspmedian
		for i=1:NStim
			if stimsumnum(i)>0
				spmeddrr(:,:,i)=median(spvaluesdrr(:,:,i,1:stimsumnum(i)),4);
			end
		end
    end
 
    
%% SP for SPcond by Hisashi on 081003 and 090125
    if flagspcond==1
        for j=NStim+1:NStim+NSPcond
            condlist1=getfield(cell2struct(SPcond(j-NStim, 2), 'junk'), 'junk');
            condlist2=getfield(cell2struct(SPcond(j-NStim, 3), 'junk'), 'junk');
            SPcondtemp=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
            SPcondsum=zeros(framenum,masknum);
            SPcondsumnum=0;
            if exist('flagSPblanksubtraction','var') == 0
                flagSPblanksubtraction=0;
            end
            if flagSPblanksubtraction
                condlistBSub1=getfield(cell2struct(SPcondBSub(j-NStim, 2), 'junk'), 'junk');
                condlistBSub2=getfield(cell2struct(SPcondBSub(j-NStim, 3), 'junk'), 'junk');
            end

			if flagSPblanksubtraction == 2 %150220 HT
				SPcondcounterBSub1=0;
				SPcondcounterBSub2=0;
				SPcondtempBSub1=zeros(framenum,masknum);
				SPcondtempBSub2=zeros(framenum,masknum);
				for k=blockselect				
					for m=condlistBSub1
						if  goodstim(k, m)==1
							SPcondcounterBSub1=SPcondcounterBSub1+1;
							cumgs=0;
							for l=blockselect
								if l<=k && goodstim(l,m)==1
									cumgs=cumgs+1;
								end
							end
							SPcondtempBSub1=SPcondtempBSub1+squeeze(spvaluesdrr(:,:,m,cumgs));
						end
					end
					for m=condlistBSub2
						if  goodstim(k, m)==1
							SPcondcounterBSub2=SPcondcounterBSub2+1;
							cumgs=0;
							for l=blockselect
								if l<=k && goodstim(l,m)==1
									cumgs=cumgs+1;
								end
							end
							SPcondtempBSub2=SPcondtempBSub2+squeeze(spvaluesdrr(:,:,m,cumgs));
						end
					end
				end % for k=blockselect
				spavgBsub1=SPcondtempBSub1/SPcondcounterBSub1;
				spavgBsub2=SPcondtempBSub2/SPcondcounterBSub2;
			end % if flagSPblanksubtraction == 2

            if sumtype == 3 % 150224 HT
				if flagSPblanksubtraction == 1 %110301 HT
					for k=blockselect  % start block by block processing                	
						SPcondcounterBSub2=0;
						SPcondtempBSub2=zeros(framenum,masknum);
						spavgBsub1=SPcondtempBSub1/SPcondcounterBSub1;
						for m=condlistBSub2
							if  goodstim(k, m)==1
								SPcondcounterBSub2=SPcondcounterBSub2+1;
								cumgs=0;
								for l=blockselect
									if l<=k && goodstim(l,m)==1
										cumgs=cumgs+1;
									end
								end
								SPcondtempBSub2=SPcondtempBSub2+squeeze(spvaluesdrr(:,:,m,cumgs));
							end
						end
						spavgBsub2=SPcondtempBSub2/SPcondcounterBSub2;
					end % for k=blockselect
				end % if flagSPblanksubtraction == 1
				
                SPcondtemp2=zeros(framenum,masknum); %100419 HT
                SPcondcounter2=0;
                for m=condlist2 %100419 HT                    
                    if  goodstim(k, m)==1
                        SPcondcounter2=SPcondcounter2+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
                        SPcondtemp2=SPcondtemp2+squeeze(spvaluesdrr(:,:,m,cumgs));
                        if flagSPblanksubtraction %110301 HT
                        	SPcondtemp2=SPcondtemp2-spavgBsub2;
                        end       
                    end
				end
			end % if sumtype == 3

            for k=blockselect  % start block by block processing
                if flagSPblanksubtraction == 1 %110301 HT
                	SPcondcounterBSub1=0;
                	SPcondcounterBSub2=0;
                	SPcondtempBSub1=zeros(framenum,masknum);
                	SPcondtempBSub2=zeros(framenum,masknum);
                	
                    for m=condlistBSub1
                        if  goodstim(k, m)==1
                            SPcondcounterBSub1=SPcondcounterBSub1+1;
                            cumgs=0;
                            for l=blockselect
                                if l<=k && goodstim(l,m)==1
                                    cumgs=cumgs+1;
                                end
                            end
                            SPcondtempBSub1=SPcondtempBSub1+squeeze(spvaluesdrr(:,:,m,cumgs));
                        end
                    end
                    spavgBsub1=SPcondtempBSub1/SPcondcounterBSub1;

					if sumtyupe ~= 3
						SPcondcounter2=0;
						SPcondtemp2=zeros(framenum,masknum); %100419 HT
						for m=condlistBSub2
							if  goodstim(k, m)==1
								SPcondcounterBSub2=SPcondcounterBSub2+1;
								cumgs=0;
								for l=blockselect
									if l<=k && goodstim(l,m)==1
										cumgs=cumgs+1;
									end
								end
								SPcondtempBSub2=SPcondtempBSub2+squeeze(spvaluesdrr(:,:,m,cumgs));
							end
						end
						spavgBsub2=SPcondtempBSub2/SPcondcounterBSub2;
					end
                end % if flagSPblanksubtraction == 1
                
                SPcondcounter1=0;
                SPcondtemp1=zeros(framenum,masknum,size(condlist1,2));
                for m=condlist1
                    if  goodstim(k, m)==1
                        SPcondcounter1=SPcondcounter1+1;
                        cumgs=0;
                        for l=blockselect
                            if l<=k && goodstim(l,m)==1
                                cumgs=cumgs+1;
                            end
                        end
                        SPcondtemp1(:,:,SPcondcounter1)=squeeze(spvaluesdrr(:,:,m,cumgs));
                        if flagSPblanksubtraction %110301 HT
                        	SPcondtemp1(:,:,SPcondcounter1)=SPcondtemp1(:,:,SPcondcounter1)-spavgBsub1;
                        end                        
                    end
                end
				if sumtype ~=3
					SPcondcounter2=0;
					SPcondtemp2=zeros(framenum,masknum); %100419 HT
					for m=condlist2 %100419 HT
						if  goodstim(k, m)==1
							SPcondcounter2=SPcondcounter2+1;
							cumgs=0;
							for l=blockselect
								if l<=k && goodstim(l,m)==1
									cumgs=cumgs+1;
								end
							end
							SPcondtemp2=SPcondtemp2+squeeze(spvaluesdrr(:,:,m,cumgs));
							if flagSPblanksubtraction %110301 HT
								SPcondtemp2=SPcondtemp2-spavgBsub2;
							end       
						end
					end
				end
                if sumtype==1 || sumtype==3
                    if flagSPblanksubtraction == 0 %110301 HT
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
                        end
                    else
						if SPcondcounter1>0 && SPcondcounterBSub1>0
							if SPcondcounter2>0 && SPcondcounterBSub2>0
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
						end
					end
                elseif sumtype==2
                    if flagSPblanksubtraction == 0 %110301 HT
                        if SPcondcounter1==size(condlist1,2)
                            if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0
                                for n=1:SPcondcounter1
                                    SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2; %need to be improved 100419 HT
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
                        end
                    else
						if SPcondcounter1==size(condlist1,2) && SPcondcounterBSub1 > 0
							if SPcondcounter2==size(condlist2,2) && SPcondcounter2~=0 && SPcondcounterBSub2 >0
								for n=1:SPcondcounter1
									SPcondtemp(:,:,SPcondsumnum+n)=SPcondtemp1(:,:,n)-SPcondtemp2/SPcondcounter2; %need to be improved 100419 HT
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
						end
					end
                else
                    fprintf('\n  Error: specify sumtype correctly!\n');
                end
%                 sizecondlist2=size(condlist2,2)

            end % for k=blockselect
            stimsumnum(j)=SPcondsumnum;
            if SPcondsumnum>0
                spavgdrr(:,:,j)=SPcondsum/SPcondsumnum;
                spstddrr(:,:,j)=std(SPcondtemp(:,:,1:SPcondsumnum), 0, 3);
                if flagspmedian
                	spmeddrr(:,:,j)=median(SPcondtemp(:,:,1:SPcondsumnum), 3); %100419 HT
				end					
            end
            clear SPcondtemp;
        end
    end

    %% Output
    spavgdrr=reshape(spavgdrr, [framenum, masknum, NStim+NSPcond]);
    spstddrr=reshape(spstddrr, [framenum, masknum, NStim+NSPcond]);

    spfidb1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_timecourse\', 'superpixel_dRR1', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfidb2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_timecourse\', 'superpixel_dRR2', '.txt')), 'w');

    % superpix1 (same mask are grouped together)
    fprintf(spfidb1, 'dR/R value\t');
    fprintf(spfidb1, '%f\t', xlabel);
    fprintf(spfidb1, '\r');
    for (j=1:masknum)
        fprintf(spfidb1, '%s\t', spmaskname(j, :));
        fprintf(spfidb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb1, 'fframe\t');
            fprintf(spfidb1, 'sumrange\t');
        end
        fprintf(spfidb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfidb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfidb1, '%3.8f\t', spavgdrr(f, j, i));			% original value
            end
            fprintf(spfidb1, '\r');
        end
        fprintf(spfidb1, '\r');
    end
    fprintf(spfidb1, '\r\r_____Median_____\r\t');
    fprintf(spfidb1, '%f\t', xlabel);
    fprintf(spfidb1, '\r');
    for (j=1:masknum)
        fprintf(spfidb1, '%s\t', spmaskname(j, :));
        fprintf(spfidb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb1, 'fframe\t');
            fprintf(spfidb1, 'sumrange\t');
        end
        fprintf(spfidb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfidb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfidb1, '%3.8f\t', spmeddrr(f, j, i));
                end
            end
            fprintf(spfidb1, '\r');
        end
        fprintf(spfidb1, '\r');
    end
    fprintf(spfidb1, '\r\r_____SD_____\r\t');
    fprintf(spfidb1, '%f\t', xlabel);
    fprintf(spfidb1, '\r');
    for (j=1:masknum)
        fprintf(spfidb1, '%s\t', spmaskname(j, :));
        fprintf(spfidb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb1, 'fframe\t');
            fprintf(spfidb1, 'sumrange\t');
        end
        fprintf(spfidb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfidb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfidb1, '%3.8f\t', spstddrr(f, j, i));
                end
            end
            fprintf(spfidb1, '\r');
        end
        fprintf(spfidb1, '\r');
    end
    fprintf(spfidb1, '\r\r_____Num_____\r\t');
    fprintf(spfidb1, '%f\t', xlabel);
    fprintf(spfidb1, '\r');
    for (j=1:masknum)
        fprintf(spfidb1, '%s\t', spmaskname(j, :));
        fprintf(spfidb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb1, 'fframe\t');
            fprintf(spfidb1, 'sumrange\t');
        end
        fprintf(spfidb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfidb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                fprintf(spfidb1, '%d\t', stimsumnum(i));
            end
            fprintf(spfidb1, '\r');
        end
        fprintf(spfidb1, '\r');
    end
    fprintf(spfidb1, '\r\r_____SEM_____\r\t');
    fprintf(spfidb1, '%f\t', xlabel);
    fprintf(spfidb1, '\r');
    for (j=1:masknum)
        fprintf(spfidb1, '%s\t', spmaskname(j, :));
        fprintf(spfidb1, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb1, 'fframe\t');
            fprintf(spfidb1, 'sumrange\t');
        end
        fprintf(spfidb1, '\r');
        for (i=1:NStim+NSPcond)
            if i<=NStim
                fprintf(spfidb1, 'stim%d\t', i);
            elseif i<=NStim+NSPcond
                fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%             elseif i<=NStim+NSPcond
%                 fprintf(spfidb1, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
            end
            for (f=1:framenum)
                if sqrt(stimsumnum(i))>0
                    fprintf(spfidb1, '%3.8f\t', spstddrr(f, j, i)/sqrt(stimsumnum(i)));
                end
            end
            fprintf(spfidb1, '\r');
        end
        fprintf(spfidb1, '\r');
    end

    % superpix2 is just another format for sp1 (same stim are grouped together)
    fprintf(spfidb2, 'dR/R value\t');
    fprintf(spfidb2, '%f\t', xlabel);
    fprintf(spfidb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfidb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfidb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb2, 'fframe\t');
            fprintf(spfidb2, 'sumrange\t');
        end
        fprintf(spfidb2, '\r');
        for (j=1:masknum)
            fprintf(spfidb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfidb2, '%3.8f\t', spavgdrr(f, j, i));			% original value
            end
            fprintf(spfidb2, '\r');
        end
        fprintf(spfidb2, '\r');
    end
    fprintf(spfidb2, '\r\r_____Median_____\r\t');
    fprintf(spfidb2, '%f\t', xlabel);
    fprintf(spfidb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfidb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfidb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb2, 'fframe\t');
            fprintf(spfidb2, 'sumrange\t');
        end
        fprintf(spfidb2, '\r');
        for (j=1:masknum)
            fprintf(spfidb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfidb2, '%3.8f\t', spmeddrr(f, j, i));
                end
            end
            fprintf(spfidb2, '\r');
        end
        fprintf(spfidb2, '\r');
    end
    fprintf(spfidb2, '\r\r_____SD_____\r\t');
    fprintf(spfidb2, '%f\t', xlabel);
    fprintf(spfidb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfidb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfidb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb2, 'fframe\t');
            fprintf(spfidb2, 'sumrange\t');
        end
        fprintf(spfidb2, '\r');
        for (j=1:masknum)
            fprintf(spfidb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfidb2, '%3.8f\t', spstddrr(f, j, i));
                end
            end
            fprintf(spfidb2, '\r');
        end
        fprintf(spfidb2, '\r');
    end
    fprintf(spfidb2, '\r\r_____Num_____\r\t');
    fprintf(spfidb2, '%f\t', xlabel);
    fprintf(spfidb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfidb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfidb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb2, 'fframe\t');
            fprintf(spfidb2, 'sumrange\t');
        end
        fprintf(spfidb2, '\r');
        for (j=1:masknum)
            fprintf(spfidb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    fprintf(spfidb2, '%d\t', stimsumnum(i));
                end
            end
            fprintf(spfidb2, '\r');
        end
        fprintf(spfidb2, '\r');
    end
    fprintf(spfidb2, '\r\r_____SEM_____\r\t');
    fprintf(spfidb2, '%f\t', xlabel);
    fprintf(spfidb2, '\r');
    for (i=1:NStim+NSPcond)
        if i<=NStim
            fprintf(spfidb2, 'stim%d\t', i);
        elseif i<=NStim+NSPcond
            fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim, 1), 'junk'), 'junk'));
%         elseif i<=NStim+NSPcond
%             fprintf(spfidb2, '%s\t', getfield(cell2struct(SPcond(i-NStim-NSPcond, 1), 'junk'), 'junk'));
        end
        fprintf(spfidb2, 'f%d\t', 1:FramesPerStim);
        if flagsaveframeave
            fprintf(spfidb2, 'fframe\t');
            fprintf(spfidb2, 'sumrange\t');
        end
        fprintf(spfidb2, '\r');
        for (j=1:masknum)
            fprintf(spfidb2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                if blockselectnum>1
                    if sqrt(stimsumnum(i))>0
                        fprintf(spfidb2, '%3.8f\t', spstddrr(f, j, i)/sqrt(stimsumnum(i)));
                    end
                end
            end
            fprintf(spfidb2, '\r');
        end
        fprintf(spfidb2, '\r');
    end
    fclose(spfidb1);
    fclose(spfidb2);
end

clear spavg;
clear spstd;
clear spavgdrr;
clear spstddrr; 
if flagspmedian
	clear spmed;
	clear spmeddrr;
end





