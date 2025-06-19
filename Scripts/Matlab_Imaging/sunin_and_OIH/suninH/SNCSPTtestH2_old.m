% This file is a section of Sunin code for OI data process
% Purpose: T-test for Super Pixel time course
% 
% 081125: Separated from Core program by Hisashi
% Originally written by Hisashi


if ~isdir([runfolder, SPresultname])
    mkdir(runfolder, SPresultname);    
end

if ~isdir([runfolder, SPresultname, 'superpixel_ttest'])
    mkdir([runfolder, SPresultname], 'superpixel_ttest');    
end

NSPttestcond=size(SPttestcond, 1);

%% raw
if spoperation==1 || spoperation==3
    for j=1:NSPttestcond
        ttestlist1=getfield(cell2struct(SPttestcond(j, 2), 'junk'), 'junk');
        ttestlist2=getfield(cell2struct(SPttestcond(j, 3), 'junk'), 'junk');
        SPttest1=zeros(framenum,masknum,blockselectnum*size(ttestlist1,2));
        SPttest2=zeros(framenum,masknum,blockselectnum*size(ttestlist2,2));
        SPttestsumnum1=0;
        SPttestsumnum2=0;
%         n=0;
%         m=0;
        if exist('flagSPTtestblanksubtraction','var') == 0
			flagSPTtestblanksubtraction=0;
        end
		if flagSPTtestblanksubtraction
			ttestlistBSub1=getfield(cell2struct(SPttestcondBSub(j, 2), 'junk'), 'junk');
			ttestlistBSub2=getfield(cell2struct(SPttestcondBSub(j, 3), 'junk'), 'junk');
		end
        for k=blockselect  % start block by block processing
			
			if flagSPTtestblanksubtraction %110301 HT
				SPttestcounterBSub1=0;
				SPttestcounterBSub2=0;
				SPttesttempBSub1=zeros(framenum,masknum);
				SPttesttempBSub2=zeros(framenum,masknum);
				
				for m=ttestlistBSub1
					if  goodstim(k, m)==1
						SPttestcounterBSub1=SPttestcounterBSub1+1;
						cumgs=0;
						for l=blockselect
							if l<=k && goodstim(l,m)==1
								cumgs=cumgs+1;
							end
						end
						SPttesttempBSub1=SPttesttempBSub1+squeeze(spvalues(:,:,m,cumgs));
					end
				end
				spttestavgBsub1=SPttesttempBSub1/SPttestcounterBSub1;
				for m=ttestlistBSub2
					if  goodstim(k, m)==1
						SPttestcounterBSub2=SPttestcounterBSub2+1;
						cumgs=0;
						for l=blockselect
							if l<=k && goodstim(l,m)==1
								cumgs=cumgs+1;
							end
						end
						SPttesttempBSub2=SPttesttempBSub2+squeeze(spvalues(:,:,m,cumgs));
					end
				end
				spttestavgBsub2=SPttesttempBSub2/SPttestcounterBSub2;
			end
			SPttestcounter1=0;
			SPttestcounter2=0;
			SPttesttemp1=zeros(framenum,masknum,size(ttestlist1,2));
			SPttesttemp2=zeros(framenum,masknum,size(ttestlist2,2));
            for m=ttestlist1
				if  goodstim(k, m)==1
					SPttestcounter1=SPttestcounter1+1;
					cumgs=0;
					for l=blockselect
						if l<=k && goodstim(l,m)==1
							cumgs=cumgs+1;
						end
					end
					SPttesttemp1(:,:,SPttestcounter1)=squeeze(spvalues(:,:,m,cumgs));
					if flagSPTtestblanksubtraction %110301 HT
						SPttesttemp1(:,:,SPttestcounter1)=SPttesttemp1(:,:,SPttestcounter1)-spttestavgBsub1;
					end                        
				end
			end
			for m=ttestlist2
				if  goodstim(k, m)==1
					SPttestcounter2=SPttestcounter2+1;
					cumgs=0;
					for l=blockselect
						if l<=k && goodstim(l,m)==1
							cumgs=cumgs+1;
						end
					end
					SPttesttemp2(:,:,SPttestcounter2)=squeeze(spvalues(:,:,m,cumgs));
					if flagSPTtestblanksubtraction %110301 HT
						SPttesttemp2(:,:,SPttestcounter2)=SPttesttemp2(:,:,SPttestcounter2)-spttestavgBsub2;
					end       
				end
            end
            if flagSPTtestblanksubtraction == 0 %110301 HT
                if SPttestcounter1>0
                    for n=1:SPttestcounter1
                        SPttest1(:,:,SPttestsumnum1+n)=SPttesttemp1(:,:,n);
                    end
                    SPttestsumnum1=SPttestsumnum1+SPttestcounter1;
                end
                if SPttestcounter2>0
                    for n=1:SPttestcounter2
                        SPttest2(:,:,SPttestsumnum2+n)=SPttesttemp2(:,:,n);
                    end 
                    SPttestsumnum2=SPttestsumnum2+SPttestcounter2;
                end
            else
                if SPttestcounter1>0 && SPttestcounterBSub1 >0
                    for n=1:SPttestcounter1
                        SPttest1(:,:,SPttestsumnum1+n)=SPttesttemp1(:,:,n);
                    end
                    SPttestsumnum1=SPttestsumnum1+SPttestcounter1;
                end
                if SPttestcounter2>0 && SPttestcounterBSub2 >0
                    for n=1:SPttestcounter2
                        SPttest2(:,:,SPttestsumnum2+n)=SPttesttemp2(:,:,n);
                    end 
                    SPttestsumnum2=SPttestsumnum2+SPttestcounter2;
                end
            end
        end
        if SPttestsumnum1 && SPttestsumnum2
            [H,p,CI,STATS] = ttest2(SPttest1(:,:,1:SPttestsumnum1), SPttest2(:,:,1:SPttestsumnum2), 0.05, tail, 'equal', 3);
            p2(:,:,j)=p;
        else
            p2(:,:,j)=1;
        end
    end


    spfidt1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_ttest\', 'superpixel_ttest1_raw1', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfidt2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_ttest\', 'superpixel_ttest1_raw2', '.txt')), 'w');    

    % superpix_original1
    fprintf(spfidt1, 'original value\t');
    fprintf(spfidt1, '%f\t', xlabel);
    fprintf(spfidt1, 'fframe\t');
    fprintf(spfidt1, 'sumrange\t');
    fprintf(spfidt1, '\r');
    for (j=1:masknum)
        fprintf(spfidt1, '%s\t', spmaskname(j, :));
        fprintf(spfidt1, 'f%d\t', 1:FramesPerStim);
        fprintf(spfidt1, '\r');
        for (i=1:NSPttestcond)
            fprintf(spfidt1, '%s\t', getfield(cell2struct(SPttestcond(i, 1), 'junk'), 'junk'));
            for (f=1:framenum)
                fprintf(spfidt1, '%1.8f\t', p2(f, j, i));			% original value
            end
            fprintf(spfidt1, '\r');
        end
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfidt2, 'original value\r\t');
    fprintf(spfidt2, '%f\t', xlabel);
    fprintf(spfidt2, 'fframe\t');
    fprintf(spfidt2, 'sumrange\t');
    fprintf(spfidt2, '\r');
    for (i=1:NSPttestcond)
        fprintf(spfidt2, '%s\t', getfield(cell2struct(SPttestcond(i, 1), 'junk'), 'junk'));
        fprintf(spfidt2, 'f%d\t', 1:FramesPerStim);
        fprintf(spfidt2, '\r');
        for (j=1:masknum)
            fprintf(spfidt2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfidt2, '%1.8f\t', p2(f, j, i));			% original value
            end
            fprintf(spfidt2, '\r');
        end
    end
end

%% dRR 	
if spoperation==2 || spoperation==3
    for j=1:NSPttestcond
        ttestlist1=getfield(cell2struct(SPttestcond(j, 2), 'junk'), 'junk');
        ttestlist2=getfield(cell2struct(SPttestcond(j, 3), 'junk'), 'junk');
        SPttest1=zeros(framenum,masknum,blockselectnum*size(ttestlist1,2));
        SPttest2=zeros(framenum,masknum,blockselectnum*size(ttestlist2,2));
        SPttestsumnum1=0;
        SPttestsumnum2=0;
%         n=0;
%         m=0;
        if exist('flagSPTtestblanksubtraction','var') == 0
			flagSPTtestblanksubtraction=0;
        end
		if flagSPTtestblanksubtraction
			ttestlistBSub1=getfield(cell2struct(SPttestcondBSub(j, 2), 'junk'), 'junk');
			ttestlistBSub2=getfield(cell2struct(SPttestcondBSub(j, 3), 'junk'), 'junk');
		end
        for k=blockselect  % start block by block processing
			
			if flagSPTtestblanksubtraction %110301 HT
				SPttestcounterBSub1=0;
				SPttestcounterBSub2=0;
				SPttesttempBSub1=zeros(framenum,masknum);
				SPttesttempBSub2=zeros(framenum,masknum);
				
				for m=ttestlistBSub1
					if  goodstim(k, m)==1
						SPttestcounterBSub1=SPttestcounterBSub1+1;
						cumgs=0;
						for l=blockselect
							if l<=k && goodstim(l,m)==1
								cumgs=cumgs+1;
							end
						end
						SPttesttempBSub1=SPttesttempBSub1+squeeze(spvaluesdrr(:,:,m,cumgs));
					end
				end
				spttestavgBsub1=SPttesttempBSub1/SPttestcounterBSub1;
				for m=ttestlistBSub2
					if  goodstim(k, m)==1
						SPttestcounterBSub2=SPttestcounterBSub2+1;
						cumgs=0;
						for l=blockselect
							if l<=k && goodstim(l,m)==1
								cumgs=cumgs+1;
							end
						end
						SPttesttempBSub2=SPttesttempBSub2+squeeze(spvaluesdrr(:,:,m,cumgs));
					end
				end
				spttestavgBsub2=SPttesttempBSub2/SPttestcounterBSub2;
			end
			SPttestcounter1=0;
			SPttestcounter2=0;
			SPttesttemp1=zeros(framenum,masknum,size(ttestlist1,2));
			SPttesttemp2=zeros(framenum,masknum,size(ttestlist2,2));
            for m=ttestlist1
				if  goodstim(k, m)==1
					SPttestcounter1=SPttestcounter1+1;
					cumgs=0;
					for l=blockselect
						if l<=k && goodstim(l,m)==1
							cumgs=cumgs+1;
						end
					end
					SPttesttemp1(:,:,SPttestcounter1)=squeeze(spvaluesdrr(:,:,m,cumgs));
					if flagSPTtestblanksubtraction %110301 HT
						SPttesttemp1(:,:,SPttestcounter1)=SPttesttemp1(:,:,SPttestcounter1)-spttestavgBsub1;
					end                        
				end
			end
			for m=ttestlist2
				if  goodstim(k, m)==1
					SPttestcounter2=SPttestcounter2+1;
					cumgs=0;
					for l=blockselect
						if l<=k && goodstim(l,m)==1
							cumgs=cumgs+1;
						end
					end
					SPttesttemp2(:,:,SPttestcounter2)=squeeze(spvaluesdrr(:,:,m,cumgs));
					if flagSPTtestblanksubtraction %110301 HT
						SPttesttemp2(:,:,SPttestcounter2)=SPttesttemp2(:,:,SPttestcounter2)-spttestavgBsub2;
					end       
				end
            end
            if flagSPTtestblanksubtraction == 0 %110301 HT
                if SPttestcounter1>0
                    for n=1:SPttestcounter1
                        SPttest1(:,:,SPttestsumnum1+n)=SPttesttemp1(:,:,n);
                    end
                    SPttestsumnum1=SPttestsumnum1+SPttestcounter1;
                end
                if SPttestcounter2>0
                    for n=1:SPttestcounter2
                        SPttest2(:,:,SPttestsumnum2+n)=SPttesttemp2(:,:,n);
                    end 
                    SPttestsumnum2=SPttestsumnum2+SPttestcounter2;
                end
            else
                if SPttestcounter1>0 && SPttestcounterBSub1 >0
                    for n=1:SPttestcounter1
                        SPttest1(:,:,SPttestsumnum1+n)=SPttesttemp1(:,:,n);
                    end
                    SPttestsumnum1=SPttestsumnum1+SPttestcounter1;
                end
                if SPttestcounter2>0 && SPttestcounterBSub2 >0
                    for n=1:SPttestcounter2
                        SPttest2(:,:,SPttestsumnum2+n)=SPttesttemp2(:,:,n);
                    end 
                    SPttestsumnum2=SPttestsumnum2+SPttestcounter2;
                end
            end
        end
        if SPttestsumnum1 && SPttestsumnum2
            [H,p,CI,STATS] = ttest2(SPttest1(:,:,1:SPttestsumnum1), SPttest2(:,:,1:SPttestsumnum2), 0.05, tail, 'equal', 3);
            p2(:,:,j)=p;
        else
            p2(:,:,j)=1;
        end
    end

    spfidu1= fopen(strcat(runfolder, SPresultname, strcat('superpixel_ttest\', 'superpixel_ttest1_dRR1', '.txt')), 'w');  % for sp avg & stdev output, 1, 2 are two different formats
    spfidu2= fopen(strcat(runfolder, SPresultname, strcat('superpixel_ttest\', 'superpixel_ttest1_dRR2', '.txt')), 'w');

    % superpix2 is just another format for sp1 (same stim are grouped together)
    fprintf(spfidu1, 'dR/R value\t');
    fprintf(spfidu1, '%f\t', xlabel);
    fprintf(spfidu1, 'fframe\t');
    fprintf(spfidu1, 'sumrange\t');
    fprintf(spfidu1, '\r');
    for (j=1:masknum)
        fprintf(spfidu1, '%s\t', spmaskname(j, :));
        fprintf(spfidu1, 'f%d\t', 1:FramesPerStim);
        fprintf(spfidu1, '\r');
        for (i=1:NSPttestcond)
            fprintf(spfidu1, '%s\t', getfield(cell2struct(SPttestcond(i, 1), 'junk'), 'junk'));
            for (f=1:framenum)
                fprintf(spfidu1, '%1.8f\t', p2(f, j, i));			% original value
            end
            fprintf(spfidu1, '\r');
        end
    end

    % superpix_original2 is just another format for sp1 (same stim are grouped together) %added by Hisashi 
    fprintf(spfidu2, 'dR/R value\t');
    fprintf(spfidu2, '%f\t', xlabel);
    fprintf(spfidu2, 'fframe\t');
    fprintf(spfidu2, 'sumrange\t');
    fprintf(spfidu2, '\r');
    for (i=1:NSPttestcond)
        fprintf(spfidu2, '%s\t', getfield(cell2struct(SPttestcond(i, 1), 'junk'), 'junk'));
        fprintf(spfidu2, 'f%d\t', 1:FramesPerStim);
        fprintf(spfidu2, '\r');
        for (j=1:masknum)
            fprintf(spfidu2, '%s\t', spmaskname(j, :));
            for (f=1:framenum)
                fprintf(spfidu2, '%1.8f\t', p2(f, j, i));			% original value
            end
            fprintf(spfidu2, '\r');
        end
    end
end


