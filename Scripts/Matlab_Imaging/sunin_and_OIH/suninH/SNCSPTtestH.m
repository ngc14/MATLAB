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
        condlist1=getfield(cell2struct(SPttestcond(j, 2), 'junk'), 'junk');
        condlist2=getfield(cell2struct(SPttestcond(j, 3), 'junk'), 'junk');
        SPttesttemp1=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
        SPttesttemp2=zeros(framenum,masknum,blockselectnum*size(condlist2,2));
        n=1;
        m=1;
        if exist('flagblanksubtraction','var') == 0
			flagblanksubtraction=0;
        end
		if flagblanksubtraction
			condlistBSub1=getfield(cell2struct(SPttestcondBSub(j-NStim, 2), 'junk'), 'junk');
			condlistBSub2=getfield(cell2struct(SPttestcondBSub(j-NStim, 3), 'junk'), 'junk');
			SPttesttempBSub1=zeros(framenum,masknum);
			SPttesttempBSub2=zeros(framenum,masknum);			
			for i=condlist1
				SPttesttempBSub1(:,:,n:n+stimsumnum(i)-1)=squeeze(spvalues(:,:,i,1:stimsumnum(i)));
				n=n+stimsumnum(i);
			end
			for i=condlist2
				SPttesttempBSub2(:,:,m:m+stimsumnum(i)-1)=squeeze(spvalues(:,:,i,1:stimsumnum(i)));
				m=m+stimsumnum(i);
			end
			% not completed (121114)
        end
	
        for i=condlist1
            SPttesttemp1(:,:,n:n+stimsumnum(i)-1)=squeeze(spvalues(:,:,i,1:stimsumnum(i)));
            n=n+stimsumnum(i);
        end
        for i=condlist2
            SPttesttemp2(:,:,m:m+stimsumnum(i)-1)=squeeze(spvalues(:,:,i,1:stimsumnum(i)));
            m=m+stimsumnum(i);
        end
        [H,p,CI,STATS] = ttest2(SPttesttemp1(:,:,1:n-1), SPttesttemp2(:,:,1:m-1), 0.05, tail, 'equal', 3);
        p2(:,:,j)=p;
        clear SPttesttemp1;
        clear SPttesttemp2;
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
        condlist1=getfield(cell2struct(SPttestcond(j, 2), 'junk'), 'junk');
        condlist2=getfield(cell2struct(SPttestcond(j, 3), 'junk'), 'junk');
        SPttesttemp1=zeros(framenum,masknum,blockselectnum*size(condlist1,2));
        SPttesttemp2=zeros(framenum,masknum,blockselectnum*size(condlist2,2));
        n=1;
        m=1;
        for i=condlist1
            SPttesttemp1(:,:,n:n+stimsumnum(i)-1)=squeeze(spvaluesdrr(:,:,i,1:stimsumnum(i)));
            n=n+stimsumnum(i);
        end
        for i=condlist2
            SPttesttemp2(:,:,m:m+stimsumnum(i)-1)=squeeze(spvaluesdrr(:,:,i,1:stimsumnum(i)));
            m=m+stimsumnum(i);
        end
        [H,p,CI,STATS] = ttest2(SPttesttemp1(:,:,1:n-1), SPttesttemp2(:,:,1:m-1), 0.05, tail, 'equal', 3);
        p2(:,:,j)=p;   
        clear SPttesttemp1;
        clear SPttesttemp2;
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


