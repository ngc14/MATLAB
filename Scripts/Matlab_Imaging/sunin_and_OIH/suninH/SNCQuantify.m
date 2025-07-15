% This file is a section of Sunin code for OI data process
% Purpose: Save average values across domains
% 
% 081125: Separated from Core program by Hisashi
% Originally written by HDL

% output detailed (including every block) averaged masked values for each domain type (masks) and stimulus. 
if (flagquantify==2)
    fundfid1=fopen(strcat(resultfolder, 'fun1detail.txt'), 'w');
    for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map
        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
        fprintf(fundfid1, '%s\t', mapname);
        for k=1:blockselectnum
            fprintf(fundfid1,'blk%d\t', k);
        end
        fprintf(fundfid1, '\r');
        for j=1:masknum
            fprintf(fundfid1, '%s\t', spmaskname(j, :));
            for k=1:blockselectnum
                fprintf(fundfid1, '%f\t', domainvalues(j, i, k));
            end
            fprintf(fundfid1, '\r');
        end
        fprintf(fundfid1, '\r');
    end
    fclose(fundfid1);
    fundfid2=fopen(strcat(resultfolder, 'fun1detail2.txt'), 'w');
    for j=1:masknum
        fprintf(fundfid1, '%s\t', spmaskname(j, :));
        for k=1:blockselectnum
            fprintf(fundfid2,'%d\t', k-1);
        end
        fprintf(fundfid2, '\r');
        for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map
            mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
            fprintf(fundfid2, '%s\t', mapname);
            for k=1:blockselectnum
                fprintf(fundfid2, '%f\t', domainvalues(j, i, k));
            end
            fprintf(fundfid2, '\r');
        end
        fprintf(fundfid2, '\r');
    end
    fclose(fundfid2);
end

% save averaged map value covered by the mask. calculate SD 
if (flagspmask>0 & SDtype~=1)  % calculate mean and SD across blocks
    domainavg=zeros(masknum, mapnum);
    domainstd=zeros(masknum, mapnum);
    if (flaggoodstim==1)    % note: for selected stim presentation analysis, no analysis map quantatification since there is no such map for every block
        for j=1:masknum
            newdomainvalues=reshape(domainvalues(j,1:NStim,:), NStim, blockselectnum);
            domainavg(j,1:NStim)=(sum(newdomainvalues.*goodstim', 2)./sum(goodstim', 2))';
            % sd?
        end
    else
        domainavg=mean(domainvalues, 3); % mean across blocks (3rd dimension)
        domainstd=std(domainvalues, 0, 3); 
    end
%    [0, p, 95 sig level] = 
    %[hypothesis, pvalue, percentile, siglevel]=ttest2(domainvalues(1,7,:),domainvalues(2,7,:), 0.05, 0);
    %pvalue
elseif (flagspmask==2&SDtype==1) % calculate mean & SD across dots (use 'maps' matrix instead 'domainvalues')
    for i=1:mapnum %NStim   % mapnum for calculate all map, 'NStim' only for single-condition map     
        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
        for j=1:masknum
            dotvalue=zeros(domaindotnum(j),1);
            for k=1:domaindotnum(j)
                L=1;
                pixelvalue=0;
                for x=round(xx(k, j)-domainradius(j)):round(xx(k, j)+domainradius(j))
                    for y=round(yy(k, j)-domainradius(j)):round(yy(k, j)+domainradius(j))
                        if ((x-xx(k,j))^2+(y-yy(k,j))^2)<=domainradius(j)^2
                            pixelvalue=pixelvalue+maps(y,x,i);
%                            if (i==1)  % for testing
%                                testmap(y,x,j)=255;
%                            end
                            L=L+1;
                        end
                    end
                end
                dotvalue(k)=pixelvalue/(L-1);
            end
            domainavg(j, i)=mean(dotvalue);     % mask in column and condition in row
            domainstd(j, i)=std(dotvalue);
            if (flagsavedetailfun==1)
                dotindextemp=dotdetail(:,1)*100+dotdetail(:,2);     % just for easy finding
                if ~isempty(find(dotindextemp == (i*100+j)))    % in dotdetail, first colum is mapnum, second column is masknum
                    dotfid=fopen(strcat(resultfolder, mapname, '-', spmaskname(j,:), '.txt'),'w');
                    for ii=1:domaindotnum(j)
                        fprintf(dotfid, '%f\r',dotvalue(ii));
                    end
                    fclose(dotfid);
                end
            end
            clear dotvalue;    
        end
    end
else
%    fprintf ('\r!!!Error001, check flag "flagspmask" & "SDtype"\r');
    %return;
end
% output


funfid1= fopen(strcat(resultfolder, strcat('fun1', '.txt')), 'a');  % for domain avg & stdev output, 1, 2 & 3 are three different formats
funfid2= fopen(strcat(resultfolder, strcat('fun2', '.txt')), 'a');
funfid3= fopen(strcat(resultfolder, strcat('fun3', '.txt')), 'a');
stdfid1 =fopen(strcat(resultfolder, strcat('std1', '.txt')), 'a');
stdfid2 =fopen(strcat(resultfolder, strcat('std2', '.txt')), 'a');
stdfid3 =fopen(strcat(resultfolder, strcat('std3', '.txt')), 'a');
if flagafunout==1
    outnum=mapnum;  %output all map values
else
    outnum=NStim;   %output only stim map values
end
if (flagspmask>0)
    fprintf(funfid1, '\t');
    fprintf(stdfid1, '\t');
    for i=1:outnum      %output to fun1, group by domains after all blocks
        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
        fprintf(funfid1, '\t%s', mapname);
        fprintf(stdfid1, '\t%s', mapname);
    end
    fprintf(funfid1, '\r');
    fprintf(stdfid1, '\r');
    for i=1:masknum       
        fprintf(funfid1, '\t%s\t', spmaskname(i,:));
        fprintf(stdfid1, '\t%s\t', spmaskname(i,:));
        for j=1:outnum
            fprintf(funfid1, '%10.7f\t', domainavg(i, j));
            fprintf(stdfid1, '%10.7f\t', domainstd(i, j));
        end
%        fprintf(funfid1, '\r');    %put these two line on for single block output
%        fprintf(stdfid1, '\r');
    end
    fprintf(funfid1, '\r');
    fprintf(stdfid1, '\r');

    fprintf(funfid2, '\t');
    fprintf(stdfid2, '\t');
    for i=1:mapnum           %output to fun2, group by SF
        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
        fprintf(funfid2, '\t%s', mapname);
        fprintf(stdfid2, '\t%s', mapname);
    end
    fprintf(funfid2, '\r');
    fprintf(stdfid2, '\r');
    for i=1:masknum     
        fprintf(funfid2, '\t%s\t', spmaskname(i,:));
        fprintf(stdfid2, '\t%s\t', spmaskname(i,:));
        for j=1:outnum
            fprintf(funfid2, '%10.7f\t', domainavg(i, j));
            fprintf(stdfid2, '%10.7f\t', domainstd(i, j));
        end
        fprintf(funfid2, '\r');
        fprintf(stdfid2, '\r');
    end
    fprintf(funfid2, '\r');
    fprintf(stdfid2, '\r');

%    for i=1:masknum           %output to fun3, group by contrast
%        fprintf(funfid3, '\t%s', spmaskname(i,:));
%        fprintf(stdfid3, '\t%s', spmaskname(i,:));
%    end
    fprintf(funfid3, '\r');
    fprintf(stdfid3, '\r');
    for i=1:outnum     
        mapname=getfield(cell2struct(mapnames(i), 'junk'), 'junk');
        fprintf(funfid3, '\t%s\t', mapname);
        fprintf(stdfid3, '\t%s\t', mapname);
        for j=1:masknum
            fprintf(funfid3, '%10.7f\t', domainavg(j, i));
            fprintf(stdfid3, '%10.7f\t', domainstd(j, i));
        end
    end
    fprintf(funfid3, '\r');
    fprintf(stdfid3, '\r');
end
fclose(funfid1);
fclose(stdfid1);
fclose(funfid2);
fclose(stdfid2);
fclose(funfid3);
fclose(stdfid3);


