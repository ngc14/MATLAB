function OIAshift2goodstim3H(shiftfile3, corrthreshold, shiftthreshold, NStim, NFrame, framerange)

% For processing within stim shifts (shift3), discard large shift stims. % 
% input: 
%   shiftfile3: output file from OIAlign (text file, 11 columns, Nblock*NStim rows)
%   shiftthreshold:  % pixels: if mean pixel shift for framerange exceeds this value, the trial is marked as '0'
%   percentage: % percentage of trials to be inluded, this percentage is for absolute correlation value only. 
%   NStim: number of stimulus, for output format
%   NFrame: number of frames per stim
%   framerange: frames range for sum frames.
% output a file with 'goodstim' info, default filename "withingoodstim + input filename"

% Method: find out those frames has large shift or poor correlation. 

manual=0;   % for testing
if manual
    shiftfile3='_shiftlog2006-11-14-19-31..3.txt';
    shiftthreshold=1;   
    percentage=0.9;     % percentage of '1's
    NStim=9;
    NFrame=16;
    framerange=[5:16];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputfile=strcat(shiftfile3(1:end-4), '_within_gs.txt');

[filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj]=OIAReadShiftLog(shiftfile3);
NBlock=size(filename, 1)/NStim/NFrame;
A=zeros(size(filename,1), 4);
A(:,1)=dx;
A(:,2)=dy;
A(:,3)=bestobj;
% A(:,4)=currentobj;
A(:,4)=bestobj; % modified by Hisashi on 2/22/2010
A=reshape(A', [4, NFrame, NStim, NBlock]);

%%% commented out by Hisashi
% index=ceil(NBlock*NStim*(1-percentage));
% tempobj=min(A(4, framerange, :, :), [], 2);
% tempobj=reshape(tempobj, [NStim*NBlock, 1]);
% currentobj=sort(tempobj);
% corrthreshold=currentobj(index);


C=ones(4, NStim, NBlock);      % goodstim
for i=1:NBlock
    for j=1:NStim
%         mean1=mean(abs(A(1, framerange, j, i)), 2);
%         mean2=mean(abs(A(2, framerange, j, i)), 2);
        mean1=max(max(abs(A(1, framerange, j, i)), 2));
        mean2=max(max(abs(A(2, framerange, j, i)), 2));

% 		mean3=mean((A(3, framerange, j, i)), 2);
% 		mean4=mean((A(4, framerange, j, i)), 2);
        if mean1>shiftthreshold 
			C(1, j, i)=0;
        end
        if mean2>shiftthreshold 
			C(2, j, i)=0;
        end
% 		if min(A(4, :, j, i))<corrthreshold
        if min(A(4, framerange, j, i))<corrthreshold
			C(4, j, i)=0;
        end
    end
end
fid=fopen(outputfile, 'w');
for i=1:NBlock
    for j=1:NStim
        fprintf(fid, '%d \t', C(1,j,i)*C(2,j,i)*C(4,j,i));
    end
    fprintf(fid, '\r\n');
end
fclose(fid);
fprintf('\rGoodstim3H:\rBased on shifts(threshold=%1.1f), goodstim=%2.1f%%\r', shiftthreshold, sum(sum(C(1,:,:).*C(2,:,:)))*100/NBlock/NStim);
fprintf('Based on the original correlation(threshold=%0.5f), goodstim=%2.1f%%\r', corrthreshold, sum(sum(C(4,:,:)))*100/NBlock/NStim);
fprintf('Based on both, goodstim=%2.1f%%\r', sum(sum(C(1,:,:).*C(2,:,:).*C(4,:,:)))*100/NBlock/NStim);
return
