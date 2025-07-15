function threshold=OIAshift2goodstim2(shiftfile2, threshold, NStim, percentage)

% fframe alignment usually contains come frames that has very low correlation coeff, this process is to discard these frames.
% input: 
%   shiftfile2: output file from OIAlign (text file, 11 columns, Nblock*NStim rows)
%   threshold:  threshold correlation-coefficient-difference, if one frame has larger than 'threshold' difference than neighbouring frames, this frame will be discard, not used if percentage is used
%   NStim: number of stimulus, for output format
%   percentage: percentage of frames that pass this discrimination, note: since using relative measurement (compare neighbours), it's meaningless to have percentage < 50%, usually 90%
% output a file with 'goodstim' info, default filename "goodstim2 + input filename"
% Method: find those frames that has big cange in correlation, used smooth line method. 

manual=0;   % for testing
if manual
    shiftfile2='_shiftlog2006-4-20-14-30..2.txt';
    threshold=0.0001;   % corelation threshold usually 0.0001
    NStim=5;
    percentage=0.9;     % usually 0.9
end
kernel=5;           % how many trials on each side to average to get a smooth curve, should be an odd number, usually 5
outputfile=strcat(shiftfile2(1:end-4), '_ff_gs.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if percentage>0.99 | percentage<0.5
    fprintf('Note: percentage is too small or too big in "OIAShift2goodstim2"!\r');
end

[filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj]=OIAReadShiftLog(shiftfile2);
N=size(filename, 1);    % Nblock*Nstim
A=zeros(N, 3);
A(:,1)=dx;
A(:,2)=dy;
A(:,3)=bestobj;

B=zeros(N,1);   %smooth line
C=ones(N,1);      % goodstim
for i=1:N
    mincorr=1;
    n=0;
    for j=max(i-kernel, 1):min(i+kernel, N) 
        B(i)=B(i)+A(j,3);
        n=n+1;
        if A(i,3)<mincorr
            mincorr=A(j,3);
        end
    end
    B(i)=(B(i)-mincorr)/(n-1);  % minus the lowest correlation to reduce noise
end

C=(B-A(:,3))>threshold;
percenttemp=sum(C(:))/N;

sd=std(A(:,3));  % these 5 lines are for fast covergence. 
oldsign=0;  
regulation=1;
tolerance=0.001;
count=0;

% find a threshold level that gives 'percentage' good stim.
if percentage~=0
    while abs(percenttemp-percentage)>tolerance
        if percenttemp>percentage   
            sign=-1;    % need decrease threshold (make difference smaller)
        else 
            sign=1;
        end
        if oldsign==0
            oldsign=sign;
        end
        if oldsign~=sign    % there is a reverse, so make the change of threshold smaller (fine tuning)
            regulation=regulation+10;   
            oldsign=sign;
            count=count+1;
            if count>=100
                tolerance=tolerance+tolerance;
                count=0;
            end
        end
        threshold=threshold+sign*sd/regulation;
%        fprintf('%d\t   %f\t   %f\r', regulation, threshold, percenttemp);
        C=B-A(:,3)<threshold;
        percenttemp=sum(C(:))/N;
    end
end

C=reshape(C, [NStim, N/NStim]);
C=C';

% output goodstim
fid=fopen(outputfile, 'w');
for i=1:N/NStim
    for j=1:NStim
        fprintf(fid, '%d \t', C(i, j));
    end
    fprintf(fid, '\r\n');
end
fclose(fid);
fprintf('\rGoodstim2: threshold=%f,  goodstim=%f%%\r\n', threshold, 100*sum(C(:))/N);

return
