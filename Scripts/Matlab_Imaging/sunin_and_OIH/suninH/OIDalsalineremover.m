function imresult=OIDalsalineremover(imin)

% remove artefact lines in images captured by Dalsa 1M60P. Those lines are made junctions between CCD panels
%  (121118) Ver 1.0    made by Hisashi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imresult=imin;

% imresult(:,128)=0;
% imresult(:,256)=0;
% imresult(:,384)=0;

% imresult(:,128)=(imin(:,127)*2+imin(:,130))/3;
% imresult(:,129)=(imin(:,127)+imin(:,130)*2)/3;
% 
% imresult(:,256)=(imin(:,255)+imin(:,259))/2;
% imresult(:,258)=(imin(:,257)+imin(:,259))/2;
% imresult(:,257)=(imresult(:,256)+imresult(:,258))/2;

% imresult(:,128)=(imin(:,127)+imin(:,129))/2;
% imresult(:,130)=(imin(:,129)+imin(:,131))/2;
% imresult(:,129)=(imresult(:,128)+imresult(:,130))/2;

imresult(:,128)=(imin(:,127)*2+imin(:,130))/3;
imresult(:,129)=(imin(:,127)+imin(:,130)*2)/3;
% imresult(:,128)=(imresult(:,127)+imresult(:,128)+imresult(:,129))/3; %temp


imresult(:,256)=(imin(:,255)+imin(:,257))/2;
imresult(:,258)=(imin(:,257)+imin(:,259))/2;
imresult(:,257)=(imresult(:,256)+imresult(:,258))/2;

imresult(:,384)=(imin(:,383)+imin(:,385))/2;
imresult(:,386)=(imin(:,385)+imin(:,387))/2;
imresult(:,385)=(imresult(:,384)+imresult(:,386))/2;



return