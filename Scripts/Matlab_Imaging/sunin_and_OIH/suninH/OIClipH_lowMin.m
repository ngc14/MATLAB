%  [imresult, frameMedian, LowClip, HighClip]=OIClipH_lowMin(imin, method, value, immask, value2)
%  imin: input image (double)
%  method:
%         0: no clipping; 
%         1: clipping at median+-SD (value);  
%         2: clipping at median+-SD (value) with a mask (masks\clipmask\default.bmp);  
%         3: clipping using the window specified in 'clipvalue' (in this case it's a 2x2 matrix x1, y1; x2, y2) 
%         4: clipping at median+-intensity change (value), Ex. 0.0005 or 0.0008;  
%         5: clipping at 0+-SD (value);
%         6: clipping at 0+-intensity change (value);
%     
%  value: for clipping: sd value
%  immask: mask image, usually blood vessel pattern (0 & 1 ), 0 is bloodvessel
%
%  imin and immask must be same size
%
%  090116 lowClip is always set the to minimum value in the image. This is for vector analysis. Hisashi
%  071023 Case 5&6 was added by Hisashi
%  070413 Case 4 was added by Hisashi
%  040417 last modified: Haidong
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imresult, frameMedian, lowClip, highClip]=OIClipH_lowMin(imin, method, value, immask);

switch method  %see above function description
case 0    %0: no clipping
    imresult=imin;
    lowClip=999;
    highClip=999;
    frameMedian=999;
    return;
case 1    % 1: clipping to median +-SD (value) 
    temp=reshape(imin, prod(size(imin)), 1);    %transform image to 1-D array for calculation of median and stdev
    frameMedian=median(temp);    
    frameStd=std(temp);
%     lowClip=frameMedian-value*frameStd;
    lowClip=min(temp);
    highClip=frameMedian+value*frameStd;
case 2    %2: clipping to median +-SD (value) with mask (immask)
    maskpixelnum=round(sum(sum(immask))/255); %modified by Hisashi
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    frameMedian=median(usefulpixels);
    frameStd=std(usefulpixels);
%     lowClip=frameMedian-value*frameStd;
    lowClip=min(usefulpixels);
    highClip=frameMedian+value*frameStd;
case 3  %3: clipping to +-SD using the window specified in value (in this case it's a 2x3 matrix: x1, y1; x2, y2; and value(3, 1)=sd)
    immasked=imin(value(1,1):value(2,1), value(1,2):value(2,2));    %x1:x2, y1:y2
    temp=reshape(immasked, prod(size(immasked)), 1);    
    frameMedian=median(temp);    
    frameStd=std(temp);
%     lowClip=frameMedian-value(3, 1)*frameStd;
    lowClip=min(temp);
    highClip=frameMedian+value(3, 1)*frameStd;
case 4  % 4: clipping to median +-certain value (value) 
    temp=reshape(imin, prod(size(imin)), 1);    %transform image to 1-D array for calculation of median and stdev
    frameMedian=median(temp);    
    frameStd=std(temp);
%     lowClip=frameMedian-value;
    lowClip=min(temp);
    highClip=frameMedian+value;
case 5    % 5: clipping to 0 +-SD (value) 
    temp=reshape(imin, prod(size(imin)), 1);    %transform image to 1-D array for calculation of median and stdev
    frameMedian=0; %median(temp);    
    frameStd=std(temp);
    lowClip=min(temp);
    highClip=value*frameStd;
case 6  % 6: clipping to 0 +-certain value (value) 
    temp=reshape(imin, prod(size(imin)), 1);    %transform image to 1-D array for calculation of median and stdev
    frameMedian=0;    
    frameStd=std(temp);
    lowClip=min(temp);
    highClip=value;
case 7
    % need implement
    imresult=frame;   
    lowClip=999;
    highClip=999;
    frameMedian;
    return;
otherwise
    fprintf(' Error: need specify cliping method in OIClip()!');
end

% added by Hisashi for test
% fprintf('frameStd = %g\n', frameStd)

%clipping operations
Templ= imin > lowClip; %logical operation 0 or 1; locations of lower value
Tempu= imin < highClip;%logical operation 0 or 1; locations of higher value
Tempul= (Templ.*Tempu).*imin;% bewteen low and high clips
frame=Tempul + (highClip*(~Tempu));%clip all avules higher than highClip
frame=frame + (lowClip*(~Templ));%clip all avules lower than lowClip
imresult=frame;       
return; 
