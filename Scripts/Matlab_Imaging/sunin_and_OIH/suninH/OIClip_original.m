%  [imresult, frameMedian, LowClip, HighClip]=OIClip(imin, method, value, immask)
%  imin: input image (double)
%  method:  0: no clipping
%           1: clipping to +-SD (value) 
%           2: clipping to +-SD (value) with mask (immask)
%           3: clipping to +-SD using the window specified in value (in this case it's a 2x3 matrix: x1, y1; x2, y2; and value(3, 1)=sd)
%           4: clipping by percentile (i.e. x% pixels above and below median value)     --need add
%  value: for clipping: sd value
%  immask: mask image, usually blood vessel pattern (0 & 1 ), 0 is bloodvessel
%
%  imin and immask must be same size
%  last modified: 040417
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imresult, frameMedian, lowClip, highClip]=OIClip(imin, method, value, immask);

switch method  %see above function description
case 0    %0: no clipping
    imresult=frame;
    lowClip=999;
    highClip=999;
    frameMedian=999;
    return;
case 1    % 1: clipping to +-SD (value) 
    temp=reshape(imin, prod(size(imin)), 1);    %transform image to 1-D array for calculation of median and stdev
    frameMedian=median(temp);    
    frameStd=std(temp);
    lowClip=frameMedian-value*frameStd;
    highClip=frameMedian+value*frameStd;
case 2    %2: clipping to +-SD (value) with mask (immask)
    maskpixelnum=sum(sum(immask));
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    frameMedian=median(usefulpixels);
    frameStd=std(usefulpixels);
    lowClip=frameMedian-value*frameStd;
    highClip=frameMedian+value*frameStd;
case 3  %3: clipping to +-SD using the window specified in value (in this case it's a 2x3 matrix: x1, y1; x2, y2; and value(3, 1)=sd)
    immasked=imin(value(1,1):value(2,1), value(1,2):value(2,2));    %x1:x2, y1:y2
    temp=reshape(immasked, prod(size(immasked)), 1);    
    frameMedian=median(temp);    
    frameStd=std(temp);
    lowClip=frameMedian-value(3, 1)*frameStd;
    highClip=frameMedian+value(3, 1)*frameStd;
case 4
    % need implement
    imresult=frame;   
    lowClip=999;
    highClip=999;
    frameMedian
    return;
otherwise
    fprintf(' Error: need specify cliping method in OIClip()!');
end
%clipping operations 
Templ= imin > lowClip; %logical operation 0 or 1; locations of lower value
Tempu= imin < highClip;%logical operation 0 or 1; locations of higher value
Tempul= (Templ.*Tempu).*imin;% bewteen low and high clips
frame=Tempul + (highClip*(~Tempu));%clip all avules higher than highClip
frame=frame + (lowClip*(~Templ));%clip all avules lower than lowClip
imresult=frame;       
return; 
