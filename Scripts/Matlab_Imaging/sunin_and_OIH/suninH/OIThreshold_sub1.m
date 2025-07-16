%  [imresult, LowClip, HighClip]=OIThreshold_sub1(imin, method, clipvalue, immask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imresult, lowClip, highClip]=OIThreshold_sub1(imin, method, clipvalue, immask);
if immask==''
    immask=ones(504,504)*255;
end

switch method  %see above function description
case 0    % 0: no thresholding
    imresult=imin;
    lowClip=999;
    highClip=999;
    frameMedian=999;
case 1    % 1: Threshold at a certain percent from the max and min, centered at 0, within mask (thresholdmask)
    maskpixelnum=round(sum(sum(immask))/255); 
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    % frameMedian=median(usefulpixels);
    minvalue=min(usefulpixels);
    maxvalue=max(usefulpixels);
    lowClip=minvalue*(1-clipvalue);
    highClip=maxvalue*(1-clipvalue);
case 2    % 2: Threshold at a certain percent from the max and min, centered at median,  within mask (thresholdmask)
    maskpixelnum=round(sum(sum(immask))/255); 
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    frameMedian=median(usefulpixels);
    minvalue=min(usefulpixels);
    maxvalue=max(usefulpixels);
    lowClip=frameMedian-(frameMedian-minvalue)*(1-clipvalue);
    highClip=frameMedian+(maxvalue-frameMedian)*(1-clipvalue);
case 3    % 3: Threshold at a certain percentile from the max and min within mask (thresholdmask)
    maskpixelnum=round(sum(sum(immask))/255); 
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    usefulpixels_sort=sortrows(usefulpixels,1);
    lowClip=usefulpixels_sort(round(maskpixelnum*clipvalue))
    highClip=usefulpixels_sort(round(maskpixelnum*(1-clipvalue)))
case 4    % 4: Threshold at a certain percent from the max and min within mask (thresholdmask)
    maskpixelnum=round(sum(sum(immask))/255); 
    combined=[reshape(imin, prod(size(imin)), 1), reshape(double(immask), prod(size(imin)), 1)];    %following 3 lines for pick up useful pixels
    combined=flipud(sortrows(combined,2));
    usefulpixels=combined(1:maskpixelnum, 1);
    % usefulpixels_sort=sortrows(usefulpixels,1);
	minvalue=min(usefulpixels);
    maxvalue=max(usefulpixels);
    range=max-min;
    lowClip=minvalue+(range*clipvalue);
    highClip=maxvalue-(range*clipvalue);
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
