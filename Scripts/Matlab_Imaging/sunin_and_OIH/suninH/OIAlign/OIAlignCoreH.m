function [imresult, pixshift, bestobj, currentobj] = OIAlignCore(imbase, imshift, range, objfun)

% function [imresult, pixshift, bestobj] = OIAlignCore(imbase, imshift, objfun)
%
% Ver H1: (150221:Hisashi) change variable name 'pixoff' to 'pixshift'
% Ver 1: (060417:HDL) modified from OIAlign2.m
%
% Align 'imshift' with 'imbase'
% Doesn't include preprocess.
%
% input
% 	imbase: pedestal image to be compared, usually is the very first frame
% 	imshift: image to be shifted & aligned with 'imbase', size should be the same as 'imbase'
% 	range(4): (x1, x2, y1, y2) for (left, right, top, bottom) max shift range
%   objfun: Objective function
%       1: fast correlation, use 'normxcorr2.m', precision=1 pixel
%       2: slow correlation, a combination of fast correlation (1) and sub-pixel correlation, default precision=0.1 pixel (change 'precision) to modify this default.  
%       3: use stdev of difference image (default precision=1 pixel)
%       4: min-difference score (default precision=1 pixel)
%       5: use mutual info to align (note: only uint8 image can be use here, precision=1 pixel)
% return: 
%   imresult: shifted image aligned with 'imbase' (now block this function to improve speed, 060417)
%   pixshift(2): (dx, dy) values: shift 'imshift' by (dx, dy) will achieve best alignment with 'imbase'
%   bestobj: value of objective function at best alignment. 
%   currentobj: value of objective function at current (non-shift) position
%
% Note: 'imshift' should be the size as 'imbase'
%       Datatype: 'imshift' and 'imbase' should be 2D grayscale image, any data type (but when objfun=5, they will be transformed to 'uint8).
% improvement: use range() to replace croppixel, will improve some speed.
%

%% Initial check
basesize = size (imbase);
imshiftsize = size (imshift);
imresult=imshift;       % note, 'imresult' is not imresult
if size(basesize)>2 | size(imshiftsize)>2
    fprintf('Error in OIAlignCore: only 2D map can be aligned (no color!)\r\r');
    return;
end
if size(basesize)~=size(imshiftsize)
    fprintf('Error in OIAlignCore: sizes do not match (imbase, imshift)!\r\r');
    return;
end
croppixel= ceil(max(max(abs(range)))+1); % how many pixel will corp off along the edges, for obtain a smaller area for comparison, since shifting creates 'blank' margins
% basesize
% croppixel
if min(basesize)<2*croppixel
%     fprintf('Error in OIAlignCore: min(basesize)<2*croppixel !\r\r');
    fprintf('Error in OIAlignCore: min(clip region)<2*max(shiftrange) !\r\r');    
    return;
end

pixshift=[0 0];
if objfun==1 | objfun==2    
	imshiftcrop=imshift(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
	cc=normxcorr2(imshiftcrop, imbase);	
	[bestobj, maxlocation] = max(cc(:));	%abs(cc)?
	[yoffset, xoffset] = ind2sub(size(cc), maxlocation(1));
	pixshift(1)=-(size(imbase,2)-croppixel-xoffset);	% note: reverse sign to be comparable with score 1-3 (shift back pixel values)
	pixshift(2)=-(size(imbase,1)-croppixel-yoffset);
    currentobj=cc(size(imbase,1)-croppixel, size(imbase,2)-croppixel);
    if currentobj>bestobj
        fprintf('strange!\r');
    end
%	imresult=OIShift(imshift, pixshift(1), pixshift(2));
end
if objfun==1    % fast correlation, no further alignment
    return;
end             % slow correlation, do fine-tuning alignment based on fast correlation results (normcorr2)

if objfun==2    % Do fine-tuning shift 
    precision=0.1;  % step size in shifting an image
    croppixel=croppixel+1;  
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    pixshifttemp=pixshift;
    for xstep=pixshift(1)-0.6:precision:pixshift(1)+0.6    % x shifting         
        for ystep=pixshift(2)-0.6:precision:pixshift(2)+0.6   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            cortemp=corr2(imbasesm, imtempsm);
            if (cortemp>bestobj)
                bestobj=cortemp;
                pixshifttemp=[xstep ystep];
            end
            if (cortemp*-1>bestobj) % 150221 HT
                bestobj=cortemp*-1;
                pixshifttemp=[xstep ystep];
            end
        end
    end
%    currentobj=corr2(imbasesm, imtempsm);  % no more calculation of currentobj, use the old one (in last coarse shift)
    pixshift=pixshifttemp;
%	imresult=OIShift(imshift, pixshift(1), pixshift(2));
    return;
end

if objfun==3    % stdev of difference
    bestobj=9999999999;
    precision=1;    % modify here if need higher precision
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    for xstep=range(1):precision:range(2)    % x shifting
        for ystep=range(3):precision:range(4)   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            diffimg=imbasesm-imtempsm;
            objtemp=std2(diffimg);
            if (objtemp<bestobj)
                bestobj=objtemp;
                pixshift=[xstep ystep];
            end
            if xstep==0 & ystep==0
                currentobj=cortemp;
            end
        end
    end
%    imresult=OIShift(imshift, pixshift(1), pixshift(2));
    return;
end
         
if objfun==4    % minimum difference 
    precision=1;    % modify here if need improve precision
    maxobj=prod(basesize)*99999;   %estimate a large number
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    for xstep=range(1):precision:range(2)    % x shifting
        for ystep=range(3):precision:range(4)   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            diffimg=imbasesm-imtempsm;
            objtemp=sum(sum(abs(diffimg)));
            if (objtemp<bestobj)
                bestobj=objtemp;
                pixshift=[xstep ystep];
            end
            if xstep==0 & ystep==0
                currentobj=cortemp;
            end
        end
    end
%    imresult=OIShift(imshift, pixshift(1), pixshift(2));
    return;
end

% note: following method need check
if objfun==5	% use "image_registr_MI.m"
   %[h,im_matched, theta,I,J]=image_registr_MI(image1, image2, angle, step,crop);   % larger image2(imshift) is register to smaller image1(imbase)
    [entropy,im_matched, theta,yoffset,xoffset]=image_registr_MI(nc(imbasesm(range(4)+1:end+range(3), range(2)+1:end+range(1))), nc(imshift), 0, step1, 0); % angle 0, no crop
    pixshift(1)=xoffset-range(2)-1;
    pixshift(2)=yoffset-range(4)-1;
    bestobj=max(entropy(:));
    currentobj=entropy(range(2), range(4));
%    imresult=OIShift(imshift, pixshift(1), pixshift(2));
	return    
end

printf('Please check shift objfun!\r');
return;