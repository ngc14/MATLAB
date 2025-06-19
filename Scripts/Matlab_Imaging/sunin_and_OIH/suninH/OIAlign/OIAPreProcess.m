function B=OIAPreProcess(A, illubk, LPKernel, shiftxy, croprange)

% Used in OI Align, some pre-process before alignment
% B=OIAPreProcess(A, bk, LFP, croprange, shiftxy)

% input
%   A: original image
%   illubk: illumination bkground to subtract
%   LPKernel: low-filter kernel, e.g. 2
%   croprange: x1 y1 x2 y2 of croprange (points included)
%   shiftxy: shift dx dy

if illubk~=0
    A=A-illubk;
end
if LPKernel~=0
    B=conv2(A, fspecial('Gaussian', LPKernel, LPKernel), 'same'); %lpk is dameter
end
if shiftxy(1)~=0 | shiftxy(2)~=0
    B=OIShift(B, shiftxy(1), shiftxy(2));
end
if croprange~=[0 0 0 0]
    B=B(croprange(2):croprange(4), croprange(1):croprange(3));
end

return;