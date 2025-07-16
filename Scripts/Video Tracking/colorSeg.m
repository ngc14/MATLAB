function [seg, cMap] = colorSeg(im, cMark, colors)
hScale = 1.4;
sScale = .1;
vScale = 0;
% larger std scale -> more lenient
stdDistScale = .2;

% converting to hsv space
im = rgb2hsv(im);
h=im(:,:,1); s=im(:,:,2); v=im(:,:,3);

cMark = cMark(colors,:);
numColor = size(cMark,1);
% classify every pixel using nearest neighbor rule
distance = zeros([size(im,1), size(im,2), numColor]);
cMap = zeros([numColor, 3]);
for i=1:numColor
    dx=double(h)-cMark(i,1);  dy=double(s)-cMark(i,2); dz=double(v)-cMark(i,3);
    % tune parameters based on color of markers
    distance(:,:,i) = sqrt(hScale*dx.^2 + sScale*dy.^2 + vScale*dz.^2);
    cMap(i,:) = [cMark(i,1); cMark(i,2); cMark(i,3)];
end
cLabels=1:numColor;
[value,index]=min(distance,[],3);
label=cLabels(index);

% get colormap from image for new labels
cMap(end+1,:) = [0, 0, 0];
cMap = circshift(cMap,1);

% top % of nearest neighbor distances from minimum distance value
mask = value > min(min(value)) + stdDistScale*std2(value);
label(ind2sub(size(value),find(mask))) = 0;
bwcc = label > 0;

filled = imfill(bwcc, 'holes');
%smooth = smoothInterp(filled,label);
cleaned = bwareaopen(filled, 5);

seg = zeros(size(im,1), size(im,2),3);
% convert segmented image to correct colors
for r = 1:size(im,1)
    for c = 1:size(im,2)
        seg(r,c,:) = cMap(cleaned(r,c)+1,:);
    end
end
seg = hsv2rgb(seg) .* 255;
end