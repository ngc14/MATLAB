function processed = processIm(im, avg)
grey = rgb2gray(uint8(im));
grey = imfilter(grey, ones(3,3)/9);
%% background subtraction
motionA = abs(grey-avg);
motionB = abs(avg-grey);
combined = imfuse(motionA, motionB, 'blend');
% bigger scaling -> more selective
thresh = combined > .1*std2(combined);
%% cleaning
dilated = bwmorph(thresh, 'dilate');
% get rid of areas smaller than 1/100th of the entire area (NOT USED)
cleaned = bwareaopen(dilated, floor((size(im,1) * size(im,2))/100), 4);
% fill holes
fill = imfill(cleaned, 'holes');
mask = repmat(fill, [1,1,3]);
maskedIm = im .* mask;
processed = medfilt3(maskedIm, [7, 7, 1]);
end