% normalize image I to [0,255] then convert it to uint8
function J = norm_to_uint8(I)

M = max(I(:));         % find the maximum intensity
m = min(I(:));         % find the minimum intensity

I = floor(255*(I-m)/(M - m)); % normalize to [0-255] 
J = uint8(I);                 % convert to uint8
return
