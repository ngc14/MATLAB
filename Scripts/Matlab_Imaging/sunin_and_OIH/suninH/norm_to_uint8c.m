% normalize image I to [0,255] then convert it to uint8
% different from 'norm_to_uint8', this one user define the clipsd (median still from map)

function J = norm_to_uint8(I, clipsd)

Imedian=median(I(:));
m = Imedian-clipsd;
M = Imedian+clipsd;

I = floor(255*(I-m)/(M - m)); % normalize to [0-255] 
J = uint8(I);                 % convert to uint8
return
