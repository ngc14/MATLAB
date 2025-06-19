function J = norm_to_01(I)
% normalize image I to [0, 1]

M = max(max(max(I)));         % find the maximum intensity
m = min(min(min(I)));         % find the minimum intensity
J = (I-m)/(M - m); % normalize to [0-1] 
return
