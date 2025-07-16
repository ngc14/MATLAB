% This file is a section of Sunin code for OI data process
% Purpose: Constructing SMap from single condition maps.
% 
% HDL 060324
% HT 101111 

for j=1:size(Smap, 1)
    Sum1=zeros(FrameHeight, FrameWidth);
    Sum2=zeros(FrameHeight, FrameWidth);
    condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
    for n=condlist1
        Sum1=Sum1+maps(:,:,n);
    end
    Sum1=Sum1/size(condlist1,2);
    condlist2=getfield(cell2struct(Smap(j, 3), 'junk'), 'junk');
    for n=condlist2
        Sum2=Sum2+maps(:,:,n);
    end
    if size(condlist2,2)~=0
        Sum2=Sum2/size(condlist2,2);
    end
%     maps(:,:,NStim+j)=Sum1-Sum2; % HT 10111
    maps(:,:,NStim+j)=Sum1./Sum2; % HT 101111
end