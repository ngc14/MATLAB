% This file is a section of Sunin code for OI data process
% Purpose: Constructing SMap from single condition maps.
% 
% HDL 060324

for j=1:smapnum
    Sum1=zeros(FrameHeight, FrameWidth);
    Sum2=zeros(FrameHeight, FrameWidth);

	Sum1=squeeze(NewSum_Smap1(:,:,j)/counter_Smap1(j));

	Sum2=squeeze(NewSum_Smap2(:,:,j)/counter_Smap2(j));

    maps(:,:,NStim+j)=Sum1-Sum2;
        
end