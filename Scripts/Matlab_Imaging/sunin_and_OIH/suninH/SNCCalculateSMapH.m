% This file is a section of Sunin code for OI data process
% Purpose: Constructing SMap from single condition maps.
% 
% HDL 060324

for j=1:size(Smap, 1)
    Sum1=zeros(FrameHeight, FrameWidth);
    Sum2=zeros(FrameHeight, FrameWidth);
    %
    condlist1=getfield(cell2struct(Smap(j, 2), 'junk'), 'junk');
    for n=condlist1
        Sum1=Sum1+maps(:,:,n);
    end
    Sum1=Sum1/size(condlist1,2);
    %
    condlist2=getfield(cell2struct(Smap(j, 3), 'junk'), 'junk');
    for n=condlist2
        Sum2=Sum2+maps(:,:,n);
    end
    if size(condlist2,2)~=0
        Sum2=Sum2/size(condlist2,2);
    end
    %
    if exist('flagblanksubtraction','var') == 0
        flagblanksubtraction=0;
    end
    if flagblanksubtraction
        Sum1BSub=zeros(FrameHeight, FrameWidth);
        Sum2BSub=zeros(FrameHeight, FrameWidth);
        condlistBSub1=getfield(cell2struct(SmapBSub(j, 2), 'junk'), 'junk');
        condlistBSub2=getfield(cell2struct(SmapBSub(j, 3), 'junk'), 'junk');
        for n=condlistBSub1
            Sum1BSub=Sum1BSub+maps(:,:,n);
        end
        Sum1BSub=Sum1BSub/size(condlistBSub1,2);        
        for n=condlistBSub2
            Sum2BSub=Sum2BSub+maps(:,:,n);
        end
        if size(condlistBSub2,2)~=0
            Sum2BSub=Sum2BSub/size(condlistBSub2,2);
        end
        
        Sum1=Sum1-Sum1BSub;
        Sum2=Sum2-Sum2BSub;
    end
    
    maps(:,:,NStim+j)=Sum1-Sum2;
        
    clear condlist1;
    clear condlist2;
    if flagblanksubtraction
        clear condlistBSub1;
        clear condlistBSub2;
    end    
end