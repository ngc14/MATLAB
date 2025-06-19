function imresult=OIShift(imin, dx, dy);

%  imresult = OIShift(imin, dx, dy)
%  imin: input image
%  dx, dy: shift distance in x and y dimension 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using linear interpretation for non-integer shifts 

imsize=size(imin);
if size(imsize,2)>2
    fprintf('OIShift not apply to >2 dimension figure\n');
    return
end
%x1=zeros(imsize); 
%x2=zeros(imsize); 
x1=zeros(imsize)+ median(median(imin));  % change to this one on 8/09/04
x2=zeros(imsize)+ median(median(imin)); 

if dx~=0
    imresult=imin;
    if (dx>0)
        intstep = floor(dx);    % e.g. dx=3.3, intstep=3
        xresidue = dx - intstep; % residual = 0.3
        x1(:,(intstep+1):end)=imin(:,1:(end-intstep)); %shift by 3 
        x2(:,(intstep+2):end)=imin(:,1:(end-intstep-1)); %shift by 4 
        imresult=x1*(1-xresidue)+x2*xresidue; % 0.7x1 + 0.3x2
    else %(dx<0)
        dx=abs(dx);
        intstep = floor(dx);
        xresidue = dx - intstep;
        x1(:,1:(end-intstep))=imin(:,(intstep+1):end);    
        x2(:,1:(end-intstep-1))=imin(:,(intstep+2):end);  
        imresult=x1*(1-xresidue)+x2*xresidue; 
    end
    imin=imresult;
end
x1=zeros(imsize)+ median(median(imin));  
x2=zeros(imsize)+ median(median(imin)); 
if (dy>0)
    intstep = floor(dy);
    yresidue = dy - intstep;
    x1((intstep+1):end, :)=imin(1:(end-intstep), :); 
    x2((intstep+2):end, :)=imin(1:(end-intstep-1), :); 
    imresult=x1*(1-yresidue)+x2*yresidue; 
elseif (dy<0)
    dy=abs(dy);
    intstep = floor(dy);
    yresidue = dy - intstep;
    x1(1:(end-intstep), :)=imin((intstep+1):end, :); 
    x2(1:(end-intstep-1), :)=imin((intstep+2):end, :); 
    imresult=x1*(1-yresidue)+x2*yresidue; 
end
if (dx==0)&(dy==0)
    imresult=imin;
end
return