function imresult=OIShift3(imin, dx, dy);

% faster version of oishift2
%  same as OIShift but only shift by pixels 
%  imresult = OIShift(imin, dx, dy)
%  imin: input image
%  dx, dy: shift distance in x and y dimension 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using linear interpretation for non-integer shifts 

imresult=zeros(504, 504);

% calculate the location of up-left corner of the shifted image in both matrix
if dx<0
	ax=-dx+1;
	bx=1;
else
	ax=1;	
	bx=dx+1;
end
if dy<0
	ay=-dy+1;
	by=1;
else
	ay=1;
	by=dy+1;
end	
imresult(by:by+dim2-abs(dy)-1, bx:bx+dim1-abs(dx)-1)=imin(ay:ay+dim2-abs(dy)-1, ax:ax+dim1-abs(dx)-1);
return