function imresult=OIShift2(imin, dx, dy);

%  same as OIShift but only shift by pixels 
%  imresult = OIShift(imin, dx, dy)
%  imin: input image
%  dx, dy: shift distance in x and y dimension 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using linear interpretation for non-integer shifts 

dim1=size(imin, 1);
dim2=size(imin, 2);
imresult=ones(dim1, dim2)*255;

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