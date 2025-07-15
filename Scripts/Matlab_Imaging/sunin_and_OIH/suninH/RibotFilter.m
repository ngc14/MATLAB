function approx_poly=RibotFilter(z, P)

% function B=RibotFilter(z, P)
% From Ribot (Ribot et al. 2005 JNM) 

% Original name: approx_poly = polyfit2d_corrected(x,y,z,P)
% Original description: 
% For the analysis of optical imaging signals, we are using matlab in our 
% laboratory. There is a very easy way with matlab for approximating the 
% recorded pictures by a polynomial functions using the last squares 
% method. You should use the command polyfit2d.
% http://web.ccr.jussieu.fr/ccr/Documentation/Calcul/matlab5v11/docs/ftp.mathworks.com/pub/contrib/v4/optim/polyfit2d.m 
% 
% You can define x and y as follows [x,y]=meshgrid(1:ydim,1:xdim) where 
% xdim and ydim are the dimensions of your pictures. These two matrices x 
% and y will represent the spatial coordinates from which the polynomial 
% functions will be calculated.
% 
% z represents the image you search to approximate.
% 
% n and m are the maximal coefficients of the polynomial functions in x 
% and y.
% 
% In our study we chose only one polynomial order P. And the order of the 
% polynomial functions P is such that, for each term x^m' .y^n' , 
% m'+n'<=P. That requires a few changes in the program.
% 
% I attach the corrected version of the program. The result of this 
% function approx_poly will be the polynomial approximation, as a vector. 
% In practice, P=2 or 3 was enough to obtain a satisfactory approximation; 
% the calculations are extremely quick.
% 
% In the hope that will help,
% Regards,
% 
% Jerome Ribot
% Laboratory for Visual Neurocomputing
% Brain Science Institute
% Riken

[height, width]=size(z);

[x,y]=meshgrid(1:height,1:width);

if any((size(x) ~= size(y)) | (size(z) ~= size(y)))
    error('X, Y,and Z vectors must be the same size')
end

x = x(:); y = y(:); z= z(:);
P = P + 1; 

%Calculation of the polynomial coefficents
clear polycoeff;
a=[];
index=1;
for i1=1:P
    j1=1;
    while(i1+j1<=P+1)
        a(:,index) = (x.^(j1-1)).*(y.^(i1-1));
        index=index+1;
        j1=j1+1;
    end
end
polycoeff = (a\z);

%Calculation of the approximated function 
approx_poly=zeros(size(x));
for index=1:length(polycoeff)
    approx_poly=approx_poly+polycoeff(index).*a(:,index);
end
approx_poly=reshape(approx_poly, [height, width]);
