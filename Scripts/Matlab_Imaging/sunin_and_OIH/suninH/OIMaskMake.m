% maskmake2.m
% modified from maskmake.m
% make domain masks (dots with different radius, crosses) based on xy accordinates
% use paste method to accelete computing
% modified so it can be used as a sub-function.

function maskmap=maskmake(manual, xs, ys, pixradius, dotnum, xdim, ydim)
% manual = 0, when be called, =1 manually input parameters.

if manual~=0
pixradius=5.0;
xdim=504;
ydim=504;
drivername='c:\';
foldername='_matlab\FrameView\';
mapname='thin';
fidx=fopen(strcat(drivername, foldername, mapname, 'x.txt'), 'r');
fidy=fopen(strcat(drivername, foldername, mapname, 'y.txt'), 'r');
xs=fscanf(fidx, '%f');    % read in x and y accordinates
ys=fscanf(fidy, '%f');
dotnum=size(xs,1);
fprintf('total dot number: %d\r', dotnum);
fclose(fidx);
fclose(fidy);
end

maskmap=zeros(ydim, xdim);

% make a circular patch
patch=zeros(2*pixradius+1, 2*pixradius+1);
for x=1:2*pixradius+1
    for y=1:2*pixradius+1
        if ((x-pixradius-1)*(x-pixradius-1)+(y-pixradius-1)*(y-pixradius-1))<=(pixradius*pixradius)
        	patch(x, y)=1;
        end
    end
end

% make a cross patch
%	    if ((abs(i-y(k))<=pixradius2)&(abs(j-x(k))<=1))|((abs(j-x(k))<=pixradius2)&(abs(i-y(k))<=1))
%    		map4(i, j)=255;
%	    end

%paste patchs on map
for k=1:dotnum
    x0=xs(k)-pixradius;
    y0=ys(k)-pixradius;
    if (x0>0 & y0>0 & (x0+2*pixradius)<=xdim & (y0+2*pixradius)<=ydim)
        maskmap(y0:y0+2*pixradius, x0:x0+2*pixradius)=maskmap(y0:y0+2*pixradius, x0:x0+2*pixradius)|patch(:,:);
    else
        fprintf('\rcoordinates #%d(%d, %d)+-%d exceed map dimention... continue anyway\r\r', k, xs(k), ys(k), pixradius);
    end
end
if manual==1
    maskmap=maskmap*255;
    imwrite(uint8(maskmap), strcat(drivername, foldername, mapname, '_1.bmp'), 'bmp');
end

return