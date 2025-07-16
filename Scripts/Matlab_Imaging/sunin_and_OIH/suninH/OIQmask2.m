function values=OIQmask2(data, xarray, yarray, diameter, mask)
% Get pixel values from matrix 'data' using coordinates (x,y) 
% Each location select a circular patch with specified 'diameter'
% Check if location is within mask befor calculating
% Return averaged pixel value for each individual domain.

% Different from OIQmask1: OIQmask1 returns pixel value

% input
%   xarray/yarray: nx1 arrays
%   diameter: 1 3 5 7

% output:
%   values: 1 x n array, n=size(x)

%Need add: deal with x/y outside range

manual=0;
if manual
	data=OIReadIVF('08_V8.0.ivf');
	xtemp=textread('blobx.txt', '%d');
	ytemp=textread('bloby.txt', '%d');
	diameter=10;
	mask=imread('v1mask.bmp');
else
    xtemp=xarray;
    ytemp=yarray;
    if nargin<5
        mask=ones(504);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndomain=size(xarray, 1);
radius=ceil(diameter/2);
imsize=size(data);

% check location
n=0;
for i=1:Ndomain
    if (xtemp(i)-radius)>=1 & (xtemp(i)+radius)<=imsize(2) & (ytemp(i)-radius)>=1 & (ytemp(i)+radius)<=imsize(1)  % check if xy is with in frame, this is necessory, otherwise next check will have problem
        if mask(ytemp(i), (xtemp(i)-radius))==1 & mask(ytemp(i), (xtemp(i)+radius))==1 & mask((ytemp(i)-radius), xtemp(i))==1 & mask((ytemp(i)+radius), xtemp(i))==1 % if spot is within v1
            n=n+1;
            xarray(n,1)=xtemp(i);
            yarray(n,1)=ytemp(i);
        end
    end
end
if manual
    fprintf('%d out of %d domain within range', n, Ndomain);
end
Ndomain=n;

% create patch
center.x=floor(radius)+1;
center.y=floor(radius)+1;
patch=zeros(diameter);
for x=1:diameter
    for y=1:diameter
        if (x-center.x)^2+(y-center.x)^2<radius^2
            patch(y,x)=1;
        end
    end
end

% get value
pixnum=sum(sum(patch));
for i=1:Ndomain
    x0=xarray(i)-center.x;
    y0=yarray(i)-center.y;
    temp=data(y0:y0+diameter-1, x0:x0+diameter-1).*patch;
    values(i)=sum(sum(temp))/pixnum;
end

return
    