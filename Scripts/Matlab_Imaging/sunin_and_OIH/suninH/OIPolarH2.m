function Polarmap=OIPolarH2(singlemap, lut, vectmask, pclipmethod, pclipvalue, LPfiltermethod, lowpass, HPfiltermethod, highpass)

% (100427 Hisashi) Bug fix in FFT filtering
% (091230) fixed in 'vectormask' function by Hisashi
% (071021) Modified by Hisashi
% (051007) Modified by HDL
% Polarmap=OIPolar(singlemap, lut, mask, clipsd, lowpass, highpass)
% Calculates polar map from a set of single/difference maps
% input: singlemap: single condition maps (height, width, anglenum)
%        lut: color lookuptable, if not defined, no color for polar
%		 mask: a mask used for normalization
% output: polarmap.ang: angle map (value: 1-256)
%         polarmap.mag: magnitude map (unit: same as input map)
%         polarmap.sum: a simple summation of all conditions (like allo or alld)
% Modified from Xiangming Xu 04 'Ivf2PolarMapsFIN4HD.m'


% try modify following flags to see the effect, which may vary from case to
% case
flagVectNorm=1; % whether or not use Bosking normalization, default 1, no use for 040113GarRun2
flagSubMap=0;   % whether or not use subtraction map (H V A O -> HV AO), default 0, no use for 040113GarRun2

%--------------------------------------
n=size(singlemap, 3);
[r, c] = size(singlemap(:, :, 1));
if nargin==1	
	lut=textread('bwfix.lut'); %% if no lookup table provided, use standard 
	vectmask=ones([r,c]);	%if no mask provided, use default: all map area. 
	pclipmethod=1;
    pclipvalue=0;
	lowpass=0;
	highpass=0;
end	
vectmaskpixnum=sum(sum(vectmask));

Polarmap.ang=zeros([r,c]);	% angle
Polarmap.mag=zeros([r,c]);	% magnitude
% Polarmap.sum=mean(singlemap, 3);	% simple average of all conditions (like lumo)

% Get difference map, according to Bosking 97, may not necessory
% skip this when n is a odd number
if ~mod(n, 2) && flagSubMap
    tempmap=singlemap;
    for i=1:n
        if i<=n/2
            singlemap(:,:,i)=tempmap(:,:,i)-tempmap(:,:,i+n/2);
        else
            singlemap(:,:,i)=tempmap(:,:,i)-tempmap(:,:,i-n/2);
        end	
    end
    clear tempmap;
end
fprintf('\r');

for k = 1:n;
	% clipping
	% low & high pass filtering	(note: this procedure take long time)
%     LPfiltermethod
%     lowpass
%     HPfiltermethod
%     highpass
%     flagmaskfilter

	singlemap(:,:,k)=OIEasyFilterH2(singlemap(:,:,k), LPfiltermethod, lowpass, HPfiltermethod, highpass);
	
	if strcmp(LPfiltermethod, 'fft') || strcmp(HPfiltermethod, 'fft')
		if r == c
			fftfilter=OIMakFIRfilterH(r,LPfiltermethod, lowpass, HPfiltermethod, highpass);
		else
			fprintf('Error: For FFT filtering, maps should be square.\n');
		end            
		singlemap(:,:,k)=ifft2(fft2(singlemap(:,:,k)).*fftfilter);
	end        

    if pclipvalue~=0
% 		singlemap(:,:,k)=OIClipH(singlemap(:,:,k), 1, clipsd); %Commeted out on 090118 by Hisashi
% 		singlemap(:,:,k)=OIClipH(singlemap(:,:,k), pclipmethod, pclipvalue); %by Hisashi on 090116; Commeted out on 090118
% 		singlemap(:,:,k)=OIClipH(singlemap(:,:,k), pclipmethod, pclipvalue);
	end
end

Polarmap.sum=mean(singlemap, 3);	% simple average of all conditions (like lumo)

% normalization	(using Bosking 97 method)
if flagVectNorm==1
    imgmean=zeros(n, 1);	% image mean for each image
    imgMAD=zeros(n, 1);	% Mean Absolute Deviation from the mean for each image
    OADmax=-999;		% OverAll max/min of the deviation, (across all maps)
    OADmin=999;
    for i=1:n
        vectmaskarray=reshape(vectmask, prod(size(vectmask)),1);
        [junk, vectmaskindex]=sort(vectmaskarray);
        imgarray=reshape(singlemap(:,:,i), prod(size(singlemap(:,:,i))), 1);
        b=imgarray(vectmaskindex(end-vectmaskpixnum+1:end));
        imgmean(i)=mean(b);
        imgMAD(i)=mad(b);
        cmax=max((b-imgmean(i))/imgMAD(i));
        cmin=min((b-imgmean(i))/imgMAD(i));
        if cmax>OADmax
            OADmax=cmax;
        end
        if cmin<OADmin
            OADmin=cmin;
        end	
    end
    imgscale=63.0/(OADmax-OADmin);
    for i=1:n
        b=(singlemap(:,:,i)-imgmean(i))/imgMAD(i);
        singlemap(:,:,i)=(b-OADmin)*imgscale;
    end
end

% calculate polar map
% fprintf('Creating polar image, vector calculating...\n');
for k = 1: n
    angle = (360/n)*(k-1);  %%double angle here  2* (angle/pi)*180
    xcomponent(k) = cos(angle*pi/180);
    ycomponent(k) = sin(angle*pi/180);
end
maxtotalmag =0;
magbuff = zeros(r,c);
preferor = zeros(r,c);
magor = zeros(r,c);
count=0;	% for online display calculating progress
for i = 1:r
	if (i/r)>count	% for online display calculating progress
%	    fprintf('%2.0f%%  ', i/r*100);
	    count=count+0.05;
	end
    for j = 1:c
         xmag = 0;
         ymag = 0;
         tmpsum=0;      
         for k = 1:n
             response = singlemap(i, j, k);
             xmag = xmag + response*xcomponent(k);
             ymag = ymag + response*ycomponent(k);
             tmpsum = tmpsum + response;
         end
         totalmag = sqrt((xmag*xmag) + (ymag*ymag));
%         if (totalmag > maxtotalmag)
%             maxtotalmag = totalmag;
%         end
         Polarmap.mag(i, j) = totalmag;
         tanangle = ymag/xmag;
         angle = atan(tanangle)*180/pi;
         if (xmag <0)
            angle = angle + 180; %%correction for anglequantrant misplacement
         end
         angle = angle + 180;
         if (angle >= 360)
            angle = angle -360;
         end
         Polarmap.ang(i, j) = floor((angle/360)*256);  %% dvided by 360 to half the angle
         Polarmap.ang(i, j) = angle;
    end  %%j
end  %%i
return;

% Note:
%Although input image OR is a range of 180 not 360 deg, atan calculation
%results are in the range of 360.
%So it needs to double each anlge to get rid of directionality before 
%compiling OR maps like, Batshelet 1981; Drgao et al., 2000; 
%Worgotter & Eysel, 1991; cos(2*ang) sin(2*ang);  resultant strength r=
%sqrt((sum(xcos 2*ang))^2 + (sum(ycos 2*ang))^2).  In the end, 
%resultant vector angle: atan((sum(ycos 2*ang))^2/(sum(xcos 2*ang))^2)/2  (doubletheta/2)