function J = OIMeanFilt(I, kernelradius, skippix);

% J = OIMeanFilt(I, kernelradius, skippix);
% Returns same size matrix with average value within a kernel size (smoothing)
% I:  Input image (Float point format), 
% kernelradius: integer to decide the radius size, e.g. 3 X 3, 5 X 5.  
% J: result image (Float point format)
% skippix: for fasten calculation, =2 means skip every other pixel,

% modified from Xiangming Xu 04 'FloatMeanfiltering'

if nargin<2
	fprintf('Too few input arguments\n');
elseif nargin<3
	skippix=1;
end
[r,c] = size(I);
J = zeros(r,c);
radius = floor(kernelradius/2);
pixelnum=(2*radius+1)^2;
for i = 1:r
	if (~mod(i, 10))
%		fprintf('%d  ', i);
	end
	if (~mod(i, 100))
%		fprintf('\r');
	end
    for j = 1:c
    	%if (i-radius)<1 | (j-radius)<1 |(i+radius)>r |(j+radius)>r
	        numpix=0;
    	    sumpixval=0;
        	for padr = (i - radius):skippix:(i + radius)
            	for  padc = (j -radius):skippix:(j + radius)
                	if ~((padr <1) | (padc <1) | (padr > r) | (padc > c) )
                    	numpix = numpix + 1;
	                    sumpixval = sumpixval + I(padr, padc);
    	            end
        	    end
        	end	
	        J(i, j) = sumpixval/numpix;
        %else
%			temparray=reshape(I(i-radius:i+radius, j-radius:j+radius), pixelnum, 1);
        %	J(i, j) = mean2(I(i-radius:i+radius, j-radius:j+radius));
        %end
    end  %for j = 1:c
%	fprintf('\r');
end % for i = 1:r
return;