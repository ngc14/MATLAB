% This file is a section of Sunin code for OI data process
% Purpose: Perform filtering on final block-averaged maps (subtraction maps (and single condition maps)). 
% 
% HDL 060324

if flaglpfilter
	for k=NStim+1:mapnum
		fprintf('\rlow-pass filtering (%s) map: %d', lpfmethod, k-NStim);
        switch lpfmethod 
            case 'fastmean'
        		maps(:,:,k) = conv2(maps(:,:,k), fspecial('disk', floor(lpkernelsize/2)), 'same');
            case 'slowmean'
        		maps(:,:,k) = OIMeanFilt(maps(:,:,k), lpkernelsize);
            case 'gaussian'
        		maps(:,:,k) = conv2(maps(:,:,k), fspecial('gaussian', lpkernelsize, floor(lpkernelsize)), 'same');
            case 'fastmedian'
        		maps(:,:,k) = medfilt2(maps(:,:,k), [lpkernelsize, lpkernelsize], 'symmetric');
            otherwise
                fprintf('Error: please specify filting method, no filtering performed\r');
        end
    end
end
if flaghpfilter
	for k=NStim+1:mapnum
		fprintf('\rhigh-pass filtering (%s) map: %d', hpfmethod, k-NStim);
        switch hpfmethod 
            case 'fastmean'
        		maps(:,:,k) = maps(:,:,k)-conv2(maps(:,:,k), fspecial('disk', floor(hpkernelsize/2)), 'same');
            case 'slowmean'
        		maps(:,:,k) = maps(:,:,k)-OIMeanFilt(maps(:,:,k), hpkernelsize);
            case 'gaussian'
        		maps(:,:,k) = maps(:,:,k)-conv2(maps(:,:,k), fspecial('gaussian', hpkernelsize, floor(hpkernelsize)), 'same');
            case 'fastmedian'
        		maps(:,:,k) = maps(:,:,k)-medfilt2(maps(:,:,k), [hpkernelsize, hpkernelsize], 'symmetric');
            case 'ribot'
        		maps(:,:,k) = maps(:,:,k)-RibotFilter(maps(:,:,k), hpkernelsize);   % note hpkernelsize here is actually order of polynomial function, usually 2 or 3
            otherwise
                fprintf('Error: please specify filting method, no filtering performed\r');
        end
    end
end

