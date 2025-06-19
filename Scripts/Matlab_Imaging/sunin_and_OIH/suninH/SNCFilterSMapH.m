% This file is a section of Sunin code for OI data process
% Purpose: Perform filtering on final block-averaged maps (subtraction maps (and single condition maps)). 
% 
%   (080229) Hisashi: startmap
%   (080129) Hisashi: Modified for mask filtering
%   HDL 060324
if flagmap==1 % by Hisashi on 080229
    startmap=NStim+1;
else
    startmap=1;
end
if flagmaskfilter==0
    if flaglpfilter
        for k=startmap:mapnum
            fprintf('\rlow-pass filtering (%s) map: %d', lpfmethod, k-startmap);
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
        for k=startmap:mapnum
            fprintf('\rhigh-pass filtering (%s) map: %d', hpfmethod, k-startmap);
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
else
    A_mask=maps(:,:,k);
    filtermaskname = strcat(resultdriver, 'expt1\', expname, 'masks\filtermask\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2
    A_mask = (double(imread(filtermaskname, 'bmp'))/255-1);
    fprintf('max:%g min: %g', max(max(A_mask)), min(min(A_mask)));
    fprintf('max:%g min: %g', max(max((1+conv2(A_mask, fspecial('gaussian', hpkernelsize, floor(hpkernelsize)), 'same')))), min(min((1+conv2(A_mask, fspecial('gaussian', hpkernelsize, floor(hpkernelsize)), 'same')))));
    
    if flaglpfilter
        for k=startmap:mapnum
            fprintf('\rlow-pass filtering (%s) map with mask: %d', lpfmethod, k-startmap);
            switch lpfmethod 
                case 'fastmean'
                    maps(:,:,k) = conv2(maps(:,:,k), fspecial('disk', floor(lpkernelsize/2)), 'same').*(1+conv2(A_mask, fspecial('disk', floor(lpkernelsize/2)), 'same'));
                case 'slowmean'
                    maps(:,:,k) = OIMeanFilt(maps(:,:,k), lpkernelsize).*(1+OIMeanFilt(A_mask, lpkernelsize));
                case 'gaussian'
                    maps(:,:,k) = conv2(maps(:,:,k), fspecial('gaussian', lpkernelsize, floor(lpkernelsize)), 'same').*(1+conv2(A_mask, fspecial('gaussian', lpkernelsize, floor(lpkernelsize)), 'same'));
                case 'fastmedian'
                    maps(:,:,k) = medfilt2(maps(:,:,k), [lpkernelsize, lpkernelsize], 'symmetric').*(1+medfilt2(A_mask, [lpkernelsize, lpkernelsize], 'symmetric'));
                otherwise
                    fprintf('Error: please specify filting method, no filtering performed\r');
            end
        end
    end
    if flaghpfilter
        for k=startmap:mapnum
            fprintf('\rhigh-pass filtering (%s) map with mask: %d', hpfmethod, k-startmap);
            switch hpfmethod 
                case 'fastmean'
                    maps(:,:,k) = maps(:,:,k)-conv2(maps(:,:,k), fspecial('disk', floor(hpkernelsize/2)), 'same').*(1+conv2(A_mask, fspecial('disk', floor(hpkernelsize/2)), 'same'));
                case 'slowmean'
                    maps(:,:,k) = maps(:,:,k)-OIMeanFilt(maps(:,:,k), hpkernelsize).*(1+OIMeanFilt(A_mask, hpkernelsize));
                case 'gaussian'
                    maps(:,:,k) = maps(:,:,k)-conv2(maps(:,:,k), fspecial('gaussian', hpkernelsize, floor(hpkernelsize)), 'same').*(1+conv2(A_mask, fspecial('gaussian', hpkernelsize, floor(hpkernelsize)), 'same'));
                case 'fastmedian'
                    maps(:,:,k) = maps(:,:,k)-medfilt2(maps(:,:,k), [hpkernelsize, hpkernelsize], 'symmetric')-maps(:,:,k).*(1+medfilt2(A_mask, [hpkernelsize, hpkernelsize], 'symmetric'));
                case 'ribot'
                    maps(:,:,k) = maps(:,:,k)-RibotFilter(maps(:,:,k), hpkernelsize)-maps(:,:,k).*(1+RibotFilter(A_mask, hpkernelsize));   % note hpkernelsize here is actually order of polynomial function, usually 2 or 3
                otherwise
                    fprintf('Error: please specify filting method, no filtering performed\r');
            end
        end
    end
end
