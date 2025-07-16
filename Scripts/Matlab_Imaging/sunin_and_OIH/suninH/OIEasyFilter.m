function B=OIEasyFilter(A, LPMethod, LPFKernel, HPMethod, HPFKernel)
% function OIEasyFilter(LPFName, LPFKernel, HPFName, HPFKernel)
% kernel=0 means no filter

B=A;
if LPFKernel~=0
    switch LPMethod
        case 'fastmean'
            B = conv2(A, fspecial('disk', floor(LPFKernel/2)), 'same');
        case 'slowmean'
            B = OIMeanFilt(A, LPFKernel);
        case 'gaussian'
            B = conv2(A, fspecial('gaussian', LPFKernel, floor(LPFKernel)), 'same');
        case 'fastmedian'
            B = medfilt2(A, [LPFKernel, LPFKernel], 'symmetric');
        otherwise
            fprintf('Error: please specify filting method, no filtering performed\r');
    end
end
if HPFKernel~=0
    switch HPMethod
        case 'fastmean'
            B = B-conv2(A, fspecial('disk', floor(HPFKernel/2)), 'same');
        case 'slowmean'
            B = B-OIMeanFilt(A, HPFKernel);
        case 'gaussian'
            B = B-conv2(A, fspecial('gaussian', HPFKernel, floor(HPFKernel)), 'same');
        case 'fastmedian'
            B = B-medfilt2(A, [HPFKernel, HPFKernel], 'symmetric');
        case 'ribot'
            B = B-RibotFilter(A, HPFKernel);   % note HPFKernel here is actually order of polynomial function, usually 2 or 3
        otherwise
            fprintf('Error: please specify filting method, no filtering performed\r');
    end
end
return