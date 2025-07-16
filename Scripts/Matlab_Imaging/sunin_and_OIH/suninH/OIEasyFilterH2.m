function B=OIEasyFilterH2(A, LPMethod, LPKernel, HPMethod, HPKernel)
% function OIEasyFilter(LPFName, LPKernel, HPFName, HPKernel)
% kernel=0 means no filter
%
% 091231 Hisashi: 'fastmeanH' was added.
% 090115 'medfilt2' was replaced by 'medfilt2H' with 'replicate'.
% 090115 'conv2' was replaced by 'imfilter' with 'replicate'.
% 080118 maskfiltering function was added by Hisashi

B=A;


if LPKernel~=0
    switch LPMethod
        case 'fastmean'
            B = imfilter(A, fspecial('disk', floor(LPKernel/2)), 'replicate'); % By Hisashi 090116
        case 'fastmeanH'
            C = imresize(A, [size(A)]*0.25, 'bilinear');
            C = imfilter(C, fspecial('disk', floor(LPKernel/2*0.25)), 'replicate'); % By Hisashi 090116
            B = imresize(C, [size(A)], 'bilinear');
        case 'slowmean'
            B = OIMeanFilt(A, LPKernel);
        case 'gaussian'
            B = imfilter(A, fspecial('gaussian', LPKernel, floor(LPKernel/2)), 'replicate'); %By Hisashi 090116
        case 'median'
%             B = medfilt2H(A, [LPKernel, LPKernel], 'replicate'); %chnaged by Hisashi on 091009
            B = medfilt2(A, [LPKernel, LPKernel], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit
        case {'fastmedian','fastmedianH'}
            C=imresize(A, [size(A)]*0.25, 'bilinear'); %[126, 126]
%             C=medfilt2H(C, [round(LPKernel*0.25), round(LPKernel*0.25)], 'replicate'); %chnaged by Hisashi on 091009
            C=medfilt2(C, [round(LPKernel*0.25), round(LPKernel*0.25)], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit
            B = imresize(C, [size(A)], 'bilinear');
        case 'ribot' % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM)
            B = RibotFilter(A, LPKernel);   % note LPKernel here is actually order of polynomial function, usually 2 or 3
        case 'log'
            B = imfilter(A, fspecial('log', LPKernel*3, floor(LPKernel/2)), 'replicate'); % By Hisashi 090116
        case 'fft'
        otherwise
            fprintf('Error: please specify LP filting method, no filtering performed\r');
    end
end
if HPKernel~=0
    switch HPMethod
        case 'fastmean' % 'fastmean' is a 'disk-like' mean filter, disk diameter is 2*floor(size/2)+1
            B = B-imfilter(A, fspecial('disk', floor(HPKernel/2)), 'replicate'); %By Hisashi 090116
        case 'fastmeanH'
            C = imresize(A, [size(A)]*0.25, 'bilinear'); %[126, 126]
            C = imfilter(C, fspecial('disk', floor(HPKernel/2*0.25)), 'replicate');
            B = B - imresize(C, [size(A)], 'bilinear');
        case 'slowmean'  % 'slowmean' 'disk-like' but better at edge
            B = B-OIMeanFilt(A, HPKernel);
        case 'gaussian' % 'gaussian' Gaussian filter with half sd
            B = B-imfilter(B, fspecial('gaussian', HPKernel, floor(HPKernel/2)), 'replicate'); %By Hisashi 090116
        case 'median'
%             B = B-medfilt2H(A, [HPKernel, HPKernel], 'replicate'); %chnaged by Hisashi on 091009
            B = B-medfilt2(A, [HPKernel, HPKernel], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit            
        case {'fastmedian','fastmedianH'}
            C = imresize(A, [size(A)]*0.25, 'bilinear'); %[126, 126]
%             C = medfilt2H(C, [round(HPKernel*0.25), round(HPKernel*0.25)], 'replicate'); %chnaged by Hisashi on 091009
            C = medfilt2(C, [round(HPKernel*0.25), round(HPKernel*0.25)], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit
            B = B - imresize(C, [size(A)], 'bilinear');
        case 'ribot' % 'ribot', fast fitting with polynomial surface (see Ribot et al. 2005 JNM)
            B = B - RibotFilter(A, HPKernel);   % note HPKernel here is actually order of polynomial function, usually 2 or 3
        case 'wholemedian'
            temp=reshape(B, prod(size(A)), 1);    %transform image to 1-D array for calculation of median and stdev
            frameMedian=double(median(temp));
            B = B-frameMedian;
        case 'fft'
        otherwise
            fprintf('Error: please specify HP filting method, no filtering performed\r');
    end
end
return