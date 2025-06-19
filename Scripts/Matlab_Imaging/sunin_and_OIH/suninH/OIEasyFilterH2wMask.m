function B=OIEasyFilterH2wMask(A, LPMethod, LPKernel, HPMethod, HPKernel, filtermask)
% function OIEasyFilter(LPFName, LPKernel, HPFName, HPKernel)
% kernel=0 means no filter
%
% 091231 Hisashi: 'fastmeanH' was added.
% 080118 maskfiltering function was added by Hisashi

B=A;
% filtermask=A;
%     filtermask = strcat(resultdriver, 'expt1\', expname, 'masks\filtermask\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2
% filtermaskname = strcat(resultdriver, outputfolder, expname, 'masks\filtermask\default.bmp');     % can use a blood vessel map for clip masking, only for clipmethod==2
% filtermask = (double(imread(filtermaskname, 'bmp'))/255-1);
% filtermasktemp = (double(imread(filtermaskname, 'bmp'))/255);
% if size(filtermasktemp,3)==1
%     fprintf('size 1');
%     filtermask = filtermasktemp;
% else
%     fprintf('size 3');
%     filtermask = filtermasktemp(:,:,1);
% end
% clear filtermasktemp;
% filtermask = filtermask*(-1)+1;
% A_mask = filtermask;

% fprintf('\rfiltermask: max:%g min: %g', max(max(filtermask)), min(min(filtermask)));
% fprintf('\rA: max:%g min: %g', max(max(imfilter(A, fspecial('disk', floor(LPKernel/2)), 'replicate'))), min(min(imfilter(A, fspecial('disk', floor(LPKernel/2)), 'replicate'))));

% filtermask = double(filtermask)/255-1;
if LPKernel~=0
    switch LPMethod
        case 'fastmean'
            filtermask2 = imfilter(filtermask, fspecial('disk', round(LPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);
            A_masked = imfilter(A.*filtermask, fspecial('disk', round(LPKernel/2))); 
            A_masked = A_masked./filtermask3;
            B = A_masked.*double(filtermask2 ~= 0)+A.*double(filtermask2 == 0);
        case 'fastmeanH'
            filtermask2 = imfilter(filtermask, fspecial('disk', round(LPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);

            C = imresize(A, [size(A)]*0.25, 'bilinear');
            sfiltermask  = imresize(filtermask, [size(A)]*0.25, 'bilinear');
            sfiltermask2 = imresize(filtermask2, [size(A)]*0.25, 'bilinear');
            sfiltermask3 = imresize(filtermask3, [size(A)]*0.25, 'bilinear');
            
            C_masked = imfilter(C.*sfiltermask, fspecial('disk', round(LPKernel/2*0.25))); 
            C_masked = C_masked./sfiltermask3;
            C = C_masked.*double(sfiltermask2 ~= 0)+C.*double(sfiltermask2 == 0);
            
            B = imresize(C, [size(A)], 'bilinear');
        case 'gaussian'
            filtermask2 = imfilter(filtermask, fspecial('gaussian', LPKernel, round(LPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);
            A_masked = imfilter(A.*filtermask, fspecial('gaussian', LPKernel, round(LPKernel/2))); 
            A_masked = A_masked./filtermask3;
            B = A_masked.*double(filtermask2 ~= 0)+A.*double(filtermask2 == 0);
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
        case 'fastmean'
            filtermask2 = imfilter(filtermask, fspecial('disk', round(HPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);
            A_masked = imfilter(A.*filtermask, fspecial('disk', round(HPKernel/2))); 
            A_masked = A_masked./filtermask3;
            B = B - A_masked.*double(filtermask2 ~= 0)+A.*double(filtermask2 == 0);
        case 'fastmeanH'
            filtermask2 = imfilter(filtermask, fspecial('disk', round(HPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);

            C = imresize(A, [size(A)]*0.25, 'bilinear');
            sfiltermask  = imresize(filtermask, [size(A)]*0.25, 'bilinear');
            sfiltermask2 = imresize(filtermask2, [size(A)]*0.25, 'bilinear');
            sfiltermask3 = imresize(filtermask3, [size(A)]*0.25, 'bilinear');
            
            C_masked = imfilter(C.*sfiltermask, fspecial('disk', round(HPKernel/2*0.25))); 
            C_masked = C_masked./sfiltermask3;
            C = C_masked.*double(sfiltermask2 ~= 0)+C.*double(sfiltermask2 == 0);
            
            B = B - imresize(C, [size(A)], 'bilinear');            
        case 'gaussian'
            filtermask2 = imfilter(filtermask, fspecial('gaussian', HPKernel, round(HPKernel/2)));
            filtermask3 = filtermask2+double(filtermask2 == 0);
            A_masked = imfilter(A.*filtermask, fspecial('gaussian', HPKernel, round(HPKernel/2))); 
            A_masked = A_masked./filtermask3;
            B = B - A_masked.*double(filtermask2 ~= 0)+A.*double(filtermask2 == 0);
        case 'median'
%             B = B-medfilt2H(A, [HPKernel, HPKernel], 'replicate'); %chnaged by Hisashi on 091009
            B = B-medfilt2(A, [HPKernel, HPKernel], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit            
        case {'fastmedian','fastmedianH'}
            C = imresize(A, [size(A)]*0.25, 'bilinear'); %[126, 126]
%             B = B-medfilt2H(A, [HPKernel, HPKernel], 'replicate'); %chnaged by Hisashi on 091009
            B = B-medfilt2(A, [HPKernel, HPKernel], 'symmetric'); %chnaged by Hisashi on 110921 for 64bit            
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