function B=OIMakFIRfilterH(newsize, LPMethod, LPcutoff, HPMethod, HPcutoff)
%
% 091221 Make a band-pass FIR filter for FFT. Hisashi

% [m,n]=size(A);
% if m~=n;
%     error('The matrix should be square.');
% end
% N=m;

[f1,f2] = freqspace(newsize,'meshgrid');
filter = ones(newsize);
r = sqrt(f1.^2 + f2.^2);

if strcmp(LPMethod, 'fft')
    filter((r>(1/newsize*2)*(newsize/LPcutoff))) = 0;
end
if strcmp(HPMethod, 'fft')
    filter((r<(1/newsize*2)*(newsize/HPcutoff))) = 0;
%     pause
end
%         colormap(jet(newsize))
%         mesh(f1,f2,filter)
%         pause;
% win = fspecial('gaussian',newsize,3);
% win = win ./ max(win(:));  % Make the maximum window value be 1.
% mesh(win)
% pause;
h = fwind2(filter,hamming(newsize));
% h = fwind2(filter,win);
i = freqz2(h,[newsize newsize]);
% freqz2(h)
% pause;
filter=fftshift(i);

B=filter;

return





% if LPcutoff~=0
%     switch LPMethod
%         case 'fft'
%     end
% end
% if HPcutoff~=0
%     switch HPMethod
%         case 'fastmean'
%             B = B-imfilter(A.*filtermask, fspecial('disk', round(HPcutoff/2)))./imfilter(filtermask, fspecial('disk', round(HPcutoff/2))); %By Hisashi 090116
%         case 'gaussian'
%             B = B-imfilter(A.*filtermask, fspecial('gaussian', HPcutoff, round(HPcutoff/2)))./imfilter(filtermask, fspecial('gaussian', HPcutoff, round(HPcutoff/2))); %By Hisashi 090116
%         case 'fft'
%         otherwise
%             fprintf('Error: Mask filtering is available only for "fastmean" and "gaussian" now'.\r');
%     end
% end
% return