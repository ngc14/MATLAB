% This file is a section of Sunin code for OI data process
% Purpose: Perform filtering on final block-averaged maps (subtraction maps (and single condition maps)). 
%
%   (091222) Hisashi: 'fft' filtering function was added.
%   (090204) Hisashi: Use OIEasyfilterH2
%   (080229) Hisashi: startmap
%   HDL 060324

if flagmap==1 % by Hisashi on 080229
    startmap=NStim+1;
else
    startmap=1;
end

for k=startmap:mapnum
    mapname=getfield(cell2struct(mapnames(k), 'junk'), 'junk');
	fprintf('spatial filtering: %s\r', mapname);
	maps(:,:,k)=OIEasyFilterH2(maps(:,:,k), lpfmethod, lpkernelsize, hpfmethod, hpkernelsize);

	if strcmp(lpfmethod, 'fft') || strcmp(hpfmethod, 'fft')
		maps(:,:,k)=ifft2(fft2(maps(:,:,k)).*fftfilter);
	end
end

