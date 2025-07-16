function pixcourse=OIReadPix(anapar, foldername, filename, stim, y, x, kernel);

% pixcourse=OIReadPix(anapar, foldername, filename, stim, y, x, kernel)
% Read one pixel from one stim, (include all frames) from all blocks (used for ttest/tmap)
% input: anapar: block info, get from OIHeadRead
%		 foldername: the folder where blocks locates
%		 filename: array of block file names including folder
%		 stim: stim # 
%		 y, x: pixel location
% 		 kernel: usually 1, 3, 5, 7..., optional, low pass filtering for that pixel, use convolution, need improve
% output: pixcourse(frames, block) contains timecourse of the pixel (y, x)

if nargin==5
	kernel=[1];
else
	kernel=ones(kernel)/kernel/kernel;
end
blockfilenum=size(filename, 1); 

Width = anapar.FrameWidth;
Height = anapar.FrameHeight;
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
DataType = anapar.DataType;
system=anapar.system;

for k=1:blockfilenum
	cfilename=[foldername, getfield(cell2struct(filename(k), 'junk'), 'junk')];
	if kernel==1
	 	fid=fopen(cfilename,'rb','ieee-le'); % open block file
        switch system
        case 'r'    % trace by trace style
		    skip= (Width*Height*NFrames+8*NFrames+Width*Height+8)*(stim-1)+(y-1)*Width*NFrames+(x-1)*NFrames;
		    fseek(fid,5120+skip*2,-1);
		    pixcourse(1:NFrames, k)=fread(fid,NFrames,'uint16=>double'); 
        case 'v'
		 	for j=1:NFrames
			    skip=(stim-1)*Width*Height*NFrames + (j-1)*Width*Height + (y-1)*Width + x-1;
		    	switch DataType
			    case 12
				    fseek(fid,1716+skip*2,-1);
				    pixcourse(j, k)=fread(fid,1,'uint16=>double'); 
			    case 13
					fseek(fid,1716+skip*4,-1);
			    	pixcourse(j, k)=fread(fid,1,'uint32=>double'); 
			    case 14
				    fseek(fid,1716+skip*4,-1);
				    pixcourse(j, k)=fread(fid,1,'float32=>double'); 
			    otherwise
			    	fprintf('Error: data type must be 12, 13 or 14!\n');
			    end
		    end
        end
		fclose(fid);
	else
%		OneStim=OIReadStim(cfilename, stim, 'v');
	 	for j=1:NFrames
			OneFrame=OIReadFrame(cfilename, 'v', stim, j);
			ksize=floor(kernel/2);
			pixcourse(j,k)=sum(sum(OneFrame(y-ksize:y+ksize, x-ksize:x+ksize).*kernel));
%			pixcourse(j,k)=sum(sum(OneStim(y-ksize:y+ksize, x-ksize:x+ksize, j).*kernel));
	    end
	end
end