function AllFrames=OIReadStim(blkfilename, stim, system)

% For OI data processing, read one stim from one block into a 3-D array
% AllFrames=OIReadStim(blkfilename, stim, system)
% Input: 'blkfilename' including the path, 'system' is either 'v' or 'r' ('r' need implement), stim is the stim number wish to read
% Return: 'AllFrames', a 3-D array, in order of [Height, Width, Nframe/stim]

anapar=OIHeadRead(blkfilename, system); % read header info
Width = anapar.FrameWidth;
Height = anapar.FrameHeight;
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
DataType = anapar.DataType;

if stim>NConds
    fprintf('\rstim wish to read exceed cond num in OIReadStim!\r');
    return;
end
fid=fopen(blkfilename,'rb','ieee-le'); % open block file
if system=='v'
    fseek(fid,1716,-1); % skip 1716 byte header for VDAQ system
    switch DataType
    case 12
	    fseek(fid, Width*Height*NFrames*(stim-1)*2, 'cof');
	    AllFrames=fread(fid,[Width*Height*NFrames],'uint16=>double'); 
    case 13
	    fseek(fid, Width*Height*NFrames*(stim-1)*4, 'cof');
	    AllFrames=fread(fid,[Width*Height*NFrames],'uint32=>double'); 
    case 14
	    fseek(fid, Width*Height*NFrames*(stim-1)*4, 'cof');
	    AllFrames=fread(fid,[Width*Height*NFrames],'float32=>double'); 
    otherwise
    	fprintf('Error: data type must be 12, 13 or 14!\n');
    end
    % reshape vector into 4D array
    AllFrames=reshape(AllFrames,Width,Height,NFrames);
    AllFrames=permute(AllFrames, [2,1,3]);  
elseif system=='r'
    fseek(fid,5120,-1); % skip 5120 byte header for RedShirt system
    fseek(fid, 2*(Width*Height*NFrames+8*NFrames+Width*Height+8)*(stim-1), 'cof');
    AllFrames=fread(fid,[Width*Height*NFrames],'uint16=>double'); 
    %tBNCFrames=fread(fid,[8*FramesPerStim], 'uint16');   %NBNC=8 Read 8 BNCs
    %tDarkFrames=fread(fid,[FrameWidth*FrameHeight],'uint16');    % read dark frames
    %tBNCDark=fread(fid, [8], 'uint16');  % read BNC dark 
    AllFrames=reshape(AllFrames, NFrames, Width, Height); % reshape into 3-D matrix
    AllFrames=permute(AllFrames, [3,2,1]);    % so the dimension order is (height, width, frame)
else 
    fprintf('\rspecify system please\r');
end
fclose(fid);
return