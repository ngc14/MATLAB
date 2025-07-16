function AllFrames=OIReadBlk(blkfilename, system)

% For OI data processing, read one block into a 4-D array
% AllFrames=OIReadBlk(blkfilename, system)

% Input: 'blkfilename' including the path, 'system' is either 'v' or 'r' ('r' need implement)
% Return: 'AllFrames', a 4-D array, in order of [Height, Width, Nframe/stim, Nstim]

anapar=OIHeadRead(blkfilename, system); % read header info
Width = anapar.FrameWidth;
Height = anapar.FrameHeight;
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
DataType = anapar.DataType;

fid=fopen(blkfilename,'rb','ieee-le'); % open block file
if system=='v'
    fseek(fid,1716,-1); % skip 1716 byte header for VDAQ system
    % read whole file from raw block file, need modify if memory size becomes a problem
    switch DataType
    case 12
    	AllFrames=fread(fid,[Width*Height*NFrames*NConds],'uint16=>double'); 
    case 13
    	AllFrames=fread(fid,[Width*Height*NFrames*NConds],'uint32=>double'); 
    case 14
    	AllFrames=fread(fid,[Width*Height*NFrames*NConds],'float32=>double'); 
    otherwise
    	fprintf('Error: data type must be 12, 13 or 14!\n');
    end
    % reshape vector into 4D array
    AllFrames=reshape(AllFrames,Width,Height,NFrames,NConds);
    AllFrames=permute(AllFrames, [2,1,3,4]);  
elseif system=='r'
    fseek(fid,5120,-1); % skip 5120 byte header for RedShirt system
    %% For RedShirt file:
    % <need add>
    fprintf('\this function not implemented yet\r');
else 
    fprintf('\rspecify system please\r');
end
fclose(fid);
return