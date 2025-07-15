function OneFrame=OIReadFrame(blkfilename, system, stim, frame)

% OneFrame=OIReadFrame(blkfilename, system, stim, frame)
% For OI data processing, read one frame into a 2-D array 
% Input: 'blkfilename' including the path, 'system' is either 'v'(VDAQ), 'r'(RedShirt) or 'd'(Dan Ts'o)
% stim is stim condition number, frame is the frame number want to read
% Return: 'AllFrames', a 4-D array, in order of [Height, Width, Nframe/stim, Nstim]

anapar=OIHeadRead(blkfilename, system); % read header info
Width = anapar.FrameWidth;
Height = anapar.FrameHeight;
NFrames = anapar.FramesPerStim;
NConds = anapar.NStim;
DataType = anapar.DataType;

fid=fopen(blkfilename,'rb','ieee-le'); % open block file
if stim>NConds | frame>NFrames | stim<=0 | frame<=0
    fprintf('check ''stim'' or ''frame'' in reading the very first frame for shifting!\r');
    return
end
if system=='d'
    skipframes=(stim-1)*NFrames+frame-1;
    fseek(fid,skipframes*Width*Height*2,-1); 
    OneFrame=fread(fid,[Width Height],'uint16=>double'); 
elseif system=='v'
    skipframes=(stim-1)*NFrames+frame-1;
    switch DataType
    case 12
	    fseek(fid,1716+skipframes*Width*Height*2,-1); % skip 1716 byte header for VDAQ system
	    OneFrame=fread(fid,[Width Height],'uint16=>double'); 
    case 13
	    fseek(fid,1716+skipframes*Width*Height*4,-1); % skip 1716 byte header for VDAQ system
	    OneFrame=fread(fid,[Width Height],'uint32=>double'); 
    case 14
	    fseek(fid,1716+skipframes*Width*Height*4,-1); % skip 1716 byte header for VDAQ system
	    OneFrame=fread(fid,[Width Height],'float32=>double'); 
    otherwise
    	fprintf('Error: data type must be 12, 13 or 14!\n');
    end
    OneFrame=permute(OneFrame, [2,1]);  
elseif system=='r'      % for redshirt file, has to read the whole stim then select one frame (since trace-by-trace style)
    fseek(fid,5120,-1); % skip 5120 byte header for RedShirt system
    fseek(fid, 2*(Width*Height*NFrames+8*NFrames+Width*Height+8)*(stim-1), 'cof');
    AllFrames=fread(fid,[Width*Height*NFrames],'uint16=>double'); 
    %tBNCFrames=fread(fid,[8*FramesPerStim], 'uint16');   %NBNC=8 Read 8 BNCs
    %tDarkFrames=fread(fid,[FrameWidth*FrameHeight],'uint16');    % read dark frames
    %tBNCDark=fread(fid, [8], 'uint16');  % read BNC dark 
    AllFrames=reshape(AllFrames, NFrames, Width, Height); % reshape into 3-D matrix
    AllFrames=permute(AllFrames, [3,2,1]);    % so the dimension order is (height, width, frame)
    OneFrame=AllFrames(:,:,frame);
else 
    fprintf('\rspecify system please\r');
end
fclose(fid);
return

%use these two lines to save the frame
%OneFrame=OIReadFrame('h:\expt\050120Ruf\run6\ruf_e6b00.blk', 'v', 1, 1);
%imwrite(norm_to_uint8(OneFrame), 'c:\_AAA\red8read.bmp',  'bmp');