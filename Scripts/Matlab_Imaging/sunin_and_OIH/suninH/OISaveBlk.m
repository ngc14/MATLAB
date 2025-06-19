%   Save Optical Imaging block file     HDL 050521
%   OISaveBlk(data, header, filename, system)

%   data: matrix(width, height, FramePerStim, NStim)
%   header: 1716 for vdaq, uint8 type
%   filename: include path
%   system: 'v' for vdaq

% Following 3 lines for testing only
% fidlastblock=fopen('Bha_E1B000.BLK', 'r');   % simply copy an existing file header since they are the same, except data type
% header=fread(fidlastblock, 1716, 'uint8');
% fclose(fidlastblock);

header(29)=uint8(14);   % force to float type
fidnewblock=fopen(filename, 'wb');
fwrite(fidnewblock, header, 'uint8');
data=permute(data, [2,1,3,4]);
fwrite(fidnewblock, avgblock, 'float32');
fclose(fidnewblock); 
    
return    