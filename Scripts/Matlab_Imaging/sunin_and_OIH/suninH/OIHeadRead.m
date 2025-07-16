function anapar=OIHeadRead(oifile, system)
% Read header information from optical imaging data file
% anapar=OIHeadRead(oifile, system)

% (most usefule ones are:
%	Frame width
%	Frame heigh
%	Frame number per stimulus
%	Stimulus number

verbose=0;
anapar.system=system;   % just transfer
fid = fopen(oifile,'r');
if fid==-1
	fprintf('Data file not exist\n');
end
if (system=='d')	% For Dan Ts'o old system (no header info, need modify following lines)
    anapar.FramesPerStim=10;
    anapar.FrameWidth=192;
    anapar.FrameHeight=144;
    anapar.NStim=33
    anapar.DataType=12; % Unsigned short
	anapar.HeaderLength=0;		
	fileinfo=dir(oifile);
	if fileinfo.bytes~=anapar.FramesPerStim*anapar.FrameWidth*anapar.FrameHeight*anapar.NStim*2
		fprintf('size does not match, check setting in "OIHeadRead.m"\n');
	end
elseif (system=='r')	% for redshirt data
    fseek(fid, (5-1)*2, 'bof');
    anapar.FramesPerStim=fread(fid, 1, 'uint16=>double');
    fseek(fid, (385-1)*2, 'bof');
    anapar.FrameWidth=fread(fid, 1, 'uint16=>double');
    fseek(fid, (386-1)*2, 'bof');
    anapar.FrameHeight=fread(fid, 1, 'uint16=>double');
    fseek(fid, (389-1)*2, 'bof');
    FrameInterval=fread(fid, 1, 'uint16=>double')/1000.0;
    fseek(fid, (392-1)*2, 'bof');
    AcquisitionRatio=fread(fid, 1, 'uint16=>double');
    fseek(fid, (2100)*2, 'bof');
    anapar.NStim=fread(fid, 1, 'uint16=>double'); % in old data (<050621) these bytes were not used, so NStim=0, need manual set NStim below
    %anapar.NStim =7;		% Modify this line if stim number is not available from data header
	if anapar.NStim<1	
	    fprintf('Note, old data (<050621) need manual enter stim number) \n');
	    fprintf('Please modify "OIReadHeader.m" at line 27 before processing.\n');
	    return;
	end
    fseek(fid, (2101-1)*2, 'bof');	% stim id
    StimID=fread(fid, anapar.NStim, 'uint16=>double'); % in old data (<050621) these bytes were not used.
    NBNC=8;
    anapar.DataType=12; % Unsigned short

elseif (system=='v') % for VDAQ data
    lFileSize=fread(fid,1,'int32=>double');         %*
    lCheckSum_Header=fread(fid,1,'int32=>double');  %*
    lCheckSum_Data=fread(fid,1,'int32=>double');    %*
    lLenHeader=fread(fid,1,'int32=>double');        %*
    lVersionID=fread(fid,1,'int32=>double');        %*
    lFileType=fread(fid,1,'int32=>double');         %*
    lFileSubtype=fread(fid,1,'int32=>double');      %*
    lDataType=fread(fid,1,'int32=>double');         %*
    lSizeOf=fread(fid,1,'int32=>double');           %??
    lFrameWidth=fread(fid,1,'int32=>double');       %**
    lFrameHeight=fread(fid,1,'int32=>double');      %**
    lNFramesPerStim=fread(fid,1,'int32=>double');   %**
    lNStimuli=fread(fid,1,'int32=>double');         %**
    lInitialXBinFactor=fread(fid,1,'int32=>double');
    lInitialYBinFactor=fread(fid,1,'int32=>double');
    lXBinFactor=fread(fid,1,'int32=>double');
    lYBinFactor=fread(fid,1,'int32=>double');
    acUserName=fread(fid,32,'uchar=>char');
    acRecordingDate=fread(fid,16,'uchar=>char');
    lX1ROI=fread(fid,1,'int32=>double');
    lY1ROI=fread(fid,1,'int32=>double');
    lX2ROI=fread(fid,1,'int32=>double');
    lY2ROI=fread(fid,1,'int32=>double');

    lStimOffs=fread(fid,1,'int32=>double');
    lStimSize=fread(fid,1,'int32=>double');
    lFrameOffs=fread(fid,1,'int32=>double');
    lFrameSize=fread(fid,1,'int32=>double');
    lRefOffs=fread(fid,1,'int32=>double');
    lRefSize=fread(fid,1,'int32=>double');          %**
    lRefWidth=fread(fid,1,'int32=>double');
    lRefHeight=fread(fid,1,'int32=>double');
    aushWhichBlocks=fread(fid,32,'uchar=>char');
    aushWhichFrames=fread(fid,32,'uchar=>char');

    fLoClip=fread(fid,1,'float32=>double');       %*  
    fHiClip=fread(fid,1,'float32=>double');         %*
    lLoPass=fread(fid,1,'int32=>double');
    lHiPass=fread(fid,1,'int32=>double');
    acOperationsPerformed=fread(fid,64,'uchar=>char');

    fMagnification=fread(fid,1,'float32=>double');
    ushGain=fread(fid,1,'uint16=>double');
    ushWavelength=fread(fid,1,'uint16=>double');
    lExposureTime=fread(fid,1,'int32=>double');
    lNRepetitions=fread(fid,1,'int32=>double');
    lAcquisitionDelay=fread(fid,1,'int32=>double');
    lInterStimInterval=fread(fid,1,'int32=>double');
    acCreationDate=fread(fid,16,'uchar=>char');
    acDataFilename=fread(fid,64,'uchar=>char');         %*
    acOraReserved=fread(fid,256,'uchar=>char');

    lIncludesRefFrame=fread(fid,1,'int32=>double');
    acListOfStimuli=fread(fid,256,'uchar=>char');
    lNVideoFramesPerDataFrame=fread(fid,1,'int32=>double');
    lNTrials=fread(fid,1,'int32=>double');
    lScaleFactor=fread(fid,1,'int32=>double');
    fMeanAmpGain=fread(fid,1,'float32=>double');        %*
    fMeanAmpDC=fread(fid,1,'float32=>double');          %*
    ucBegBaselineFrameNo=fread(fid,1,'uchar=>double');
    ucEndBaselineFrameNo=fread(fid,1,'uchar=>double');
    ucBegActivityFrameNo=fread(fid,1,'uchar=>double');
    ucEndActivityFrameNo=fread(fid,1,'uchar=>double');
    acVdaqReserved=fread(fid,252,'uchar=>char');
    acUser=fread(fid,256,'uchar=>char');
    acComment=fread(fid,256,'uchar=>char');

    anapar.FrameWidth=lFrameWidth;
    anapar.FrameHeight=lFrameHeight;
    anapar.NStim=lNStimuli;
    anapar.FramesPerStim=lNFramesPerStim;
    anapar.DataType=lDataType;      % 11: unsigned Char, 12: unsigned short (16bit?), 13: unsigned  long (32bit), 14: float(32bit)
else
    fprintf('has to be v or r');
end
if verbose
	fprintf('\nFile Info: Width (%d), Height(%d), Stim (%d), Frame (%d), DataType (%d)', ...
    	anapar.FrameWidth, anapar.FrameHeight, anapar.NStim, anapar.FramesPerStim, anapar.DataType);
end
status = fclose(fid);
return;
