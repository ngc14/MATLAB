function meanBG = loadAvgVoltage(hFile,dataChannelID)
%% get analog info contains things like range and sampling rate
[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, dataChannelID(1));

%% extract data and data time
TimeStamps = hFile.FileInfo(hFile.Entity(dataChannelID(1)).FileType).TimeStamps;
TimeStamps(:,1) = cumsum(TimeStamps(:,1));
numSamples = sum(TimeStamps(end,:));

analogInputData = NaN(length(dataChannelID),1);
ds = 5000;
fileSzs = [hFile.FileInfo.FileSize];
if(fileSzs(strcmp({hFile.FileInfo.Type},'ns5'))/(1024^3) > 900)
        
end
if(length(dataChannelID)>1)
    [ns_AnalogReadResult,analogInputData] = ns_GetAnalogDataBlock(hFile, dataChannelID, TimeStamps(1,1), TimeStamps(size(TimeStamps,1),size(TimeStamps,2)));
else
    [~, ~, tempData] = ns_GetAnalogData(hFile, dataChannelID,TimeStamps(1,1), TimeStamps(size(TimeStamps,1),size(TimeStamps,2)));
    analogInputData(:,end+1:end+length(tempData)) = tempData';
    clear tempData;
end

%% filter out noise spikes (abnormally large data values)
% detect noise spike and replace it with NaN
analogInputData(analogInputData > 3E38)=NaN;
% make sure data and data time matrix are an even number in length
if rem(length(analogInputData),2)==1
    analogInputData = analogInputData(1:end-1,:);
end
%% calculate avg trace value
analogInputData = downsample(analogInputData,ds);
meanBG = nanmean(highpass(analogInputData,300,ds),1);
%clear analogInputData;
end