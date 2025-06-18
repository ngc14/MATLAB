function [times, ids, fileHeader] = loadSpikeData(hFile,entityID)
%% Extract channel entity info
[ns_RESULT, fileHeader] = ns_GetFileInfo(hFile);

%create fileInfo Structure
fileInfo = hFile.FileInfo(hFile.Entity(entityID).FileType);
PacketIndex = find(fileInfo.MemoryMap.Data.PacketID == ...
    hFile.Entity(entityID).ElectrodeID);
times = double(fileInfo.MemoryMap.Data.TimeStamp(PacketIndex)).*fileHeader.TimeStampResolution;
ids = fileInfo.MemoryMap.Data.Class(PacketIndex);
end

