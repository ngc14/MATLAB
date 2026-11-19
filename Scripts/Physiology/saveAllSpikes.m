function saveAllSpikes(sortedDir)
[~,mostRecentSort] = max(cellfun(@(n,d) d.*(contains(n,wildcardPattern+("S"|"s")+"ort")./...
    contains(n,wildcardPattern+("S"|"s")+"ort")),{sortedDir.name},{sortedDir.datenum}));
sortedFile = [sortedDir(1).folder,'\',sortedDir(mostRecentSort).name];
resDir = dir([sortedDir(1).folder,'\Results\*.mat']);
[~,sortRes]= natsort({resDir.name});
resDir = resDir(sortRes);
if(strcmp(sortedFile(end-3:end),'.plx'))
    hFile = plx_info(sortedFile,1);
elseif(strcmp(sortedFile(end-3:end),'.nev'))
    [~,hFile] = ns_OpenFile(sortedFile);
end
for r = 1:length(resDir)
    mFile = matfile([resDir(r).folder,'\',resDir(r).name],'Writable',true);
    if(strcmp(sortedFile(end-3:end),'.plx'))
        dataTime = nonzeros(hFile(2:end,r+1));
        ids = find(hFile(2:end,r+1));
    elseif(strcmp(sortedFile(end-3:end),'.nev'))
        [dataTime, ids] = loadSpikeData(hFile,r);
    end
    ss=mFile.sortedSpikeData;
    if(~isempty(dataTime))
        sortedIDs = unique(ids);
        sortedIDs = sortedIDs(sortedIDs>0 & sortedIDs<255);
        allUnits = cell(1,length(sortedIDs));
        if(size(ss.SpikeTimes,1)~=length(sortedIDs))
            disp(['Ch ', num2str(r), ': Mismatch unit labels']);
        end
        for u = 1:length(sortedIDs)
            if(strcmp(sortedFile(end-3:end),'.nev'))
                allUnits{u} = dataTime(ids==sortedIDs(u));
            elseif(strcmp(sortedFile(end-3:end),'.plx'))
                [~,allUnits{u}] = plx_ts(sortedFile,r, sortedIDs(u));
            end
        end
    else
        allUnits{u} = {};
    end
    mFile.allUnits = allUnits;
end
ns_CloseFile(hFile);
clear mFile;
end