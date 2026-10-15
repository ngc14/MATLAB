function chMap = loadChannelMap(sessionDir)
pathInds = regexp(sessionDir,'\');
hFilePath = dir([sessionDir(1:pathInds(end)),'*.nev']);
hFilePath = hFilePath(cellfun(@(s) ~contains(s,"withdraw",'IgnoreCase',true) & ...
    ~contains(s,'stim','IgnoreCase',true) & ~contains(s,'sort','IgnoreCase',true) & ...
    contains(s,["Gilligan","Skipper"]), {hFilePath.name}));
dirPath = dir(hFilePath(1).folder);
noteFiles = dirPath(find(cellfun(@(f) contains(f,extract(hFilePath.name,wildcardPattern+lookAheadBoundary(characterListPattern("_")+...
    digitsPattern+characterListPattern(".")))+"_Note"),{dirPath.name}),1)).name;
notesLines = readlines(hFilePath.folder+"\"+noteFiles);
infoLines = notesLines(contains(notesLines,"Channel")|contains(notesLines,"Electrode"));
[res,hFile] = ns_OpenFile([hFilePath.folder,'\',hFilePath.name],'single');
if(strcmp(res, 'ns_OK'))
    oldCI = [];
    chMapR = [hFile.Entity.Label];
    chMapR = chMapR(~cellfun(@isempty,chMapR));
    chs = cellfun(@(r) cell2mat(regexp(r,'(\d+)(?!.*\d)','match')),chMapR,'Uniformoutput',false);
    chs = chs(~cellfun(@isempty,chs));
    if(isempty(chMapR))
        if(any(contains({dirPath.name},'map')))
            chMapR = cellstr(readlines(dirPath(1).folder+"\"+dirPath(contains({dirPath.name},'map')).name));
            chMapR = cellfun(@(r) regexp(r,';','split'),chMapR(cellfun(@(h) ~contains(h,"#") &  length(h)>1,chMapR)),'UniformOutput',false);
            chMapR = vertcat(chMapR{:});
            chs = cellfun(@(t) t{end},cellfun(@(s) regexp(s,'[.](\w*)','tokens'),chMapR(:,1)','UniformOutput',false));
            chMapR = cellfun(@strtrim,chMapR(:,2)','UniformOutput',false);
        else
            chMapR = [hFile.Entity.ElectrodeID];
            chMapR = arrayfun(@(e) "e"+sprintf('%02d',e), chMapR(chMapR~=0));
            chs= arrayfun(@num2str,[hFile.Entity.ElectrodeID],'UniformOutput',false);
            chs = chs(cellfun(@(c) ~strcmp(c,"0"), chs));
        end
    end
    [chs,ci,~] = unique(cellfun(@str2double,chs),'sorted');
    chs = chs-(min(chs)-1);
    if(any(contains(infoLines,"Microprobes",'IgnoreCase',true)) && ...
            any(contains(infoLines, "Channel 1"+wildcardPattern+"deep")) ||...
            any(contains(infoLines, "Channel 32"+wildcardPattern+"superficial")))
        ci = flipud(ci);
    elseif(any(contains(infoLines,"VProbe","IgnoreCase",true)))
        oldCI = ci;
        maxCh = str2double(cellfun(@string,regexp(reshape(char(infoLines)',1,[]),'(\d*){1} contacts','tokens')));
        ci = [1:2:maxCh,2:2:maxCh];
        missingChs = find(~ismember(1:maxCh,chs));
        if(~isempty(missingChs))
            ci = ci(ci~=missingChs);
            ci(ci>=min(missingChs)) = ci(ci>=min(missingChs))-length(missingChs);
        end
        ci(ci>maxCh) = NaN;
    end
    chMap{1,1} = chMapR(ci);
    if(~isempty(oldCI))
        chMap{1,2} = oldCI;
    end
    chMap{2} = chs(ci);
end
ns_CloseFile(hFile);
end