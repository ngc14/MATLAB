function monkeySiteMask = getVoronoiMask(datesTable,mm,refMask,plotReps)
allReps = cellfun(@unique,datesTable.SiteRep,'UniformOutput',false);
allReps = unique([allReps{:}]);
excludeReps = allReps(~ismember(allReps,plotReps));
repThreshInd = cellfun(@(r)  ~(ismember(r,excludeReps)), datesTable.SiteRep,'UniformOutput',false);
datesTable.SiteRep = cellfun(@(s,r) s(r), datesTable.SiteRep, repThreshInd, 'UniformOutput',false);
datesTable.Thresh = cellfun(@(t,r) t(r), datesTable.Thresh, repThreshInd, 'UniformOutput',false);
datesTable = datesTable(~cellfun(@isempty,datesTable.SiteRep),:);
siteLocation = [datesTable.x,datesTable.y];
[verticies, vCells] = voronoin(fliplr([siteLocation; [0 size(refMask,2); ...
    size(refMask,1) 0; 0 0; size(refMask,1) size(refMask,2)]]));
siteMasks = [];
for i = 1:length(siteLocation)
    currSite = siteLocation(i,:);
    tempCircle = zeros(size(refMask,[1 2]) + 2*mm.tileBuffer);
    tempCircle((currSite(2)-mm.siteRadius+mm.tileBuffer):...
        (currSite(2)+mm.siteRadius+mm.tileBuffer),...
        (currSite(1)-mm.siteRadius+mm.tileBuffer):...
        (currSite(1)+mm.siteRadius+mm.tileBuffer)) = mm.poolCircle;
    tempCircle = tempCircle(mm.tileBuffer:end-(mm.tileBuffer+1),...
        mm.tileBuffer:end-(mm.tileBuffer+1));
    siteMasks(:,:,i) = tempCircle & poly2mask(verticies(vCells{i},2),...
        verticies(vCells{i},1),size(tempCircle,1),size(tempCircle,2));
end
monkeySiteMask = any(siteMasks(:,:,:),3) & refMask;
end