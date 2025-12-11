function unitInds = mapSites2Units(unit2site,currInds)
if(iscell(currInds))
    unit2site = num2cell(unit2site);
    unitInds = cellfun(@(sr,si) repmat({sr},si,1),currInds,...
        unit2site,'UniformOutput', false);
else
    unitInds = arrayfun(@(sr,si) repmat(sr,si,1),currInds,...
        unit2site,'UniformOutput', false);
end
unitInds = vertcat(unitInds{:});
end