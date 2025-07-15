function countAUC(reps,unitIndsIn,AUCVals,sessionToUnitIndsIn,monkeyVMap,vMask,...
    activityInds,saveDirPath,saveNameIn,phaseNames,plotMaps,countRange,AUCRange)
AUCRangeD = [0 5];
countRangeD = [1 10];
numBins = 4;
if(~exist('AUCRange', 'var'))
    AUCRange = AUCRangeD;
end
if(~exist('countRange', 'var'))
    countRange = countRangeD;
end
activityInd = activityInds;
unitInds = unitIndsIn;

sessionToUnitInds = sessionToUnitIndsIn;
saveName = saveNameIn;
% count units
unitCountFig = modulatedUnitsPerRep(reps,unitInds,saveName,"");
cellfun(@(f)saveFigures(f, strjoin([saveDirPath,"Count\"],''),saveName,[]),unitCountFig);
if(plotMaps)
    siteUnitMods = mapUnitInds2Sites(sessionToUnitInds, unitInds);
    emptyInds = ~cellfun(@any,siteUnitMods);
    countMapFig = mapUnitVals(vMask,monkeyVMap,cellfun(@nansum,siteUnitMods),...
        emptyInds,true,numBins,countRange);
    title(strjoin([saveName,"Counts"]));
    saveFigures(countMapFig, strjoin([saveDirPath,"Maps\Count\"],''),saveName,[]);
    countMapFig = mapUnitVals(vMask,monkeyVMap,cellfun(@(s) 100*(nansum(s)/length(s)),...
        siteUnitMods),emptyInds,true,numBins,[5 65]);
    title(strjoin([saveName,"Percs"]));
    saveFigures(countMapFig, strjoin([saveDirPath,"Maps\Perc\"],''),saveName,[]);
end

unitTypeInds = regexp(saveName,'_');
if(any(unitTypeInds))
    cSaveName = char(saveName);
    unitTypeName = strjoin({cSaveName(unitTypeInds+1:end),'\'},'');
    saveName = cSaveName(1:unitTypeInds-1);
else
    unitTypeName = "";
end
% unit AUC
for p = 1:length(phaseNames)
    if(~isempty(AUCRange))
        AUCColor = AUCRange{p};
    else
        AUCColor = AUCRangeD;
    end
    currAUC = AUCVals{p};
    FRFig = modulatedUnitsPerRep(reps(unitInds),currAUC(unitInds),strjoin(["FR:",saveName,...
        phaseNames(p)]),strjoin([saveName,"Units"]));
    ylim([AUCColor(1),AUCColor(end)*1.5]);
    cellfun(@(f) saveFigures(f, strjoin([saveDirPath,"FR\",string(unitTypeName)],''),...
        saveName+phaseNames(p),[]),FRFig);
    if(plotMaps)
        siteAUC = mapUnitInds2Sites(sessionToUnitInds,currAUC.*(unitInds./unitInds));
        emptyInds = ~cellfun(@any,siteAUC);
        FRMapFig = mapUnitVals(vMask,monkeyVMap,cellfun(@nanmean,siteAUC),emptyInds,false,numBins,AUCColor);
        title(strjoin([saveName,phaseNames(p)]));
        saveFigures(FRMapFig, strjoin([saveDirPath,"Maps\FR\",string(unitTypeName)],''),...
            saveName+phaseNames(p),[]);
    end
end
end

function sessionInds = mapUnitInds2Sites(session2Units, currInds)
sessionInds = arrayfun(@(a) currInds(session2Units==a),...
    min(session2Units):max(session2Units),'UniformOutput', false);
end
