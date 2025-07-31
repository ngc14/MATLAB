function heatmap_distribution_plots(repsIn,unitIndsIn,siteUnitModsIn,trialSegs,bins,PSTHS,...
    maxFRs,unitLocationIn,alignLim,saveDirPath,saveName,FRRange)
HISTALIGNMENT = 1;
HOLDIND = 5;
condInd = 0;
if(contains(saveName,"photocell",'IgnoreCase',1))
    condInd = 1;
elseif(contains(saveName,"rest",'IgnoreCase',1))
    condInd = 3;
end
rp = fieldnames(MotorMapping.repColors);
if(exist('FRRange', 'var'))
   FRLim = FRRange;
else
    FRLim = [5 50];
end
%% unit PSTHS
reps = repsIn;
unitInds = unitIndsIn;
unitLocation = unitLocationIn;
sessionInds = ~isnan(siteUnitModsIn) & unitInds;
siteUnitMods = siteUnitModsIn;
siteUnitMods(~sessionInds) = 0;
[siteIndsN,siteInds] = unique(siteUnitMods);
siteInds = siteInds(siteIndsN>0);
avgSegTimes = findBins(vertcat(trialSegs{:}),bins);

unitPSTHS = cellfun(@(p) (unitInds./unitInds)'.*p, PSTHS,'UniformOutput',false);
sortVals = [abs(diff(horzcat(unitPSTHS{HISTALIGNMENT}),1,2)),zeros(size(maxFRs,1),1)];
sortValsNaN = cellfun(@(s,a) mean(s(a(1):a(end))),num2cell(sortVals,2),...
    num2cell([avgSegTimes(:,1),avgSegTimes(:,HOLDIND-condInd)]+...
    int64(alignLim{HISTALIGNMENT}./mode(diff(bins))),2));
sortValsNaN(~unitInds,:) = NaN;

unitNormPSTHS =  cellfun(@(u) u./maxFRs, unitPSTHS,'UniformOutput', false);
peakTimeHistograms(linspace(min(sortValsNaN),max(sortValsNaN),15),reps,sortValsNaN,...
    {trialSegs{HISTALIGNMENT}},siteInds,strcat(saveDirPath,"Histograms\"),saveName);
phaseHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,sortValsNaN,{trialSegs},siteInds,[],alignLim);
labelUnitPSTHS(phaseHMFig);
cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "Peak\",rh,"\"),...
    strcat(saveName,"_Heatmaps"),[]),phaseHMFig,[rp(ismember(string(rp),...
    (reps())))',"All"]);

mlSort = unitLocation(:,2);
mlSort(~unitInds,:) = NaN;
rcSort = unitLocation(:,1);
rcSort(~unitInds) = NaN;
mlHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,mlSort,{trialSegs},siteInds,[],alignLim);
labelUnitPSTHS(mlHMFig);

cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "MedLat\",rh,"\"),...
    strcat(saveName,"_Heatmaps"),[]),mlHMFig,[rp(ismember(string(rp),...
    unique(reps())))',"All"]);
%%
rcHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,rcSort,{trialSegs},siteInds,[],alignLim);
labelUnitPSTHS(rcHMFig);
cellfun(@(f) camroll(gca(f),90),rcHMFig,'UniformOutput',false);
cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "RostCaud\",rh,"\"),...
    strcat(saveName,"_Heatmaps"),[]),rcHMFig,[rp(ismember(string(rp),...
    (reps())))',"All"]);
pause(1);
end
function labelUnitPSTHS (cellFigs)
ax1 = cellfun(@gca, cellFigs,'UniformOutput',false);
ax2 = cellfun(@(x1) axes('Position', x1.Position), ax1,'UniformOutput',false);
cellfun(@(a1,a2) linkprop([a1,a2], {'XLim','YLim','Position','view'}), ax1,ax2,'UniformOutput',false);
cellfun(@(y1,y2) ylim(y2,[y1.YLim]), ax1,ax2,'UniformOutput',false);
cellfun(@(y1,y2)  yticks(y2,y1.YTick), ax1,ax2,'UniformOutput',false);
cellfun(@(a2) set(a2,'Visible', 'off'), ax2,'UniformOutput', false);
cellfun(@(a1,a2) text(a1,repmat(max(xlim(a1))*1.01,size(a2.YTick)), ...
    a2.YTick,arrayfun(@(s) sprintf('%d ',fix(s)), a2.YTick,'UniformOutput',false),...
    'HorizontalAlignment','left','VerticalAlignment','middle'),...
    ax1,ax2,'UniformOutput', false);
end