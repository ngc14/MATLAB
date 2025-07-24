function heatmap_distribution_plots(repsIn,unitIndsIn,siteUnitModsIn,trialSegs,bins,PSTHS,...
    maxFRs,unitLocationIn,alignLim,saveDirPath,saveName,FRRange)
HISTALIGNMENT = 1;
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
unitPSTHS = cellfun(@(p) (unitInds./unitInds)'.*p, PSTHS,'UniformOutput',false);

unitNormPSTHS = cellfun(@(u) u./maxFRs, unitPSTHS,'UniformOutput', false);
avgSegTimes = mean(cell2mat(vertcat(trialSegs{:})),1);
if(contains(saveName,"Sphere","IgnoreCase",true))
    avgSegTimes = [avgSegTimes(1),avgSegTimes(5)];
elseif(contains(saveName,"Photocell","IgnoreCase",true))
    avgSegTimes = [avgSegTimes(1),avgSegTimes(4)];
else
    avgSegTimes = [avgSegTimes(1),avgSegTimes(2)];
end
avgSegTimes = findBins(avgSegTimes+alignLim{HISTALIGNMENT},bins);
sortVals = [abs(diff(horzcat(unitNormPSTHS{HISTALIGNMENT}),1,2)),zeros(size(maxFRs,1),1)];
sortValsNaN = mean(abs(sortVals(:,avgSegTimes(1):avgSegTimes(end))),2,'omitnan');
sortValsNaN(~unitInds) = NaN;
unitNormPSTHS =  {[abs(diff(horzcat(unitPSTHS{HISTALIGNMENT}),1,2)),zeros(size(maxFRs,1),1)]};

peakTimeHistograms(linspace(min(sortValsNaN),max(sortValsNaN),15),reps,sortValsNaN,...
    trialSegs{HISTALIGNMENT},siteInds,strcat(saveDirPath,"Histograms\"),saveName);
phaseHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,sortValsNaN,trialSegs,siteInds,[],alignLim);
labelUnitPSTHS(phaseHMFig);
cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "Peak\",rh,"\"),...
    strcat(saveName,"_Heatmaps"),[]),phaseHMFig,[rp(ismember(string(rp),...
    (reps())))',"All"]);

mlSort = unitLocation(:,2);
mlSort(~unitInds,:) = NaN;
rcSort = unitLocation(:,1);
rcSort(~unitInds) = NaN;
mlHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,mlSort,trialSegs,siteInds,[],alignLim);
labelUnitPSTHS(mlHMFig);

cellfun(@(ph,rh) saveFigures(ph, strcat(saveDirPath, "MedLat\",rh,"\"),...
    strcat(saveName,"_Heatmaps"),[]),mlHMFig,[rp(ismember(string(rp),...
    unique(reps())))',"All"]);
%%
rcHMFig = unitJointPSTH(bins,unitNormPSTHS,reps,rcSort,trialSegs,siteInds,[],alignLim);
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