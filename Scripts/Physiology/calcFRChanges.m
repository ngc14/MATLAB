function [FRChanges, phaseFR, baselineFR] = calcFRChanges(PSTHS, baselineBin, phaseBins)
if(~isempty(PSTHS) & ~isempty(baselineBin) & ~isempty(phaseBins))
validInds = ~any(isnan(phaseBins),2) & ~any(isnan(baselineBin),2) ;
baselineBin = baselineBin(validInds,:);
phaseBins = phaseBins(validInds,:);
PSTHS = PSTHS(validInds,:,:);
%PSTHS organized as (TRIALS x BINS x UNITS)
baselineFR = squeeze(nanmean(PSTHS(:,baselineBin(:,1):baselineBin(:,2),:),2));
% replace 0 baseline FR with minimum report baseline FR for each unit
baselineFR(baselineFR==0) = 1;
phaseFR  = squeeze(nanmean(PSTHS(:,phaseBins(:,1):phaseBins(:,2),:),2));
FRChanges = phaseFR./baselineFR;

% minBasePerUnit = min(baselineFR,[],1);
% minBasePerUnit(isinf(minBasePerUnit)) = (1/10)/
% [zR, zC] = find(isinf(baselineFR));
% baselineFR(zR,zC) = [ones(size(zR))*minBasePerUnit(zC)]';

% =(phase mean firing rate) - (baseline mean firing rate)/
%  (phase mean firing rate + baseline mean firing rate)
% Cohen’s d, defined as follows: score=(phase mean firing rate -
% baseline mean firing rate)/SD
% Extent of Single-Neuron Activity Modulation by
% Hippocampal Interictal Discharges Predicts Declarative
% Memory Disruption in Humans
% Chrystal M. Reed, Clayton P. Mosher, and Chandravadia,
% Jeffrey M. Chung, Adam N. Mamelak, and Ueli Rutishauser1
% FRVals{p}= (taskModulatedFR-baselineFR)./(taskModulatedFR+baselineFR);

else
    [FRChanges, phaseFR, baselineFR] = deal(NaN);
end
end

