function medianCorr = correlateSession(muscleModel, candidateFile)
%% get model derivatives
derivativesModel = load(['\\130.49.229.252\gharbawie\Lab\Gilligan\Sorted EMG Data\Results\Models\'...
    muscleModel]);
derivativesModel = derivativesModel.derivativesModel;
%% load data from directory
l = load(candidateFile);
derivativesCandidate = EMG_plot_derivative(l.sortedEMGData,0);
%% correlate
corrVals = [];
conds = fieldnames(derivativesModel);
for c=1:length(conds)
    tempCond = conds{c};
    tempCond = tempCond(~isspace(tempCond));
    condDerivativesM = derivativesModel.(tempCond);
    condDerivativesC = derivativesCandidate.(tempCond);
    for s = 1:size(condDerivativesM,2)
        averageDeriv = nanmean(cell2mat(condDerivativesM(:,s)));
        % photocell
        if(~isempty(averageDeriv) & iscell(condDerivativesC))
            corrVals(c,s) = corr(averageDeriv', nanmean(cell2mat(condDerivativesC(:,s)),1)');
        end
    end
end
medianCorr = median(median(corrVals));
end