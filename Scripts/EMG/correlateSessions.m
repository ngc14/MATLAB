directoryModel = dir('\\130.49.229.252\gharbawie\Lab\Gilligan\All Data');
dateRange = {'04_01_2018', '08_10_2018'};
muscleModel = 'ProximalRadialFlexor';
exclude = {'08_02_2018', '08_09_2018','08_10_2018'};
dates = [];
directoryModel = directoryModel(3:end);
for d = 1:size(directoryModel,1)
    underscores = regexp(directoryModel(d).name, '_');
    fileDate = directoryModel(d).name(underscores(1)+1:end);
    dates(end+1) = datenum(fileDate);
end
filesToProcess = find(dates>=datenum(dateRange{1}) & dates<=datenum(dateRange{2}));
for e = 1:size(exclude,2)
    excludeInd = find(dates==datenum(exclude{e}));
    if(~isempty(excludeInd))
        filesToProcess = filesToProcess(filesToProcess~=excludeInd);
    end
end
directoryModel = directoryModel(filesToProcess);

sessionDate = '08_01_2018';
muscleCandidate = 'DistalDigitFlexor';
fileCandidate = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
    sessionDate, '\EMG\Results_New\Gilligan_', sessionDate, '_', muscleCandidate, '.mat'];

showPlot = 1;

eventData = [];
segmentTimes = [];
arduinoTrials = [];
Fs = [];
conds = [];
condSegs = [];
numTrials = 1;
for f = 1:size(directoryModel,1)
    resultFolder = [directoryModel(f).folder, '\',directoryModel(f).name,'\EMG\Results_New\',...
        directoryModel(f).name, '_', muscleModel, '.mat'];
    if(exist(resultFolder, 'file'))
        l = load(resultFolder);
        conds = [conds, l.sortedEMGData.Conditions];
        condSegs = [condSegs, l.sortedEMGData.ConditionSegments];
        Fs = [Fs, l.sortedEMGData.SampleRate];
        eventData = [eventData, l.sortedEMGData.EMGData;];
        arduinoTrials = [arduinoTrials; l.sortedEMGData.ArduinoData];
        segmentTimes = [segmentTimes, l.sortedEMGData.SegTimes];
        numTrials(end+1,:) = length(l.sortedEMGData.ArduinoData);
    end
end
conds = conds(1:4);
condSegs = condSegs(1:4);
Fs = unique(Fs);

sortedEMGData.EMGData= eventData;
sortedEMGData.SegTimes = segmentTimes;
sortedEMGData.SampleRate=Fs;
sortedEMGData.ArduinoData = arduinoTrials;
sortedEMGData.Conditions = conds;
sortedEMGData.ConditionSegments = condSegs;
sortedEMGData.Muscle = muscleModel;
%% get derivatives
disp('Obtaining model derivatives');
derivativesModel = EMG_plot_derivative(sortedEMGData,showPlot);
%%
% load data from directory
disp('Obtaining candidate derivatives');
l = load(fileCandidate);
derivativesCandidate = EMG_plot_derivative(l.sortedEMGData,showPlot);
%% correlate
correlationCondSeg = [];
for c=1:length(conds)
    tempCond = conds{c};
    tempCond = tempCond(~isspace(tempCond));
    condDerivativesM = derivativesModel.(tempCond);
    condDerivativesC = derivativesCandidate.(tempCond);
    for s = 1:size(condDerivativesM,2)
        averageDeriv = nanmean(cell2mat(condDerivativesM(:,s)));
        % photocell
        %         correlationTrial = [];
        if(~isempty(averageDeriv))
            %             for t = 1:size(condDerivativesC,1)
            %                 correlationTrial(t) = corr(averageDeriv',condDerivativesC{t,s}');
            %             end
            %             correlationCondSeg(c,s) = median(correlationTrial);
            correlationCondSeg(c,s) = corr(averageDeriv', nanmean(cell2mat(condDerivativesC(:,s)))');
        end
    end
end
disp(median(median(correlationCondSeg)));