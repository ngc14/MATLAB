clear all;
close all;
monkey = 'Gilligan';
directory = dir(['S:\Lab\',monkey,'\All Data']);
dateRange = {'07_01_2018', '01_01_2021'};
modelNames = {'Deltoid','Biceps', 'Triceps', 'Extensor', 'ProximalRadialFlexor',...
    'D4D5Extensor','DistalFlexor'};
modelNames = {'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor',...
    'Digit Extensor','Digit Flexor'};
segNames = {'RT', 'R', 'L', 'H', 'W'};
alignSegs = {'StartReach'};
saveVals = 1;
showFigs = 1;
xLims = {[-1 3]};
gap = .25;
Fs = 2000;
winSize = .1; %in seconds
winSize = winSize*Fs;
xLims = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);
exclude = {};

dates = [];
directory = directory(cellfun(@(c) contains(c, monkey), {directory.name}) & [directory.isdir]);
for d = 1:size(directory,1)
    underscores = regexp(directory(d).name, '_');
    fileDate = directory(d).name(underscores(1)+1:end);
    dates(end+1) = datenum(fileDate);
end
filesToProcess = find(dates>=datenum(dateRange{1}) & dates<=datenum(dateRange{2}));
for e = 1:size(exclude,2)
    excludeInd = find(dates==datenum(exclude{e}));
    if(~isempty(excludeInd))
        filesToProcess = filesToProcess(filesToProcess~=excludeInd);
    end
end
directory = directory(filesToProcess);

for mm = 1:length(modelNames)
    eventData = [];
    segmentTimes = [];
    arduinoTrials = [];
    Fs = [];
    conds = [];
    numTrials= [1];
    fName = [];
    condSegs = [];
    numSessions = 0;
    muscle = modelNames{mm};

    for f = 1:size(directory,1)
        if(strcmp(monkey,'Gilligan'))
            if(strcmp(muscle,'Digit Extensor'))
                resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                    directory(f).name, '_D4D5Extensor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                        directory(f).name, '_Digit Extensor.mat'];
                end
            elseif(strcmp(muscle,'Digit Flexor'))
                resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                    directory(f).name, '_DistalFlexor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                        directory(f).name, '_Digit Flexor.mat'];
                end
            elseif(strcmp(muscle,'Wrist Extensor'))
                resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                    directory(f).name, '_Extensor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                        directory(f).name, '_Wrist Extensor.mat'];
                end
            elseif(strcmp(muscle,'Wrist Flexor'))
                resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                    directory(f).name, '_ProximalRadialFlexor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                        directory(f).name, '_Wrist Flexor.mat'];
                end
            else
                resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                    directory(f).name, '_', muscle, '.mat'];
            end
        else
            resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
                directory(f).name, '_', muscle, '.mat'];
        end
        if(exist(resultFolder, 'file'))
            l = load(resultFolder);
            fName = [fName; {l.Date}];
            if(length(l.Conditions)==size(conds,1))
                conds = [conds, l.Conditions];
                condSegs = [condSegs, l.ConditionSegments];
            end
            Fs = [Fs, l.Fs];
            eventData = [eventData, l.EMGData;];
            arduinoTrials = [arduinoTrials; l.ArduinoData];
            segmentTimes = [segmentTimes, l.SegTimes];
            numTrials(end+1,:) = length(l.ArduinoData);
            numSessions = numSessions + 1;
        end
    end


    [conds,uniqueInds,~] = unique(conds, 'stable');
    condSegs = condSegs(uniqueInds);
    Fs = unique(Fs);

    badTrial = zeros(size(eventData));
    for t=1:size(eventData,2)
        trialSig = eventData{t};
        rectified = movmean(abs(trialSig),20);
        baseline = nanmean(rectified(1:segmentTimes{t}(1)));
        Y = fft(trialSig);
        P2 = abs(Y/length(trialSig));
        P1 = P2(1:length(trialSig)/2 + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        freq = Fs*(0:(length(trialSig)/2))/length(trialSig);
        if(any(find(freq<15 & P1>7)))
            badTrial(t) = 1;
        elseif(range(abs(trialSig)) < 5)
            badTrial(t) = 2;
        elseif(baseline>80)
            badTrial(t) = 3;
        elseif(P1>10 & mod(freq,notch)~=0 & mod(freq,notch)<notch+.5 & mod(freq,notch)>notch-.5)
            badTrial(t) = 4;
        end
    end
    %%
    maxTrialLength = max(cellfun(@length, eventData));
    paddedEventData = cellfun(@(a) padArrays(abs(a),maxTrialLength), eventData, 'UniformOutput', false);
    paddedEventData = cellfun(@(a) conv(a,gausswin(winSize)/sum(gausswin(winSize)),'same'),  paddedEventData, 'UniformOutput', false);
    %paddedEventData = cellfun(@(a) filter(1/100*ones(1,100),1,a,[],2), paddedEventData, 'UniformOutput', false);
    %paddedEventData = cellfun(@(a) smoothdata(a,2,'movmean',[100,0]), paddedEventData, 'UniformOutput', false);
    for d = 2:size(numTrials,1)
        sessionTrials = zeros(size(eventData));
        sessionTrials(sum(numTrials(1:d-1)):sum(numTrials(1:d))-1) = 1;
        sessionVals = paddedEventData(sessionTrials & ~badTrial);
        sessionVals = cell2mat(sessionVals');
        sessionVals(sessionVals==0) =NaN;
        maxVal = max(nanmean(sessionVals));
        minVal = min(nanmean(sessionVals));

        if(isempty(maxVal))
            maxVal = NaN;
        end
        if(isempty(minVal))
            minVal = NaN;
        end
        [allMaxVals{sum(numTrials(1:d-1)):sum(numTrials(1:d))-1}] = deal(maxVal);
        [allMinVals{sum(numTrials(1:d-1)):sum(numTrials(1:d))-1}] = deal(minVal);
    end
    paddedEventData = reshape(cell2mat(paddedEventData),length(paddedEventData{1}),length(paddedEventData))';

    %%
    for c = 1:length(conds)
        f=1;

        condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
        badCondTrials = badTrial(condInd);
        currTrials = paddedEventData(condInd,:);
        currSegs = segmentTimes(condInd);
        currMinVals= allMinVals(condInd);
        currMaxVals= allMaxVals(condInd);

        meanData = nanmean(paddedEventData);
        stdData = nanstd(paddedEventData);

        hold on;
        numSegs = mode(cellfun(@length,currSegs));

        badSegTrials = cellfun(@(a) length(a)~=numSegs, currSegs);
        currTrials = currTrials(~badSegTrials & ~badCondTrials,:);
        currSegs = currSegs(~badSegTrials & ~badCondTrials);
        currMinVals= currMinVals(~badSegTrials & ~badCondTrials);
        currMaxVals= currMaxVals(~badSegTrials & ~badCondTrials);
        % for session differentiation
        trialChange = find([0 diff(cell2mat(currMinVals))]);
        plotColors = cell(1,length(currSegs));
        for t = 1:length(trialChange)
            if(t==1)
                [plotColors{1:trialChange(t)}] = deal(t);
                if(t==length(trialChange))
                    [plotColors{trialChange(t):end}] = deal(t);
                end
            elseif(t==length(trialChange))
                [plotColors{trialChange(t-1):end}] = deal(t);
            else
                [plotColors{trialChange(t-1):trialChange(t)}] = deal(t);
            end
        end
        if(sum(~badSegTrials)~=length(badSegTrials))
            disp(fName{f});
            disp([num2str(sum(badSegTrials==1)), ' trials ignored']);
        end

        %normTrials = (currTrials-cell2mat(currMinVals)') ./ (cell2mat(currMaxVals)-cell2mat(currMinVals))';
        normTrials = currTrials;

        currTrials = cell(1,size(normTrials,1));
        for z = 1:size(normTrials,1)
            currTrials{z} = normTrials(z,:);
        end

        averageSegs = mean(reshape(cell2mat(currSegs),numSegs,size(currSegs,2))');
        allSegs = reshape(cell2mat(currSegs),numSegs,size(currSegs,2))';
        plotStart = 0;
        xTicks = [];
        labels = [];
        alignedEMGSigs = [];
        for align = 1:length(alignSegs)
            if(align==1)
                plotStart = 0;
            else
                plotStart = plotStart + (gap*Fs) + xLims{align-1}(2) + abs(xLims{align}(1));
            end
            alignSegInd = find(strcmp(condSegs{c},alignSegs{align}))-1;
            if(isempty(alignSegInd) & strcmp(conds{c}, 'Rest'))
                alignSegInd = 1;
            end
            alignedEMGSigs = [alignedEMGSigs, cellfun(@(a,b) a(b(alignSegInd)+xLims{align}(1):...
                b(alignSegInd)+xLims{align}(2)),currTrials, currSegs,'UniformOutput', false)'];
            if(showFigs)
                subplot(2,2,c);
                hold on;
                colors = colormap('colorcube');
                cellfun(@(a,b) plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2),...
                    a, 'Color',colors(b,:)), alignedEMGSigs(:,align), plotColors');
                plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2), ...
                    nanmean(reshape([alignedEMGSigs{:,align}],size(alignedEMGSigs{1,align},2),size(currTrials,2))'),...
                    'k', 'LineWidth', 2);
                beforeSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:-1:1)));
                afterSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:end)));
                avgSegs = [beforeSegs plotStart afterSegs];
                for s = 1:length(avgSegs)
                    if(avgSegs(s)>=plotStart+xLims{align}(1) && avgSegs(s)<=plotStart+xLims{align}(2))
                        plot([avgSegs(s) avgSegs(s)],[0 max([alignedEMGSigs{:,align}])],'r--')
                    end
                end
                xTicks = [xTicks, plotStart+[xLims{align}(1) 0 xLims{align}(2)]];
                labels = [labels, [xLims{align}(1) 0 xLims{align}(2)]];
                plotLims = [xTicks(1),xTicks(end)];
            end
        end
        saveEMG = alignedEMGSigs;
        rawTrials = normTrials;
        if(saveVals)
            save(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\',muscle,'_',conds{c},'_Trials'], 'trialChange');
            save(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\',muscle,'_',conds{c}], 'rawTrials');
            save(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\',muscle,'_',conds{c},'_Times'], 'allSegs');
        end
        if(showFigs)
            set(gca, 'XTick', xTicks);
            set(gca, 'XTickLabel', arrayfun(@num2str, labels./Fs, 'UniformOutput', false));
            set(gca, 'XLim', plotLims);
            set(gca, 'YLim', [0 150]);
            title(sprintf('%s\n\t%s%d', conds{c}, 'n=', length(alignedEMGSigs)));
        end
    end
    saveFigures(gcf,['S:\Lab\', monkey,'\Sorted EMG Data\Averages_Raw\'], muscle,[]);

    clearvars -except modelNames directory xLims alignSegs gap segNames Fs showFigs saveVals winSize monkey
    close all;
end
function saveFigures(fHandle,saveDir,saveName,imgTrace)
if(~exist(saveDir, 'dir'))
    mkdir(saveDir)
end
if(strcmp(get(fHandle,'type'),'figure'))
    set(fHandle,'Renderer','painters');
end
saveas(fHandle,[saveDir,saveName,'.epsc']);
saveas(fHandle,[saveDir,saveName,'.png']);
saveas(fHandle,[saveDir,saveName,'.fig']);
if(~isempty(imgTrace))
    boundaries = bwboundaries(imgTrace);
    fHandle2 = copyobj(gca(fHandle),figure()); hold on;
    cellfun(@(b) plot(b(:,2), b(:,1), 'g', 'LineWidth',lineWidth), boundaries);
    saveFigures(fHandle2,[saveDir,'Imaging\'],saveName,[]);
end
end
function padded = padArrays(arr, length)
padded = arr;
padded(end+1:length)  = 0;
end