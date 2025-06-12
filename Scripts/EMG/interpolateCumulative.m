directory = dir('S:\Lab\Gilligan\All Data');
dateRange = {'04_01_2018', '08_10_2019'};
muscle = 'Deltoid';
segNames = {'Rx', 'R', 'L', 'H', 'W'};

exclude = {};
dates = [];

eventData = [];
segmentTimes = [];
arduinoTrials = [];
Fs = [];
conds = [];
numTrials= [1];
fName = [];
directory = directory(3:end);
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

for f = 1:size(directory,1)
    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\',...
        directory(f).name, '_', muscle, '.mat'];
    if(exist(resultFolder, 'file'))
        l = load(resultFolder);
        if(isfield(l,"sortedEMGData"))

        fName = [fName; {l.sortedEMGData.Date}];
        %conds = [conds, l.sortedEMGData.Conditions];
        Fs = [Fs, l.sortedEMGData.SampleRate];
        eventData = [eventData, l.sortedEMGData.EMGData;];
        arduinoTrials = [arduinoTrials; l.sortedEMGData.ArduinoData];
        segmentTimes = [segmentTimes, l.sortedEMGData.SegTimes];
        numTrials(end+1,:) = length(l.sortedEMGData.ArduinoData);
        end
    end
end
conds = conds(1:4);
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
for d = 2:size(numTrials,1)
    sessionTrials = zeros(size(eventData));
    sessionTrials(sum(numTrials(1:d-1)):sum(numTrials(1:d))-1) = 1;
    maxVal = nanmedian(cellfun(@(a) max(filter(1/100*ones(1,100),1,abs(a))), eventData(sessionTrials & ~badTrial)));
    minVal = nanmedian(cellfun(@(a) min(filter(1/100*ones(1,100),1,abs(a))), eventData(sessionTrials & ~badTrial)));
    if(isempty(maxVal))
        maxVal = NaN;
    end
    if(isempty(minVal))
        minVal = NaN;
    end
    [allMaxVals{sum(numTrials(1:d-1)):sum(numTrials(1:d))-1}] = deal(maxVal);
    [allMinVals{sum(numTrials(1:d-1)):sum(numTrials(1:d))-1}] = deal(minVal);
end
%%
figure('units','normalized','outerposition',[0 0 1 1]);
axes('Position',[0 0 1 1],'Visible','off');
text(.5,.98, muscle, 'FontSize', 28, 'HorizontalAlignment', 'Center');
text(.5,.02,'Time (ms)', 'FontSize',26, 'HorizontalAlignment', 'Center');
tx = text(.02,0.75,'Voltage (uV)', 'FontSize', 26, 'HorizontalAlignment', 'Center');
set(tx, 'Rotation', 90);
tx = text(.02,0.25,'Voltage (uV)', 'FontSize', 26, 'HorizontalAlignment', 'Center');
set(tx, 'Rotation', 90);
hx = tight_subplot(2,2,[.1 .01],[.1 .1],[.07 .01]);
colors = colormap('colorcube');
printName = {'Small Sphere', 'Large Sphere', 'Extra Small Sphere', 'Photocell'};
for c = 1:length(conds)
    axes(hx(c));
    set(gca, 'FontSize', 24);
    
    f=1;
    condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
    badCondTrials = badTrial(condInd);
    currTrials = eventData(condInd);
    currSegs = segmentTimes(condInd);
    currMinVals= allMinVals(condInd);
    currMaxVals= allMaxVals(condInd);
    
    interpEMG = {};
    interpSegEMG = {};
    hold on;
    maxTrialLength=max(cellfun(@length,currTrials));
    numSegs = mode(cellfun(@length,currSegs));
    badSegTrials = cellfun(@(a) length(a)==numSegs, currSegs);
    badCondTrials = badCondTrials(badSegTrials);
    currTrials = currTrials(badSegTrials);
    currSegs = currSegs(badSegTrials);
    currMinVals= currMinVals(badSegTrials);
    currMaxVals= currMaxVals(badSegTrials);
    trialChange = diff(cell2mat(currMinVals));
    trialChange = [0 trialChange];
    if(sum(badSegTrials)~=length(badSegTrials))
        disp(fName{f});
        disp([num2str(sum(badSegTrials==0)), ' trials ignored']);
    end
    allSegLengths = reshape(cell2mat(currSegs),numSegs, size(currSegs,2))';
    maxSegLengths = mean(diff(allSegLengths,1,2));
    maxSegLengths = [mean(allSegLengths(:,1)), maxSegLengths];
    maxSegLengths = cumsum(maxSegLengths);
    tInd = 0;
    color = 1;
    for t = 1:length(currTrials)
        if(trialChange(t)~=0 & ~isnan(trialChange(t)))
            color = color+3;
            f = f + 1;
        end
        if(~badCondTrials(t))
            tInd = tInd + 1;
            sig = filter(1/100*ones(1,100),1,abs(currTrials{t}));
            %sig = (sig-allMinVals{t})/(allMaxVals{t}-allMinVals{t});
            trialLength = length(sig);
            newX = linspace(1,trialLength,trialLength);
            interpEMG{tInd} = interp1(1:trialLength,sig,newX);
            interpEMG{tInd}(end:end+(maxTrialLength-trialLength)) = nan(1,1+maxTrialLength-trialLength);
            %interpEMG{tInd}=smoothdata(interp1(1:trialLength,sig,newX), 'movmean', 100);
            trialSegs = currSegs{t};
            segStart = 1;
            maxStart = 1;
%             for s = 1:length(trialSegs)
%                 segEnd = trialSegs(s);
%                 maxEnd = maxSegLengths(s);
%                 interpSegEMG{tInd,s} = (interp1(segStart:segEnd,sig(segStart:segEnd),...
%                     linspace(segStart,segEnd,maxEnd-maxStart)));...
%                     %-allMinVals{t})/(allMaxVals{t}-allMinVals{t});
%                 plot(maxStart:maxEnd-1, interpSegEMG{tInd,s}, 'Color', colors(color,:));
%                 segStart = trialSegs(s);
%                 maxStart = maxSegLengths(s);
%             end
        end
    end
    error = [];
    currentTrial = {};
    start = maxSegLengths(2);
    finish = maxSegLengths(4);
    %interpEMG = cellfun(@(a) movmean(a, 100), interpEMG, 'UniformOutput', false);
%     for t = 1:length(interpEMG)
%         [currentTrial{1:length(interpEMG)}] = deal(interpEMG{t});
%         error(t) = sum(cellfun(@(a,b) nansum((a(start:finish)-b(start:finish)).^2), interpEMG, currentTrial));
%     end
    [~,avgTrial] = min(error);
    %plot([interpSegEMG{avgTrial,:}],'k','LineWidth',2);
    plot(nanmean(reshape(cell2mat(interpEMG'),length(interpEMG),maxTrialLength)), 'k' , 'LineWidth', 2);
    saveEMG = reshape(cell2mat(interpEMG'),length(interpEMG),maxTrialLength);
    %save(['S:\Lab\Gilligan\Sorted EMG Data\Results\Models\Averages_Norm\',muscle,'_',conds{c}], 'saveEMG');
    %save(['S:\Lab\Gilligan\Sorted EMG Data\Results\Models\Averages_Norm\',muscle,'_',conds{c},'_Times'], 'maxSegLengths');
    ax = gca;
    if(strcmp(conds{c},'Photocell'))
        set(ax,'XLim', [maxSegLengths(1)-2000/4 maxSegLengths(4)+2000/2]);
        titleLoc = (maxSegLengths(3) + maxSegLengths(4))/2;
    else
        set(ax,'XLim', [maxSegLengths(1)-2000/8 maxSegLengths(5)+2000]);
        titleLoc = (maxSegLengths(4) + maxSegLengths(5))/2;
        
    end
    text(titleLoc, 175, sprintf('%s\n\t%s%d', printName{c}, 'n=', length(interpEMG)), 'HorizontalAlignment', 'center', 'FontSize', 26);
    set(ax, 'YLim', [0 1]);
    ax.XTickLabel = arrayfun(@(a) num2str(a), ax.XTick*.5, 'UniformOutput', false);
    
    
    yl = ylim;
    for s = 1:length(maxSegLengths)
        plot([maxSegLengths(s) maxSegLengths(s)],yl,'r--')
        if(s>1 && s<length(currSegs{t})-1)
            xLoc = (maxSegLengths(s-1) + maxSegLengths(s))/2;
            segT = text(xLoc, yl(2)+12, segNames{s-1}, 'Color', 'r', 'HorizontalAlignment', 'center', 'FontSize', 26);
        end
    end
    if(c==2 || c==4)
        ax.YTickLabel = {};
    end
    ax.XTickMode = 'manual';
end
