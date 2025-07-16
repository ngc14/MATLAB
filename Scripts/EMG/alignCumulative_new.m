clear all;
close all;
monkey = 'Gilligan';
directory = dir(['S:\Lab\',monkey,'\All Data']);
modelNames = {'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor',...
    'Digit Extensor','Digit Flexor'};
segNames = {'RT', 'R', 'L', 'H', 'W'};
alignCondSegs = containers.Map({'Extra Small Sphere','Large Sphere', 'Photocell','Rest'},...
    {{'GoSignal','StartReach', 'StartHold'},{'GoSignal','StartReach', 'StartHold'},...
    {'GoSignal','StartReach', 'StartHold'},{'GoSignal','StartReplaceSuccess', 'StartReward'}});
saveVals = 1;
saveFigs = 0;
xLims = {[-.6, .15], [-.15, .6],[-.15,.6]};
gap = .15;
Fs = 2000;
smoothKernel = .15; %in seconds
smoothKernel = smoothKernel*Fs;
xLims = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);
colors = colorcube;
exclude = string([datetime('04/01/2018','InputFormat','MM/dd/yyyy'):datetime('12/30/2018'),...
    datetime('12/31/2019','InputFormat','MM/dd/yyyy'):datetime('today')],'dd-MMM-yyyy');
directory = directory(arrayfun(@(c) isfolder(strjoin([c.folder,c.name,"EMG"],'\')),directory));
pDates = strings(1,size(directory,1));
for d = 1:size(directory,1)
    underscores = regexp(directory(d).name, '_');
    fileDate = directory(d).name(underscores(1)+1:end);
    if(all(diff(underscores)==3))
        pDates(d) = string(datetime(fileDate,'InputFormat','MM_dd_yyyy'),'dd-MMM-yyyy');
    else
        pDates(d) = string(datetime(fileDate,'InputFormat','yyyy_MM_dd'),'dd-MMM-yyyy');
    end
end
filesToProcess = ~contains(pDates,exclude);
%"04_15_2019"
%allDates = ["04_16_2019","04_17_2019","04_18_2019","04_19_2019","04_22_2019","04_23_2019","04_26_2019","04_29_2019","05_03_2019","05_17_2019","05_29_2019","05_30_2019","06_02_2020","06_03_2020","06_05_2020","06_08_2020","06_10_2020","06_11_2020","06_12_2019","06_12_2020","07_02_2020","07_03_2020","07_07_2020","07_08_2019","07_09_2019","07_09_2020","07_10_2019","07_11_2019","07_12_2019","07_14_2020","07_15_2020","07_16_2019","07_16_2020","07_17_2019","07_19_2019","07_23_2019","07_25_2019","08_01_2019","08_05_2019","08_07_2019","08_09_2019","09_05_2019","09_09_2019","09_16_2019","09_23_2019","09_26_2019","11_12_2019","11_19_2019","11_22_2019","11_26_2019","12_03_2019","12_06_2019"];
% filesToProcess = arrayfun(@(d) find(dates==datenum(d)), allDates);

directory = directory(filesToProcess);
for mm =1:length(modelNames)
    eventData = [];
    segmentTimes = [];
    arduinoTrials = [];
    conds = [];
    fName = [];
    condSegs = [];
    muscle = modelNames{mm};
    numSessions = 0;
    alignedEMG = repmat({{}},4,1);
    allSegsCond = cell(4,1);
    rawTrialsCond = cell(4,1);
    allDates = cell(4,1);
    locationCond = cell(4,1);
    
    for f = 1:size(directory,1)
        rootFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\'];
        if(isfolder([rootFolder,'Close']))
            rootFolder = [rootFolder,'Close\'];
            disp(["Close ",directory(f).name, " ",modelNames(mm)])
        end
        if(strcmp(monkey,'Gilligan'))
            if(strcmp(muscle,'Digit Extensor'))
                resultFolder = [rootFolder,directory(f).name, '_D4D5Extensor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [rootFolder,directory(f).name, '_Digit Extensor.mat'];
                end
            elseif(strcmp(muscle,'Digit Flexor'))
                resultFolder = [rootFolder,directory(f).name, '_DistalFlexor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [rootFolder,directory(f).name, '_Digit Flexor.mat'];
                end
            elseif(strcmp(muscle,'Wrist Extensor'))
                resultFolder = [rootFolder,directory(f).name, '_Extensor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [rootFolder,directory(f).name, '_Wrist Extensor.mat'];
                end
            elseif(strcmp(muscle,'Wrist Flexor'))
                resultFolder = [rootFolder,directory(f).name, '_ProximalRadialFlexor.mat'];
                if(~exist(resultFolder,'file'))
                    resultFolder = [rootFolder,directory(f).name, '_Wrist Flexor.mat'];
                end
            else
                resultFolder = [rootFolder, directory(f).name, '_', muscle, '.mat'];
            end
        else
            resultFolder = [rootFolder, directory(f).name, '_', muscle, '.mat'];
        end
        if(exist(resultFolder, 'file'))
            l = load(resultFolder);
            if(isfield(l,'sortedEMGData'))
                eventData = {};
                segmentTimes = {};
                fName = {l.sortedEMGData.Date};
                conds = l.sortedEMGData.Conditions;
                condSegs = l.sortedEMGData.ConditionSegments;
                condTrials = l.sortedEMGData.ArduinoData(:,1);
                allEMG = l.sortedEMGData.EMGData;
                condSegTimes = cellfun(@str2num, l.sortedEMGData.ArduinoData(:,2:end-2),'UniformOutput',false);
                condSegTimes(cellfun(@isempty,condSegTimes)) = {NaN};
                condSegTimes = cellfun(@(t) cumsum([10002 t]), num2cell(cell2mat(condSegTimes).*Fs/1000,2), 'UniformOutput', false);
                for c = 1:length(conds)
                    if(any(contains(condTrials,conds{c})))
                        eventData{c} = allEMG(cellfun(@(s) strcmp(s,conds{c}),condTrials));
                        segmentTimes{c} = l.sortedEMGData.SegTimes(cellfun(@(s) strcmp(s,conds{c}),condTrials));
                    end
                end
            else
                fName = {directory(f).name(10:end)};
                conds = l.Conditions;
                condSegs = l.ConditionSegs;
                eventData = l.rawTrials;
                segmentTimes = l.segs;
                
            end
                badConds = ~contains(conds,alignCondSegs.keys);
                conds = conds(~badConds);
                condSegs = condSegs(~badConds);
                eventData = eventData(~badConds);
                segmentTimes = segmentTimes(~badConds);
            %numTrials = length(l.sortedEMGData.ArduinoData);
            numSessions = numSessions+1;

                maxTrialLength = max(cellfun(@(a)max(cellfun(@(l) max(0,length(l)), a)), eventData));
                %%
            paddedEventData = cellfun(@(a) cellfun(@(b) padArrays(abs(b),maxTrialLength), ...
                a, 'UniformOutput', false), eventData, 'UniformOutput', false);
            eventData = cellfun(@(a) cellfun(@(b) conv(b,gausswin(smoothKernel)/...
                sum(gausswin(smoothKernel)),'same'),a,'UniformOutput', false),  paddedEventData, 'UniformOutput', false);
            %%
            mismatch = cellfun(@(e) find(cellfun(@length,e)~=mode(cellfun(@length,e))), eventData, 'UniformOutput',false);
            if(any(cellfun(@any,mismatch)))
                disp('wrong')
            end
            for c = 1:length(conds)
                alignSegs = values(alignCondSegs, conds(c));
                alignSegs = alignSegs{1};
                currTrials = eventData{c};
                currSegs = segmentTimes{c};
                currEventSegs = condSegs{c};                
                numSegs = length(currEventSegs);
                alignSession = {};
                plotted = false(1,length(currEventSegs));
                badCondTrials = cellfun(@(n) all(isnan(n)), currTrials);
                badSegTrials = cellfun(@(s) any(isnan(s(cellfun(@(e) ...
                    find(strcmp(currEventSegs,e)),alignSegs)))), currSegs);
                goodTrials = ~(badCondTrials | badSegTrials)';
                if(any(cellfun(@(t) ~all(isnan(t)),currTrials)))
                    currSegs = cellfun(@(s) [s(1:2), NaN(1,numSegs-length(s)), s(3:end)], currSegs, 'UniformOutput',false);
                    allSegs = reshape(cell2mat(currSegs),numSegs,size(currSegs,2))';
                    for align = 1:length(alignSegs)
                        alignSegInd = find(strcmp(currEventSegs,alignSegs{align}));
                        alignSession(:,align) = cellfun(@(a,b) a(max(1,b(alignSegInd)+xLims{align}(1)):...
                            min(length(a),b(alignSegInd)+xLims{align}(end))),currTrials, currSegs,'UniformOutput', false)';
                        if(saveFigs)
                            subplot(2,2,c); hold on;
                            if(align==1)
                                plotStart = 0;
                            else
                                plotStart = plotStart + (gap*Fs) + xLims{align-1}(end) + abs(xLims{align}(1));
                            end
                            cellfun(@(a) plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2),...
                                a, 'Color',colors(numSessions,:)),alignSession(goodTrials,align));
                            plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2), ...
                                nanmean(reshape([alignSession{goodTrials,align}],...
                                size(alignSession{find(goodTrials,1),align},2),...
                                sum(goodTrials))'),'k', 'LineWidth', 3);
                            averageSegs = nanmean(allSegs,1);
                            beforeSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:-1:1)));
                            afterSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:end)));
                            avgSegs = [beforeSegs plotStart afterSegs];
                            for s = 1:length(avgSegs)
                                if(avgSegs(s)>=plotStart+xLims{align}(1) && ...
                                        avgSegs(s)<=plotStart+xLims{align}(end) && ...
                                        (~plotted(s) ||avgSegs(s)==plotStart))
                                    if(avgSegs(s)==plotStart)
                                        plotColor = 'b';
                                    else
                                        plotColor = 'r';
                                    end
                                    plotted(s) = true;
                                    plot([avgSegs(s) avgSegs(s)],...
                                        [0 maxTrials],...
                                        'Color',plotColor,'LineStyle','--')
                                end
                            end
                            if(align>1)
                                set(gca, 'XTick', sort([get(gca, 'XTick'),plotStart+xLims{align}(1), plotStart, plotStart+xLims{align}(end)]));
                                set(gca, 'XTickLabels', [string(get(gca,'XTickLabel')'),num2str(xLims{align}(1)/Fs), "",num2str(xLims{align}(end)/Fs)]);
                            else
                                xticks([xLims{align}(1), plotStart, xLims{align}(end)]);
                                xticklabels({num2str(xLims{align}(1)/Fs), "",num2str(xLims{align}(end)/Fs)});
                            end
                            if(align==length(alignSegs))
                                allticks = get(gca,'XTick');
                                xlim([allticks(1), allticks(end)]);
                            end
                        end
                    end
                    allDates{c}(end+1:end+sum(goodTrials)) = fName;
                    alignedEMG{c}(end+1:end+sum(goodTrials),:) = alignSession(goodTrials,:);
                    allSegsCond{c}(end+1:end+sum(goodTrials)) = currSegs(goodTrials)';
                    rawTrialsCond{c}(end+1:end+sum(goodTrials)) =currTrials(goodTrials);
                    if(c==length(conds))
                        if(saveFigs)
                            linkaxes(get(gcf,'Children'));
                            ylim([0, maxTrials])
                            saveas(gcf,[directory(f).folder,'\',directory(f).name,'\EMG\Results_New\',directory(f).name,'_',muscle,'.png']);
                            savefig(gcf,[directory(f).folder,'\',directory(f).name,'\EMG\Results_New\',directory(f).name,'_',muscle],'compact')
                        end
                        if(saveVals)
                            currSessionInds = cellfun(@(dc) cellfun(@(d) datenum(d)==datenum(fName{1}), dc), allDates,'UniformOutput', false);
                            sessionEMG.alignedSig = cellfun(@(s,i) s(i,:), alignedEMG,currSessionInds,'UniformOutput',false);
                            sessionEMG.rawTrials = cellfun(@(s,i) s(i),rawTrialsCond,currSessionInds,'UniformOutput',false);
                            sessionEMG.segs = cellfun(@(s,i) s(i),allSegsCond,currSessionInds,'UniformOutput',false);
                            sessionEMG.Fs = Fs;
                            sessionEMG.alignments = values(alignCondSegs,conds);
                            sessionEMG.alignWindows = cellfun(@(xl) xl/Fs, xLims, 'UniformOutput',false);
                            sessionEMG.Conditions = conds;
                            sessionEMG.ConditionSegs = condSegs;
                            save([directory(f).folder,'\',directory(f).name,'\EMG\Results_New\',directory(f).name,'_',muscle],...
                                '-struct','sessionEMG','-v7.3');
                        end
                        clear sessionEMG;
                    end
                end
            end
            close all;
        end
        
    end
    % all sessions
    for c = 1:length(conds)
        allSegs = allSegsCond{c};
        plotted = false(1,length(condSegs{c}));
        sessionInds = cumsum([1 diff(cellfun(@length,rawTrialsCond{c}))~=0]);
        if(saveFigs)
            subplot(2,2,c);
            hold on;
            for align = 1:length(alignSegs)
                if(align==1)
                    plotStart = 0;
                else
                    plotStart = plotStart + (gap*Fs) + xLims{align-1}(end) + abs(xLims{align}(1));
                end
                alignSegInd = find(strcmp(condSegs{c},alignSegs{align}));
                if(isempty(alignSegInd)  & strcmp(conds{c}, 'Photocell'))
                    alignSegInd = find(strcmp(condSegs{c},'StartHold'));
                end
                badCondTrials = cellfun(@length, alignedEMG{c}(:,align))==1;
                badSegTrials = cellfun(@(s) all(isnan(s)), allSegs);
                badTrials = badCondTrials | badSegTrials';
                cellfun(@(s,b) plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2),...
                    s, 'Color',colors(b,:)),alignedEMG{c}(~badTrials,...
                    align),num2cell(sessionInds(~badTrials)'), 'UniformOutput', false);
                plot(plotStart+xLims{align}(1):plotStart+xLims{align}(2),...
                    nanmean(cell2mat(alignedEMG{c}(~badTrials,align)),1), 'k','LineWidth',3);                
                averageSegs = nanmean(cell2mat(allSegs'),1);
                beforeSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:-1:1)));
                afterSegs = plotStart + cumsum(diff(averageSegs(alignSegInd:end)));
                avgSegs = [beforeSegs plotStart afterSegs];
                for s = 1:length(avgSegs)
                    if(avgSegs(s)>=plotStart+xLims{align}(1) && ...
                            avgSegs(s)<=plotStart+xLims{align}(end) && ...
                            (~plotted(s) ||avgSegs(s)==plotStart))
                        if(avgSegs(s)==plotStart)
                            plotColor = 'b';
                        else
                            plotColor = 'r';
                        end
                        plotted(s) = true;
                        plot([avgSegs(s) avgSegs(s)],[0 max(cellfun(@nanmax,rawTrialsCond{c}))],...
                            'Color', plotColor, 'LineStyle','--');
                    end
                end
            end
            title(sprintf('%s\n\t%s%d', conds{c}, 'n=', length(rawTrialsCond{c})));
        end
    end
    if(saveVals)
        sessionEMG.alignedSig = alignedEMG;
        sessionEMG.rawTrials = rawTrialsCond;
        sessionEMG.segs = allSegsCond;
        sessionEMG.dates = allDates;
        sessionEMG.Fs = Fs;
        sessionEMG.alignments = alignSegs;
        sessionEMG.alignWindows = cellfun(@(xl) xl/Fs, xLims, 'UniformOutput',false);
        sessionEMG.Conditions = conds;
        sessionEMG.ConditionSegs = condSegs;
        sessionEMG.LocationConditions = locationCond;
        save(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\',muscle],...
            '-struct','sessionEMG','-v7.3');
    end
    if(saveFigs)
        linkaxes(get(gcf,'Children'));
        arrayfun(@(ax) xticklabels(ax,xticks./Fs),get(gcf,'Children'));
        saveas(gcf,['S:\Lab\', monkey,'\Sorted EMG Data\Averages_Raw\',muscle,'.png']);
        savefig(gcf,['S:\Lab\', monkey,'\Sorted EMG Data\Averages_Raw\',muscle],'compact');
    end
    clearvars -except modelNames directory xLims alignSegs gap segNames Fs saveFigs saveVals winSize monkey colors smoothKernel alignCondSegs
    close all;
end

function padded = padArrays(arr, length)
padded = arr;
padded(:,end+1:length)  = NaN;
end