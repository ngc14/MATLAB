function EMG_derivatives = EMG_plot_derivative(sortedEMGData, showPlot)
events = {'StartReach', 'StartLift', 'StartWithdraw'};
eventData = sortedEMGData.EMGData;
segmentTimes = sortedEMGData.SegTimes;
arduinoTrials = sortedEMGData.ArduinoData;
Fs = sortedEMGData.SampleRate;
conds = sortedEMGData.Conditions;
condSegs = sortedEMGData.ConditionSegments;

badTrial = EMG_getBadTrials(eventData, segmentTimes, Fs);
if(showPlot)
    figure('Name', sortedEMGData.Muscle);
end
% split EMG data by condition
for c = 1:length(conds)
    condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
    badCondTrials = badTrial(condInd);
    currTrials = eventData(condInd);
    currSegs = segmentTimes(condInd);
    % exclude trials that have wrong segmentations
    numSegs = mode(cellfun(@length,currSegs));
    badSegTrials = cellfun(@(a) length(a)==numSegs, currSegs);
    badCondTrials = badCondTrials(badSegTrials);
    currTrials = currTrials(badSegTrials);
    currSegs = currSegs(badSegTrials);
    if(sum(badSegTrials)~=length(badSegTrials))
        disp([num2str(sum(badSegTrials==0)), ' trials ignored']);
    end
    
    centeredSig = {};
    for s = 1:3
        if(showPlot)
            subplot(4,3,3*(c-1)+s)
            hold on;
        end
        tInd = 0;
        eventInd = find(strcmp(condSegs{c}, events(s)));
        for t = 1:length(currTrials)
            % skip bad trials
            if(~badCondTrials(t))
                tInd = tInd + 1;
                windowSize = 5;
                b = (1/100)*ones(1,100);
                sig = filter(b,1,abs(currTrials{t}));
                % get smoothed derivative
                centeredSig{tInd,s} = filter(b,1,diff(sig...
                    (currSegs{t}(eventInd)-Fs:currSegs{t}(eventInd)+Fs)));
                if(showPlot)
                    plot(centeredSig{tInd,s});
                end
            end
        end
        if(showPlot)
            avgLine = mean(cell2mat(centeredSig(:,s)));
            plot(avgLine, 'k', 'LineWidth',2);
            if(c==1)
                title(events(s));
            end
            if(s==1)
                ylabel([conds{c}, ': ', num2str(tInd)]);
            end
            if(~isempty(avgLine))
                avgLine = avgLine(~isnan(avgLine));
                ylim([min(avgLine)-std(avgLine) max(avgLine)+std(avgLine)]);
            end
        end
    end
    tempCond = conds{c};
    tempCond = tempCond(~isspace(tempCond));
    if(isempty(centeredSig))
        centeredSig = 0;
    end
    EMG_derivatives.(tempCond) = centeredSig;
end
end