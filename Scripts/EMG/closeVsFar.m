sessionDate = '10_11_2018';
directory = dir(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
    sessionDate, '\EMG\Results_New\Close']);
directory = directory(~[directory.isdir]);
muscle = '';
modelNames = {'Deltoid.mat', 'Biceps.mat', 'Triceps.mat', 'WristExtensor.mat',...
    'WristFlexor.mat','DigitExtensor.mat', 'DigitFlexor.mat'};
printName = {'Deltoid', 'Biceps', 'Triceps', 'WristEx', 'WristFlex','DigitEx', 'DigitFlex',};

graspSegs = {'Rx', 'R', 'G','L', 'H', 'W'};
noGraspSegs = {'Rx', 'R', 'H', 'W'};
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fax = axes('Position',[0 0 1 1],'Visible','off');

text(.5,.02,'Times (ms)', 'FontSize', 26, 'FontWeight', 'bold','HorizontalAlignment', 'Center');
tx = text(.02,0.5,'Voltage (uV)', 'FontSize', 26, 'FontWeight', 'bold', 'HorizontalAlignment', 'Center');
set(tx, 'Rotation', 90);
hx = tight_subplot(7,2,[.05 .02],[.1 .1],[.07 .03]);
if(~isempty(muscle))
    for d = 1:size(directory,1)
        if(strcmp(directory(d).name, ['Gilligan_', sessionDate, '_',muscle,'.mat']))
            directory = [directory(d)];
            break;
        end
    end
end


for d=1:size(directory,1)
    alignSegs = {'StartReach', 'StartLift','StartWithdraw'};
    xLims = {[-.25 .2],[-.1, .3],[-.5 .5]};
    gap = .35;
    
    
    
    
    
    for f = 1:2
        alignSegs = {'StartReach', 'StartLift','StartWithdraw'};
        xLims = {[-.25 .2],[-.1, .3],[-.5 .5]};
        gap = .15;
        directory = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
            sessionDate, '\EMG\Results_New\'];
        if (f==1)
            directory = dir([directory, 'Close']);
            directory = directory(~[directory.isdir]);
            di = find(strcmp({directory.name}, ['Gilligan_', sessionDate,'_',modelNames{d}]));
            indsU = regexp(directory(di).name, '_');
            %text(fax, .15,.97, directory(d).name(indsU(end)+1:end-4), 'FontSize', 28, 'FontWeight', 'bold','HorizontalAlignment', 'Center');
            muscle = directory(di).name(indsU(end)+1:end-4);
        elseif(f==2)
            directory = dir([directory, 'Far']);
            directory = directory(~[directory.isdir]);
            di = find(contains({directory.name}, muscle));
        end
        if(~isempty(di))
            l = load([directory(di).folder, '\', directory(di).name]);
            eventData = l.sortedEMGData.EMGData;
            segmentTimes = l.sortedEMGData.SegTimes;
            arduinoTrials = l.sortedEMGData.ArduinoData;
            Fs = l.sortedEMGData.SampleRate;
            conds = l.sortedEMGData.Conditions;
            condSegs = l.sortedEMGData.ConditionSegments;
            window = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);
            badTrial = EMG_getBadTrials(eventData, segmentTimes, Fs);
            conds = {'Extra Small Sphere', 'Large Sphere'};
            for c = 1:length(conds)
                axes(hx(2*(d-1)+f))
                
                set(gca, 'FontSize', 16);
                condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
                badCondTrials = badTrial(condInd);
                currTrials = eventData(condInd);
                currSegs = segmentTimes(condInd);
                
                
                hold on;
                numSegs = mode(cellfun(@length,currSegs));
                badSegTrials = cellfun(@(a) length(a)==numSegs, currSegs);
                badCondTrials = badCondTrials(badSegTrials);
                currTrials = currTrials(badSegTrials);
                currSegs = currSegs(badSegTrials);
                
                currTrials = currTrials(~badCondTrials);
                currSegs = currSegs(~badCondTrials);
                
                currTrials = cellfun(@(a) filter(1/100*ones(1,100),1,abs(a)), currTrials, 'UniformOutput', false);
                if(sum(badSegTrials)~=length(badSegTrials))
                    disp([num2str(sum(badSegTrials==0)), ' trials ignored']);
                end
                averageSegs = mean(reshape(cell2mat(currSegs),numSegs,size(currTrials,2))');
                if(strcmp(conds{c}, 'Photocell'))
                    win = {[-.5*Fs, .35*Fs], [-.50*Fs, .5*Fs]};
                    alignSegs = {'StartReach','StartHold'};
                    segNames = noGraspSegs;
                else
                    win = window;
                    reachLiftTime = averageSegs(4) - averageSegs(2);
                    reachLiftTime = round(reachLiftTime/2);
                    win{1}(2) = reachLiftTime+5;
                    win{2}(1) = -(reachLiftTime+5);
                    segNames = graspSegs;
                    if(c==1)
                        prevSegs = win;
                    end
                end
                if(c==1)
                    cf = 'r';
                    averageSegsT{f} = averageSegs;
                    
                else
                    if(f==2)
                        averageSegs = [averageSegs(1:4) averageSegs(4) + 49 averageSegs(5:end)];
                    end
                    %averageSegsT{f} = (averageSegsT{f} + averageSegs)./2;
                    cf = 'b';
                end
                plotStart = 0;
                xTicks = [];
                labels = [];
                for align = 1:length(win)
                    if(align==1)
                        plotStart = 0;
                    else
                        if(c==2)
                            plotStart = plotStart+ (gap*Fs) + prevSegs{align-1}(2) + abs(prevSegs{align}(1));
                        else
                            plotStart = plotStart + (gap*Fs) +win{align-1}(2) + abs(win{align}(1));
                        end
                    end
                    alginSegInd = find(strcmp(condSegs{c},alignSegs{align}))-1;
                    alignedEMGSigs = cellfun(@(a,b) a(b(alginSegInd)+win{align}(1):...
                        b(alginSegInd)+win{align}(2)),currTrials, currSegs,'UniformOutput', false);
                    
                    %cellfun(@(a) plot(plotStart+win{align}(1):plotStart+win{align}(2),a, cf), alignedEMGSigs);
                    %                     plot(plotStart+win{align}(1):plotStart+win{align}(2), ...
                    %                         mean(reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))'),...
                    %                         cf, 'LineWidth', 2);
                    yvals = reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))';
                    shadedErrorBar(plotStart+win{align}(1):plotStart+win{align}(2), mean(yvals), std(yvals),'lineprops', {cf,'LineWidth',2});
                    if(c==2 || c ==1)
                        beforeSegs = plotStart + cumsum(diff(averageSegsT{f}(alginSegInd:-1:1)));
                        afterSegs = plotStart + cumsum(diff(averageSegsT{f}(alginSegInd:end)));
                        avgSegs = [beforeSegs plotStart afterSegs];
                        for s = 1:length(avgSegs)
                            if(avgSegs(s)>=plotStart+win{align}(1) && avgSegs(s)<=plotStart+win{align}(2))
                                plot([avgSegs(s) avgSegs(s)],[0 max(cellfun(@max, currTrials))],['k--']);
                                if(s>1 && s<=length(segNames)+1 && d==1)
                                    ax = gca();
                                    yl = get(ax, 'YLim');
                                    xLoc = (avgSegs(s-1) + avgSegs(s))/2;
                                    segT = text(xLoc, yl(2)+range(yl)/3 - (range(yl)/6*~(mod(s,2))), segNames{s-1}, 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');
                                    %set(segT,'Rotation',25)
                                end
                            end
                        end
                    end
                    xTicks = [xTicks, plotStart+[win{align}(1) 0 win{align}(2)]];
                    labels = [labels, [win{align}(1) 0 win{align}(2)]];
                end
                
                for l = 1:2
                    if l==1
                        lc = 'r';
                    elseif(l==2)
                        lc = 'b';
                    end
                    h(l) = plot(hx(1),NaN,NaN,'Color', lc, 'LineWidth', 2);
                end
                legend(h, {'Small Sphere', 'Large Sphere'}, 'FontSize', 14);
                if(f==1 && c==1)
                    tx = text(fax, .98, 1-.13*d,printName{d},'FontSize', 20, 'HorizontalAlignment', 'center');
                end
                set(tx, 'Rotation', -90);
                set(gca, 'XTick', xTicks);
                if(d~=length(printName))
                    set(gca, 'XTickLabel', []);
                else
                    set(gca, 'XTickLabel', arrayfun(@num2str, (labels .*.5), 'UniformOutput', false));
                end
                if(d==1)
                    title(sprintf('%s\n\t%s%d', 'Close'));
                end
            end
        end
    end
end