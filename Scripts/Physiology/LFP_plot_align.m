sessionDate = '03_01_2019';
dirName = ['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
    sessionDate, '\Physiology\Results\'];
directory = dir(dirName);
channel = [];
channelMap = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32];
segNames = {'Rx', 'R', 'G', 'L', 'H', 'W'};

if(~isempty(channel))
    for d = 1:size(directory,1)
        if(strcmp(directory(d).name, ['Gilligan_', sessionDate, '_',channel,'.mat']))
            directory = [directory(d)];
            break;
        end
    end
end
directory = directory(~[directory.isdir]);
directory = directory(find(contains({directory.name}, 'Gilligan')));
[~,sortInd] = natsort({directory.name});
directory = directory(sortInd);
for d=1:size(directory,1)
    if(mod(d,4)==1)
        figure();
        count = 0;
    end
    xLims = {[-1, 1.5]};
    gap = .1;
    l = load([directory(channelMap(d)).folder, '\', directory(channelMap(d)).name]);
    eventData = l.sortedLFPData.LFPData;
    segmentTimes = l.sortedLFPData.SegTimes;
    arduinoTrials = l.sortedLFPData.ArduinoData;
    Fs = l.sortedLFPData.SampleRate;
    conds = l.sortedLFPData.Conditions;
    %conds = conds(~cellfun(@(a) strcmp(a,'Rest'), conds));
    condSegs = l.sortedLFPData.ConditionSegments;
    window = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);
    %conds = {'Extra Small Sphere', 'Large Sphere','Photocell'};
    for c = 1:length(conds)
        subplot(4,4,(4*count)+c);
        if(mod(d,4)==1)
            title(conds(c));
        end
        if(c==1)
        ylabel(d);
        end
        %set(gca, 'FontSize', 24);
        condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
        currTrials = eventData(condInd);
        currSegs = segmentTimes(condInd);
        
        hold on;
        numSegs = mode(cellfun(@length,currSegs));
        badSegTrials = cellfun(@(a) length(a)==numSegs, currSegs);

        currTrials = currTrials(badSegTrials);
        currSegs = currSegs(badSegTrials);
        
        baselines = cellfun(@(a,b) mean(a(1:b(1))), currTrials, currSegs, 'UniformOutput', false);
        %currTrials = cellfun(@(a,b) (a-b)/b, currTrials, baselines, 'UniformOutput', false);
        %currTrials = cellfun(@(a) (a-normMin)/(normMax-normMin), currTrials, 'UniformOutput', false);
        if(sum(badSegTrials)~=length(badSegTrials))
            disp([num2str(sum(badSegTrials==0)), ' trials ignored']);
        end
       averageSegs = mean(reshape(cell2mat(currSegs),numSegs,size(currTrials,2))');
        if(strcmp(conds{c}, 'Photocell'))
            alignSegs = {'StartHold'};
        elseif(strcmp(conds{c}, 'Rest'))
            alignSegs = {'StartReach'};
        else
            alignSegs = {'StartLift'};
            reachLiftTime = averageSegs(find(strcmp(condSegs{c}, 'StartLift')))-averageSegs(find(strcmp(condSegs{c}, 'StartReach')));
            reachLiftTime = averageSegs(4) - averageSegs(2);
            reachLiftTime = round(reachLiftTime/2);
            win{1}(2) = reachLiftTime+5;
            win{2}(1) = -(reachLiftTime+5);
        end
        plotStart = 0;
        xTicks = [];
        labels = [];
        win = window;
        numBad = 0;

        for align = 1:length(win)
            if(align==1)
                plotStart = 0;
            else
                plotStart = plotStart + (gap*Fs) +win{align-1}(2) + abs(win{align}(1));
            end
            alginSegInd = find(strcmp(condSegs{c},alignSegs{align}))-1;
            alignedEMGSigs = cellfun(@(a,b) a(b(alginSegInd)+win{align}(1):...
                b(alginSegInd)+win{align}(2)),currTrials, currSegs,'UniformOutput', false);
            yvals = reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))';
            alignedDataY{c,align} = yvals;
            alignedDataX{c,align} = plotStart+win{align}(1):plotStart+win{align}(2);
            
            cellfun(@(a) plot(plotStart+win{align}(1):plotStart+win{align}(2),a, 'b'), alignedEMGSigs);
            
            plot(plotStart+win{align}(1):plotStart+win{align}(2), ...
                nanmean(reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))'),...
                'k', 'LineWidth', 2);
            beforeSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:-1:1)));
            afterSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:end)));
            avgSegs = [beforeSegs plotStart afterSegs];
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=plotStart+win{align}(1) && avgSegs(s)<=plotStart+win{align}(2))
                    plot([avgSegs(s) avgSegs(s)],[-3000 3000],'k--')
                end
            end
            segs{c,align} = win{align}./Fs;
            xTicks = [xTicks, plotStart+[win{align}(1) 0 win{align}(2)]];
            labels = [labels, [win{align}(1) 0 win{align}(2)]];
        end
            ylim([-3000, 3000]);
    end

    expind = regexp(directory(d).name,'_');
        count = count + 1;
        if(count==4)
            saveas(gcf,[dirName, '\Figures\Group', num2str(d)]);
        end
end