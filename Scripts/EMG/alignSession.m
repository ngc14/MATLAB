
sessionDate = '06_06_2019';
monkey = 'Gilligan';
dirName = ['S:\Lab\',monkey,'\All Data\',monkey,'_', ...
    sessionDate, '\EMG\All_Trials\'];
directory = dir(dirName);
muscle = '';


segNames = {'Rx', 'R', 'G', 'L', 'H', 'W'};

% if(~isempty(muscle))
%     for d = 1:size(directory,1)
%         if(strcmp(directory(d).name, [monkey,'_', sessionDate, '_',muscle,'.mat']))
%             directory = [directory(1); directory(2); directory(d)];
%             break;
%         end
%     end
% end
directory = directory(~[directory.isdir]);
directory = directory(find(contains({directory.name}, monkey)));
figure();
    hx = [subplot(1,3,1),subplot(1,3,2),subplot(1,3,3)];%tight_subplot(1,3,[.1 .01],[.1 .1],[.07 .01]);
 %%   
 colorP = distinguishable_colors(length(directory));
for d=1:size(directory,1)
    alignSegs = {'StartReach'};
    figure('units','normalized','outerposition',[0 0 1 1]);
    xLims = {[-.5 2.5]};
    gap = .25;
    fullDir = dir([directory(d).folder, '\', directory(d).name]);
    l = load([fullDir(1).folder,'\',fullDir(1).name]);
    eventData = l.sortedEMGData.EMGData;
    segmentTimes = l.sortedEMGData.SegTimes;
    arduinoTrials = l.sortedEMGData.ArduinoData;
    Fs = l.sortedEMGData.SampleRate;
    segmentTimes = cellfun(@(s) Fs*((s-(s(1)-1))./(Fs/1000)),segmentTimes,'UniformOutput',false);
    conds = l.sortedEMGData.Conditions;
    condSegs = l.sortedEMGData.ConditionSegments;
    window = cellfun(@(a) a.*Fs, xLims, 'UniformOutput', false);
    badTrial = getBadTrials(eventData, segmentTimes, Fs);
    normMax = max(cellfun(@(a) max(filter(1/100*ones(1,100),1,abs(a))), eventData(~badTrial)));
    normMin = min(cellfun(@(a) min(filter(1/100*ones(1,100),1,abs(a))), eventData(~badTrial)));

%    axes('Position',[0 0 1 1],'Visible','off');
%     indsU = regexp(directory(d).name, '_');
%     text(.1,.97, directory(d).name(indsU(end)+1:end-4), 'FontSize', 28, 'FontWeight', 'bold','HorizontalAlignment', 'Center');
%     text(.5,.02,'Times (ms)', 'FontSize', 26, 'FontWeight', 'bold','HorizontalAlignment', 'Center');
%     tx = text(.02,0.8,'Voltage (uV)', 'FontSize', 26, 'FontWeight', 'bold', 'HorizontalAlignment', 'Center');
%     set(tx, 'Rotation', 90);
%     tx = text(.02,0.2,'Voltage (uV)', 'FontSize', 26, 'FontWeight', 'bold', 'HorizontalAlignment', 'Center');
%     set(tx, 'Rotation', 90);
    
    conds = {'Extra Small Sphere', 'Large Sphere','Photocell'};
    for c = 1:length(conds)
        axes(hx(c));

        hold on;
        set(gca, 'FontSize', 18);
        condInd = find(strcmp(arduinoTrials(:,1),conds{c}));
        badCondTrials = badTrial(condInd);
        badCondTrials = false(1,length(condInd));
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
        baselines = cellfun(@(a,b) mean(a(b(1):b(2))), currTrials, currSegs, 'UniformOutput', false);
        %currTrials = cellfun(@(a,b) (a-b)/b, currTrials, baselines, 'UniformOutput', false);
        %currTrials = cellfun(@(a) (a-normMin)/(normMax-normMin), currTrials, 'UniformOutput', false);
        if(sum(badSegTrials)~=length(badSegTrials))
            disp([num2str(sum(badSegTrials==0)), ' trials ignored']);
        end
        averageSegs = mean(reshape(cell2mat(currSegs),numSegs,size(currTrials,2))');
        if(strcmp(conds{c}, 'Photocell'))
            win = {[-.25*Fs, .35*Fs], [-.5*Fs, .5*Fs]};
            alignSegs = {'StartReach','StartHold'};
        elseif(strcmp(conds{c},'Rest'))
            win = {[-.25*Fs, .5*Fs]};
            alignSegs = {'StartReplaceHold'};
        else
            win = window;
            %reachLiftTime = averageSegs(find(strcmp(condSegs{c}, 'StartLift')))-averageSegs(find(strcmp(condSegs{c}, 'StartReach')));
%             reachLiftTime = averageSegs(4) - averageSegs(2);
%             reachLiftTime = round(reachLiftTime/2);
%             win{1}(2) = reachLiftTime+5;
%             win{2}(1) = -(reachLiftTime+5);
        end
        plotStart = 0;
        xTicks = [];
        labels = [];

        for align = 1%1:length(win)
            if(align==1)
                plotStart = 0;
            else
                plotStart = plotStart + (gap*Fs) +win{align-1}(2) + abs(win{align}(1));
            end
            alginSegInd = find(strcmp(condSegs{c},alignSegs{align}));
            alignedEMGSigs = cellfun(@(a,b,c) (a(b(alginSegInd)+win{align}(1):...
                min(b(alginSegInd)+win{align}(2),length(a))))/c,currTrials, currSegs,baselines,'UniformOutput', false);
            yvals = reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))';
            yvals = yvals./cell2mat(baselines)';
            alignedDataY{c,align} = yvals;
            alignedDataX{c,align} = plotStart+win{align}(1):plotStart+win{align}(2);
            % shadedErrorBar(plotStart+win{align}(1):plotStart+win{align}(2),cell2mat(alignedEMGSigs'),...
            %     [std(cell2mat(alignedEMGSigs'),0,1)+mean(cell2mat(alignedEMGSigs'),1,'omitnan');...
            %     max(0,mean(cell2mat(alignedEMGSigs'),1,'omitnan')-std(cell2mat(alignedEMGSigs'),0,1))],...
            %     'lineProps',{['-',colorP(d)],'LineWidth', 3})
            cellfun(@(a) plot(plotStart+win{align}(1):plotStart+win{align}(2),a, 'Color',colorP(d,:)), alignedEMGSigs);
%             
            plot(plotStart+win{align}(1):plotStart+win{align}(2), ...
                nanmean(reshape([alignedEMGSigs{:}],size(alignedEMGSigs{1},2),size(currTrials,2))'),...
                'k', 'LineWidth', 2);
            beforeSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:-1:1)));
            afterSegs = plotStart + cumsum(diff(averageSegs(alginSegInd:end)));
            avgSegs = [beforeSegs plotStart afterSegs];
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=plotStart+win{align}(1) && avgSegs(s)<=plotStart+win{align}(2))
                    plot([avgSegs(s) avgSegs(s)],[0 750],'--')
                end
            end
            segs{c,align} = win{align}./Fs;
            xTicks = [xTicks, plotStart+[win{align}(1) 0 win{align}(2)]];
            labels = [labels, [win{align}(1) 0 win{align}(2)]];
        end
        set(gca, 'XTick', xTicks);
        set(gca, 'XTickLabel', arrayfun(@num2str, (labels .*.5), 'UniformOutput', false));
        set(gca, 'YLim',[0 10])
        title(sprintf('%s\n\t%s%d', conds{c}, 'n=', length(alignedEMGSigs)));
    end
    expind = regexp(directory(d).name,'_');
%     save([dirName, '\',directory(d).name(expind(end)+1:end-4), '_aligned'], 'alignedDataY', 'alignedDataX', 'conds','segs');
%     saveas(gcf,[dirName, 'Figures\', erase(directory(d).name(expind(end)+1:end-4),' ')]);
end