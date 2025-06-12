monkey = 'Gilligan';
muscles = dir(['S:\Lab\',monkey,'\Sorted EMG Data\Averages_Raw\*.mat']);
conditions = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};
gap = .15;
muscleColors = struct('Deltoid',[1 0 0], 'Biceps',[.7 .4 0], 'Triceps',...
    [1 .64 0],'WristExtensor',[0 1 0], 'WristFlexor',[0 .5 0], 'DigitExtensor',...
    [0 0 1],'DigitFlexor',[0 0 .5]);
saveDirPath = ['S:\Lab\ngc14\Figures\Physiology\Results\EMG\'];

for c = 1:length(conditions)
    cond = conditions{c};
    figure('Units','normalized', 'Position', [0 0 1 1]);
    set(gcf,'Renderer','painters')
    hold on;
    muscleSegs = {};
    allMuscAligned = {};
    for m = 1:length(muscles)
        loadFile = load([muscles(m).folder,'\',muscles(m).name]);
        conditionInd = find(strcmp(loadFile.Conditions,cond));
        alignedTrials = loadFile.alignedSig{conditionInd};
        trialSegs = loadFile.segs{conditionInd};
        numSegs = loadFile.ConditionSegs{conditionInd};
        alignWinds = loadFile.alignWindows;
        alignSegments = loadFile.alignments;
        Fs = loadFile.Fs;
        plotted = false(1,length(numSegs));
        alignWinds = cellfun(@(a) a.*Fs, alignWinds, 'UniformOutput', false);
        muscleName = muscles(m).name(1:end-4);
        %clear loadFile;
        muscleSegs{m} = reshape(cell2mat(trialSegs),length(numSegs),size(trialSegs,2))';
        for a = 1:size(alignedTrials,2)
            alignSegInd = find(strcmp(alignSegments{a},numSegs));
            if(strcmp(cond,'Photocell') & isempty(alignSegInd))
                alignSegInd = find(strcmp('StartHold',numSegs));
            end
            currAlignment = alignedTrials(:,a);
           [currAlignment{cellfun(@(n) all(isnan(n)) | length(n)~=diff(alignWinds{a})+1,...
               currAlignment)}] = deal(NaN(1,diff(alignWinds{a})+1));
            currAlignment = cell2mat(currAlignment);
            allMuscAligned{m,a} = currAlignment;
            meanAligned = nanmean(currAlignment,1);
            if(strcmp(cond,'Rest'))
                summedXLim = cell2mat(alignWinds');
                plotXLim = [min(summedXLim(:,1)),sum(summedXLim(:,end))];
                plot([0 0],segSpan,'Color',plotColor,'LineStyle','--');
            else
                if(a==1)
                    plotStart = 0;
                else
                    plotStart = plotStart + (gap*Fs) + alignWinds{a-1}(end) + abs(alignWinds{a}(1));
                end
                plotXLim = plotStart+alignWinds{a}(1):plotStart+alignWinds{a}(2);
            end
            %cellfun(@(a) plot(plotXLim,a), alignedTrials(:,a));
            %plot(plotXLim,meanAligned,'Color',musclesColors.(muscleNames{m}), 'LineWidth', 3);
            %[prctile(currAlignment,84);prctile(currAlignment,16)]
            shadedErrorBar(plotXLim,meanAligned,nanstd(currAlignment,[],1)./sqrt(size(currAlignment,1)),'patchSaturation', 0.03,'transparent',true,...
                'lineprops',{'Color',muscleColors.(muscleName(~isspace(muscleName))),'LineWidth',3});
            % plot seg lines
            if(m==length(muscles))
                if(a==1)
                    g =cellfun(@(cl) plot(NaN, NaN, 'Color', muscleColors.(cl),'LineWidth', 3),fieldnames(muscleColors));
                    legend(g,fieldnames(muscleColors),'AutoUpdate','off');
                    segSpan = get(gca,'YLim');
                end
                allAligned = vertcat(allMuscAligned{:,a});
                shadedErrorBar(plotXLim,nanmean(allAligned,1),...
                    nanstd(allAligned,[],1)./sqrt(size(allAligned,1))...
                    ,'patchSaturation', 0.03,'transparent',true,...
                'lineprops',{'Color','k','LineWidth',3});
            
                beforeSegs = cellfun(@(t) t(1:alignSegInd)-t(alignSegInd), trialSegs,'UniformOutput', false);
                afterSegs = cellfun(@(t) (t(alignSegInd+1:end) - t(alignSegInd)), trialSegs, 'UniformOutput', false);
                averageSegs = cellfun(@(b,a) plotStart + [b, a], beforeSegs, afterSegs, 'UniformOutput', false);
                avgSegs = (nanmean(cell2mat(averageSegs'),1));
                for s = 1:length(avgSegs)
                    if(avgSegs(s)>=plotXLim(1) && ...
                            avgSegs(s)<=plotXLim(end) && ...
                            (~plotted(s) ||avgSegs(s)==plotStart))
                        if(avgSegs(s)==plotStart)
                            plotColor = 'b';
                        else
                            plotColor = 'r';
                        end
                        %plotted(s) = true;
                        plot([avgSegs(s) avgSegs(s)],segSpan,'Color',plotColor,'LineStyle','--')
                    end
                end
                if(a>1)
                    set(gca, 'XTick', sort([get(gca, 'XTick'),plotXLim(1), plotStart, plotXLim(end)]));
                    set(gca, 'XTickLabels', [string(get(gca,'XTickLabel')'),num2str(alignWinds{a}(1)/Fs), "",num2str(alignWinds{a}(end)/Fs)]);
                else
                    xticks([plotXLim(1), plotStart, plotXLim(end)]);
                    xticklabels({num2str(alignWinds{a}(1)/Fs), '',num2str(alignWinds{a}(end)/Fs)});
                end
                if(a==3)
                    allticks = get(gca,'XTick');
                    xlim([allticks(1), allticks(end)]);
                end
            end
        end
    end
    title(cond);
    savefig(gcf,[saveDirPath,monkey,'_',cond,'_All Muscles'],'compact');
    saveas(gcf,[saveDirPath,monkey,'_',cond,'_All Muscles','.png']);
    saveas(gcf,[saveDirPath,monkey,'_',cond,'_All Muscles','.epsc']);
    close all;
end