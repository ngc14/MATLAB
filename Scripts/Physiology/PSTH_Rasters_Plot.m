date = '06_13_2019';
monkey = 'Gilligan';
binSize=.01; % bin size in seconds
alignLimits = {[-.5, .15], [-.25, .25],[-.15,.5]};
colors = [[224,144,38]./255; 0 1 0; 0 0 1; 1 0 1];
sigma = 5; % smoothing window in bin sizes
gapWind = 0.10;
saveFolder = ['S:\Lab\ngc14\Working\EMG_UNITS\Stim_Triggered'];
conds = ["Extra Small Sphere", "Large Sphere", "Photocell", "Rest"];
allUnits = cellfun(@(c) cell(1,length(alignLimits)),conds,'UniformOutput',false);
sessionDir = dir(['S:\Lab\',monkey,'\All Data\',monkey,'_',date,'\Physiology\Results']);
sessionDir = sessionDir(~[sessionDir.isdir]);
[~,sortInd] = natsort({sessionDir.name});
sessionDir = sessionDir(sortInd);
MIN_BLOCKS_FOR_UNIT=15;
set(0, 'DefaultFigureRenderer', 'painters');
for f = 1:length(sessionDir)+1
    goodUnitsOnChannel = [];
    if(f<=length(sessionDir))
        m = matfile_m([sessionDir(f).folder, '\', sessionDir(f).name]);
        m = m.sortedSpikeData;
        if(size(m.SegTimes,1)>0)
            load([sessionDir(f).folder,'\',sessionDir(f).name]);
            units = size(sortedSpikeData.SpikeTimes,1);
            conds = sortedSpikeData.Conditions;
            locs = sortedSpikeData.Locations;
            allGoodTrials = ~(cellfun(@(a,b) length(a) <=(b(end)-b(1)) ...
                | length(a)>200*(b(end)-b(1)) | any(isnan(a)), sortedSpikeData.SpikeTimes, sortedSpikeData.SegTimes));
            blockInds = cumsum(mod(1:length(sortedSpikeData.ArduinoData),length(conds))==1);
            unitTrials = num2cell(allGoodTrials.*blockInds,2);
            goodUnitsOnChannel = find(cellfun(@(tc) length(unique(tc))-1, unitTrials)>MIN_BLOCKS_FOR_UNIT);
        end
    else
        goodUnitsOnChannel = 1;
    end
    for u = 1:length(goodUnitsOnChannel)
        figure('units','normalized','outerposition',[0 0 1 1])
        xtickAlign = cell(1,length(conds));
        maxFR = 0;
        ax = {};
        axR = {};
        for c = 1:length(conds)
            currCond = conds{c};
            alignSegs =  alignSegsConds{c};
            currTrialInds = find(cellfun(@(a) contains(a, conds{c}),sortedSpikeData.ArduinoData(:,1)));
            if(f<=length(sessionDir))
                currTrialInds = currTrialInds(cellfun(@(a) ~any(isnan(a)), sortedSpikeData.SpikeTimes(goodUnitsOnChannel(u),currTrialInds)));
                currSpikes = sortedSpikeData.SpikeTimes(goodUnitsOnChannel(u),currTrialInds);
                currSegs = sortedSpikeData.SegTimes(goodUnitsOnChannel(u),currTrialInds);
            else
                currSegs = sortedSpikeData.SegTimes(:,currTrialInds);
            end
            averageSegs = currSegs;
            alignLims = alignLimits;
            averageSegs(cellfun(@length,averageSegs)==1) = [];
            if(c==1||c==2)
                segLabs = {'Rn', 'R', 'G', 'L','H', 'W', 'RH', 'SR', 'R'};
                averageSegs(cellfun(@length,averageSegs)==8) = [];
            elseif (c==length(conds))
                alignSegs = {'GoSignal'};
                alignLims = {[-.5 .5]};
                segLabs = {'H','RH', 'SR', 'R'};
            else
                segLabs = {'Rn','R', 'G', 'H', 'W','RH', 'SR', 'R'};
            end
            numSegs = mode(cellfun(@length,currSegs));
            averageSegs = averageSegs(cellfun(@length,averageSegs)==numSegs);
            plotStart = 0;
            if(length(currTrialInds)>10)
                averageSegs = mean(reshape(cell2mat(averageSegs),numSegs,size(averageSegs,2))',1);
                ax{c} = subplot(2, length(conds), c);
                set(ax{c},'PositionConstraint','outerposition');
                set(ax{c},'FontSize',16);
                title(sortedSpikeData.Conditions{c}, 'FontSize', 18, 'FontWeight', 'bold');
                hold on;
                axR{c} = subplot(2, length(conds), length(conds)+c);
                set(axR{c},'PositionConstraint','outerposition');
                set(axR{c},'FontSize',14)
                hold on;
                %%
                for a = 1:length(alignSegs)
                    secondsBeforeEvent = alignLims{a}(1);
                    secondsAfterEvent = alignLims{a}(end);
                    alignInd = find(strcmp(sortedSpikeData.ConditionSegments{c}, alignSegs(a)));
                    alignTimes = cellfun(@(a) getAlignedTimes(a,alignInd), currSegs, 'UniformOutput', false);
                    if(f<=length(sessionDir))
                        alignedSpikes = cellfun(@minus, currSpikes, alignTimes, 'UniformOutput', false);
                        trialHists = cellfun(@(a) histcounts(a,...
                            secondsBeforeEvent:binSize:secondsAfterEvent)./binSize, alignedSpikes,'UniformOutput', false);
                        trialHistsSmooth = cellfun(@(a) conv(a,gausswin(sigma)/sum(gausswin(sigma)),'same'), trialHists, 'UniformOutput', false);
                        totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
                    else
                        totalHists = cell2mat(cellfun(@(p) mean(p,1,'omitnan'),allUnits{c}{a},'UniformOutput',false)');
                    end
                    % totalHists = siteTrialPSTHS{1}{1}{a};
                    % [~,is] = (min(abs(params.bins-alignLimits{a}(1))))
                    % [~,ie] = (min(abs(params.bins-alignLimits{a}(2))))
                    % totalHists= totalHists(:,is:ie);
                    p(c) = plot(ax{c},plotStart+(secondsBeforeEvent:binSize:...
                        secondsAfterEvent-binSize),nanmean(totalHists,1),...
                        'Color',colors(c,:), 'LineWidth', 2);
                    %errorbar
                    mainLineColor=get(p(c), 'color');
                    edgeColor=mainLineColor+(1-mainLineColor)*0.55;
                    uE=mean(totalHists,1)+std(totalHists,1);
                    lE=mean(totalHists,1)-std(totalHists,1);
                    %Make the patch
                    yP=[lE,fliplr(uE)];
                    xP=[plotStart+(secondsBeforeEvent:binSize:secondsAfterEvent-...
                        binSize),plotStart+fliplr(secondsBeforeEvent:binSize:secondsAfterEvent-binSize)];
                    xP(isnan(yP))=[];
                    yP(isnan(yP))=[];
                    d = patch(ax{c},xP,yP,colors(c,:));
                    set(d,'facecolor',colors(c,:),'edgecolor','none','facealpha',.5)
                    maxFR = max(maxFR,max(yP));
                    PSTHY = 250;
                    %% Rasters
                    yPos = 1;
                    for t = 1:length(alignedSpikes)
                        rasterSpikes = alignedSpikes{t};
                        rasterSpikes = rasterSpikes(rasterSpikes < secondsAfterEvent & rasterSpikes > secondsBeforeEvent);
                        arrayfun(@(a) line(axR{c},[a,a]+plotStart, [yPos-1 yPos+1], 'Color', colors(c,:), 'LineWidth', 1), rasterSpikes);
                        yPos = yPos + 1;
                    end
                    averageSegs = averageSegs - averageSegs(alignInd);
                    beforeSegs = plotStart + cumsum(diff(averageSegs(alignInd:-1:1)));
                    afterSegs = plotStart + cumsum(diff(averageSegs(alignInd:end)));
                    avgSegs = [beforeSegs plotStart afterSegs];
                    plot(ax{c}, [plotStart plotStart],...
                        [1 PSTHY],'Color', [0 0 0], 'LineStyle', '-');
                    plot(axR{c}, [plotStart plotStart],...
                        [1 yPos],'Color',  [0 0 0], 'LineStyle', '-');
                    text(axR{c},plotStart, yPos+1, segLabs{alignInd},...
                        'EdgeColor','none','Color', [0 0 0], 'FontSize',...
                        16, 'FontWeight', 'bold','VerticalAlignment','bottom')
                    text(ax{c},plotStart, 1, segLabs{alignInd},...
                        'EdgeColor','none','Color', [0 0 0], 'FontSize',...
                        16, 'FontWeight', 'bold','VerticalAlignment','bottom');
                    plotSegInds = (avgSegs-plotStart<=secondsBeforeEvent | ...
                        avgSegs-plotStart>=secondsAfterEvent);
                    plotSegInds(alignInd) = true;
                    xtickAlign{c}{a} = [secondsBeforeEvent+plotStart, ...
                        plotStart, secondsAfterEvent+plotStart];
                    plotLabs = segLabs;
                    plotLabs(plotSegInds) = [];
                    avgSegs(plotSegInds) = [];
                    for s = 1:length(avgSegs)
                        plot(ax{c}, [avgSegs(s) avgSegs(s)],[1 PSTHY],'Color', [.5 .5 .5], 'LineStyle', '--');
                        plot(axR{c}, [avgSegs(s) avgSegs(s)],[1 yPos],'Color',  [.5 .5 .5], 'LineStyle', '--');
                    end
                    allUnits{c}{a}{end+1} = totalHists;
                    plotStart = plotStart + abs(secondsBeforeEvent) + gapWind + ...
                        secondsAfterEvent;
                    if(c==1)
                        ylabel(ax{c}, 'Firing Rate (Hz)', 'FontSize', 16, 'FontWeight', 'bold');
                        ylabel(axR{c}, 'Trials', 'FontSize', 16, 'FontWeight', 'bold');
                    end
                end
                xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
                xticks(ax{c},cell2mat(xtickAlign{c}));
                xticks(axR{c},cell2mat(xtickAlign{c}));
                xlabs = cellfun(@(a) [num2str(a(1)) "" num2str(a(end))],...
                    alignLims, 'UniformOutput',false);
                xticklabels(ax{c},[xlabs{:}]);
                xticklabels(axR{c},[xlabs{:}]);
                xlim(ax{c},[min(cell2mat(xtickAlign{c})), max(cell2mat(xtickAlign{c}))]);
                xlim(axR{c},[min(cell2mat(xtickAlign{c})), max(cell2mat(xtickAlign{c}))]);
            end
        end
        if(~isempty(ax))
            ylim([ax{:}],[-.5,maxFR]);
        end
        if(~isempty(axR))
            ylim([axR{:}],[0,yPos+1]);
        end
        saveDirDate = [saveFolder,'\',date];
        if(~exist(saveDirDate,'dir'))
            mkdir(saveDirDate);
        end
        saveas(gcf,[saveDirDate,'\Channel',num2str(f),'_Unit',num2str(goodUnitsOnChannel(u))], 'fig');
        saveas(gcf,[saveDirDate,'\Channel',num2str(f),'_Unit',num2str(goodUnitsOnChannel(u))], 'png')
        saveas(gcf,[saveDirDate,'\Channel',num2str(f),'_Unit',num2str(goodUnitsOnChannel(u))], 'epsc')
        close all;
    end
end



function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end