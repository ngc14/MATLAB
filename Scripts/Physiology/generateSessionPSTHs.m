clear all;
close all;

set(0, 'DefaultFigureRenderer', 'painters');
dateStart ='06/17/2019';
dateEnd = '06/18/2019';
monkey = 'Gilligan';
singleOrAll = 'Single';
alignmentPhaseName = 'StartReach';
savePSTHS = 1;

binSize=.01; % bin size in seconds
sigma = .05; % smoothing window in seconds
sigma = round(sigma/binSize);
secondsBeforePSTHAlignmentPoint = 6; % time (in seconds) before alignement point to create PSTH
secondsAfterPSTHAlignmentPoint = 7; % time (in seconds) after alignement point to create PSTH
PSTHDisplayLimits = [-1, 2]; % seconds before and afer reach onset to display the saved PSTHs
bins = -secondsBeforePSTHAlignmentPoint:binSize:secondsAfterPSTHAlignmentPoint; % histogram bin values

segMap = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'},...
    {{'RT', 'R', 'G', 'L','H', 'W', 'RH', 'SR', 'R'},{'RT', 'R', 'G', 'L','H', 'W', 'RH', 'SR', 'R'}, ...
    {'RT','R', 'G', 'H', 'W','RH', 'SR', 'R'},{'H','RH', 'SR', 'R'}});
colorMap = containers.Map({'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'},...
    {'r', 'g', 'b', 'y'});

excelFileNameRead = ['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx'];
headingRow = 1;
if(strcmp(monkey, 'Gilligan'))
    dateFormat = 'mm_dd_yyyy';
else
    dateFormat = 'yyyy_mm_dd';
end
%%
[~,~,rawRef]=xlsread(excelFileNameRead);
rawRef = cell2mat(rawRef(headingRow+1:end,strcmp(rawRef(headingRow,:),'Site #')));
dateArray = [datetime(dateStart):datetime(dateEnd)]';
[~,~,dateToPNValue] = xlsread(['S:\Lab\', monkey,'\', monkey, 'Session.xlsx']);
[~, dateInd] = find(strcmp(dateToPNValue, 'Date'));
[~,domainInd] = find(strcmp(dateToPNValue, 'Domain'));
[~,PNInd] = find(strcmp(dateToPNValue, 'Site #'));
sessionDates = dateToPNValue(:,dateInd);
validInds = find(cellfun(@(a) any(isstrprop(a,'digit')), sessionDates));
inRangeDates = cellfun(@(a) datetime(a)>= dateStart & datetime(a) <=dateEnd, sessionDates(validInds));
sessionDates = sessionDates(validInds(inRangeDates));
sessionDomains = dateToPNValue(validInds(inRangeDates),domainInd);
sessionSites = dateToPNValue(validInds(inRangeDates),PNInd);
clear dateToPNValue;
sessionDates = cellfun(@(a) datetime(a), sessionDates, 'UniformOutput', false);
%%
cd(['S:\Lab\',monkey,'\All Data\']);
for d = 1:length(dateArray)
    currName = [monkey,'_', datestr(dateArray(d), dateFormat)];
    if(exist([currName, '\Physiology\Results\',currName,'_1.mat'], 'file'))
        disp(dateArray(d))
        sessionInd = find(cellfun(@(a) isequal(dateArray(d), a), sessionDates));
        PN = sessionSites{sessionInd};
        savePath = ['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\', singleOrAll, '\', num2str(PN), '\'];
        if(~exist(savePath, 'dir'))
            mkdir(savePath);
        end
        if(strcmp(sessionDomains(sessionInd), 'M1'))
            [allSpikes,alignedTrialSegs,weights,sessionTrials,sessionConds,channel,~,~] = ...
                getSessionInfo(['S:\Lab\',monkey,'\All Data\',currName], singleOrAll, alignmentPhaseName);
            if(~isempty(allSpikes))
                condAbbreviations = cellfun(@(a) arrayfun(@(b) a(b+1), [0,find(isspace(a))]), sessionConds, 'UniformOutput', false);
                graphNames = cellfun(@(a) cellfun(@(c) [a, ' - ', c], ...
                    {'WUnits', 'Trials'}, 'UniformOutput', false), condAbbreviations, 'UniformOutput', false);
                segLabs = values(segMap,sessionConds);
                colors = values(colorMap, sessionConds);
                trialHists = cellfun(@(a) histcounts(a,bins)./binSize, allSpikes,'UniformOutput', false);
                trialHists = cellfun(@(a) conv(a,gausswin(sigma)/sum(gausswin(sigma)),'full'), trialHists, 'UniformOutput', false);
                trialHists = cellfun(@(a) a(1:end-(sigma-2)), trialHists, 'UniformOutput', false);
                figure('units', 'normalize', 'position', [0 0 1 1], 'Visible', 'on');
                for c = 1:length(sessionConds)
                    condInds = cellfun(@(a) strcmp(a, sessionConds{c}), sessionTrials(:,1))';
                    condWeightInd = cellfun(@(a) strcmp(a, sessionConds{c}), sessionConds);
                    condPSTHS=reshape(cell2mat(trialHists(:,condInds)'), [sum(condInds), length(bins),size(trialHists,1)]);
                    p = plot(bins,nanmean(nanmean(bsxfun(@times,condPSTHS,...
                        reshape(weights(:,condWeightInd),1,1,[])),3),1),[colors{c},'-'], 'LineWidth', 1.5);
                    % p = plot(bins,nanmean(weights(:, condWeightInd).*...
                    %  squeeze(nanmean(allPSTHS(condInds,:,:),1))',1),[colors{c},'--'], 'LineWidth', 1.5);
                    hold on;
                    % plot(bins, nanmean(nanmean(allPSTHS(condInds,:,:),3),1),[colors{c},'-'], 'LineWidth', 1.5);
                    plot(bins, nanmean(nanmean(condPSTHS,1),3),[colors{c},'--'], 'LineWidth', 1.5);
                    totalHists = weights(:,condWeightInd).*squeeze(nanmean(condPSTHS,1))';
                    mainLineColor=get(p, 'color');
                    if(all(size(totalHists)>1))
                        totalHists = nanmean(totalHists,1);
                    else
                        totalHists = totalHists';
                    end
                    edgeColor=mainLineColor+(1-mainLineColor)*0.55;
                    uE=mean(totalHists,1)+std(totalHists,1);
                    lE=mean(totalHists,1)-std(totalHists,1);
                    yP=[lE,fliplr(uE)];
                    xP=[bins,fliplr(bins)];
                    xP(isnan(yP))=[];
                    yP(isnan(yP))=[];
                    pat = patch(xP,yP,1);
                    set(pat,'facecolor',colors{c},'edgecolor','none','facealpha',.5);
                end
                pt = [];
                ylimit = get(gca,'YLim');
                ylimit(1) = 0;
                for c = 1:length(sessionConds)
                    pt{c}(1) = plot(NaN, NaN, [colors{c}, '-'], 'LineWidth', 2);
                    pt{c}(2) = plot(NaN, NaN, [colors{c}, '--'], 'LineWidth', 2);
                    condInds = cellfun(@(a) strcmp(a, sessionConds{c}), sessionTrials(:,1))';
                    segBins =  nanmean(cell2mat(alignedTrialSegs(condInds)'),1);
                    segLabsLocs = segBins + [diff(segBins)/2, 0];
                    segLabsLocs(end) = segLabsLocs(end) + (PSTHDisplayLimits(end) + segLabsLocs(end))/2;
                    segLabsLocs = segLabsLocs(segBins> -secondsBeforePSTHAlignmentPoint & segBins<secondsAfterPSTHAlignmentPoint);
                    segBins = segBins(segBins> -secondsBeforePSTHAlignmentPoint & segBins<secondsAfterPSTHAlignmentPoint);
                    segLabsToPlot = segLabs{c}(segBins> -secondsBeforePSTHAlignmentPoint & segBins<secondsAfterPSTHAlignmentPoint);
                    cellfun(@(a,b) text(a,ylimit(mod(c,2)+1), b, 'Color', colors{c}, ...
                        'FontSize', 16,'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'Units', 'data'),...
                        num2cell(segLabsLocs), segLabsToPlot);
                    arrayfun(@(a) line([a a], ylimit, 'Color', colors{c}, 'LineStyle','--'),segBins);
                end
                xlabel('Time in relation to reach onset (s)');
                ylabel('Firing Rate');
                legend([pt{:}], [graphNames{:}],'Location', 'northwest','AutoUpdate','off');
                xlim(PSTHDisplayLimits);
                if(savePSTHS)
                    saveas(gcf,[savePath,'',num2str(PN),'.epsc' ]);
                    saveas(gcf,[savePath,'',num2str(PN),'.fig' ]);
                    pause(0.1);
                end
                close all;
            end
        else
            disp('Not a M1 recording');
        end
    end
end

function updatedispstat(a,b)
persistent count total;
if(a==0)
    count = 1;
    total = b;
else
    dispstat(['Processing session... ',num2str(100*round(count/total,2)),'%']);
    count = count + 1;
end
end