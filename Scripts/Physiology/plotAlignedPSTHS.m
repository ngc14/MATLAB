function plotAlignedPSTHS(bins,alignLimits,PSTHunits,segs,saveDir)
zeroBinInd = find(bins==0);
alignmentGap = .25;
binSize = mode(diff(bins));
alignmentGap = alignmentGap/binSize;

for u = 1:size(PSTHunits{1},2)
    
    PSTH = cellfun(@(a,w) a(:,(fix(w(1)/binSize)+zeroBinInd):...
        fix((w(end)/binSize)+zeroBinInd),:),PSTHunits, alignLimits, 'UniformOutput', false);
    allSegs = segs;
    figure(); hold on;
    plotStart = 0;
    plotted = false(1,size(allSegs,2));
    xAlignTicks = {};
    for a = 1:length(PSTH)
        xAlignTicks{a} = plotStart+[1:size(PSTH{a},2)];
        plot(xAlignTicks{a},nanmean(squeeze(PSTH{a}(u,:,:)),2), 'LineWidth',2,'Color', 'b')
        ue = nanmean(squeeze(PSTH{a}(u,:,:)),2) + (nanstd(squeeze(PSTH{a}(u,:,:)),0,2)/sqrt(size(PSTH{a},3)));
        le = nanmean(squeeze(PSTH{a}(u,:,:)),2) - (nanstd(squeeze(PSTH{a}(u,:,:)),0,2)/sqrt(size(PSTH{a},3)));
        yP = [le',fliplr(ue')];
        xP = [xAlignTicks{a},fliplr(xAlignTicks{a})];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        d = patch(xP,yP,1);
        set(d,'edgecolor','none','facealpha',.5,'facecolor','r');
        avgSegs = nanmean(allSegs{a}{u},1);
        for s = 1:length(avgSegs)
            if(avgSegs(s)>=alignLimits{a}(1) && ...
                    avgSegs(s)<=alignLimits{a}(end) && ...
                    (~plotted(s) ||avgSegs(s)==plotStart))
                if(avgSegs(s)==0)
                    plotColor = 'k';
                else
                    plotColor = 'r';
                end
                %plotted(s) = true;
                pSeg = find(isalmost(alignLimits{a}(1):binSize:...
                    alignLimits{a}(end),avgSegs(s),binSize/1.99),1);
                plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 FRLim(end)],'Color',plotColor,'LineStyle','--')
            end
        end
        if(a==length(PSTH))
            allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs([pd(1):binSize:pd(end)])...
                ==min(abs([pd(1):binSize:pd(end)]))), ta(end)], xAlignTicks, ...
                alignLimits, 'UniformOutput', false);
            xticks([allXTicks{:}]);
            allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                "0"; num2str(pd(end),'%.2f')],alignLimits, 'UniformOutput', false);
            xticklabels([allLabels{:}]);
        end
        plotStart = plotStart + size(PSTH{a},2) + alignmentGap;
    end
    saveFigures(gcf,saveDir, num2str(u),[]);
    if(u==size(PSTH{1},1))
        
        figure(); hold on;
        plotStart = 0;
        plotted = false(1,9);
        xAlignTicks = {};
        for a = 1:length(PSTH)
            xAlignTicks{a} = plotStart+[1:size(PSTH{a},2)];
            plot(xAlignTicks{a},nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1), 'LineWidth',2,'Color', 'b')
            ue = squeeze(nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1))+  squeeze(nanmean(nanstd(squeeze(PSTH{a}(:,:,:)),0,3)/sqrt(size(PSTH{a},3)),1));
            le = squeeze(nanmean(nanmean(squeeze(PSTH{a}(:,:,:)),3),1)) - squeeze(nanmean(nanstd(squeeze(PSTH{a}(:,:,:)),0,3)/sqrt(size(PSTH{a},3)),1));
            yP = [le,fliplr(ue)];
            xP = [xAlignTicks{a},fliplr(xAlignTicks{a})];
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            d = patch(xP,yP,1);
            set(d,'edgecolor','none','facealpha',.5,'facecolor','r');
            avgSegs = nanmean(allSegs{a}{u},1);
            for s = 1:length(avgSegs)
                if(avgSegs(s)>=alignLimits{a}(1) && ...
                        avgSegs(s)<=alignLimits{a}(end) && ...
                        (~plotted(s) ||avgSegs(s)==plotStart))
                    if(avgSegs(s)==0)
                        plotColor = 'k';
                    else
                        plotColor = 'r';
                    end
                    %plotted(s) = true;
                    pSeg = find(isalmost(alignLimits{a}(1):binSize:...
                        alignLimits{a}(end),avgSegs(s),binSize/1.99),1);
                    plot([xAlignTicks{a}(pSeg) xAlignTicks{a}(pSeg)],[0 FRLim(end)],'Color',plotColor,'LineStyle','--')
                end
            end
            if(a==length(PSTH))
                allXTicks = cellfun(@(ta,pd) [ta(1),ta(1)+ find(abs([pd(1):binSize:pd(end)])...
                    ==min(abs([pd(1):binSize:pd(end)]))), ta(end)], xAlignTicks, ...
                    alignLimits, 'UniformOutput', false);
                xticks([allXTicks{:}]);
                allLabels = cellfun(@(pd)[num2str(pd(1),'%.2f'); ...
                    "0"; num2str(pd(end),'%.2f')],alignLimits, 'UniformOutput', false);
                xticklabels([allLabels{:}]);
            end
            plotStart = plotStart + size(PSTH{a},2) + alignmentGap;
        end
    end
    saveFigures(gcf,saveDir, 'Session',[]);
    
end
end