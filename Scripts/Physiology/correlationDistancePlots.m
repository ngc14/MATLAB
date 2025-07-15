function correlationDistancePlots(saveFig, monkey, singleOrAll, unitOrSession, ...
    masterPSTH, masterBinLims, masterJoints, masterSites, activityInd, masterChannels)
pixelToMM = 0.018;
laminarToMM = 1;
distBins = [0, 0:5];
lamBins = 0:31;
corrLim = [-1, 1];
distType = '';
corrType = '';
fileNameEnd = '_Max_All';
if(~isempty(distType))
    fileNameEnd = [fileNameEnd, '_', distType];
end
if(~isempty(corrType))
    fileNameEnd = [fileNameEnd, '_', corrType];
    corrLim(1) = 0;
end
prefix = {'Ina', 'A'};

jointName = sort(unique(masterJoints));
jointName = jointName([4 1 3 2]);
newSessionInd = find(diff([masterChannels{:}])<0);
activityInd = logical(activityInd);
activityVec = arrayfun(@(a) strcat(prefix{a+1}, 'ctive'),...
    activityInd, 'UniformOutput', false);
activityMat = [categorical(activityVec).*categorical(activityVec')]';

distMatrix = [];
distX = [];
distY = [];
corrMatrix = [];
jointComps = {};
parfor p = 1:length(masterPSTH)
    [distRow, corrRow, distXRow, distYRow] = deal(NaN(length(masterPSTH),1));
    jointRow = {};
    taskPSTH1 = masterPSTH{p}(masterBinLims{p}(1):masterBinLims{p}(2));
    for pp = 1:length(masterPSTH)
        condName = [masterJoints(p); masterJoints(pp)];
        condName = [condName{:}];
        jointRow{pp} = condName;
        taskPSTH2 = masterPSTH{pp}(masterBinLims{pp}(1):masterBinLims{pp}(2));
        if(length(taskPSTH2)>length(taskPSTH1))
            shortNaNTask = NaN(1,length(taskPSTH2));
            longNaNTask = taskPSTH2;
            startInd = (masterBinLims{pp}(1)-masterBinLims{p}(1))+1;
            shortNaNTask(startInd:startInd+length(taskPSTH1)-1) = taskPSTH1;
        elseif(length(taskPSTH2)<length(taskPSTH1))
            shortNaNTask = NaN(1,length(taskPSTH1));
            longNaNTask = taskPSTH1;
            startInd = (masterBinLims{p}(1)-masterBinLims{pp}(1))+1;
            shortNaNTask(startInd:startInd+length(taskPSTH2)-1) = taskPSTH2;
        else
            shortNaNTask = taskPSTH1;
            longNaNTask = taskPSTH2;
        end
        rMat = corrcoef(shortNaNTask,longNaNTask,'Rows', 'pairwise');
        distRow(pp) = pdist([masterSites{p}; masterSites{pp}]).*pixelToMM;
        distXRow(pp) = pixelToMM*abs(masterSites{p}(1) - masterSites{pp}(1));
        distYRow(pp) = pixelToMM*abs(masterSites{p}(2) - masterSites{pp}(2));
        if(p==pp)
            corrRow(pp) = NaN;
        else
            if(strcmp(corrType, 'Clip'))
                corrRow(pp) = max(0,rMat(2,1));
            elseif(strcmp(corrType, 'Abs'))
                corrRow(pp) = abs(rMat(2,1));
            else
                corrRow(pp) = rMat(2,1);
            end
        end
    end
    distMatrix(p,:) = distRow;
    corrMatrix(p,:) = corrRow;
    distX(p,:) = distXRow;
    distY(p,:) = distYRow;
    jointComps(p,:) = jointRow;
end
limitingDist =  distMatrix(cellfun(@(a) strcmp(a, 'HandHand'), jointComps));
limitingDist = nanmax(limitingDist);
corrMatrix(distMatrix > limitingDist) = NaN;
if(~ismissing(distType))
    if(strcmp(distType, 'Weighted'))
    normDistance = distMatrix./max(distMatrix(:));
    end
end
if(~ismissing(corrType))
    if(strcmp(corrType, 'Clip')| strcmp(corrType,'Abs'))
        corrMatrix = 1-((1 - corrMatrix + normDistance)./2);
    else
        corrMatrix = (corrMatrix + 1)./2;
        corrMatrix = -2.*((1 - corrMatrix + normDistance)./2)+1;
    end
end
fInd = cellfun(@(a) strcmp(a,'ArmArm'), jointComps) | ...
    cellfun(@(a) strcmp(a,'ArmHand'), jointComps)  | cellfun(@(a) strcmp(a,'HandArm'), jointComps)...
    | cellfun(@(a) strcmp(a,'HandHand'), jointComps)
fInd = tril(fInd);
corrVec = corrMatrix(fInd);
corrVec(isnan(corrVec)) = 1;
distVec = distMatrix(fInd);
b = regress(corrVec,[distVec, ones(size(distVec))]);
y_regressed = corrVec - [distVec,ones(size(distVec))]*b;

if(strcmp(unitOrSession, 'Unit'))
    channelDist = [];
    corrLam = [];
    activeLamSites = [];
    repLam = {};
    startChannelInd = 1;
    for s = 1:length(newSessionInd)
        currUnits = startChannelInd:newSessionInd(s);
        startChannelInd = newSessionInd(s)+1;
        currChannels = [masterChannels{currUnits}];
        if(length(currChannels)>1)
            [rCh, cCh] =  find(tril(ones(length(currUnits)),-1));
            currChannelDist = laminarToMM*(squareform(abs(bsxfun(@minus,currChannels',currChannels))));
            channelDist(end+1:end+length(rCh)) = currChannelDist;
            corrLam(end+1:end+length(rCh)) =arrayfun(@(a,b) ...
                corrMatrix(a,b), rCh, cCh);
            repLam(end+1:end+length(rCh)) = repmat(sort(unique([...
                masterJoints(currUnits(rCh)), masterJoints(currUnits(cCh))])),...
                1,length(rCh));
            activeLamSites(end+1:end+length(rCh)) = repmat(mode(activityInd(currUnits(rCh))),...
                1,length(rCh));
        else
            channelDist(end+1) = 0;
            corrLam(end+1) = NaN;
            repLam(end+1) = masterJoints(currUnits);
            activeLamSites(end+1) = activityInd(currUnits);
        end
    end
    channelDist(channelDist>=lamBins(end)) = lamBins(end);
    laminarVals = discretize(channelDist, lamBins, 'categorical', 'IncludedEdge', 'right');
    activeLamSites = categorical(arrayfun(@(a) strcat(prefix{a+1}, 'ctive'),...
        activeLamSites, 'UniformOutput', false));
    if(strcmp(distType, 'Weighted'))
        normChan = channelDist./max(channelDist(:));
        if(strcmp(corrType, 'Clip')| strcmp(corrType,'Abs'))
            corrLam = 1-((1 - corrLam + normChan)./2);
        else
            corrLam = (corrLam + 1)./2;
            corrLam = -2.*((1 - corrLam + normChan)./2)+1;
        end
    end
    clear channelDist;
end

sortInds = reshape(1:numel(distMatrix), size(distMatrix));
[~, sortedDistInds] = sortrows(sortrows(distMatrix, 'ascend')');
sortedDistInds = sortInds(sortedDistInds,sortedDistInds);
[~, jointInds] = ismember(masterJoints,jointName);
jointMat = [categorical(jointInds) .* categorical(jointInds')]';
remapJoints = jointMat(sortedDistInds);
[~, sortedJointInds] = sortrows(sortrows(remapJoints, 'ascend')');
sortedInds = sortedDistInds(sortedJointInds, sortedJointInds);
organizedMat = corrMatrix(sortedInds);
activityMat = activityMat(sortedInds);
jointComps = jointComps(sortedInds);
activityCat = mergecats(activityMat, {'Inactive Active', 'Active Inactive'});

[distMatrix(distMatrix>=distBins(end)), distX(distX>=distBins(end)), ...
    distY(distY>=distBins(end))] = deal(distBins(end));
distVals = discretize(distMatrix, distBins, 'categorical', 'IncludedEdge', 'right');
clear distMatrix
xVals = discretize(distX, distBins, 'categorical', 'IncludedEdge', 'right');
clear distX
yVals = discretize(distY, distBins, 'categorical', 'IncludedEdge', 'right');
clear distY

jointCompsOrdered = cellstr(categorical(jointName).*categorical(jointName'));
jointCompsOrdered = cellfun(@(a) a(~isspace(a)), jointCompsOrdered, 'UniformOutput', false);
jointCompsOrdered = jointCompsOrdered(:);
for a = 1:2
    if(a==1)
        colorBox = 'b';
        jointOrderCond = jointCompsOrdered;
        jointValsCond = jointComps;
        if(strcmp(unitOrSession, 'Unit'))
            laminarValsCond = string(laminarVals(:));
        end
    elseif(a==2)
        colorBox = 'brg';
        if(strcmp(unitOrSession, 'Unit'))
            laminarValsCond = string(laminarVals.*activeLamSites);
        end
        distVals = distVals.*activityCat;
        xVals = xVals.*activityCat;
        yVals = yVals.*activityCat;
        jointValsCond = jointComps.*activityCat;
        jointOrderCond = cellfun(@(a) arrayfun(@(n) strcat(a,{' '}, n), categories(activityCat)), jointCompsOrdered, 'UniformOutput', false);
        jointOrderCond = [jointOrderCond{:}];
        jointOrderCond = jointOrderCond(:);
        fileNameEnd = [fileNameEnd, '_Acitivity'];
    end
    distValsCond = string(distVals(:));
    xValsCond = string(xVals(:));
    yValsCond = string(yVals(:));
    jointValsCond = string(jointValsCond(:));
    
    boxFig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
    set(0, 'CurrentFigure', boxFig);
    
    subplot(1,3,1);
    boxplot(corrMatrix(:), distValsCond(:),'GroupOrder', sortCategories(distValsCond),  'Colors', colorBox);
    hold on;
    title('Correlation by total euclidean distance');
    xtickangle(45);
    ylim(corrLim);
    
    subplot(1,3,2);
    boxplot(corrMatrix(:), xValsCond(:),'GroupOrder',sortCategories(xValsCond),'Colors', colorBox);
    hold on;
    title('Correlation by rostral/caudal distance');
    xtickangle(45);
    ylim(corrLim);
    
    subplot(1,3,3);
    boxplot(corrMatrix(:), yValsCond(:),'GroupOrder',sortCategories(yValsCond), 'Colors', colorBox);
    hold on;
    title('Correlation by medial/lateral distance');
    xtickangle(45);
    ylim(corrLim);
    
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'k');
    
    if(strcmp(unitOrSession, 'Unit'))
        lamFig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
        set(0, 'CurrentFigure', lamFig);
        subplot(2,1,1);
        if(a==1)
            jointOrder = jointName;
        elseif(a==2)
            jointOrder =cellfun(@(a) {strcat(a, ' Active'); strcat(a, ' Inactive')}, jointName, 'UniformOutput', false);
            jointOrder = [jointOrder{:}];
            jointOrder = jointOrder(:);
            repLam = repLam.*activeLamSites;
        end
        boxplot(corrLam(:), string(repLam(:)),'GroupOrder', jointOrder,  'Colors', colorBox(1:end-1));
        hold on;
        title('Correlation by representation');
        ylim(corrLim)
        subplot(2,1,2);
        boxplot(corrLam(:), laminarValsCond(:),'GroupOrder', sortCategories(laminarValsCond),  'Colors', colorBox(1:end-1));
        hold on;
        title('Correlation by total laminar distance');
        ylim(corrLim)
        xtickangle(45);
        xlamLab = get(gca, 'XTickLabel');
        xlamLab{1} = '(0, 1]';
        set(gca, 'XTickLabel', xlamLab);
        lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
        set(lines, 'Color', 'k');
    end
    
    corrFig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
    set(0, 'CurrentFigure', corrFig);
    allJointCounts = {};
    currJointInds = {};
    if(a==1)
        subplot(2,4,[1,2,5,6]);
        imagesc(organizedMat);
        hold on;
        colormap('jet');
        caxis(corrLim);
        xlim([1,size(organizedMat,2)]);
        ylim([1,size(organizedMat,1)]);
        currJointInds = {reshape(jointValsCond',size(jointComps))};
        [jointCountsX, jointCountsY] = deal({cumsum(cellfun(@(c) sqrt(sum(sum(strcmp(jointComps,[c,c])))), jointName))});
        subplotRightOrder = [3,4,7,8];
        addedConditional = {true(1, length(masterPSTH))};
    elseif(a==2)
        categoryConds = categories(activityMat);
        for c = 1:length(categoryConds)
            currCat = categoryConds{c};
            mismatchCond = ~strcmp(currCat(1:find(isspace(currCat))-1),currCat(find(isspace(currCat))+1:end));
            if(mismatchCond)
                subplot(3,4,[5,6]+((c-2).*[2,2]));
            else
                subplot(3,4,[1,2]+(double(c==length(categoryConds)).*[8,8]));
            end
            hold on;
            set(gca,'YDir', 'reverse');
            title(string(categoryConds{c}), 'FontSize', 20, 'Color', colorBox(c-double((c>length(categoryConds)/2))));
            condInd = cellfun(@(a) strcmp(a, currCat), cellstr(activityMat));
            currMat = organizedMat(condInd);
            allJointCond = string(jointComps(condInd));
            if(mismatchCond)
                if(startsWith(categoryConds{c}, 'Active'))
                    currMat = reshape(currMat,sum(activityInd), sum(~activityInd));
                    allJointCond = reshape(allJointCond,sum(activityInd), sum(~activityInd));
                else
                    currMat = reshape(currMat,sum(~activityInd), sum(activityInd));
                    allJointCond = reshape(allJointCond,sum(~activityInd), sum(activityInd));
                end
            else
                currMat = reshape(currMat,sqrt(length(currMat)), sqrt(length(currMat)));
                allJointCond = reshape(allJointCond,sqrt(length(allJointCond)), sqrt(length(allJointCond)));
            end
            imagesc(currMat);
            colormap('jet');
            caxis(corrLim);
            xlim([1,size(currMat,2)]);
            ylim([1,size(currMat,1)]);
            
            currJointInds(end+1) = {allJointCond};
            if(mismatchCond)
                allJointCounts(end+1) = {[NaN]};
            else
                allJointCounts(end+1) = {cumsum(cellfun(@(j) sqrt(sum(sum(strcmp(...
                    allJointCond,[j,j])))), jointName))};
            end
        end
        subplotRightOrder = [3,4,7,8];
        jointCountsCond = cellfun(@max, allJointCounts, 'UniformOutput', false);
        jointCountsCond = cellfun(@(a,b) arrayfun(@(s) find(s==[jointCountsCond{:}]),size(a)), currJointInds, 'UniformOutput', false);
        
        addedConditional = {activityInd; true(1, length(masterPSTH)); ...
            true(1, length(masterPSTH)); ~activityInd};
        jointCountsX = cellfun(@(xi) allJointCounts(xi(2)), jointCountsCond)';
        jointCountsY = cellfun(@(yi) allJointCounts(yi(1)), jointCountsCond)';
        currJointInds = currJointInds';
    end
    ax = get(gcf, 'Children');
    axsInds = arrayfun(@(c) arrayfun(@(x) strcmp(string(get([x{:}], 'Type')), 'image'),...
        arrayfun(@(h) h.Children, c, 'UniformOutput', false), 'UniformOutput', false),...
        ax, 'Uniformoutput', false);
    ax = ax(cellfun(@(a) any([a{:}]), axsInds));
    [~,axsOrder] = sort(arrayfun(@(a) a.Title.String, ax, 'UniformOutput', false));
    ax = ax(axsOrder);
    axYPos = cellfun(@(y) (-2-y(1))/abs(y(2)-y(1)), arrayfun(@(ay) get(ay,'YLim'), ax, 'UniformOutput', false));
    axXPos = cellfun(@(x) (-2-x(1))/abs(x(2)-x(1)), arrayfun(@(axx) get(axx,'YLim'), ax, 'UniformOutput', false));
    
    
    arrayfun(@(cx, v) arrayfun(@(x) line(cx, [x,x], [get(cx,'YLim')], 'Color', 'k',...
        'LineWidth', 1.5),[v{:}]), ax, jointCountsX);
    arrayfun(@(cx, v) arrayfun(@(y) line(cx, [get(cx, 'XLim')],[y,y], 'Color', 'k',...
        'LineWidth', 1.5),[v{:}]), ax, jointCountsY);
    jointCountsX = cellfun(@(jx) [1.5 zeros(1,length(jx)-1)]+...
        (([-jx(1) jx(1:end-1)] + jx)./2), jointCountsX, 'UniformOutput', false);
    jointCountsY = cellfun(@(jy) ([-jy(1) jy(1:end-1)] + jy)./2, jointCountsY, 'UniformOutput', false);
    for j = 1:length(jointName)
        textCondX = cellfun(@(x) x(j), jointCountsX, 'UniformOutput', true);
        arrayfun(@(a,x,p) text(a, x, p, jointName{j}, 'FontSize', ...
            14, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom',...
            'HorizontalAlignment', 'center'),ax, textCondX, axYPos);
        for jj = 1:length(jointName)
            jointUnitCond1 = contains(masterJoints, jointName{j});
            jointUnitCond2 = contains(masterJoints, jointName{jj});
            textCondY = cellfun(@(y) y(jj), jointCountsY, 'UniformOutput', true);
            avgVal = cellfun(@(a) corrcoef(nanmean(reshape([masterPSTH{...
                jointUnitCond1 & a}], sum(jointUnitCond1 & a), length(...
                masterPSTH{1})),1), nanmean(reshape([masterPSTH{...
                jointUnitCond2 & a}], sum(jointUnitCond2 & a), length(...
                masterPSTH{1})),1)), addedConditional, 'UniformOutput', false);
            avgVal = cellfun(@(a) a(2), avgVal, 'UniformOutput', true);
            currMats = arrayfun(@(a) a.Children(arrayfun(@(c) strcmp(string(get([c(:)], 'Type')), 'image'),a.Children)),ax);
            vecVals = arrayfun(@(cm,in) cm.CData(strcmp(in{1}, [jointName{j},jointName{jj}])),...
                currMats,currJointInds, 'UniformOutput', false);
            medianVals = cellfun(@(v) nanmedian(v(:)), vecVals);
            meanVals = cellfun(@(v) nanmean(v(:)), vecVals);
            if(j==1)
                arrayfun(@(a,x,p) text(a,p,x, jointName{jj}, 'FontSize', ...
                    14, 'FontWeight', 'bold', 'HorizontalAlignment', 'right'),...
                    ax,textCondY,axXPos);
            end
            if(j==1)
                hAlignment = 'left';
            else
                hAlignment = 'center';
            end
            if(mod(j,2) && jj~=1)
                vAlignment = 'bottom';
            else
                vAlignment = 'top';
            end
            if(j<=jj)
                arrayfun(@(a,x,y,m,u,v) text(a,x,y,['M: ',num2str(m,2),...
                    ', MU: ',num2str(u,2),', PSTH:',num2str(v,2)],'FontSize', 12,'FontWeight','bold',...
                    'Color', 'Black','HorizontalAlignment',hAlignment,...
                    'VerticalAlignment',vAlignment, 'BackgroundColor','white'),...
                    ax, textCondX, textCondY,medianVals,meanVals,avgVal);
            end
        end
    end
    for sj = 1:length(jointName)
        subplot(2,4,subplotRightOrder(sj));
        sjInds = cellfun(@(a) startsWith(a, jointName{sj}), jointValsCond);
        boxplot(organizedMat(sjInds), jointValsCond(sjInds), 'GroupOrder', ...
            jointOrderCond(cellfun(@(a) startsWith(a, jointName{sj}), ...
            jointOrderCond) & ismember(jointOrderCond,...
            unique(jointValsCond(sjInds)))), 'Colors', colorBox);
        hold on;
        title([jointName{sj}, ' correlations by representation']);
        xtickangle(45);
        ylim(corrLim);
    end
    
    if(saveFig)
        saveas(corrFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
            singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,...
            fileNameEnd,'.fig']);
        saveas(corrFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
            singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,...
            fileNameEnd,'.png']);
        saveas(corrFig, ['G:\My Drive\Xcorr', fileNameEnd,'.png']);
        saveas(boxFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
            singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,'_Distance',...
            fileNameEnd,'.fig']);
        saveas(boxFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
            singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,'_Distance',...
            fileNameEnd,'.png']);
        saveas(boxFig, ['G:\My Drive\Dist', fileNameEnd, '.png'])
        if(strcmp(unitOrSession, 'Unit'))
            saveas(lamFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
                singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,'_Laminar',...
                fileNameEnd,'.fig']);
            saveas(lamFig,['S:\Lab\', monkey, '\Mapping\Encoding Maps\PSTHs\',...
                singleOrAll,'\FRs\Quantification\XCorr_',unitOrSession,'_Laminar',...
                fileNameEnd,'.png']);
            saveas(lamFig, ['G:\My Drive\lam', fileNameEnd, '.png'])
        end
    end
end
%     plotInds = reshape(1:length(axesHandles),length(unitRow), length(phaseCol)+1)';
%     plotTotals = plotInds(:,end);
%     plotInds = plotInds(:,1:end-1);
%     linkaxes([axesHandles{plotInds(:)}],'y');
%     linkaxes([axesHandles{plotTotals(1:end-1)}],'y');
%     yVals = cellfun(@(a) arrayfun(@(b) b(1).YData + b(1).YPositiveDelta, ...
%         findobj(a,'type', 'errorbar'), 'UniformOutput', false), axesHandles(plotInds(:)));
%     yVals = [yVals{:}];
%     set(axesHandles{plotInds(1)}, 'YLim', [0 nanmean(yVals) + nanstd(yVals)/3]);
% if(a==1)
%

%     textCondX = (jointUnitsCount(j) +  jointUnitsCount(j+1))/2;
%     arrayfun(@(a,b) text(a, b, -2, jointName{j}, 'FontSize', ...
%         14, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'center'),axs,textCondX);
%     textCondY = (jointUnitsCount(jj) +  jointUnitsCount(jj+1))/2;
%     addedConditional = {ones(size(masterJoints))};
% end

end
function sortedOrder = sortCategories(labeledCategories)
removedChar = cellfun(@(a) {a(1), a(2:end)}, cellstr(unique(labeledCategories)), 'UniformOutput', false);
sortedOrder = cellfun(@(b) b{2},removedChar,'UniformOutput', false);
[sortedOrder, sortedInd] = natsort(sortedOrder);
sortedOrder = cellfun(@(a,b) strcat(a{1}, b), removedChar(sortedInd), sortedOrder, 'UniformOutput', false);
end