function figureMap = mapUnitVals(vMask,siteMask,valsToGraph,emptyInds,countZeroIn,numStepsIn,rangeIn)
zeroColor = [.7 .7 .7];
txColor = [0 .5 0];
fontSize = 8;

countZero = countZeroIn;
numSteps = numStepsIn;
imDims = size(vMask);
numSites = length(siteMask);

if(~isempty(rangeIn))
    countRange = linspace(rangeIn(1), rangeIn(2),numSteps);
else
    countRange = linspace(max(min(valsToGraph),nanmedian(valsToGraph)- ...
        nanstd(valsToGraph)),min(nanmedian(valsToGraph) +...
        nanstd(valsToGraph),max(valsToGraph)),numSteps);
end
if(all(round(valsToGraph)==valsToGraph))
    countRange = round(countRange);
end
positiveMap = flipud([linspace(.7,1,numSteps)',linspace(.1,.6,numSteps)',linspace(.2,0,numSteps)']);
positiveMap = flipud(colormap(summer(numSteps)));
negativeMap = flipud([linspace(0,.7,numSteps)', linspace(0,.7,numSteps)',linspace(.5,1,numSteps)']);
%positiveMap = flipud(hot(numSteps*2));
negativeMap = (cool(numSteps));
%positiveMap = positiveMap(all(positiveMap<=[1 .7 0],2),:);
positiveMap = positiveMap(1:numSteps,:);
%%
cMap = positiveMap;
if(countZero)
    if(~any(countRange==0))
        zeroInd = find(countRange<0,1,'last');
        if(isempty(zeroInd))
            cMap = [zeroColor; cMap];
            zeroInd = 1;
        end
        %countRange = [countRange(1:zeroInd), 0, countRange(zeroInd+1:end)];
    else
        zeroInd = find(countRange==0);
    end
    cMap(zeroInd,:) = zeroColor;
end
cMap(1:sum(countRange<0)-countZero,:) = flipud(negativeMap(1:sum(countRange<0)-countZero,:));
cMap = [cMap;1 1 1;];
tickNames = arrayfun(@(t) num2str(t,'%0.2f'), countRange, 'UniformOutput', false);
%%
MM_image = Inf(imDims(1),imDims(2));
binUnits = discretize(valsToGraph,countRange);
binUnits(valsToGraph<=countRange(1)) = 1;
binUnits(valsToGraph>=countRange(end)) = numSteps;
if(countZero)
    binUnits(binUnits>=zeroInd) = binUnits(binUnits>=zeroInd)+1;
    binUnits(valsToGraph==0) = zeroInd;
    tickNames = [{'0.00'},tickNames];
end
for i = 1:numSites
    % currProps = regionprops(siteMask{i}, 'Centroid');
    % currMask = images.roi.Circle;
    % currMask.Center = currProps.Centroid;
    % currMask.Radius = 10;
    % siteMask{i} = createMask(currMask, image(siteMask{i}));
    MM_image(siteMask{i}) = binUnits(i);
end
if(any(isnan(MM_image(:))))
    MM_image = MM_image + 1;
    MM_image(isnan(MM_image)) = 1;
    cMap = [.85 .85 .85; cMap];
    tickNames = [{'NaN'}, tickNames];
end
% mask is black, nans are white
% MM_image(vMask==0) = length(cMap)-1;
MM_image(isinf(MM_image) & vMask==1) = length(cMap);
figureMap = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
h = imshow(ind2rgb(MM_image,cMap),'XData', 1:size(MM_image,2),...
    'YData',1:size(MM_image,1));hold on;
bnds = arrayfun(@(b) bwboundaries(siteMask{b}), find(emptyInds), 'UniformOutput',false);
cellfun(@(tr) plot(tr{1}(:,2), tr{1}(:,1), 'LineStyle','--', 'LineWidth',.5,...
    'Color','k'),bnds,'UniformOutput',false);
set(h,'alphadata',~isinf(MM_image));
h2 = imshow(ind2rgb(vMask,[.3 .3 .3; 1 1 1]));
set(h2, 'alphadata',double(vMask~=1));
colormap(h.Parent,cMap(1:end-1,:));
cb = colorbar(h.Parent);
cb.Ticks = linspace(cb.Limits(1),cb.Limits(end),size(cMap,1));
cb.TickLabels = [tickNames(:);' ';]';
for i = 1:numSites
    currVal = valsToGraph(i);
        % text(gca(figureMap), vCentroids(i,1), vCentroids(i,2),...
        %     num2str(currVal,'%0.2f'),'FontWeight', 'bold', ...
        %     'HorizontalAlignment', 'center','Color', txColor, 'FontSize', fontSize);
end