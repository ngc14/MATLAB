%% User Input
segments = {'Reach', 'Grasp'};
conditions = {'Large Sphere', 'Small Sphere', 'Photocell', 'Small Cube'};

outputdriver = 'S:\';
outputfolder = 'Lab\ngc14\';

date = '12_07_2017';
dataname = ['Track_Videos\', date, '\'];
%%
prompt = 'What is the color of the distal thumb? 1 for yellow, 0 for orange: ';
yellowThumb = input(prompt);

resultFolder = [outputdriver, outputfolder, dataname, 'Results\'];
condFig = cell(size(conditions));
for s = 1:size(segments,2)
    for c = 1:size(conditions,2)
        condFold = [resultFolder, conditions{c}, '\'];
        segFold = [condFold, segments{s}, '\'];
        % for each condition and each segment, get the corresponding
        % summary plot
        synthesizeGraphs(segFold, yellowThumb);
        condFig{c} = openfig([segFold, regexprep(conditions{c},' ','_'),...
            '_Summary_Plot.fig'], 'invisible');
    end
    
    finalFig = figure();
    yLims = [Inf -Inf];
    xLims = [Inf -Inf];
    doAvgFig = ~exist([resultFolder, 'Reach_Avg_Graphs.fig'], 'file') ...
        && strcmp(segments{s}, 'Reach');
    if(doAvgFig)
        avgFig = figure();
    end
    for c = 1:size(conditions,2)
        plotPos = sprintf('22%d',c);
        p = subplot(plotPos, axes(finalFig));
        % black background
        set(p, 'color', [0 0 0]);
        info = get(condFig{c}, 'CurrentAxes');
        t = get(info, 'title');
        title(p, t.String);
        % copy current condition plot to subplot
        copyobj(allchild(info), p);
        if(strcmp(segments{s}, 'Reach') && doAvgFig)
            data = get(info, 'Children');
            lineInds = arrayfun(@(x) strcmp(x.Type, 'line'), data);
            lines = data(lineInds);
            % get superimposed background image
            im = data(~lineInds);
            % get number of markers (identified by different number of
            % colors)
            allColors = arrayfun(@(x) x.Color, lines, 'UniformOutput', false);
            colors = [];
            hashed = [];
            for a = 1:size(allColors,1)
                % convert color to single unique rgb value
                val = rgb2gray(allColors{a});
                % only add unique colors to colors array
                if(~ismember(val(1), hashed))
                    colors(end+1,:) = allColors{a};
                    hashed(end+1) = val(1);
                end
            end
            currAvgAx = subplot(plotPos, axes(avgFig));
            xIm = [0 size(im.CData,2); 0 size(im.CData,2)];
            yIm = [0 0; size(im.CData,1) size(im.CData,1)];
            % superimpose background image onto subplot
            surf(currAvgAx, xIm, yIm, [0 0; 0 0;], 'cData', im.CData, ...
                'FaceColor', 'texturemap')
            view(2);
            set(currAvgAx, 'XLim', [0 size(im.CData,2)]);
            set(currAvgAx, 'YLim', [0 size(im.CData,1)]);
            title(currAvgAx, t.String);
            
            hold(currAvgAx, 'on');
            for m = 1:size(colors,1)
                currColor = colors(m,:);
                % get indicies of all data points with the same color (same
                % marker)
                inds = arrayfun(@(x) sum(x.Color==currColor)==3, lines);
                allCurrFingerPoints = lines(inds);
                yPositions = get(allCurrFingerPoints, 'YData');
                xPositions = get(allCurrFingerPoints, 'XData');
                largestTrial = max(cellfun(@(x) size(x,2), xPositions));
                [xPlot, yPlot] = deal(NaN(size(xPositions,1),largestTrial));
                % interpolate each line plot to the longest trial using
                % spline interpolation
                for n = 1:size(xPositions,1)
                    space = (largestTrial-1)/(size(xPositions{n},2)-1);
                    sp = 1:(1/space):size(xPositions{n},2);
                    xPlot(n,:) = interp1(1:size(xPositions{n},2), xPositions{n}, sp, 'spline');
                    yPlot(n,:) = interp1(1:size(xPositions{n},2), yPositions{n}, sp, 'spline');
                end
                plot(currAvgAx, nanmean(xPlot), nanmean(yPlot), 'LineWidth',...
                    2, 'Color', currColor);
            end
        end
        
        % set all subplots to the smallest x and y axis limits that
        % includes all data points
        currLimsY = get(info, 'YLim');
        currLimsX = get(info, 'XLim');
        if(currLimsY(1) < yLims(1))
            yLims(1) = currLimsY(1);
        end
        if(currLimsY(2) > yLims(2))
            yLims(2) = currLimsY(2);
        end
        if(currLimsX(1) < xLims(1))
            xLims(1) = currLimsX(1);
        end
        if(currLimsX(2) > xLims(2))
            xLims(2) = currLimsX(2);
        end
    end
    
    % formatting
    axInd = arrayfun(@(x) strcmp(x.Type, 'axes'), finalFig.Children);
    allAx = finalFig.Children(axInd);
    arrayfun(@(x) set(x, 'YLim', yLims), allAx);
    arrayfun(@(x) set(x, 'XLim', xLims), allAx);
    arrayfun(@(x) set(x, 'YDir', 'reverse'), allAx);
    saveas(finalFig, [resultFolder, segments{s}, '_Graphs'], 'fig');
    
    % formatting for average reach plot
    if(doAvgFig)
        allAx = avgFig.Children;
        axesAvg =  arrayfun(@(x) strcmp(x.Type, 'axes'), avgFig.Children);
        allAx = allAx(axesAvg);
        arrayfun(@(x) set(x, 'YLim', yLims), allAx);
        arrayfun(@(x) set(x, 'XLim', xLims), allAx);
        arrayfun(@(x) set(x, 'YDir', 'reverse'), allAx);
        saveas(avgFig, [resultFolder, segments{s}, '_Avg_Graphs'], 'fig');
    end
end
