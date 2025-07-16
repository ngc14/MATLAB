function synthesizeGraphs(currDir, yellowThumb)
foldersOnPath = regexp(currDir, '[\\]');
segment = currDir(foldersOnPath(end-1)+1:end-1);
cond = currDir(foldersOnPath(end-2)+1:foldersOnPath(end-1)-1);
figDir = dir([currDir, 'Track\']);
pathToReferenceImages = 'S:\Lab\ngc14\Track_Videos\Reference_Images\';

% formatting
figHandle = figure();
axHandle = gca(figHandle);
set(axHandle, 'color', [0 0 0]);
hold(axHandle, 'on');
title(axHandle, cond);
set(axHandle, 'YDir', 'reverse');

cond = regexprep(cond,' ','_');
if(strcmp(segment, 'Reach'))
    % superimpose background image onto reach plots
    im = imread([pathToReferenceImages, cond, '.png']);
    xIm = [0 size(im,2); 0 size(im,2)];
    yIm = [0 0; size(im,1) size(im,1)];
    zIm = [0 0; 0 0];
    surf(axHandle, xIm, yIm, zIm, 'cData', uint8(im),...
        'FaceColor', 'texturemap');
    set(axHandle, 'xLim', [0 size(im,2)]);
    set(axHandle, 'yLim', [0 size(im,1)]);
end
disp(cond);
% colors used for the 5 different markers
colorMap = {'r', 'g', 'b', 'y', 'w'};
for i = 1:size(figDir,1)
    if(~figDir(i).isdir)
        
        currFig = openfig([figDir(i).folder, '\', figDir(i).name], 'invisible');
        axesObjs = gca(currFig);
        dataObjs = get(axesObjs, 'Children');
        
        scatterInds = arrayfun(@(x) strcmp(x.Type, 'scatter'), dataObjs);
        % get all plotted data points
        scatterInfo = dataObjs(scatterInds);
        droppedInds = arrayfun(@(x) strcmp(x.Marker, 'x'), scatterInfo);
        % get number of different colored markers (yellow and orange -> 2)
        colors = unique(arrayfun(@(x) x.MarkerFaceColor(2), scatterInfo(~droppedInds)));
        
        if(size(find(~droppedInds),1)==11)
            % x direction
            ydata = get(dataObjs, 'YData');
            % y direction
            zdata = get(dataObjs, 'ZData');
            % frame number
            xdata = get(dataObjs, 'XData');
            
            % only look at data points that were not dropped
            ydata = ydata(scatterInds);
            ydata = ydata(~droppedInds);
            zdata = zdata(scatterInds);
            zdata = zdata(~droppedInds);
            xdata = xdata(scatterInds);
            xdata = xdata(~droppedInds);
            
            % get position of all markers in first frame
            firstFrameY = cellfun(@(x) x(1), zdata);
            firstFrameX = cellfun(@(x) x(1), ydata);
            
            % get yellow markers (s value of hsv space is larger for yellow
            % than orange)
            yellowInds = arrayfun(@(x) x.MarkerFaceColor(2) == max(colors),...
                scatterInfo(~droppedInds));
            
            % sort all markers from right most to left most
            [~, xmap] = sort(firstFrameX, 'descend');
            % sort all markers from bottom most to top most
            [~, ymap] = sort(firstFrameY, 'descend');
            
            % get potential thumb markers based on the thumb color given
            if(yellowThumb)
                thumbCands = ymap(yellowInds(ymap));
            else
                thumbCands = ymap(~yellowInds(ymap));
            end
            
            % the four most right markers on the most distal digits (except
            % for in large sphere case)
            midMarkers = xmap(1:4);
            % sort the distal digits from bottom most to top most
            [~, rankedMids] = sort(firstFrameY(midMarkers), 'descend');
            % defined ordering used for plotting
            midOrdering = midMarkers(rankedMids)';
            
            if(strcmp(cond, 'Large_Sphere') && strcmp(segment, 'Grasp'))
                % sort thumb candidates from right most to left most
                if(yellowThumb)
                    checkThumbCands = xmap(yellowInds(xmap));
                    % get rid of wrist value
                    checkThumbCands = checkThumbCands(1:end-1);
                else
                    checkThumbCands = xmap(~yellowInds(xmap));
                end
                % if the thumb candidate is not also within the range of the
                % left most markers (leeway given bc proximal knuckles may 
                % be further back than thumb in this conformation) then
                % get the next candidate
                while(abs(firstFrameX(thumbCands(1)) - ...
                        firstFrameX(checkThumbCands(end))) > 10)
                    thumbCands = thumbCands(2:end);
                end
                
                % get markers furthest to the right and the bottom
                [~, bottomRight] = sort((firstFrameX + firstFrameY), 'descend');
                % do not include the thumb candidate
                bottomRight = bottomRight(bottomRight~=thumbCands(1));
                midMarkers = bottomRight(1:4);
                % split into yellow and orange markers
                yellowMids = midMarkers(ismember(midMarkers, xmap(yellowInds(xmap))));
                orangeMids = midMarkers(ismember(midMarkers, xmap(~yellowInds(xmap))));
                % sort from left to right and top to bottom
                [yellowX, sortYellowX] = sort(firstFrameX(yellowMids), 'descend');
                [orangeX, sortOrangeX] = sort(firstFrameX(orangeMids), 'descend');
                [yellowY, sortYellowY] = sort(firstFrameY(yellowMids), 'ascend');
                [orangeY, sortOrangeY] = sort(firstFrameY(orangeMids), 'ascend');
                % whichever sorting is more exclusive to define the
                % different markers, use that sorting
                if(abs(diff(yellowX)) < abs(diff(yellowY)))
                    sortedYellow = sortYellowY;
                else
                    sortedYellow = sortYellowX;
                end
                if(abs(diff(orangeX)) < abs(diff(orangeY)))
                    sortedOrange = sortOrangeY;
                else
                    sortedOrange = sortOrangeX;
                end

                if(min(firstFrameY(yellowMids)) < min(firstFrameY(orangeMids)))
                    % ordering if yellow is the top most marker (digit 5)
                    midOrdering = [orangeMids(sortedOrange(2)),yellowMids(sortedYellow(2)), ...
                        orangeMids(sortedOrange(1)), yellowMids(sortedYellow(1))];
                else
                    % ordering if orange is the top most marker (digit 5)
                    midOrdering = [yellowMids(sortedYellow(2)), orangeMids(sortedOrange(2)),...
                       yellowMids(sortedYellow(1)), orangeMids(sortedOrange(1))];
                end
            end
            
            ordered = [thumbCands(1) midOrdering];
            
            if(strcmp(segment, 'Grasp'))
                framesNotDropped = xdata(ordered);
                if(cellfun(@(x) x(1) == 1, framesNotDropped))
                    for p = 1:size(ordered,2)
                        % plot individual points for each marker in grasp
                        pinfo = plot(axHandle, firstFrameX(ordered(p)), firstFrameY(ordered(p)),...
                            'Marker', 'o', 'MarkerFaceColor', colorMap{p}, ...
                            'MarkerEdgeColor', colorMap{p}, 'MarkerSize', 3);
                        % used for callback function
                        pinfo.UserData = figDir(i).name;
                        %plot(axHandle, xVals, yVals, ':m');
                    end
                end
            elseif(strcmp(segment, 'Reach'))
                xVals = ydata(ordered);
                yVals = zdata(ordered);
                for p = 1:size(xVals,1)
                    % plot trajectory for each marker in reach
                    pinfo = plot(axHandle, xVals{p}, yVals{p}, 'Color',colorMap{p},...
                        'LineStyle', '-');
                    pinfo.UserData = figDir(i).name;
                end
            end
        end
    end
end
% add callback function when clicking line or point to identify
% corresponding trial and video
dcm = datacursormode(figHandle);
datacursormode on;
set(dcm, 'updatefcn', @videonamecallback);
saveas(figHandle, [currDir, cond,'_Summary_Plot'], 'fig');
end
function output_txt = videonamecallback(~,event_obj)
% Display the position of the data cursor
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
t = get(event_obj, 'Target');
output_txt = {['X: ', num2str(pos(1),4)],...
    ['Y: ', num2str(pos(2),4)],...
    ['Vid: ', get(t, 'UserData')]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
end
