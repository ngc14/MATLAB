function track_2(vidInfo, start, finish, currFold, saveName)
manual = false;
debug = false;
% from Scale.avi video 65 pixels = 1 inch
pixelToInch = 65;
videoPath = [vidInfo.folder, '\', vidInfo.name];
vidObj = VideoReader(videoPath);

folders = regexp(vidInfo.folder, '[\\]');

datadriver = vidInfo.folder(1:folders(1));
datafolder = vidInfo.folder(folders(1)+1:folders(3));
expname = vidInfo.folder(folders(3)+1:folders(end));

% create average image for background subtraction or load it from dir
if(~exist([datadriver, datafolder,expname, 'average.mat'], 'file'))
    % image average of first x frames when there is no movement
    avgObj = VideoReader([datadriver, datafolder, expname, 'average.avi']);
    avg = imfilter(double(read(avgObj,1)), ones(3,3)/9);
    for i = 2:avgObj.NumberOfFrames
        avg = avg + imfilter(double(read(avgObj,i)), ones(3,3)/9);
    end
    avg = avg ./ (avgObj.NumberOfFrames);
    avg = rgb2gray(uint8(avg));
    save([datadriver, datafolder, expname,'average.mat'], 'avg');
else
    avg = load([datadriver, datafolder, expname, 'average.mat']);
    avg = avg.avg;
end

props = {};
firstFrame = double(read(vidObj, start));
% manual color selection
if(manual)
    prompt = 'Enter number of colors to track: ';
    numColor = input(prompt);
    % select point(s) to track
    handle = imshow(rgb2hsv(uint8(firstFrame)));
    % processing
    firstFrameProcessed = medfilt3(uint8(firstFrame), [7, 7, 1]);
    firstFrameProcessed = imfilter(firstFrameProcessed, ones(3,3)/9);
    hsv = rgb2hsv(firstFrameProcessed);
    h=hsv(:,:,1); s=hsv(:,:,2); v=hsv(:,:,3);
    imshow(hsv);
    % predefined regions for 3 different colors present
    for i = 1:numColor
        mask(:,:,i) = createMask(imellipse(gca));
        %imellipse, imrect
        %mask(:,:,i) = roipoly(firstFrameProcessed);
    end
    % color markers (region average values of h, s, and v)
    for i = 1:numColor
        cMark(i,1) = mean2(h(mask(:,:,i)));
        cMark(i,2) = mean2(s(mask(:,:,i)));
        cMark(i,3) = mean2(v(mask(:,:,i)));
    end
else
    % load saved color template for matching
    cMark = load([datadriver, datafolder,'Track_Videos\cMark.mat']);
    cMark = cMark.cMark;
    cMark(:,3) = .7;
    firstFrameProcessed = processIm(double(read(vidObj, start)), avg);
end

for c = 1:size(cMark,1)
    % segment first frame for each different color from template
    seg = colorSeg(uint8(firstFrameProcessed), cMark, c);
    % save information for each segmentation
    currProps = regionprops(logical(seg(:,:,1)), 'BoundingBox', ...
        'Centroid', 'Area', 'MinorAxisLength', 'EquivDiameter', 'Orientation');
    % angles NOT USED
    angles = [];
    for i = 1:size(currProps,1)
        % in first frame, marker must have an equivdiameter > 4 and an area
        % > 25
        if(currProps(i).EquivDiameter > 4 && currProps(i).Area > 25)
            for z = 1:size(currProps,1)
                if(currProps(z).Area> 25 && currProps(z).EquivDiameter>4)
                    angle = atan2(currProps(i).Centroid(1)-currProps(z).Centroid(1), ...
                        currProps(i).Centroid(2)-currProps(z).Centroid(2));
                    angles(end+1) = angle * (180/pi);
                end
            end
            % initalize props of each marker
            props(end+1).Centroid = currProps(i).Centroid;
            props(end).Area = currProps(i).Area;
            props(end).MinorAxisLength = currProps(i).MinorAxisLength;
            props(end).Angles = sort(angles);
            props(end).Orientation = currProps(i).Orientation;
            
            % static features of props that will not change
            props(end).ColorIndex = c;
            props(end).Dropped = 0;
            props(end).Offset = [0 0];
            props(end).Marker = size(props,2)+(6*floor(c/2));
            props(end).Prev = currProps(i);
        end
    end
end

numTrack = size(props,2);

% get all centroids
cents = cat(1, props.Centroid);
[~, val] = min(cents(:,1));
% sort them from left to right (except for the reference marker)
xPos = sort(cents(1:end ~= val,1));
% calculate the mean xposition to serve as the dividing line between
% markers
divLine = mean(xPos);
xs = cents(:,1);
ys = cents(:,2);
for c = 1:size(cMark,1)
    % get markers that are of the same color (excluding reference marker)
    yMap{c} = ys([props.ColorIndex]==c & 1:end ~= val);
    xMap{c} = xs([props.ColorIndex]==c & 1:end ~= val);
end
for i = 1:numTrack
    % get the position of all markers with the same color as the current
    % marker
    subset = yMap{props(i).ColorIndex};
    % identify as medial or proximal marker
    if(props(i).Centroid(1) > divLine)
        subset = subset(xMap{props(i).ColorIndex} > divLine);
    else
        subset = subset(xMap{props(i).ColorIndex} < divLine);
    end
    % sort from bottom to top to get relative position of marker
    yInd = find(sort(subset, 'descend')==props(i).Centroid(2));
    if(isempty(yInd))
        yInd = 0;
    end
    % initailize ordering
    props(i).Ordering = [props(i).Centroid(1)>divLine, yInd];
end
% colormap for track and speed plots
cm = colorcube.*255;
ind = (size(cm,1)/numTrack)*[1:numTrack]-(size(cm,1)/numTrack);
ct = cm(round(ind+1),:);

iterTrack = zeros(numTrack,2);
% initalize the plotting data
for i = 1:numTrack
    tracked(1, i, :) = [props(i).Centroid];
end

%for each consecutive frame
for i = start:finish
    prevProcessed = processIm(double(read(vidObj, i)), avg);
    nextProcessed = processIm(double(read(vidObj, i+1)), avg);
    % for culmination of individual color segmentations
    prevAdd = zeros(size(prevProcessed));
    nextAdd = zeros(size(nextProcessed));
    for c = 1:size(cMark,1)
        % get individual color segmentations for each previous and current
        % frame
        prevs{c} = colorSeg(uint8(prevProcessed), cMark, c);
        nexts{c} = colorSeg(uint8(nextProcessed), cMark, c);
        nextP{c} = regionprops(logical(nexts{c}(:,:,1)), 'BoundingBox',...
            'Centroid', 'Area', 'MinorAxisLength', 'EquivDiameter', 'Orientation');
        prevAdd = prevAdd + prevs{c};
        nextAdd = nextAdd + nexts{c};
    end
    % get reference marker offsest
    offset = getOffset(prevAdd, nextAdd);
    for j = 1:numTrack
        props(j).Offset = offset;
    end
    
    if(debug)
        fprintf('%d: [%f, %f]\n', i, offset(1), offset(2));
    end
    % assign markers and calculate error
    [errors, props] = calcError(nextP, props, size(firstFrame), []);
    checkRepeat = ones(1,numTrack);
    % check each marker to see if there are any repeat assignments
    while(any(checkRepeat))
        map = find(checkRepeat);
        for m = 1:size(map,2)
            j = map(m);
            checkRepeat(j) = 0;
            nextBest = 1;
            exclude = {};
            % stop looking for next best marker if you have looked through
            % more than the half of the available markers
            while(nextBest && size(exclude,2) < floor(numTrack/2))
                nextBest = 0;
                for a=1:numTrack
                    % do not compare marker against itself
                    aInd = 1+mod((a-1)+j+numTrack-1,numTrack);
                    % if markers are assigned to the same position
                    if(norm(props(j).Centroid-props(aInd).Centroid,2)<.5 && j~=aInd)
                        % if the current marker has the larger error
                        if(errors(j) > errors(aInd))
                            % reset it to its previous position and find
                            % its next best match (exclude the current
                            % marker when matching)
                            centroid = props(j).Prev.Centroid;
                            props(j) = props(j).Prev;
                            props(j).Centroid = centroid;
                            
                            exclude(end+1) = {props(aInd).Centroid};
                            nextBest = 1;
                        % the other marker has a larger error so check that
                        % marker again
                        elseif(errors(j) < errors(aInd))
                            checkRepeat(aInd) = 1;
                        end
                    end
                end
                if(nextBest)
                    % finding the next best match for the marker
                    [errors(j), props(j)] = calcError(nextP, props(j),...
                        size(firstFrame), exclude);
                end
            end
        end
    end
    % all markers that have errors smaller than .5 but were dropped in the
    % previous frame
    [props(errors <= .5 & [props.Dropped]).Dropped] = deal(.5);
    for j = 1:numTrack
        if(debug)
            map = {'Tracked', 'Other', 'Dropped'};
            fprintf('%d: %f - %s\n', props(j).Marker, 100*errors(j), map{ceil(1.1*props(j).Dropped)+1});
        end
        % update plotting data
        if(props(j).Dropped)
            iterTrack(j,:) = [NaN, NaN];
        else
            iterTrack(j,:) = props(j).Centroid;
        end
        props(j).Prev.Centroid = props(j).Centroid;
    end
    tracked((i+2-start), :, :) = iterTrack;
    
    if(debug)
        cents = cat(1, props.Centroid);
        if(i+2-start > 2)
            subplot(221);
            subimage(insertShape(uint8(prevProcessed), 'circle', ...
                [cents(:,1), cents(:,2), repmat(5, [numTrack,1])], 'color', ct, 'LineWidth', 3));
            text(cents(:,1)-20, cents(:,2), arrayfun(@num2str, ...
                [props.Marker], 'UniformOutput', false), 'color', 'white', 'FontSize', 7)
        end
        subplot(222);
        subimage(insertShape(uint8(nextProcessed), 'circle', [cents(:,1), ...
            cents(:,2), repmat(5, [numTrack,1])], 'color', ct,'LineWidth', 3));
        subplot(223);
        subimage(uint8(prevAdd));
        ax = subplot(224);
        subimage(uint8(nextAdd));
        %pause();
    end
end
%% FIGURE
cMap = hsv2rgb(cMark(:,1:3));
tr = figure();
ax2 = gca;
vl = figure();
ax3 = gca;
hold(ax2, 'on');
hold(ax3, 'on');
% plot frame vs movement
for i = 1:numTrack
    nans = isnan(tracked(:,i,1));
    if(any(nans))
        % get the frames in which the marker was found
        timesFound = find(diff(nans)==-1);
        % get the frames in which the marker was dropped
        timesDropped = find(diff(nans)==1);
        for inter = 1:size(timesFound,1)
            interval = timesDropped(inter)+1:timesFound(inter);
            midp = ceil(length(interval)/2);
            % if marker was only dropped for one/two frame(s)
            if(midp==1)
                offset = 0;
            else
                offset = 1;
            end
            % all frames prior to middle of dropped/found interval are
            % equal to the marker position when dropped
            tracked(interval(1):interval(midp-offset), i, :) = ...
                repmat(tracked(interval(1)-1, i, :),...
                [length(interval(1):interval(midp-offset)),1]);
            % middle of dropped/found interval should equal the midpoint of
            % the last detected position and next found position
            tracked(interval(midp), i, :) = (tracked(interval(1)-1,i,:)...
                + tracked(interval(end)+1,i,:))/2;
            % all frames after the middle of dropped/found interval are
            % equal to the marker position when found
            tracked(interval(midp+1):interval(end), i, :) = ...
                repmat(tracked(interval(end)+1, i, :),...
                [length(interval(midp+1):interval(end)),1]);
        end
        % if marker is still dropped at the end of the video
        if(nans(end))
            tracked(timesDropped(end)+1:end, i ,:)  = ...
                repmat(tracked(timesDropped(end), i, :),[size(tracked,1)-timesDropped(end),1]);
        end
        % plot dropped markers with red X
        dr = scatter3(ax2, find(nans)', tracked(nans,i,1),...
            tracked(nans,i,2), 100, 'x');
        dr.MarkerEdgeColor = [1 0 0];
    end
    % velocity values
    v = diff(tracked(:,i,1));
    w = diff(tracked(:,i,2));
    velocity(i,:) = abs(v+w)./pixelToInch;
    [v(end+1), w(end+1)] = deal(0);
    % vector plots
    quiv = quiver3(ax2,(1:size(tracked,1))', tracked(:,i,1), ...
        tracked(:,i,2), ones(size(tracked,1),1),v,w,0);
    % point plots
    pt = scatter3(ax2, find(~nans)', tracked(~nans,i,1), ...
        tracked(~nans,i,2), 16, 'MarkerEdgeColor', 'black');
    % coloring
    line = plot(ax3, 1:size(tracked,1)-1, velocity(i,:));
    pt.MarkerFaceColor= deal(cMap(props(i).ColorIndex,:));
    [line.Color, quiv.Color] = deal(ct(i,:)./255);
end
% formatting
firstFrame = double(read(vidObj, start));
xlabel(ax3, 'Frame Number');
ylabel(ax3, 'Inches');
hold(ax3, 'off');
xlim(ax2, [0 size(tracked,1)]);
ylim(ax2, [0 size(firstFrame, 2)]);
zlim(ax2, [0 size(firstFrame, 1)]);
set(ax2,'YDir','reverse');
set(ax2,'ZDir','reverse');
xlabel(ax2, 'Frame Number');
ylabel(ax2, 'Forward/Backward');
zlabel(ax2, 'Left/Right');
xIm = [size(tracked,1), size(tracked,1); size(tracked,1), size(tracked,1)];
yIm = [0 size(firstFrame,2); 0 size(firstFrame,2)];
zIm = [0 0; size(firstFrame,1) size(firstFrame,1)];
surf(ax2, xIm, yIm, zIm, 'cData', uint8(firstFrame),...
    'FaceColor', 'texturemap');
hold(ax2, 'off');
view(ax2, [-90 0 0]);
hgsave(tr, [currFold, 'Track\', saveName, '.fig']);
hgsave(vl, [currFold, 'Speed\', saveName, '.fig']);
clear tracked velocity;
close(tr); close(vl);
end
