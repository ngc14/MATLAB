function [error, updatedProps] = calcError(allCurrProps, props, sz, exclude)
error = zeros(1,size(props,2));
% get all centroid candidates in current frame
t = allCurrProps{1};
for c = 2:size(allCurrProps,2)
    t = cat(1,t,allCurrProps{c});
end
allCents = cat(1, t.Centroid);

wristVal = 0;
% identify wrist value if it is within a certain range
if(any(allCents(:,1) < quantile(allCents(:,1), .1)))
    wristVal = min(allCents(:,1));
end
% exclude centroids that exceed 2stds in all directions
inds = allCents(:,1)~= wristVal & allCents(:,1) < quantile(allCents(:,1),.9)...
    & allCents(:,2) < mean(allCents(:,2)) + 2*std(allCents(:,2)) & ...
    allCents(:,2) > mean(allCents(:,2)) - 2*std(allCents(:,2));
% using the likely tracked markers, calculate the mean and std
divLine = mean(allCents(inds,:));
stdC = std(allCents(inds,:));

for j = 1:size(props,2)
    % only look at markers of the same color
    currProps = allCurrProps{props(j).ColorIndex};
    centroids = cat(1, currProps.Centroid);
    % if the reference marker is moving substantially in the
    % forward/backward direction, exclude markers in the current frame 
    % that did not move enough in the direction for the current marker to 
    % match.
    if(props(j).Offset(1) < -7)
        xLim = centroids(:,1)' > props(j).Centroid(1);
    elseif(props(j).Offset(1) > 7)
        xLim = centroids(:,1)' < props(j).Centroid(1);
    else
        xLim = 1:size(currProps,1);
    end
    % do not limit matching if all centroids are past offset
    if(size(find(xLim),1) < 1)
        xLim = 1:size(currProps,1);
    end
    
    % do not exclude any markers for matching unless specificed by user
    notExcluded = ones(1,size(currProps,1));
    for ex = 1:size(exclude,2)
        vals = cellfun(@(x) norm(x - [exclude{ex}],2), {currProps.Centroid}, 'UniformOutput',false);
        notExcluded = notExcluded & cellfun(@(x) x>.5, vals);
    end
    % exclude current centroids that are not within a certain range of the
    % mean centroid of all current centroids
    notExcluded = notExcluded & centroids(:,1)' < divLine(1)+2.5*stdC(1) & ...
        centroids(:,2)' < divLine(2) + 3.5*stdC(2) & ...
        centroids(:,2)' > divLine(2) - 3.5*stdC(2) & xLim;
    xs = centroids(centroids(:,1) ~= wristVal & notExcluded',1);
    ys = centroids(centroids(:,1) ~= wristVal & notExcluded',2);
    
    % feature vector for detected markers
    [feats, err] = deal(zeros(size(currProps,1),6));
    % array of angles from one centroid to all other centroid (NOT USED)
    angles = zeros(size(currProps,1));
    % array specifying relative x and y ordering of the marker
    ordering = zeros(size(currProps,1),2);
    for n = 1:size(currProps,1)
        for z = 1:size(currProps,1)
            % angle between centroid and all other centroids
            angle = atan2(currProps(n).Centroid(1)-currProps(z).Centroid(1), ...
                currProps(n).Centroid(2)-currProps(z).Centroid(2));
            angles(n,z) = angle * (180/pi);
        end
        % decide if a marker is medial marker or a proximal marker based on
        % whether its x position is more than or less than the mean
        % centroid x position
        med = currProps(n).Centroid(1) > divLine(1);
        if(med)
            % get all medial markers
            subset = ys(xs > divLine(1));
        else
            % get all proximal markers
            subset = ys(xs < divLine(1));
        end
        % sort subset of markers from bottom to top (digit 1 -> 5)
        yPos = sort(subset, 'descend');
        yInd = find(yPos==currProps(n).Centroid(2));
        if(isempty(yInd))
            yInd = 0;
        end
        % cap ordering in case of more centroids than markers
        if(med)
            yInd = min(yInd,2);
        else
            yInd = min(yInd,3);
        end
        ordering(n,:) = [med, yInd(1)];
        
        % features used for calculating error from marker in previous frame
        % to all markers in current frame
        % 1: offset distance
        % 2: minor axis length of marker
        % 3: area of marker
        % 4: orientation of marker (NOT USED)
        % 5: medial or proximal marker
        % 6: digit 1->5
        feats(n,:) = [norm(currProps(n).Centroid-props(j).Centroid+props(j).Offset), ...
            abs(currProps(n).MinorAxisLength - props(j).MinorAxisLength),...
            abs(currProps(n).Area - props(j).Area), ...
            abs(currProps(n).Orientation - props(j).Orientation),...
            abs(ordering(n,1) - props(j).Ordering(1))...
            abs(ordering(n,2) - props(j).Ordering(2))];
    end
    
    % angle differences (NOT USED)
    angles = sort(angles,2);
    propAngles = repmat(props(j).Angles, [size(currProps,1) 1]);
    if(~isequal(size(propAngles), size(angles)))
        errAngles = sum(zeros(size(angles)),2);
    else
        diffAngles = abs(angles - propAngles);
        errAngles = sum(diffAngles,2);
        errAngles = errAngles./max(errAngles);
    end
    
    % normalize all features to maximum area except for distance metric
    normalized = feats(:,2:end)./max(feats(:,2:end));
    normalized(isnan(normalized)) = 0;
    
    % error weighting
    err(:,1) = .1*feats(:,1);
    err(:,2) = .7*normalized(:,1);
    err(:,3) = .7*normalized(:,2);
    err(:,4) = .05*normalized(:,3);
    err(:,5) = .7*normalized(:,4);
    err(:,6) = .8*normalized(:,5);
    err = sum(err,2);
    
    % if all markers have been excluded -> marker is dropped
    if(~any(notExcluded))
        % error is 100%
        error(j) = 1;
    else
        % assign marker to centroid with minimum error
        [minErr,closest] = min(err(notExcluded));
        
        error(j) = minErr/max(err);
        propsLim = currProps(find(notExcluded));
        matchedProps = propsLim(closest);
    end
    
    % marker is considered dropped if error is above 50%
    if(error(j) > .5)
        updatedProps(j).Dropped = 1;
        % move centroid to its previous location plus the offset from the
        % reference marker
        updatedProps(j).Centroid = props(j).Prev.Centroid - props(j).Offset;
        
        % if marker is out of frame
        if(updatedProps(j).Centroid(1) > sz(2))
            updatedProps(j).Centroid(1) = sz(2);
        end
        if(updatedProps(j).Centroid(2) > sz(1))
            updatedProps(j).Centroid(2) = sz(1);
        end
        
        % keep features from marker when it was not dropped
        updatedProps(j).Area = props(j).Area;
        updatedProps(j).MinorAxisLength = props(j).MinorAxisLength;
        updatedProps(j).Angles = props(j).Angles;
        updatedProps(j).Orientation = props(j).Orientation;
        updatedProps(j).Prev = props(j).Prev;
    else
        % if it is no longer dropped in 'current' frame, do not detect it
        % until it is visible in 'current' frame and 'next' frame
        if(props(j).Dropped == 1)
            updatedProps(j) = props(j);
        else
            updatedProps(j).Dropped = 0;
        end
        % update features
        updatedProps(j).Centroid = matchedProps.Centroid;
        updatedProps(j).Area = matchedProps.Area;
        updatedProps(j).MinorAxisLength = matchedProps.MinorAxisLength;
        updatedProps(j).Angles = angles(closest,:);
        updatedProps(j).Orientation = matchedProps.Orientation;
        updatedProps(j).Prev = props(j);
        updatedProps(j).Ordering = ordering(closest,:);
    end
    
    % static features of markers
    updatedProps(j).Marker = props(j).Marker;
    updatedProps(j).ColorIndex = props(j).ColorIndex;
    updatedProps(j).Offset = [0 0];
    updatedProps(j).Ordering = props(j).Ordering;
end

end