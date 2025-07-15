monkey = 'Skipper';
date = '02_07_2020';
view1 = 'Arm 1';
view2 = 'Arm 2';
dir1 = dir(['U:\Lab\',monkey,'\All Data\', monkey,'_', date,'\Videos\Camera_',view1,'\Renamed']);
dir2 = dir(['U:\Lab\',monkey,'\All Data\', monkey,'_', date,'\Videos\Camera_',view2,'\Renamed']);
dir1 = dir1(~[dir1.isdir]);
dir2 = dir2(~[dir2.isdir]);
cameraParams1 = load(['U:\Lab\ngc14\Camera Calibration\',view1]);
cameraParams2 = load(['U:\Lab\ngc14\Camera Calibration\',view2]);
cameraParams1 = cameraParams1.cameraParams;
cameraParams2 = cameraParams2.cameraParams;
goodTransform = 'n';
iter = 1;
numFramesToMatch = 10;

while(~strcmp(goodTransform, 'y'))
    vid1 = VideoReader([dir1(iter).folder,'\',dir1(iter).name]);
    vid2 = VideoReader([dir2(iter).folder,'\',dir2(iter).name]);
    
    frames1 = read(vid1);
    frames2 = read(vid2);
    
    if(strcmp(view1,'Arm 1'))
        for i = 1:size(frames1,4)
            frames1(:,:,:,i) = flip(flip(frames1(:,:,:,i),2),1);
        end
    end
    frames = randperm(min(size(frames1,4),size(frames2,4)));

    for f = 1:numFramesToMatch
        f1 = undistortImage(frames1(:,:,:,frames(f)),cameraParams1);
        f2 = undistortImage(frames2(:,:,:,frames(f)),cameraParams2);
        p1=detectKAZEFeatures(f1);
        p2=detectKAZEFeatures(f2);
        [features1, valid_points1] = extractFeatures(f1, p1);
        [features2, valid_points2] = extractFeatures(f2, p2);
        indexPairs = matchFeatures(features1, features2, 'Unique', true);
        matchedPoints1 = valid_points1(indexPairs(:,1),:);
        matchedPoints2 = valid_points2(indexPairs(:,2),:);
        [~, epipolarInliers] = estimateEssentialMatrix(...
            matchedPoints1, matchedPoints2, cameraParams1, cameraParams2);

        matches(f) = sum(epipolarInliers);
    end
    [~,bestMatch] = max(matches);
    
    f1 = undistortImage(frames1(:,:,:,frames(bestMatch)),cameraParams1);
    f2 = undistortImage(frames2(:,:,:,frames(bestMatch)),cameraParams2);
    p1=detectKAZEFeatures(f1);
    p2=detectKAZEFeatures(f2);
    [features1, valid_points1] = extractFeatures(f1, p1);
    [features2, valid_points2] = extractFeatures(f2, p2);
    indexPairs = matchFeatures(features1, features2, 'Unique', true);
    matchedPoints1 = valid_points1(indexPairs(:,1),:);
    matchedPoints2 = valid_points2(indexPairs(:,2),:);
    [E, epipolarInliers] = estimateEssentialMatrix(...
            matchedPoints1, matchedPoints2, cameraParams1, cameraParams2);
        
    inlierPoints1 = matchedPoints1(epipolarInliers, :);
    inlierPoints2 = matchedPoints2(epipolarInliers, :);
    showMatchedFeatures(f1, f2, inlierPoints1, inlierPoints2);

    [orient, loc] = relativeCameraPose(E, cameraParams1, inlierPoints1, inlierPoints2);
    [R, t] = cameraPoseToExtrinsics(orient, loc);
    camMatrix1 = cameraMatrix(cameraParams1, eye(3), [0 0 0]);
    camMatrix2 = cameraMatrix(cameraParams2, R, t);
    
    points3D = triangulate(matchedPoints1, matchedPoints2, camMatrix1, camMatrix2);
    cameraSize = 0.3;
    figure
    [orient, loc] = estimateWorldCameraPose(matchedPoints2.Location, points3D, cameraParams2);

    plotCamera('Size', cameraSize, 'Color', 'r', 'Label', '1', 'Opacity', 0);
    hold on
    pcshow(pointCloud(points3D), 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 45);
    plotCamera('Location', loc, 'Orientation', orient, 'Size', cameraSize, ...
        'Color', 'b', 'Label', '2', 'Opacity', 0);
    camorbit(0, -30);
    camzoom(1.5);

    goodTransform = input('Transform good? (y/n)');
    iter = iter+1;
end