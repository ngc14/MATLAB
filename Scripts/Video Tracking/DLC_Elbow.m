monkey = 'Skipper';
date = '02_28_2020';
cameras = {'Arm 1', 'Arm 2'};
conditions = {'Extra Small Sphere', 'Large Sphere', 'Photocell'};
numFramesToMatch = 15;
alignSegs = {'Reaction_time', 'Grasp Duration', 'Withdrawal Duration'};
before_time = -.3;
after_time = .3;
FR = 120;

cameraPaths = dir(['S:\Lab\',monkey,'\All Data\', monkey,'_',date,'\Videos\Camera*']);
cameraPaths = cameraPaths(contains({cameraPaths.name}, cameras),:);
delimiter = '\t';
formatSpec = '%*s%s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

fileName = {};
coords = {};
for c = 1:length(cameras)
    cameraPath = cameraPaths(c,:);
    allTrials = dir([cameraPath.folder, '\',cameraPath.name, '\Results\*.csv']);
    
    for t = 1:length(allTrials)
        [bodyParts{c}, coords{c,t}] = loadCoordinates([allTrials(t).folder, '\', allTrials(t).name]);
        fileName{c,t} = allTrials(t).name;
    end
end
camCombs = nchoosek(1:length(cameras),2);
p3d = cell(1,length(bodyParts{1}));
fileName = cellfun(@(a) a(1:regexp(a,'D')-1), fileName,'UniformOutput', false);
[~,sortInd] = natsort(fileName(1,:));
fileName = fileName(:,sortInd);
coords = coords(:,sortInd);
for b = 1:size(camCombs,1)
    camComb = camCombs(b,:);
    cam1 = camComb(1);
    cam2 = camComb(2);
    cam1Coords = coords(cam1,:);
    cam2Coords = coords(cam2,:);
    cam1Mat = load(['S:\Lab\ngc14\Camera Calibration\', cameras{cam1},'.mat']);
    cam2Mat = load(['S:\Lab\ngc14\Camera Calibration\', cameras{cam2},'.mat']);
    cam1Mat = cam1Mat.cameraParams;
    cam2Mat = cam2Mat.cameraParams;
    
    frames1 = read(VideoReader([cameraPaths(cam1).folder, '\', ...
        cameraPaths(cam1).name,'\Renamed\',fileName{cam1,round(length(fileName)/2)},'.avi']));
    frames2 = read(VideoReader([cameraPaths(cam2).folder, '\', ...
        cameraPaths(cam2).name,'\Renamed\',fileName{cam2,round(length(fileName)/2)},'.avi']));
    frames = randperm(min(size(frames1,4),size(frames2,4)),numFramesToMatch);
    for f = 1:numFramesToMatch
        f1 = undistortImage(frames1(:,:,:,frames(f)),cam1Mat);
        f2 = undistortImage(frames2(:,:,:,frames(f)),cam2Mat);
        p1=detectKAZEFeatures(f1);
        p2=detectKAZEFeatures(f2);
        [features1, valid_points1] = extractFeatures(f1, p1);
        [features2, valid_points2] = extractFeatures(f2, p2);
        indexPairs = matchFeatures(features1, features2, 'Unique', true);
        matchedPoints1 = valid_points1(indexPairs(:,1),:);
        matchedPoints2 = valid_points2(indexPairs(:,2),:);
        [~, epipolarInliers] = estimateEssentialMatrix(matchedPoints1, matchedPoints2, cam1Mat, cam2Mat);
        matches(f) = sum(epipolarInliers);
    end
    [~,bestMatch] = max(matches);
    f1 = undistortImage(frames1(:,:,:,frames(bestMatch)),cam1Mat);
    f2 = undistortImage(frames2(:,:,:,frames(bestMatch)),cam2Mat);
    p1=detectKAZEFeatures(f1);
    p2=detectKAZEFeatures(f2);
    [features1, valid_points1] = extractFeatures(f1, p1);
    [features2, valid_points2] = extractFeatures(f2, p2);
    indexPairs = matchFeatures(features1, features2, 'Unique', true);
    matchedPoints1 = valid_points1(indexPairs(:,1),:);
    matchedPoints2 = valid_points2(indexPairs(:,2),:);
    [E, epipolarInliers] = estimateEssentialMatrix(matchedPoints1, matchedPoints2, cam1Mat, cam2Mat);
    inlierPoints1 = matchedPoints1(epipolarInliers, :);
    inlierPoints2 = matchedPoints2(epipolarInliers, :);
    [orient, loc] = relativeCameraPose(E, cam1Mat, inlierPoints1, inlierPoints2);
    camMatrix1 = cameraMatrix(cam1Mat, eye(3), [0 0 0]);
    [R, t] = cameraPoseToExtrinsics(orient, loc);
    camMatrix2 = cameraMatrix(cam2Mat, R, t);
    for d = 1:length(bodyParts{1})
        p3d{d} = cellfun(@(a,b) triangulate(a{d}(1:min(length(a{d}),length(b{d})),:), ...
            b{d}(1:min(length(a{d}),length(b{d})),:),camMatrix1,camMatrix2),...
            cam1Coords,cam2Coords,'UniformOutput', false);
    end
    
    FW3D = p3d{cellfun(@(a) strcmp(a,'ForearmW'), bodyParts{1})};
    FW3D = cellfun(@(a) squeeze(smooth3(reshape(a,[length(a),1,3]))), FW3D, 'UniformOutput', false);
    FE3D = p3d{cellfun(@(a) strcmp(a,'ForearmE'), bodyParts{1})};
    FE3D = cellfun(@(a) squeeze(smooth3(reshape(a,[length(a),1,3]))), FE3D, 'UniformOutput', false);
    FVec = cellfun(@(a,b) a-b, FW3D, FE3D, 'UniformOutput', false);
    D3D = p3d{cellfun(@(a) strcmp(a,'Deltoid'), bodyParts{1})};
    D3D = cellfun(@(a) squeeze(smooth3(reshape(a,[length(a),1,3]))), D3D, 'UniformOutput', false);
    B3D = p3d{cellfun(@(a) strcmp(a,'Biceps'), bodyParts{1})};
    B3D = cellfun(@(a) squeeze(smooth3(reshape(a,[length(a),1,3]))), B3D, 'UniformOutput', false);
    BVec = cellfun(@(a,b) a-b, B3D, D3D, 'UniformOutput', false);
    crossVecs = cellfun(@(a,b) cross(a,b), FVec, BVec, 'UniformOutput', false);
    crossVecs = cellfun(@(a) sqrt(sum(a.^2,2)), crossVecs, 'UniformOutput', false);
    dotVecs = cellfun(@(a,b) dot(a,b,2), FVec, BVec, 'UniformOutput', false);
    FE{b} = cellfun(@(a,b) atan2d(a,-b), crossVecs, dotVecs, 'UniformOutput', false);
    clear p3d;
end
%% Load trial informaiton from Arduino output text file (.txt)
% Open the text file.
fileID = fopen(['S:\Lab\',monkey,'\All Data\', monkey,'_',date, '\Arduino\', monkey, '_', date, '.txt'],'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
raw = [dataArray{1:end-1}];

% remove lines that are not task related (pause, headers)
trialInfo = raw(ismember(raw(:,3), conditions),:);

% Remove false start trials. Videos are not recorded for false starts.
falseStartIdx=strcmp(trialInfo(:,10),'False Start');
trialInfo=trialInfo(~falseStartIdx,:);

restIdx = strcmp(trialInfo(:,3), 'Rest');
trialInfo = trialInfo(~restIdx,:);
times = cell2mat(cellfun(@str2double, trialInfo, 'UniformOutput', false));

[~, alignInds] = cellfun(@(a) find(strcmp(raw, a)), alignSegs, 'UniformOutput', false);
alignInds = cellfun(@(a) a-1, alignInds, 'UniformOutput', false);
[~, startTimeInd] = find(strcmp(raw, 'Reaction_time'));
[~, stopTimeInd] = find(strcmp(raw, 'Withdrawal Duration'));
%%
for b =1:length(FE)
    for c = 1:length(conditions)
        trials = find(strcmp(trialInfo(:,3), conditions{c}) & ~isnan(str2double(trialInfo(:,10))));
        offsets = frameArduino(trials)-cellfun(@length, FE{b}(trials))';
        trials = trials(offsets<8);
        %figure();
        hold on;
        startInd = 1;
        for a = 1:length(alignInds)
            if(strcmp(alignSegs{a}, 'Reaction_time'))
                [startFrame{1:length(trials)}] = deal(1);
                [endFrame{1:length(trials)}] = deal(after_time-before_time*FR);
            else
                centerFrame = round(sum(times(trials,startTimeInd:alignInds{a}),2)*(FR/1000));
                startFrame = num2cell(centerFrame+round(FR*before_time))';
                if(strcmp(alignSegs{a}, 'Withdrawal Duration'))
                    endFrame = num2cell(centerFrame + round(FR/1000*.01))';
                else
                    endFrame = num2cell(centerFrame + round(FR*after_time))';
                end
            end
            alignedKin{c,a} = cellfun(@(a,b,c) a(b:c), FE{b}(trials),startFrame,endFrame, 'UniformOutput', false);
            %cellfun(@(a) plot(a,'b'), FE(trials));
        end
        clear startFrame endFrame
    end
end