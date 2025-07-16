username = 'ngc14';

monkey = 'Skipper';                   % animal name
referenceDate = '12_07_2018';                % experiment date
referenceDate = '10_30_2020';
dates = {'12_18_2020'};

referencePath = ['S:\Lab\',monkey,'\All Data\', monkey, '_',referenceDate, '\Imaging\'];
sizeInfo = folderSizeTree(referencePath);
runFolders = find([sizeInfo.level{:}]==1);
[~, runFolderInd] = max([sizeInfo.size{runFolders}]);
referenceGreenPath = sizeInfo.name{runFolders(runFolderInd)};
runNumber = referenceGreenPath(end-1:end);
refGreen = double(imread([referenceGreenPath, '\green', runNumber, '_edited.bmp']));

refGreen = rgb2gray(refGreen - mean(refGreen(:))); % mean-subtract and clip image to increase contrast
[refGreen, ~, clipLow, clipHigh] = OIClipH(refGreen, 1, 1.5, ones(768,768));

optimizer = registration.optimizer.OnePlusOneEvolutionary;
optimizer.InitialRadius = .0013;
optimizer.GrowthFactor = 2;
optimizer.MaximumIterations = 500;
metric = registration.metric.MattesMutualInformation;
for i = 1:length(dates)
    runPath = ['S:\Lab\',monkey,'\All Data\', monkey, '_',dates{i}, '\Imaging\'];
    
    sizeInfo = folderSizeTree(runPath);
    runFolders = find([sizeInfo.level{:}]==1);
    [~, runFolderInd] = max([sizeInfo.size{runFolders}]);
    runGreenPath = sizeInfo.name{runFolders(runFolderInd)};
    cd(runGreenPath);
    runNumber = runGreenPath(end-1:end);
    runGreen = double(imread([runGreenPath, '\green', runNumber, '_edited.bmp']));
    runGreen = rgb2gray(runGreen - mean(runGreen(:))); % mean-subtract and clip image to increase contrast
    
    [runGreen, ~, clipLow, clipHigh] = OIClipH(runGreen, 1, 1.5, ones(768,768));
    tform = imregtform(runGreen,refGreen, 'rigid', optimizer,metric);
    
    %%
    f = figure();
    imA = imshow(refGreen);
    hold on;
    title('Press "s" to save', 'FontSize', 20);
    imB = imshow(imwarp(runGreen,tform,'OutputView',imref2d(size(refGreen))));
    imB.AlphaData = 0.5;
    %%GUI
    Tinv = tform.invert.T;
    ss = Tinv(2,1);
    sc = Tinv(1,1);
    startTheta = atan2(ss,sc)*180/pi;
    startTx = tform.T(3,1);
    startTy = tform.T(3,2);
    imB.UserData = {runGreen, startTheta, startTx, startTy, tform};
    set(f, 'KeyPressFcn', @updateImage);
    waitfor(f);
    
end
function updateImage(~, evnt)
b = evnt.Key;
oi = evnt.Source.Children.Children(1).UserData{1};
theta = evnt.Source.Children.Children(1).UserData{2};
tx = evnt.Source.Children.Children(1).UserData{3};
ty = evnt.Source.Children.Children(1).UserData{4};
tform = evnt.Source.Children.Children(1).UserData{5};
switch b
    case 'rightarrow'
        tx = tx + 1;
    case 'leftarrow'
        tx = tx -1;
    case 'uparrow'
        ty = ty - 1;
    case 'downarrow'
        ty = ty +1;
    case 'comma'
        theta = theta + .5;
    case 'period'
        theta = theta - .5;
    case 'equal'
        evnt.Source.Children.Children(1).AlphaData =  evnt.Source.Children.Children(1).AlphaData + 0.05;
    case 'hyphen'
        evnt.Source.Children.Children(1).AlphaData = evnt.Source.Children.Children(1).AlphaData - 0.05;
    case 's'
        save('transformation_matrix', 'tform');
        disp('Saved!');
end
tform = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; tx ty 1;];

evnt.Source.Children.Children(1).UserData{2} = theta;
evnt.Source.Children.Children(1).UserData{3} = tx;
evnt.Source.Children.Children(1).UserData{4} = ty;
evnt.Source.Children.Children(1).UserData{5} = tform;
evnt.Source.Children.Children(1).CData = ...
imwarp(oi,affine2d(tform),'OutputView',imref2d(size(oi)));
end