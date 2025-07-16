monkey = 'Skipper';                   % animal name
%dates = {'12_03_2018', '12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
%dates = {'02_01_2019','02_06_2019'};
dates = {'10_21_2020', '10_30_2020', '11_07_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
dates = {'10_21_2020', '10_30_2020', '11_07_2020','11_09_2020'}
%runs = {[1], [0], [0], [0], [1], [0], [0], [1], [0]};
runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0],[0],[0]};

RIGID = false;
LPk = 15;
HPk = 250;
frames = '51-56';
for d = 1:length(dates)
    runFolder = ['S:\Lab\',monkey,'\All Data\', monkey, '_',dates{d}, '\Imaging\run0',num2str(runs{d})];
    %     sizeInfo = folderSizeTree(sessionPath);
    %     runFolders = find([sizeInfo.level{:}]==1);
    %     [~, runFolderInd] = max([sizeInfo.size{runFolders}]);
    %     runFolder = sizeInfo.name{runFolders(runFolderInd)};
    if(RIGID)
        tform = load([runFolder,'\transformation_matrix.mat']);
        tform = tform.tform;
    else
        tform = load([runFolder,'\tform_nonrigid.mat']);
    end
    tTestMain = dir([runFolder, '\Results\Ttest_nsc15_V3\']);
    for frameRanges = 1:1:70
        frameRange = max(1,frameRanges-2):min(70,frameRanges+2);
        frames = [num2str(frameRange(1)),'-',num2str(frameRange(end))];
        tTestPath = tTestMain(contains({tTestMain.name}, ['Frames', frames]));
        
        for f = 1:size(tTestPath,1)
            folderPath = [tTestPath(f).folder, '\',tTestPath(f).name];
            folderPath = dir([folderPath,'\LP',num2str(LPk),'HP',num2str(HPk),'\Align_Offset\FF_1']);
            folderPath = folderPath(contains({folderPath.name},'.bmp'));
            if ~isempty(folderPath)
                if ~exist([folderPath(f).folder, '\Warped\'], 'dir')
                    mkdir(folderPath(f).folder, '\Warped\');
                end
                for t = 1:size(folderPath,1)
                    image = imread([folderPath(t).folder, '\', folderPath(t).name]);
                    if(RIGID)
                        warped = imwarp(image, affine2d(tform),'OutputView',imref2d(size(image)));
                        imwrite(warped,[folderPath(t).folder, '\Warped\', folderPath(t).name]);
                    else
                        disp([frames, folderPath(t).name]);
                        warped = bspline_transform(tform.O_trans,image,tform.Spacing,3);
                        imwrite(logical(warped),[folderPath(t).folder, '\Warped\NONRIGID_', folderPath(t).name]);
                    end
                end
            end
        end
    end
end