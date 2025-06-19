animal = 'Gilligan';
date = '02_01_2019';
run = 'run00';
clipValues = [-0.30,0.30];
avgFrames = [50:55];

parentFolder = ['S:\Lab\',animal,'\All Data\',animal,'_', date,'\Imaging\',...
    run];

condDirs = dir([parentFolder,'\Results\Blockview\']);
condDirs = condDirs(3:end);

load([parentFolder,'\tform_nonrigid.mat']);

for c = 1:length(condDirs)
    frames = dir([parentFolder,'\Results\Blockview\', condDirs(c).name,...
        '\FF_1\LP5_HP550\Align_Offset\noClip\Frame*']);
    if(~exist([parentFolder,'\Results\Blockview\', condDirs(c).name,...
            '\FF_1\LP5_HP550\Align_Offset\IndFrames\AllClip\'], 'dir'))
        mkdir([parentFolder,'\Results\Blockview\', condDirs(c).name,...
            '\FF_1\LP5_HP550\Align_Offset\IndFrames\AllClip\']);
    end
     avgFrame = zeros(768,768);
    for f = 1:length(frames)
         frame = load([frames(f).folder,'\', frames(f).name]);
         tr = bspline_transform(O_trans,frame.saveFrame,Spacing,3);
        clipped = OIClipH2(tr,9,clipValues, []);
        imwrite((clipped-clipValues(1))./(clipValues(2)-clipValues(1)).*length(colormap(flipud(jet))), colormap(flipud(jet)),...
            [parentFolder, '\Results\Blockview\', condDirs(c).name,...
            '\FF_1\LP5_HP550\Align_Offset\IndFrames\AllClip\',frames(f).name(1:end-4),'_COLOR.png']);
        if(ismember(str2double(frames(f).name(7:8)),avgFrames))
            avgFrame = avgFrame + tr;
        end
    end
    avgFrame = avgFrame ./ (length(avgFrames));
    avgFrame = OIClipH2(avgFrame,9,clipValues, []);
    
    imwrite((avgFrame-clipValues(1))./(clipValues(2)-clipValues(1)).*length(colormap(flipud(jet))), colormap(flipud(jet)),...
        [parentFolder, '\Results\Blockview\', condDirs(c).name,...
        '\FF_1\LP5_HP550\Align_Offset\IndFrames\AllClip\Frame_AVG_COLOR.png']);
end