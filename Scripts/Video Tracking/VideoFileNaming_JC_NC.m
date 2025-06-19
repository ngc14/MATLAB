%{
This program renames video files created by the camera script by matching
them with the appropriate information in an Excel document. Video files
should be named 1.avi, 2.avi, etc. Renamed video files are copied into a
new folder in the video directory, original files are left untouched in
case of error.

User parameters:
infoFile: String for location of Excel file that contains desired trial info.
secondFile: Optional string for second Excel file if one is needed.
videoDirectory: String for location of folder that contains all of the
video files. This folder needs to be empty of any files besides the video
files.
%}

% NOTE: NEED TO ADD ABILITY TO USE ignoreConditions, deleteConditions,
% deleteUnsuccessful, AND deleteFalseStart
%% Parameters
% infoFile should be the Arduino output text file  (.txt)
date = '05_20_2020';
monkey = 'Skipper';
conds = {'Extra Small Sphere', 'Large Sphere', 'Photocell', 'Rest'};

% Initialize variables.
delimiter = '\t';
startRow = 3;

% Read columns of data as text:
formatSpec = '%*s%s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

if(strcmp(monkey, 'Gilligan'))
    mainDir = ['S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\Gilligan_',...
        date];
    secondDir = ['S:\Lab\Gilligan\All Data\Gilligan_', date];
else
    mainDir = ['S:\Lab\BehaveBackup\Monkey_Training\Skipper_Macaque_152_17\Skipper_',...
        date];
    secondDir = ['S:\Lab\Skipper\All Data\Skipper_', date];
end
infoFile = [secondDir,'\Arduino\', monkey, '_',date,'.txt'];

cameras = dir([secondDir,'\Videos\']);
% cameras = dir([mainDir,'\Videos\']);

cameras = cameras(contains({cameras.name},'Camera'),:);
% if(~exist(secondDir,'dir'))
%     mkdir(secondDir);
% end
if(~exist([secondDir,'\Arduino'],'dir'))
    mkdir([secondDir,'\Arduino']);
    copyfile([mainDir,'\', monkey, '_',date,'.txt'], infoFile);
end
if(~exist([secondDir,'\Videos'],'dir'))
    mkdir([secondDir,'\Videos']);
end
for i = 1:length(cameras)
    if(~exist([secondDir,'\Videos\',cameras(i).name],'dir'))
        mkdir([secondDir,'\Videos\',cameras(i).name]);
    end
    if(length(dir([secondDir,'\Videos\',cameras(i).name]))<3)
        copyfile([mainDir,'\Videos\',cameras(i).name], [secondDir,'\Videos\',cameras(i).name]);
    end
end



extension = '.avi'; % Extension of saved videos
%% Load trial informaiton from Arduino output text file (.txt)
% Open the text file.
fileID = fopen(infoFile,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
raw = [dataArray{1:end-1}];

% remove lines that are not task related (pause, headers)
raw = raw(ismember(raw(:,3), conds),:);

% Create output variable
trialInfo = raw;

% Remove false start trials. Videos are not recorded for false starts.
falseStartIdx=strcmp(trialInfo(:,10),'False Start');
trialInfo=trialInfo(~falseStartIdx,:);

for c = 1:length(cameras)
    renamedDirectory = [secondDir,'\Videos\',cameras(c).name,'\Renamed\'];
    if(~exist(renamedDirectory, 'dir'))
        mkdir(renamedDirectory);
    end
    videoDirectory = [secondDir, '\Videos\', cameras(c).name];
    videoNum = 1;
    for i = 1:length(trialInfo)
        % i directly maps to trialInfo with removed instances
        info = trialInfo(i, :);
        
        block = num2str(cell2mat(info(1)));
        trial = num2str(cell2mat(info(2)));
        condition = cell2mat(info(3));
        if(~strcmp(condition, 'Rest'))
            oldname = [videoDirectory, '\', num2str(videoNum), extension];
            newname = [renamedDirectory, block, '_', num2str(videoNum), '_',trial,'_',condition, extension];
            
            % Save renamed file in new directory
            copyfile(oldname, newname);
        end
        videoNum = videoNum + 1;
    end
end