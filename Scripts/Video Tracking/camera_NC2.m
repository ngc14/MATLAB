
% This program sets camera parameters such as exposure and generates video
% files for received data. Videos are named 1.avi, 2.avi, etc.
%
% User parameters:
% exposureTime: The amount of time spent gathering each frame, in microseconds.
% This value will directly affect frame rate.
% vidPath: The directory where generated video files will be stored. NOTE:
% need to stop proIMgram if this is not empty
% blackwhite: Set this to true to record video in black & white, false for
% color.

clear all;
imaqreset;
spinnaker = imaqhwinfo('mwspinnakerimaq');
numCameras = length(spinnaker.DeviceIDs);
date = datetime('now', 'Format', 'MM_dd_yyyy');
triggerLine = 1;

%% Parameters
frameRate = 120;
% vidPath = strcat('O:\Lab\ngc14\Track_Videos\',...

%%
numVids = 0;
videos = {};
diskLoggers = {};
vidPath = strcat('C:\Users\Omar Lab\Documents\Monkey Training\Skipper_Macaque_152_17\Skipper_',...,...
    char(datetime('now', 'Format', 'MM_dd_yyyy')),'\Videos\');
for i = 1:numCameras
    videos{end+1} = videoinput('mwspinnakerimaq',i,'Mono8');
    if(~exist([vidPath,'Camera_',num2str(i)], 'dir'))
        mkdir([vidPath,'Camera_',num2str(i)]);
    end
end

src = cellfun(@(a) getselectedsource(a), videos,'UniformOutput', false);
for i = 1:numCameras
        src{i}.DeviceLinkThroughputLimit = 250000000;
    src{i}.ExposureAuto = 'Off';
    src{i}.ExposureTime = (1/(2*frameRate))*1000000;
    src{i}.AcquisitionFrameRateEnable = 'True';
    src{i}.AcquisitionFrameRate = 100;
      
    % trigger configurations for acquiring frames
    src{i}.LineSelector = 'Line3';
    triggerconfig(videos{i}, 'manual');
    videos{i}.FramesPerTrigger = Inf;
end
disp('Ready to capture video.');
%% Session Videos
while(numVids < 1500)
    numVids = numVids + 1;
    for i = 1:numCameras
        diskLoggers{end+1} = VideoWriter([vidPath,'Camera_',num2str(i),'\',num2str(numVids)], 'Uncompressed AVI');
        diskLoggers{i}.FrameRate = src{i}.AcquisitionFrameRate;
    end
    
    start(cellfun(@(a) a, videos));
    while(1)
        if(strcmp(src{triggerLine}.LineStatus, 'True'))
            trigger(cellfun(@(a) a, videos));
            while(strcmp(src{triggerLine}.LineStatus, 'True'))
                disp(videos{1}.FramesAcquired);
            end
            stop(cellfun(@(a) a, videos));
            break;
        end
    end
    
    data = cellfun(@(a) getdata(a, a.FramesAvailable), videos, 'UniformOutput', false);
    dataSZ = cellfun(@(a) size(a), data, 'UniformOutput', false);
    open(cellfun(@(a) a, diskLogger));
    
    parfor n = 1:numCameras
        for ii = 1:size(dataSZ{n}(end))
            writeVideo(diskLogger{n}, uint8(data{n}(:,:,:,ii)));
        end
    end
    
    disp(numVids);
    close(cellfun(@(a) a, diskLogger));
end
