
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
date = datetime('now', 'Format', 'MM_dd_yyyy');

%% Parameters
exposureTime = 2500;
% vidPath = strcat('O:\Lab\ngc14\Track_Videos\',...
vidPath = strcat('\\130.49.229.252\gharbawie\Lab\ngc14\Track_Videos\',...
    char(datetime('now', 'Format', 'MM_dd_yyyy')),'\');
blackwhite = false;
%%

trialNum = 1;
numVids = 0;
if ~isdir(vidPath)
    mkdir(vidPath)
end

if ~blackwhite
    vid = videoinput('gige', 1, 'BayerRG8');
else
    vid = videoinput('gige', 1, 'Mono8');
end

src = getselectedsource(vid);
src.ExposureTimeAbs = exposureTime;
src.AcquisitionFrameRateAbs = 120;
src.PacketSize = 8000;

%% Average Video
vid.FramesPerTrigger = 240;
k = input('Press return when ready to capture average video.', 's');
while(~strcmp(k,''))
    k = input('', 's');
end
start(vid);
while(vid.FramesAvailable < vid.FramesPerTrigger)
    
end
stop(vid);
close;
data = getdata(vid, vid.FramesAvailable);

diskLogger = VideoWriter([vidPath, 'average'], 'Uncompressed AVI');
diskLogger.FrameRate = src.AcquisitionFrameRateAbs;
open(diskLogger);
for ii = 1:size(data,4)
    writeVideo(diskLogger, data(:,:,:,ii));
end
close(diskLogger);

% trigger configurations for acquiring frames
triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
src.TriggerMode = 'On';
src.TriggerSelector = 'FrameStart';
src.TriggerSource = 'Line1';
src.TriggerActivation = 'LevelHigh';
vid.TriggerRepeat = 500;
vid.FramesPerTrigger = Inf;

disp('Ready to capture video.');
%% Session Videos
while(numVids < 1500)
    numVids = numVids + 1;
    filename = [vidPath, num2str(numVids)];
    diskLogger = VideoWriter(filename, 'Uncompressed AVI');
    diskLogger.FrameRate = src.AcquisitionFrameRateAbs;
    %diskLogger.Quality = 25;
    start(vid);
    framesAcq = 0;
    lastFrame = uint64(0);
    
    while (vid.FramesAcquired < 2000)
        newFramesAcq = vid.FramesAcquired;
        if (newFramesAcq > framesAcq)
            lastFrame = tic;
            framesAcq = newFramesAcq;
            % disp(framesAcq);
        end
        % time since the last frame acquired has exceeded 2x the frame rate
        % timing, thus video has stopped acquiring frames -> trial is done,
        % save video
        if (newFramesAcq && toc(lastFrame) > 0.03)
            break;
        end
    end
    
    stop(vid);
    data = getdata(vid, vid.FramesAvailable);
    dataSZ = size(data);
    open(diskLogger);
    
    for ii = 1:dataSZ(end)
        if ~blackwhite
            writeVideo(diskLogger, data(:,:,:,ii));
        else
            writeVideo(diskLogger, data(:,ii));
        end
    end
    
    disp(numVids);
    close(diskLogger);
end
