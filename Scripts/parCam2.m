%% Logging to disk from two cameras using Parallel Computing Toolbox
% Acquisition for each camera takes place in a separate worker
date = datetime('now', 'Format', 'MM_dd_yyyy');
frameRate = 120;
videoNum = 1;
dirPath = ['X:\Monkey Training\Skipper_Macaque_152_17\Skipper_',char(date),'\Videos\'];
%dirPath = ['X:\Monkey Training\Gilligan_Macaque_110-16\Gilligan_',char(date),'\Videos\'];
spinInfo = imaqhwinfo('mwspinnakerimaq');
numCamerasFound = numel(spinInfo.DeviceIDs);

%% Disconnect from all cameras from main MATLAB process and workers
delete(imaqfind);

%% Create videoinput objects (one camera per worker)

v = videoinput('mwspinnakerimaq', i, 'Mono8');
s = v.Source;
if ~exist([dirPath,'Camera_',s.DeviceUserID],'dir')
    mkdir([dirPath,'Camera_',s.DeviceUserID]);
end

s.DeviceLinkThroughputLimit = 300000000;

% Configure properties common for both cameras
s.ExposureAuto = 'Off';
s.ExposureTime = 1000;
s.AcquisitionFrameRateEnable = 'True';
s.AcquisitionFrameRate = frameRate;
% trigger configurations for acquiring frames
s.LineSelector = 'Line3';
v.FramesPerTrigger = Inf;

%% Configure manual triggering and wait for acquisition trigger
triggerconfig(v,'manual');
disp('Ready to capture videos.');
%%
while(videoNum<1000)
    diskLogger = VideoWriter([dirPath,'Camera_',s.DeviceUserID,'\',num2str(videoNum)], 'Grayscale AVI');
    diskLogger.FrameRate = s.AcquisitionFrameRate;
    v.LoggingMode = 'disk';
    v.DiskLogger = diskLogger;
    
    start(v);
    button = 0;
    while(1)
        if((button==10))
            trigger(v);
            while(button<2480)
                button = button + 1;
            end
            stop(v);
            break;
        end
        button = button + 1;
    end
    %% Display acquisition and logging status while logging
    % Wait until acquisition is complete and specify wait timeout
    wait(v,12,'logging');
    %disp([v.FramesAcquired v.DiskLoggerFrameCount]);
    disp(videoNum);
    videoNum = videoNum + 1;
end


