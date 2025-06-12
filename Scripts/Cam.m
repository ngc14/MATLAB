%% Logging to disk from two cameras using Parallel Computing Toolbox
% Acquisition for each camera takes place in a separate worker
date = datetime('now', 'Format', 'MM_dd_yyyy');
frameRate = 120;
videoNum = 1;
dirPath = ['X:\Monkey Training\Skipper_Macaque_152_17\Skipper_',char(date),'\Videos\'];
%dirPath = ['X:\Monkey Training\Gilligan_Macaque_110-16\Gilligan_',char(date),'\Videos\'];
spinInfo = imaqhwinfo('mwspinnakerimaq');
numCamerasFound = numel(spinInfo.DeviceIDs);
startCameraName = 'Arm 1';
%% Disconnect from all cameras from main MATLAB process and workers

for c = 1:numCamerasFound
    delete(imaqfind);
end

%% Create videoinput objects (one camera per worker)
for c = 1:numCamerasFound
    
    cameraID = c;
    v{c} = videoinput('mwspinnakerimaq', cameraID, 'Mono8');
    s{c} = v{c}.Source;
    if ~exist([dirPath,'Camera_',s{c}.DeviceUserID],'dir')
        mkdir([dirPath,'Camera_',s{c}.DeviceUserID]);
    end
    s{c}.DeviceLinkThroughputLimit = 300000000;
    
    % Configure properties common for both cameras
    s{c}.ExposureAuto = 'Off';
    s{c}.ExposureTime = 1000;
    s{c}.AcquisitionFrameRateEnable = 'True';
    s{c}.AcquisitionFrameRate = frameRate;
    % trigger configurations for acquiring frames
    if(strcmp(startCameraName, s{c}.DeviceUserID))
        triggerInd = c;
        s{c}.LineSelector = 'Line3';
        s{c}.V3_3Enable = 'True';
    else
        s{c}.LineSelector = 'Line3';
    end
    v{c}.FramesPerTrigger = Inf;
    triggerconfig(v{c}, 'manual');
    
end

disp('Ready to capture videos.');
%%
while(videoNum<1000)
    for c = 1:numCamerasFound
        diskLogger{c} = VideoWriter([dirPath,'Camera_',s{c}.DeviceUserID,'\',num2str(videoNum)], 'Grayscale AVI');
        diskLogger{c}.FrameRate = s{c}.AcquisitionFrameRate;
        v{c}.LoggingMode = 'disk';
        v{c}.DiskLogger = diskLogger{c};
        start(v{c});
    end
    while(1)
        if(strcmp(s{triggerInd}.LineStatus, 'True'))
            for c = 1:numCamerasFound
                trigger(v{c});
            end
            while(strcmp(s{triggerInd}.LineStatus, 'True'))
            end
            for c = 1:numCamerasFound
                stop(v{c});
            end
            break;
        end
    end
    %% Display acquisition and logging status while logging
    % Wait until acquisition is complete and specify wait timeout
    for c = 1:numCamerasFound
        wait(v{c},12, 'logging');
    end
    %disp([v.FramesAcquired v.DiskLoggerFrameCount]);
    disp(videoNum);
    videoNum = videoNum + 1;
end
