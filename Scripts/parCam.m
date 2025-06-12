%% Logging to disk from two cameras using Parallel Computing Toolbox
% Acquisition for each camera takes place in a separate worker
date = datetime('now', 'Format', 'MM_dd_yyyy');
frameRate = 120;
videoNum = 1;
dirPath = ['X:\Monkey Training\Skipper_Macaque_152_17\Skipper_',char(date),'\Videos\'];
%dirPath = ['X:\Monkey Training\Gilligan_Macaque_110-16\Gilligan_',char(date),'\Videos\'];
spinInfo = imaqhwinfo('mwspinnakerimaq');
numCamerasFound = numel(spinInfo.DeviceIDs);
%% Create a parallel pool with two workers, one per camera
if isempty(gcp('nocreate'))
    parpool(numCamerasFound);
end

%% Disconnect from all cameras from main MATLAB process and workers

spmd(numCamerasFound)
    delete(imaqfind);
end

%% Create videoinput objects (one camera per worker)
spmd(numCamerasFound)
    % labBarrier ensures that the camera detection code is called
    % by only one worker at a time.
    for idx = 1:numlabs
        if idx == labindex
            imaqreset
            % Detect cameras
            fprintf('Worker %d detected %d cameras.\n', ...
                labindex, numCamerasFound);
        end
        labBarrier
    end
    
    cameraID = labindex;
    v = videoinput('mwspinnakerimaq', cameraID, 'Mono8');
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
    
end

disp('Ready to capture videos.');
%%

while(videoNum<1000)
    spmd(numCamerasFound)
        diskLogger = VideoWriter([dirPath,'Camera_',s.DeviceUserID,'\',num2str(videoNum)], 'Grayscale AVI');
        diskLogger.FrameRate = s.AcquisitionFrameRate;
        v.LoggingMode = 'disk';
        v.DiskLogger = diskLogger;
        
        start(v);
        labBarrier;
        while(1)
            if(strcmp(s.LineStatus,'True'))
                while(strcmp(s.LineStatus,'True'))
                end
                labBarrier;
                stop(v);
                break;
            end
        end
        %% Display acquisition and logging status while logging
        % Wait until acquisition is complete and specify wait timeout
        wait(v,12, 'logging');
        disp([v.FramesAcquired v.DiskLoggerFrameCount]);
    end
    labBarrier;
    disp(videoNum);
    videoNum = videoNum + 1;
end
