%% Logging to disk from two cameras using Parallel Computing Toolbox
% Acquisition for each camera takes place in a separate worker
date = datetime('now', 'Format', 'MM_dd_yyyy');
frameRate = 120;
%% Create a parallel pool with two workers, one per camera
if isempty(gcp('nocreate'))
    parpool(1)
end

%% Disconnect from all cameras from main MATLAB process and workers
delete(imaqfind);
spmd(1)
    delete(imaqfind);
end

%% Create videoinput objects (one camera per worker)
spmd(1)
    % labBarrier ensures that the camera detection code is called 
    % by only one worker at a time.
    for idx = 1:numlabs
        if idx == labindex
            imaqreset
            
            % Detect cameras
            spinInfo = imaqhwinfo('mwspinnakerimaq');
            numCamerasFound = numel(spinInfo.DeviceIDs);
            fprintf('Worker %d detected %d cameras.\n', ...
                labindex, numCamerasFound);
        end
        labBarrier
    end
    
    cameraID = labindex;
    v = videoinput('mwspinnakerimaq', cameraID);
    s = v.Source;
    s.DeviceUserID = num2str(cameraID);
    s.DeviceLinkThroughputLimit = 290000000;  
    
    % Configure properties common for both cameras
    s.ExposureAuto = 'Off';
    s.ExposureTime = (1/(2*frameRate))*1000000;
    s.AcquisitionFrameRateEnable = 'True';
    s.AcquisitionFrameRate = frameRate; 
    % trigger configurations for acquiring frames
    s.LineSelector = 'Line3';
    v.FramesPerTrigger = 1000;

    v.LoggingMode = 'disk';
    diskLogger = VideoWriter(['C:\Users\Omar Lab\Documents\Monkey Training\Skipper_Macaque_152_17\Skipper_',...,...
    char(date),'\Videos\Camera_',num2str(cameraID)], 'Grayscale AVI');
    v.DiskLogger = diskLogger;
    diskLogger.FrameRate = s.AcquisitionFrameRate;
    % Configure properties that are camera specific
end

%% Configure manual triggering and wait for acquisition trigger
spmd(1)
    triggerconfig(v, 'manual');
    start(v);
end

%% Trigger acquisition
spmd(1)
    while(1)
        if(strcmp(s.LineStatus, 'True'))
            trigger(v);
            disp('trug');
            while(strcmp(s.LineStatus, 'True'))
            end
            stop(v);
            disp('done');
            break;
        end
    end
end

%% Display acquisition and logging status while logging
spmd(1)
    % Display number of frames acquired and logged while acquiring
    while strcmp(v.Logging, 'on')
        disp([v.FramesAcquired , v.DiskLoggerFrameCount])
        pause(1)
    end
    
    % Wait until acquisition is complete and specify wait timeout
    wait(v, 10);

    % Wait until all frames are logged
    while (v.FramesAcquired ~= v.DiskLoggerFrameCount) 
        pause(1);
    end
    disp([v.FramesAcquired v.DiskLoggerFrameCount]);
    
end

% Clean up
spmd(1)
    delete(imaqfind);
end