monkey = 'Gilligan';
date = '01_17_2019';
frames_select = [1:30];

dataPath = ['S:\Lab\',monkey,'\All Data\',monkey,'_',date, '\Imaging\'];
runs = dir([dataPath,'run*']);

figure();
hold on;
startPlot = 1;
for r = 1:length(runs)
    currRun = [runs(r).folder, '\', runs(r).name];
    blocks = dir([currRun,'\*.BLK']);
    [~,indx] = natsort({blocks.name});
    blocks = blocks(indx);
    
    anapar = OIHeadRead([blocks(1).folder,'\',blocks(1).name], 'v');
    W = anapar.FrameWidth;
    H = anapar.FrameHeight;
    NFrames = anapar.FramesPerStim;
    
    trial_value = zeros(1,length(blocks));
    trial_std = zeros(1,length(blocks));
    parfor b = 1:length(blocks)
        all_frames = zeros(H,W,NFrames);
        all_frames(:,:,:) = OIReadStim([blocks(b).folder,'\',blocks(b).name],0,'v');
        frames = all_frames(:,:,frames_select);
        frame_value = zeros(1,length(frames_select));
        for f = 1:length(frames_select)
            currFrame = frames(:,:,f);
            frame_value(f) = median(currFrame(:));
        end
        trial_value(b) = median(frame_value);
        trial_std(b) = std(frame_value);
    end
    shadedErrorBar(startPlot:startPlot + length(blocks)-1, trial_value,trial_std);
    startPlot = startPlot + length(blocks);
end