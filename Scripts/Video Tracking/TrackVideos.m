clear all; close all;
%% User Input
datadriver = 'S:\'; % Data disk name
datafolder = 'Lab\ngc14\'; % Data folder name on data disk, results will be saved in outputfolder

outputdriver = 'S:\'; % Data disk name
outputfolder = 'Lab\ngc14\'; % Output folder name on data disk, results will be saved here

date = '12_11_2017';
dataname = ['Track_Videos\', date, '\'];
segments = {'Grasp', 'Reach'};
frameRate = 120;
%%
resultFolder = [outputdriver, outputfolder, dataname, 'Results\'];
if ~isdir(resultFolder)
    mkdir(resultFolder)
end
vidFolder = [datadriver, datafolder, dataname, 'Renamed\'];

infoFold = 'S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\';
infoFile = [infoFold, 'Gilligan_', date, '\Gilligan_', date,'.txt'];

% Read columns of data as text:
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(infoFile,'r');
% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'HeaderLines',1,...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
% get rid of inavalid lines (pauses)
validLines = cellfun(@(x) size(x,2) > 0, dataArray{:,3});
dataArray = cellfun(@(x) x(validLines), dataArray, 'UniformOutput', false);
% information from text file
trials = cellfun(@(x) str2double(x), dataArray{:,3});
conditions = dataArray{:,4};
reactionTimes = cellfun(@(x) str2double(x), dataArray{:,5}, 'UniformOutput', false);
reachDurations = cellfun(@(x) str2double(x), dataArray{:,6}, 'UniformOutput', false);
liftDurations = cellfun(@(x) str2double(x), dataArray{:,7}, 'UniformOutput', false);
holdDurations = cellfun(@(x) str2double(x), dataArray{:,8}, 'UniformOutput', false);

% parallel pooling
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = gcp();
    addAttachedFiles(poolobj,{'Track_2.m'});
end

videos = dir(vidFolder);
parfor i = 1:size(videos,1)
    underscores = regexp(videos(i).name, '[_]');
    extension = regexp(videos(i).name, '[.]');
    % only analyze .avi files in current directory
    if(~videos(i).isdir && strcmp(videos(i).name(extension:end), '.avi'))
        for s = 1:size(segments,2)
            
            trialNum = str2double(videos(i).name(underscores(2)+1:underscores(3)-1));
            cond = videos(i).name(underscores(end)+1:extension-1);
            % find trial index with corresponding information based on the
            % video trial name and its condition (if paused in the middle
            % of the session, trial name will reset to 1 creating
            % duplicates, adding condition is another filter to ensure
            % correct information is being used.
            trialIndex = find(trials == trialNum & strcmp(conditions, cond));
            
            currFold = [resultFolder,  cond, '\', segments{s}, '\'];
            if(~exist([currFold, 'Track\'], 'dir'))
                mkdir([currFold, 'Track\']);
            end
            if(~exist([currFold, 'Speed\'], 'dir'))
                mkdir([currFold, 'Speed\']);
            end
            
            vidName = videos(i).name(1:underscores(end)-1);
            
            if(~exist([currFold, 'Track\', vidName, '.fig'], 'file') ...
                    || ~exist([currFold, 'Speed\', vidName, '.fig'], 'file'))
                disp([cond '_' vidName '_' segments{s}]);
                % start and end frame calculation for reach segmentation
                if(strcmp(segments{s}, 'Reach'))
                    start = frameRate * (reactionTimes{trialIndex}/1000);
                    finish = start + (frameRate * (reachDurations{trialIndex}/1000));
                % start and end frame calculation for grasp segmentation
                % (not photocell)
                elseif(strcmp(segments{s}, 'Grasp') && ~strcmp(cond, 'Photocell'))
                    start = (frameRate * (reactionTimes{trialIndex}/1000)) + ...
                        (frameRate * (reachDurations{trialIndex}/1000)) + 1;
                    finish = start + (frameRate * (liftDurations{trialIndex}/1000));
                % start and end frame calculation for grasp segmentation of
                % photocell
                elseif(strcmp(segments{s}, 'Grasp') && strcmp(cond, 'Photocell'))
                    start = (frameRate * (reactionTimes{trialIndex}/1000)) + ...
                        (frameRate * (reachDurations{trialIndex}/1000)) + 1;
                    finish = start + (frameRate * (holdDurations{trialIndex}/1000));
                else
                    assert(0 > 1, 'Need to define case for segmentation other than Reach or Grasp');
                end
                % get track and speed plots for current video and segment
                track_2(videos(i), ceil(start), ceil(finish), currFold, vidName);
            end
        end
    end
end