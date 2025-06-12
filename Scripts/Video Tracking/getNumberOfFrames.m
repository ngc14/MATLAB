currDir = dir('S:\Lab\ngc14\Track_Videos\02_23_2018');

[numFrames, numFramesSorted] = deal([]);
index = 0;

for i = 1:size(currDir,1)
    dot = regexp(currDir(i).name, '[.]');
    underscores = regexp(currDir(i).name, '[_]');
    % only process .avi files in current directory
    if(~isdir(currDir(i).name) && strcmp(currDir(i).name(dot:end), '.avi'))
        index = index+1;
        vidObj = VideoReader([currDir(i).folder, '\' ,currDir(i).name]);
        % if underscores is empty, videos have not been renamed yet
        if(isempty(underscores))
            % video number
            numFrames(index,1) = str2double(currDir(i).name(1:dot-1));
            numFrames(index,2) = vidObj.NumberOfFrames;
        else
            % trial number
            numFrames(index,1) =  str2double(currDir(i).name(1:underscores(1)-1));
            % block number
            numFrames(index,2) =  str2double(currDir(i).name(underscores(2)+1:underscores(3)-1));
            numFrames(index,3) = vidObj.NumberOfFrames;
        end
    end
end

[~,Idx]=sort(numFrames(:,1));
numFramesSorted=numFrames(Idx,:);

save('numFrames', 'numFramesSorted');

