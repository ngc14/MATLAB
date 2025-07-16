%% Initialize variables.
function [bodyParts, coordinates] = loadCoordinates(filePath)
delimiter = ',';
startRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filePath,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
headerInfo = cellfun(@(a) a(1:2), dataArray, 'UniformOutput', false);
headerInfo = [headerInfo{:}];
numericalValues = cell2mat(cellfun(@(a) str2double(a(3:end)), dataArray, 'UniformOutput', false));
bodyParts = unique(headerInfo(1,:));
bodyParts = bodyParts(bodyParts~="" & bodyParts~="bodyparts");
fclose(fileID);
coordinates = cell(1,length(bodyParts));

for b = 1:length(bodyParts)
    xInds = find(strcmp(headerInfo(1,:),bodyParts(b)) & strcmp(headerInfo(2,:),"x"));
    yInds = find(strcmp(headerInfo(1,:),bodyParts(b)) & strcmp(headerInfo(2,:),"y"));
    coordinates{b} = [numericalValues(:,xInds), numericalValues(:,yInds)];
end

end