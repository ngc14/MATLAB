directory = dir('\\130.49.229.252\gharbawie\Lab\Gilligan\All Data');
finalT = {};
for f = 3:size(directory,1)
    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New'];
    if(exist(resultFolder, 'dir'))
        results = dir(resultFolder);
        for r = 3:size(results,1)
            matF = matfile([results(r).folder, '\', results(r).name]);
            disp([results(r).folder, '\', results(r).name]);
            data = matF.sortedEMGData;
            muscle = data.Muscle;
            muscle = muscle(~isspace(muscle));
            arduino = data.ArduinoData;
            arduino = arduino(:,1);
            conds = data.Conditions;
            badTrials = EMG_getBadTrials(data.EMGData, data.SegTimes, data.SampleRate);
            for c = 1:length(conds)
                num(c) = sum(cellfun(@(a) strcmp(a,conds{c}), arduino) & ~badTrials');
            end
            if(isempty(finalT))
                finalT(end+1,:) = {muscle, num};
            elseif(~any(cellfun(@(s) ~isempty(strfind(muscle, s)), finalT(:,1))))
                finalT(end+1,:) = {muscle, num};
            else
                fStr = cellfun(@(s) max(strfind(muscle, s),0), finalT(:,1), 'UniformOutput', false);
                matches = ~cellfun(@isempty, fStr);
                bestMatch = min([fStr{matches}]);
                [fStr{~matches}] = deal(0);
                fInd = cellfun(@(s) s==bestMatch, fStr);
                finalT(fInd,2) = {finalT{fInd,2} + num};
            end
        end
    end
end
save(['\\130.49.229.252\gharbawie\Lab\Gilligan\Sorted EMG Data\Results\count'], 'finalT');