fullDir = dir('\\130.49.229.252\gharbawie\Lab\Gilligan\All Data');
fullDir = fullDir(3:end);
corrTable = array2table(zeros(0,6), 'VariableNames', {'FileName', 'D4D5Extensor', ...
    'Extensor','DistalFlexor', 'ProximalRadialFlexor', 'ProximalUlnarFlexor'});
for f = 1:size(fullDir,1)
    resultFolder = [fullDir(f).folder, '\',fullDir(f).name,'\EMG\Results_New\'];
    if(exist(resultFolder, 'file'))
        currDir = dir(resultFolder);
        for m = 1:size(currDir,1)
            if(contains(currDir(m).name, 'Extensor'))
                disp(currDir(m).name);
                D4D5Ex = EMG_CorrelateSession('D4D5Extensor', ...
                    [currDir(m).folder, '\', currDir(m).name]);
                Ex = EMG_CorrelateSession('Extensor', ...
                    [currDir(m).folder, '\', currDir(m).name]);
                corrTable = [corrTable; {currDir(m).name, D4D5Ex, Ex, 0, 0, 0}];
            elseif(contains(currDir(m).name, 'Flexor'))
                disp(currDir(m).name);
                DistalFl = EMG_CorrelateSession('DistalFlexor', ...
                    [currDir(m).folder, '\', currDir(m).name]);
                ProxRadFl = EMG_CorrelateSession('ProximalRadialFlexor', ...
                    [currDir(m).folder, '\', currDir(m).name]);
                ProxUlnFl = EMG_CorrelateSession('ProximalUlnarFlexor', ...
                    [currDir(m).folder, '\', currDir(m).name]);
                corrTable = [corrTable; {currDir(m).name, 0, 0, DistalFl, ProxRadFl, ProxUlnFl}];
            end
        end
    end
end
save('\\130.49.229.252\gharbawie\Lab\Gilligan\Sorted EMG Data\Results\Models\CorrTable', 'corrTable');