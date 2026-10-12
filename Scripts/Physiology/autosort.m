nevFile = 'C:\Users\ngc14\Desktop\TODO-Spike_Sort\autosort\Gilligan_2022_06_23_0001.nev';
nb_before = 15;
nb_after = 37;
%% 2. Read Metadata Header using NPMK
NEV = openNEV(nevFile, 'report','nosave','noread','nomat');
numChannels = length([NEV.ElectrodesInfo.ElectrodeID]);
%% 3. Parallel Processing Loop Over Channels
fprintf('Processing %d channels in parallel...\n', numChannels); tic;
wNEV = openNEV(nevFile, 'report','read', 'nosave', 'nomat');
nsxFile = nevFile; nsxFile(end-3:end) = '.ns5';
sf = wNEV.MetaTags.TimeRes;
electrodes = wNEV.Data.Spikes.Electrode;
wfs = wNEV.Data.Spikes.Waveform;
ts = wNEV.Data.Spikes.TimeStamp;
units = cell(1,numChannels);
%%
mat_spikes = cell(1, numChannels);
mat_timestamps = cell(1, numChannels);
mat_units = cell(1, numChannels);
groupCh = 8;
parfor ch_idx = 1:numChannels
    ch_spks = find(electrodes==ch_idx);
    ch_wfs = wfs(:,ch_spks)';
    ch_ts = double(ts(ch_spks));
    [sorted_ts, sort_idx] = sort(ch_ts, 'ascend');
    try
        [~, score, ~] = pca(double(ch_wfs));
        features = score(:, 1:min(3, size(ch_wfs,1)))';
        clusterRes = isosplit5(features);

        mat_spikes{ch_idx} = double(ch_wfs(sort_idx, :));
        mat_timestamps{ch_idx} = double(sorted_ts(:) ./ sf); % Digital ticks to Seconds
        mat_units{ch_idx} = double(clusterRes(sort_idx)');
        fprintf('Channel %d complete: Sorted.\n', ch_idx);
    catch ME
        warning('Failed on channel %d: %s', ch_idx, ME.message);
    end
    if(isempty(mat_spikes{ch_idx}))
        mat_spikes{ch_idx} = zeros(0, waveform_length, 'double');
        mat_timestamps{ch_idx} = zeros(0, 1, 'double');
        mat_units{ch_idx} = zeros(0, 1, 'double');
    end
end
for ch_idx = 1:numChannels
    if(mod(ch_idx,groupCh)==0)
        gInd = [1:groupCh]+floor((ch_idx-1)/groupCh)*groupCh;
        spikes = mat_spikes(gInd)';
        timestamps = mat_timestamps(gInd)';
        units = mat_units(gInd)';
        save([nevFile(1:end-4),'-autosort',num2str(gInd),'.mat'], 'spikes', 'timestamps', 'units','-v7');
    end
end
spikes = mat_spikes';
timestamps = mat_timestamps';
units = mat_units';
save([nevFile(1:end-4),'-autosort.mat'], 'spikes', 'timestamps', 'units', '-v6');
fprintf('TOTAL PROCESSING TIME: %.2f seconds\n', toc);
