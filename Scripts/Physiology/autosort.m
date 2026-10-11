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
%%
parfor ch_idx = 1:numChannels
    ch_spks = electrodes==ch_idx;
    ch_wfs = wfs(:,ch_spks)';
    try 
        [~, score, ~] = pca(double(ch_wfs));
        features = score(:, 1:min(3, size(ch_wfs,1)))'; 
        clusterRes = uint8(isosplit5(features));
        units{ch_idx} = [find(ch_spks);clusterRes];
        fprintf('Channel %d complete: Sorted.\n', ch_idx);
    catch ME
        warning('Failed on channel %d: %s', ch_idx, ME.message);
    end
end
for ch_idx = 1:numChannels
    wNEV.Data.Spikes.Units(units{ch_idx}(1,:)) = units{ch_idx}(2,:);
end
[~, sortIdx] = sort(wNEV.Data.Spikes.TimeStamp);
wNEV.Data.Spikes.Electrode = wNEV.Data.Spikes.Electrode(sortIdx);
wNEV.Data.Spikes.TimeStamp = wNEV.Data.Spikes.TimeStamp(sortIdx);
wNEV.Data.Spikes.Unit      = wNEV.Data.Spikes.Unit(sortIdx);
wNEV.Data.Spikes.Waveform = wNEV.Data.Spikes.Waveform(:, sortIdx);

saveNEV(wNEV,[nevFile(1:end-4),'-autosort.nev']);
fprintf('TOTAL PROCESSING TIME: %.2f seconds\n', toc);