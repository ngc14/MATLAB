nevFile = 'C:\Users\ngc14\Desktop\TODO-Spike_Sort\autosort\Gilligan_2022_06_23_0001.nev'; 
matOutputFile = 'C:\Users\ngc14\Desktop\TODO-Spike_Sort\autosort\Gilligan_2022_06_23_0001-autosort.mat';
nb_before = 15;
nb_after = 37;
waveform_length = nb_before + nb_after;
%% 2. Read Metadata Header using NPMK
fprintf('Reading NEV file headers...\n');
NEV = openNEV(nevFile, 'report','nosave','noread','nomat');
channelIDs = [NEV.ElectrodesInfo.ElectrodeID];
numChannels = length(channelIDs);
sf = NEV.MetaTags.TimeRes;
%% 3. Parallel Processing Loop Over Channels
% Initialize cell array structures to hold final outputs per channel
spikes_cell = cell(1, numChannels);
timestamps_cell = cell(1, numChannels);
units_cell = cell(1, numChannels);
fprintf('Processing %d channels in parallel...\n', numChannels);
tic;
parfor ch_idx = 1:numChannels
    ch_id = channelIDs(ch_idx);
    % Open a localized independent reader handle inside the worker thread
    wNEV = openNEV(nevFile, ['c:',num2str(ch_id)], 'read', 'nosave', 'nomat');
    % Check if any spikes exist on this channel
    if isempty(wNEV.Data.Spikes.TimeStamp)
        spikes_cell{ch_idx} = zeros(0, waveform_length, 'double');
        timestamps_cell{ch_idx} = zeros(0, 1, 'double');
        units_cell{ch_idx} = zeros(0, 1, 'int32');
        continue;
    end
    % Get original timestamps and convert from digital ticks to Seconds
    raw_ts = double(wNEV.Data.Spikes.TimeStamp(:));
    sec_ts = raw_ts ./ sf;
    % Get waveforms: shape is [Samples x Num_Spikes] from openNEV
    raw_wf = double(wNEV.Data.Spikes.Waveform); 
    % Transpose to [Num_Spikes x Samples] matching WaveClus/OFS expectations
    wf_clean = raw_wf'; 
    num_spikes = size(wf_clean, 1);
    %% Run Core WaveClus Clustering Functions directly in RAM
    try
        par_struct = struct();
        par_struct.w_before = nb_before;
        par_struct.w_after = nb_after;
        par_struct.inputs = 1;
        par_struct.features = 'wav';  % Use Wavelet Transform features
        par_struct.scales = 4;
        
        % Extract Wavelet Features Output shape: [Num_Spikes x Num_Features]
        features = Get_features(wf_clean, par_struct);
        % Run Superparamagnetic Clustering (SPC) 
        [cluster_labels, ~] = Do_clustering(features, par_struct);        
        % WaveClus outputs row configurations; transpose to [Num_Spikes x 1]
        unit_assignments = int32(cluster_labels(:));
        
        spikes_cell{ch_idx} = wf_clean;
        timestamps_cell{ch_idx} = sec_ts;
        units_cell{ch_idx} = unit_assignments;
        
        fprintf('Channel %d complete: Sorted %d spikes cleanly.\n', ch_id, num_spikes);
    catch ME
        warning('Failed processing channel %d: %s', ch_id, ME.message);
        spikes_cell{ch_idx} = zeros(0, waveform_length, 'double');
        timestamps_cell{ch_idx} = zeros(0, 1, 'double');
        units_cell{ch_idx} = zeros(0, 1, 'int32');
    end
end
total_time = toc;
spikes = spikes_cell';
timestamps = timestamps_cell';
units = units_cell';
save(matOutputFile, 'spikes', 'timestamps', 'units', '-v7');
fprintf('TOTAL PROCESSING TIME: %.2f seconds\n', total_time);
