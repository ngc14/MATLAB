function result = fullsort(nsxFile, cfg)
%AUTOSPIKESORTBLACKROCK  Automatic, probe-wide spike sorting of ONE Blackrock
%   recording (.ns5, with its paired .nev if present) from an acute (non-chronic)
%   32-channel linear multielectrode, using SpikeInterface through MATLAB's
%   Python interface.
%
%   Design choices for this use case
%   - The whole probe is sorted jointly (multi-channel templates / neighbourhoods).
%     Nothing is sorted channel-by-channel; each unit gets an "extremum channel"
%     only as a label for downstream QC.
%   - Every recording is treated as independent: bad channels, noise levels,
%     reference and sorting are all estimated per recording. No cross-recording
%     unit tracking or chronic-drift assumptions are made.
%   - Global median reference (robust if some top contacts sit outside tissue).
%   - Nothing is discarded automatically. Units get a good/mua/noise label plus
%     quality metrics so you can QC them manually (e.g. in Offline Sorter).
%
%   One-time setup (before the first py.* call in this MATLAB session):
%     pyenv('Version','C:\path\to\python.exe');   % env with the packages below
%     % pip install "spikeinterface[full]" probeinterface neo mountainsort5
%     % If you hit DLL/MKL conflicts, use pyenv(...,'ExecutionMode','OutOfProcess').
%
%   Usage:
%     cfg = struct('contactPitchUm', 50, 'sorter', 'mountainsort5');
%     res = autoSpikeSortBlackrock('D:\data\rec001.ns5', cfg);
%
%     % Many acute recordings: sort each one independently
%     files = dir(fullfile(dataDir,'**','*.ns5'));
%     res = cell(numel(files),1);
%     for k = 1:numel(files)
%         res{k} = autoSpikeSortBlackrock(fullfile(files(k).folder, files(k).name), cfg);
%     end
%
%   Key cfg fields (all optional; see applyDefaults for the full list)
%     nChannels            32
%     channelIdx           1-based positions of the probe channels inside the NS5
%                          stream ([] -> 1:nChannels)
%     deviceChannelIndices 0-based; entry i = which selected channel feeds contact i,
%                          counting from contact 1 (y = 0). Use to undo headstage /
%                          adaptor remapping ([] -> 0:nChannels-1)
%     contactPitchUm       contact spacing along the shank (CHECK for your probe)
%     sorter               'mountainsort5' (CPU, default) | 'kilosort4' | 'spykingcircus2' ...
%     sorterParams         struct of overrides, e.g. struct('detect_threshold',6)
%     motionPreset         '' (off) or e.g. 'rigid_fast' if you see drift in a recording
%
%   Output: result struct (also saved to <outDir>/sorting_result.mat) with
%     units         table: unit_id, extremum channel, depth along probe (um),
%                   n_spikes, quality metrics, auto label
%     spikeSamples  {nUnits x 1} 0-based sample indices on the sorted timeline
%     spikeTimesSec spikeSamples / fs (seconds from start of the NS5 data)
%     segmentStartSamples / segmentLengthsSamples  (if the NS5 had paused segments,
%                   they were concatenated for sorting; use these to split again)

arguments
    nsxFile (1,:) char
    cfg (1,1) struct = struct()
end

cfg = applyDefaults(cfg);
chMap = loadChannelMap(nsxFile);
[~,cfg.deviceChannelIndices] = sort(chMap{2});
cfg.deviceChannelIndices = cfg.deviceChannelIndices-1;
if ~isfile(nsxFile), error('autoSpikeSort:noFile', 'File not found: %s', nsxFile); end
[folder, base] = fileparts(nsxFile);
outDir = cfg.outDir;
if isempty(outDir), outDir = fullfile(folder, [base '_si_sort']); end
if ~isfolder(outDir), mkdir(outDir); end

si = loadSpikeInterface();
siVersion = pyStr(py.getattr(si, '__version__'));
si.set_global_job_kwargs(pyargs('n_jobs', int64(cfg.nJobs), ...
    'chunk_duration', '1s', 'progress_bar', false,'pool_engine', 'loop'));
say('SpikeInterface %s | sorter: %s', siVersion, cfg.sorter);

%% 1. Read the recording and attach the probe ------------------------------------
rec_Raw = readBlackrock(si, nsxFile, cfg.streamId);
fs = rec_Raw.get_sampling_frequency();
if abs(fs - 30000) > 1
    warning('autoSpikeSort:fs', 'Sampling rate is %.1f Hz (expected 30000 for NS5).', fs);
end

allIds = cellfun(@pyStr, cell(py.list(rec_Raw.get_channel_ids().tolist())), 'UniformOutput', false);
sel = cfg.channelIdx;
if isempty(sel), sel = 1:cfg.nChannels; end
if numel(sel) ~= cfg.nChannels || max(sel) > numel(allIds)
    error('autoSpikeSort:channels', ...
        'Need %d channels; stream has %d. Set cfg.channelIdx.', cfg.nChannels, numel(allIds));
end
rec_Raw = rec_Raw.select_channels(py.list(allIds(sel)));
pif = py.importlib.import_module('probeinterface');
probe = pif.generate_linear_probe(pyargs('num_elec', int64(cfg.nChannels), ...
    'ypitch', cfg.contactPitchUm, 'contact_shapes', 'circle', ...
    'contact_shape_params', py.dict(pyargs('radius', cfg.contactRadiusUm))));
idx = cfg.deviceChannelIndices;
if isempty(idx), idx = 0:nCh-1; end
probe.set_device_channel_indices(py.numpy.array(int64(idx(:)')));
rec_Raw = rec_Raw.set_probe(probe);

% Paused NS5 recordings come in as several segments: concatenate for sorting.
nSeg = double(rec_Raw.get_num_segments());
segLens = zeros(1, nSeg);
for s = 1:nSeg
    segLens(s) = double(rec_Raw.get_num_samples(pyargs('segment_index', int64(s-1))));
end
segStarts = [0, cumsum(segLens(1:end-1))];
if nSeg > 1
    say('%d segments found; concatenating for sorting.', nSeg);
    rec_Raw = si.concatenate_recordings(si.split_recording(rec_Raw));
end

%% 2. Preprocess (per recording, no chronic assumptions) -------------------------
rec_filt = si.scale_to_uV(rec_Raw);
rec_filt = si.bandpass_filter(rec_filt, pyargs('freq_min', cfg.freqMin, 'freq_max', cfg.freqMax));
rec = rec_filt;%si.whiten(rec_filt, pyargs('mode', 'local', 'radius_um',2*cfg.contactPitchUm,'dtype', 'float32'));
[rec, badIds] = handleBadChannels(si, rec, cfg);
rec = applyReference(si, rec, cfg);
if ~isempty(cfg.motionPreset)
    say('Motion correction: %s', cfg.motionPreset);
    rec = si.correct_motion(rec, pyargs('preset', cfg.motionPreset));
end

binDir = fullfile(outDir, 'preprocessed');
if isfolder(binDir), rmdir(binDir, 's'); end
say('Caching preprocessed recording...');
%rec = rec.save(pyargs('folder', binDir, 'format', 'binary'));
usedIds = cellfun(@pyStr, cell(py.list(rec.get_channel_ids().tolist())), 'UniformOutput', false);
%% 3. Sort the whole probe jointly -----------------------------------------------
sortDir = fullfile(outDir, 'sorter');
cfg.sorterParams = si.get_default_sorter_params('spykingcircus');
sp = cfg.sorterParams.dictionary;
sp{"filter"} = false;
sp{"detect_sign"} = py.int(0);
sp{"detect_threshold"} = 5.5;
%sp{"whiten"} = false;
%sp{"max_threads_per_worker"} = py.int(8);
cfg.sorterParams = py.dict(sp);
kw = buildSorterKwargs(si, cfg.sorter, cfg.sorterParams);
%dkwargs = struct('installation_mode', 'dev','spikeinterface_folder_source','C:\ProgramData\anaconda3\envs\si\Lib\site-packages\spikeinterface');
args = [{'folder', sortDir, 'remove_existing_folder', true, 'verbose', true}, kw];
say('Running %s...', cfg.sorter);
fprintf('Step 1: Extracting and writing raw data payload natively via MATLAB...\n');

% Extract the raw matrix using the original file reader before any filters were applied
% (Replace 'si' with the variable name of your raw loaded dataset if different)
py_raw_traces = rec_Raw.get_traces(); 
mat_traces = int16(py.numpy.array(py_raw_traces));

% Save the raw binary file out to disk safely
binary_file_path = fullfile(outDir, 'raw_recording.dat');
fid = fopen(binary_file_path, 'wb');
fwrite(fid, mat_traces', 'int16'); 
fclose(fid);

fprintf('Step 2: Generating a standalone Python pipeline script...\n');

% 2. Write a Python script that applies filters AND runs the sorter inside the container environment
script_path = fullfile(outDir, 'run_pipeline.py');
fid = fopen(script_path, 'w');
fprintf(fid, 'import spikeinterface.core as si\n');
fprintf(fid, 'import spikeinterface.preprocessing as spr\n');
fprintf(fid, 'import spikeinterface.sorters as ss\n');
fprintf(fid, 'import probeinterface as pif\n');
fprintf(fid, 'if __name__ == "__main__":\n');
% Load raw binary file
fprintf(fid, '    rec = si.read_binary(r"/data/raw_recording.dat", sampling_frequency=%f, dtype="int16", num_channels=%d)\n',...
    fs,cfg.nChannels);
% Safely apply filtering/common reference inside the Linux container space
fprintf(fid, '    rec_filtered = spr.bandpass_filter(rec, freq_min=%f, freq_max=%f)\n',...
    cfg.freqMin, cfg.freqMax);
fprintf(fid, '    rec_cmr = spr.common_reference(rec_filtered, reference="global", operator="median")\n');
fprintf(fid, '    probe = pif.generate_linear_probe(num_elec=%d,ypitch=%d,contact_shapes="circle",contact_shape_params={"radius":%d})\n',...
    cfg.nChannels,cfg.contactPitchUm,cfg.contactRadiusUm);
fprintf(fid, '    probe.set_device_channel_indices([%d])\n',...
    cfg.deviceChannelIndices);
fprintf(fid, '    rec_cmr = rec_cmr.set_probe(probe)\n');
fprintf(fid, '    print("Running complete processing pipeline inside the container environment...")\n');
fprintf(fid, '    sorting = ss.run_sorter("spykingcircus", rec_cmr, output_folder="/data/spykingcircus_output", verbose=True, remove_existing_folder=True)\n');
fclose(fid);

fprintf('Step 3: Launching Docker container outside of MATLAB...\n');

% 3. Build a Docker command to run everything safely inside a SpikeInterface Linux environment
% This mounts your local output directory to the /data directory inside the container
docker_cmd = sprintf('docker run --rm -v "%s:/data" spikeinterface/spyking-circus-base:1.0.7 /bin/bash -c "pip install spikeinterface && python /data/run_pipeline.py"', outDir);
% Execute via native Windows terminal (Zero MATLAB-Python thread crash exposure)
[status, cmdout] = system(docker_cmd, '-echo');

% 4. Load the results back cleanly into MATLAB upon completion
if status == 0
    fprintf('Pipeline completed successfully! Reading sorting data back...\n');
    si_sorters = py.importlib.import_module('spikeinterface.sorters');
    sorting = si_sorters.read_sorter_folder(fullfile(outDir, 'spykingcircus_output'));
else
    error('Pipeline failed inside Docker container with error:\n%s', cmdout);
end




say('%d units found.', nUnits);

result = struct('file', nsxFile, 'fs', fs, 'siVersion', siVersion, 'sorter', cfg.sorter, ...
    'sorterKwargs', {kw}, 'channelIdsUsed', {usedIds}, 'badChannelIds', {badIds}, ...
    'segmentStartSamples', segStarts, 'segmentLengthsSamples', segLens, ...
    'units', table(), 'spikeSamples', {{}}, 'spikeTimesSec', {{}}, 'cfg', cfg);
if nUnits == 0
    warning('autoSpikeSort:none', 'Sorter returned no units.');
    save(fullfile(outDir, 'sorting_result.mat'), 'result', '-v7.3');
    return
end

%% 4. Per-unit analysis: waveforms, locations, quality metrics -------------------
anDir = fullfile(outDir, 'analyzer');
if isfolder(anDir), rmdir(anDir, 's'); end
an = si.create_sorting_analyzer(sorting, rec_filt, pyargs('format', 'binary_folder', ...
    'folder', anDir, 'sparse', true, 'method', 'radius', 'radius_um', cfg.sparsityRadiusUm));

tryCompute(an, 'random_spikes', 'method', 'uniform', 'max_spikes_per_unit', int64(cfg.maxSpikesPerUnit));
tryCompute(an, 'noise_levels');
tryCompute(an, 'waveforms', 'ms_before', .5, 'ms_after', 1.2);
tryCompute(an, 'templates');
tryCompute(an, 'spike_amplitudes');
tryCompute(an, 'unit_locations', 'method', 'center_of_mass');
tryCompute(an, 'quality_metrics', 'metric_names', py.list({'num_spikes', 'firing_rate', ...
    'presence_ratio', 'snr', 'isi_violation', 'amplitude_cutoff'}));

uidObjs = cell(py.list(sorting.get_unit_ids().tolist()));
ids = cellfun(@double, uidObjs(:));

% Spike trains (0-based sample indices) -------------------------------------------
spikeSamples = cell(nUnits, 1);
for u = 1:nUnits
    st = sorting.get_unit_spike_train(pyargs('unit_id', uidObjs{u}));
    spikeSamples{u} = double(py.array.array('d', st.astype('float64').tolist()))';
end
nSpk = cellfun(@numel, spikeSamples);

% Extremum channel (label only; sorting itself was probe-wide) ---------------------
extLabel = repmat({''}, nUnits, 1);
try
    exCh = si.get_template_extremum_channel(an, pyargs('peak_sign', 'neg', 'outputs', 'id'));
    extLabel = cellfun(@(x) pyStr(exCh.get(x)), uidObjs(:), 'UniformOutput', false);
    pr =an.recording.get_probe();
    contactIDs = an.recording.ids_to_indices(string(py.list(exCh.values)));
    contactIDs = double(contactIDs) + 1;
    all_positions = double(pr.contact_positions);
    all_positions = all_positions(contactIDs, :);
catch ME
    warning('autoSpikeSort:ext', 'Extremum channel failed: %s', ME.message);
end

% Depth along the probe axis (um from contact 1) ----------------------------------
depthUm = nan(nUnits, 1);
try
    loc = an.get_extension('unit_locations').get_data();
    nCol = double(loc.shape);
    L = reshape(double(py.array.array('d', loc.flatten().tolist())), nCol(2), [])';
    depthUm = L(:, 2);
catch ME
    warning('autoSpikeSort:loc', 'Unit locations unavailable: %s', ME.message);
end

T = table(ids, str2double(extLabel), string(extLabel), depthUm, nSpk, contactIDs', all_positions(:,2), 'VariableNames', ...
    {'unit_id', 'extremum_channel', 'extremum_channel_id', 'depth_um', 'n_spikes', 'contactID', 'position'});

try
    csvPath = fullfile(outDir, 'quality_metrics.csv');
    an.get_extension('quality_metrics').get_data().to_csv(csvPath, pyargs('index_label', 'unit_id'));
    Q = readtable(csvPath);
    T = outerjoin(T, Q, 'Keys', 'unit_id', 'Type', 'left', 'MergeKeys', true);
catch ME
    warning('autoSpikeSort:qm', 'Quality metrics unavailable: %s', ME.message);
end

% Automatic labels (nothing is removed; thresholds are in cfg) --------------------
snr = getVar(T, 'snr');
isi = getVar(T, 'isi_violations_ratio');
amp = getVar(T, 'amplitude_cutoff');
isGood = snr >= cfg.snrMin & (isnan(isi) | isi <= cfg.isiViolMax) & ...
         (isnan(amp) | amp <= cfg.ampCutoffMax) & T.n_spikes >= cfg.minSpikes;
isNoise = snr < cfg.noiseSnr;
lab = repmat("mua", height(T), 1);
lab(isNoise) = "noise";
lab(isGood) = "good";
T.auto_label = lab;

%% 5. Save -------------------------------------------------------------------------
result.units = T;
result.spikeSamples = spikeSamples;
result.spikeTimesSec = cellfun(@(x) x / fs, spikeSamples, 'UniformOutput', false);
save(fullfile(outDir, 'sorting_result.mat'), 'result', '-v7.3');
writetable(T, fullfile(outDir, 'unit_summary.csv'));
if ~cfg.keepPreprocessed && isfolder(binDir), rmdir(binDir, 's'); end
say('Done: %d good / %d mua / %d noise. Output: %s', sum(lab=="good"), sum(lab=="mua"), ...
    sum(lab=="noise"), outDir);
end

%% ------------------------------- local functions ---------------------------------

function cfg = applyDefaults(cfg)
d = struct( ...
    'nChannels', 32, ...
    'channelIdx', [], ...
    'deviceChannelIndices', [], ...
    'contactPitchUm', 100, ...          % CHECK against your probe's datasheet
    'contactRadiusUm', 15, ...
    'streamId', '5', ...               % '5' = NS5 (30 kHz)
    'freqMin', 300, 'freqMax', 6000, ...
    'badChannelMethod', 'std', 'badChannelThreshold', 5, ...
    'badChannelAction', 'none', ...  % 'remove' | 'interpolate' | 'none'
    'reference', 'global', ...         % 'global' | 'local'
    'motionPreset', '', ...
    'sorter', 'spykingcircus', ...
    'sorterParams', struct(), ...
    'sparsityRadiusUm', 150, ...
    'maxSpikesPerUnit', 150, ...
    'snrMin', 4, 'isiViolMax', 0.5, 'ampCutoffMax', 0.2, 'minSpikes', 50, 'noiseSnr', 2.5, ...
    'nJobs', 1, ...                    % >1 is safest with pyenv ExecutionMode OutOfProcess
    'outDir', '', ...
    'keepPreprocessed', true);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = d.(f{i}); end
end
end

function si = loadSpikeInterface()
pe = pyenv;
if strlength(string(pe.Version)) == 0
    error('autoSpikeSort:noPython', ...
        'No Python configured. Run pyenv(''Version'', ''<path to python>'') first.');
end
try
    si = py.importlib.import_module('spikeinterface.full');
catch ME
    error('autoSpikeSort:noSI', ['Could not import spikeinterface.full from %s.\n' ...
        'pip install "spikeinterface[full]" probeinterface neo mountainsort5\n%s'], ...
        pe.Executable, ME.message);
end
end

function rec = readBlackrock(si, f, streamId)
attempts = {{'stream_id', streamId}, {'stream_name', ['nsx' streamId]}};
err = [];
for k = 1:numel(attempts)
    a = attempts{k};
    try
        rec = si.read_blackrock(pyargs('file_path', f, a{1}, a{2}));
        return
    catch ME
        err = ME;
    end
end
rethrow(err);
end

function [rec, badIds] = handleBadChannels(si, rec, cfg)
badIds = {};
if strcmp(cfg.badChannelAction, 'none'), return; end
out = cell(si.detect_bad_channels(rec, pyargs('method', cfg.badChannelMethod, ...
    'std_mad_threshold', cfg.badChannelThreshold)));
badIds = cellfun(@pyStr, cell(py.list(out{1}.tolist())), 'UniformOutput', false);
if isempty(badIds), return; end
say('Bad channels: %s', strjoin(badIds, ', '));
if numel(badIds) > cfg.nChannels / 4
    warning('autoSpikeSort:manyBad', ...
        '%d of %d channels flagged bad - check this recording / the channel mapping.', ...
        numel(badIds), cfg.nChannels);
end
if strcmp(cfg.badChannelAction, 'interpolate')
    try
        rec = si.interpolate_bad_channels(rec, py.list(badIds));
        return
    catch ME
        warning('autoSpikeSort:interp', 'Interpolation failed (%s); removing instead.', ME.message);
    end
end
rec = rec.remove_channels(py.list(badIds));
end

function rec = applyReference(si, rec, cfg)
switch cfg.reference
    case 'global'
        rec = si.common_reference(rec, pyargs('reference', 'global', 'operator', 'median'));
    case 'local'
        r = [1.5 5] * cfg.contactPitchUm;   % annulus must exceed one pitch
        rec = si.common_reference(rec, pyargs('reference', 'local', 'operator', 'median', ...
            'local_radius', py.tuple({r(1), r(2)})));
    otherwise
        error('autoSpikeSort:ref', 'Unknown reference mode "%s".', cfg.reference);
end
end

function kw = buildSorterKwargs(si, name, user)
defaults = si.get_default_sorter_params(name);
p = struct();
switch lower(name)
    case 'kilosort4'
        p.do_CAR = false;               % already referenced
        p.nblocks = int64(0);           % no drift registration by default
    case 'spykingcircus2'
        p.apply_preprocessing = false;  % already preprocessed
end
user = user.struct;
f = fieldnames(user);
for i = 1:numel(f), p.(f{i}) = user.(f{i}); end
f = fieldnames(p);
kw = {};
for i = 1:numel(f)
    if logical(py.operator.contains(defaults, f{i}))
        kw = [kw, {f{i}, p.(f{i})}]; %#ok<AGROW>
    else
        warning('autoSpikeSort:param', ...
            '"%s" has no parameter "%s" in this SpikeInterface version; skipped.', name, f{i});
    end
end
end

function tryCompute(an, name, varargin)
try
    if isempty(varargin)
        an.compute(name);
    else
        an.compute(name, pyargs(varargin{:}));
    end
catch ME
    warning('autoSpikeSort:ext', 'Extension "%s" failed: %s', name, ME.message);
end
end

function v = getVar(T, name)
if ismember(name, T.Properties.VariableNames), v = T.(name); else, v = nan(height(T), 1); end
end

function s = pyStr(x)
s = char(py.str(x));
end

function say(fmt, varargin)
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), sprintf(fmt, varargin{:}));
end