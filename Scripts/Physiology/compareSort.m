% Load NEV file and convert waveform amplitudes to microvolts
NEV = openNEV('C:\Users\ngc14\Desktop\autosort\Gilligan_05_17_2019_0001.nev', 'uV','nomat','nosave'); 
sortedRes = load([NEV.MetaTags.FilePath,'\',NEV.MetaTags.Filename,'_si_sort\sorting_results.mat']);
% 1. Extract raw flat arrays from the NEV structure
allElectrodes = NEV.Data.Spikes.Electrode; % Channel assignment for each spike
allUnits      = NEV.Data.Spikes.Unit;      % Sorted Unit ID (0 = unsorted, 1=Unit A, 2=Unit B, etc.)
allWaveforms  = NEV.Data.Spikes.Waveform;  % Waveform data (Rows = samples, Columns = individual spikes)
allTimestamps = NEV.Data.Spikes.TimeStamp; % Timestamps of the spikes

% 2. Identify unique channel-unit combinations (ignoring unsorted unit 0)
uniqueCombinations = unique([allElectrodes; allUnits]', 'rows');
sortedUnitsMask    = uniqueCombinations(:, 2) > 0; 
validUnits         = uniqueCombinations(sortedUnitsMask, :); % Column 1: Electrode, Column 2: Unit ID

% 3. Group and extract waveforms/timestamps per sorted unit
unitData = struct();
for i = 1:size(validUnits, 1)
    elecID = validUnits(i, 1);
    unitID = validUnits(i, 2);

    % Find logical indices matching this specific sorted unit on this electrode
    spikeIdx = (allElectrodes == elecID) & (allUnits == unitID);

    % Save sorted values into a structured format
    unitData(i).Electrode  = elecID;
    unitData(i).UnitID     = unitID;
    unitData(i).Timestamps = allTimestamps(spikeIdx);
    unitData(i).Waveforms  = allWaveforms(:, spikeIdx); % Each column is a waveform snippet
end
