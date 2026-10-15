sessionName = 'Gilligan_06_05_2019';
savedDir = dir("S:\Lab\Gilligan\All Data\"+sessionName+"\Physiology\Results\*.mat");
[~,si] = natsort({savedDir.name});
savedDir = savedDir(si);
savedChs = []; savedUnits = string([]);
chMap = loadChannelMap(savedDir(1).folder);
for n = 1:length(savedDir)
    chSorted = load([savedDir(n).folder+"\"+savedDir(n).name]);
    chUnits = 1:size(chSorted.sortedSpikeData.SpikeTimes,1);
    label = chSorted.label;
    if(length(label)~=length(chUnits) & contains(fieldnames(chSorted.sortedSpikeData),'label'))
        label = chSorted.sortedSpikeData.label;
    end
    savedChs(end+1:end+length(label),:) = [repmat(n,length(label),1),(1:length(label))'];
    savedUnits(end+1:end+length(chUnits)) = chUnits;
end
%allUnitChannels = unique([savedChs,savedUnits'],'rows');
%validUnits = allUnitChannels(arrayfun(@str2num,allUnitChannels(:,2))>0,:);
NEV = openNEV(['C:\Users\ngc14\Desktop\autosort\',sessionName(10:end),'\',sessionName,'_0001-sorted.nev'],'read','nomat','nosave');
allElectrodes = NEV.Data.Spikes.Electrode; % Channel assignment for each spike
allUnits      = NEV.Data.Spikes.Unit;      % Sorted Unit ID (0 = unsorted, 1=Unit A, 2=Unit B, etc.)
allWaveforms  = NEV.Data.Spikes.Waveform;  % Waveform data (Rows = samples, Columns = individual spikes)
allTimestamps = NEV.Data.Spikes.TimeStamp; % Timestamps of the spikes

sortedRes = load([NEV.MetaTags.FilePath,'\',sessionName,'_0001_si_sort\sorting_result.mat']);
sortedTable = table(sortedRes.result.units.unit_id,sortedRes.result.units.extremum_channel-...
    min(cellfun(@str2num,sortedRes.result.channelIdsUsed))+1,sortedRes.result.units.depth_um,...
    sortedRes.result.units.num_spikes, sortedRes.result.units.auto_label,'VariableNames',{'UnitID','Electrode','Depth','Num_Spikes','Label'});


unitData = table('Size',[0,5],'VariableTypes',{'double','double','cell','cell','string'},'VariableNames',{'UnitID','Electrode','Timestamps','Waveforms','Label'});
for i = 1:size(savedDir, 1)
    elecID = i;%savedChs(i, 1);
    unitID = unique(allUnits((allElectrodes == elecID) & (allUnits > 0 & allUnits<255)));
    if(~isempty(unitID))
        spikeIdx = cell(1,length(unitID));
        for u = 1:length(unitID)
            spikeIdx{u} = find(allElectrodes == elecID & (allUnits == unitID(u)));
        end
        unitData{end+1:end+length(unitID),'UnitID'}     = unitID';
        unitData{end-length(unitID)+1:end,'Electrode'}  = repmat(elecID,length(unitID),1);
        unitData{end-length(unitID)+1:end,'Timestamps'} = cellfun(@(s) allTimestamps(s), spikeIdx, 'UniformOutput',false)';
        unitData{end-length(unitID)+1:end,'Waveforms'}  = cellfun(@(s) allWaveforms(:, s), spikeIdx, 'UniformOutput',false)';
        %unitData{end,'Label'} = allUnitChannels(i,3);
    end
end

