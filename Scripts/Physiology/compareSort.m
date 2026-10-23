si = py.importlib.import_module('spikeinterface.full');
sessionName = 'Gilligan_07_16_2019';
NEV = openNEV(['C:\Users\ngc14\Desktop\autosort\',sessionName,'\',sessionName,'_0001-sorted.nev'],'read','nomat','nosave');
sortedDir = dir(['S:\Lab\Gilligan\All Data\',sessionName,'\Physiology\Results\*.mat']);
analyzer = si.load(['C:\Users\ngc14\Desktop\autosort\',sessionName,'\',sessionName,'_0001_si_sort\analyzer']);
[~,so] = natsort({sortedDir.name});
sortedDir = sortedDir(so);
savedChs = []; savedLabels = strings(1, 0);

wf = double(analyzer.get_extension('waveforms').get_data());
autoSortedRes = load([NEV.MetaTags.FilePath,'\',sessionName,'_0001_si_sort\sorting_result.mat']);
autoSortedTable = table(autoSortedRes.result.units.unit_id,autoSortedRes.result.units.extremum_channel-...
    min(cellfun(@str2num,autoSortedRes.result.channelIdsUsed))+1,autoSortedRes.result.units.depth_um,...
    autoSortedRes.result.units.num_spikes, autoSortedRes.result.units.auto_label,...
    arrayfun(@(n) wf(1+(n*100):(n+1)*100,:,1), 0:length(double(analyzer.unit_ids))-1, 'UniformOutput',false)',...
    cellfun(@squeeze,num2cell(double(analyzer.get_extension("templates").get_data()),[2,3]),'UniformOutput',false),'VariableNames',...
    {'UnitID','Electrode','Depth','Num_Spikes','Label','Waveforms','MeanWaveforms'});
for n = 1:length(sortedDir)
    label = "";
    chSorted = load([sortedDir(n).folder+"\"+sortedDir(n).name]);
    chUnits = 1:size(chSorted.sortedSpikeData.SpikeTimes,1);
    label = chSorted.label;
    if(length(label)~=length(chUnits))
        if(contains(fieldnames(chSorted.sortedSpikeData),'label'))
        label = chSorted.sortedSpikeData.label;
        end
    end
    savedChs(end+1:end+length(chUnits),:) = [chUnits',repmat(n,length(chUnits),1)];
    savedLabels(end+1:end+length(label)) = label;
end
allElectrodes = NEV.Data.Spikes.Electrode; % Channel assignment for each spike
allUnits      = NEV.Data.Spikes.Unit;      % Sorted Unit ID (0 = unsorted, 1=Unit A, 2=Unit B, etc.)
allWaveforms  = NEV.Data.Spikes.Waveform;  % Waveform data (Rows = samples, Columns = individual spikes)
allTimestamps = NEV.Data.Spikes.TimeStamp; % Timestamps of the spikes

unitData = table('Size',[0,6],'VariableTypes',{'double','double','double','cell','cell','string'},'VariableNames',{'UnitID','Electrode','Depth','Timestamps','Waveforms','Label'});
chMap = loadChannelMap(sortedDir(1).folder);
depthMap = [0:length(chMap{2})-1].*100;
for i = 1:size(sortedDir, 1)
    elecID = i;%savedChs(i, 1);
    unitID = double(unique(allUnits((allElectrodes == elecID) & (allUnits > 0 & allUnits<255))));
    if(~isempty(unitID))
        spikeIdx = cell(1,length(unitID));
        for u = 1:length(unitID)
            spikeIdx{u} = find(allElectrodes == elecID & (allUnits == unitID(u)));
        end
        unitData{end+1:end+length(unitID),'UnitID'}     = unitID';
        unitData{end-length(unitID)+1:end,'Electrode'}  = repmat(elecID,length(unitID),1);
        unitData{end-length(unitID)+1:end,'Depth'} = depthMap(chMap{2}==elecID)+(unitID-1)';
        unitData{end-length(unitID)+1:end,'Timestamps'} = cellfun(@(s) allTimestamps(s), spikeIdx, 'UniformOutput',false)';
        unitData{end-length(unitID)+1:end,'Waveforms'}  = cellfun(@(s) allWaveforms(:, s), spikeIdx, 'UniformOutput',false)';
    end
end
unitData{:,'Label'} = savedLabels'; figure(); tiledlayout(2,1);
nexttile(); hold on; title("Manual Sort (# spikes:"+string(java.text.DecimalFormat().format(sum(cellfun(@length,unitData.Timestamps))))+")");
bar(unitData.Depth,cellfun(@length,unitData.Timestamps)); xlim([0 3200]);
nexttile(); hold on; title("Auto Sort (# spikes: "+string(java.text.DecimalFormat().format(sum(autoSortedTable.Num_Spikes)))+")");
bar(autoSortedTable.Depth,autoSortedTable.Num_Spikes); xlim([0 3200]);
colors = [1 .5 .8; 1 1 0; 1 0 0; 0 1 1; 0 .8 1; .7 1 0; 0 1 0;];
figure(); tiledlayout(6,6,'TileIndexing','rowmajor');
for h = 1:length(sortedDir)
    nexttile(h); hold on;
    title(num2str(h));
    chMap = round(unitData.Depth/100);
    plot(1:52,zeros(52,1),'Color', 'k');
    grid off
    xticks(17)
    xticklabels(0);
    xlim([0 52]);
    if(any(chMap==h))
        chUnits = find(chMap==h);
        for a = 1:length(chUnits)
            data = double(cell2mat(unitData{chUnits(a),'Waveforms'}));
            shadedErrorBar(1:52, mean(data,2)', std(data,0,2)', 'lineprops', {'color',colors(a,:)},'patchSaturation',0.3)
        end
    end
    yL = ylim;
    if(abs(yL(end) - yL(1))>50)
        yticks(yL(1):50:yL(end));
        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)),yL(1):50:yL(end)));
    else
        yticks(yL(1):25:yL(end));
        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)), yL(1):25:yL(end)));
    end
    ylim(yL)
end
figure(); tiledlayout(7,6);
for h = 1:length(sortedDir)
    nexttile(h); hold on;
    title(num2str(h));
    chMap = round(autoSortedTable.Depth/100);
    plot(1:51,zeros(51,1),'Color', 'k');
    grid off
    xticks(17)
    xticklabels(0);
    xlim([0 51]);
    if(any(chMap==h))
        chUnits = find(chMap==h);
        for a = 1:length(chUnits)
            data = double(cell2mat(autoSortedTable{chUnits(a),'Waveforms'}))';
            shadedErrorBar(1:51, mean(data,2)', std(data,0,2)', 'lineprops', {'color',colors(a,:)},'patchSaturation',0.3)
        end
    end
    yL = ylim;
    if(abs(yL(end) - yL(1))>50)
        yticks(yL(1):50:yL(end));
        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)),yL(1):50:yL(end)));
    else
        yticks(yL(1):25:yL(end));
        yticklabels(arrayfun(@(s) string(sprintf('%.0f',s)), yL(1):25:yL(end)));
    end
    ylim(yL)
end