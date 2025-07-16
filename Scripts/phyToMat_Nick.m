function [] = phyToMat_Nick(params)
date = '04_19_2019';
monkey = 'Gilligan';
phyDataDir = ['\\univ.pitt.edu\sni\Gharbawie\Lab\', monkey,'\All Data\', monkey '_', date,'\Physiology\'];
%read phy results
spikeTimes = readNPY([phyDataDir 'spike_times.npy']); %spikeTimes are uint64
spikeClusters = readNPY([phyDataDir 'spike_clusters.npy']);
spikeTemplates = readNPY([phyDataDir 'spike_templates.npy']); spikeTemplates = spikeTemplates+1;
clusterStatus =  tdfread([phyDataDir 'cluster_group.tsv'],'\t');
clusterID = clusterStatus.cluster_id; clusterStatus = clusterStatus.group;

load([phyDataDir, 'rez.mat'],'rez');

%do this if you saved U in rez, otherwise read rez.template2channel
[~, template2channel]       = max(gather(rez.U(:,:,1)), [], 1);

nSpikes = numel(spikeTemplates); spikeChannels = nan(nSpikes,1);
for iSpike = 1:nSpikes, spikeChannels(iSpike) = template2channel(spikeTemplates(iSpike));
end %find out each spikes' channel (only way to work with current code)

chansWithSpikes = unique(spikeChannels);

%% get events
%ignore this. This part sets the trial length for the spikes struct and
%reads the alignment events

% eventsFile = [dataDirServer monkey session 'x_events.mat'];
% load(eventsFile,['periOn' array(1:2) '_30k'],'isSuccess','eventsVidLat','stimPattern');
% eval(['periOn_30k = periOn' array(1:2) '_30k;']);
% 
% keepTrials = isSuccess&(stimPattern==3)'; %only do successful no stim trials
% periOn_30k = periOn_30k(keepTrials);
% nTrials = numel(periOn_30k);
msBeforePeriOn = 1500; msAfterPeriOn = 3500; %for instructed delay

spikesT = -msBeforePeriOn:1:msAfterPeriOn; nTimePoints = numel(spikesT);

%% read phy and sort phy data channel, unit, trial
nChannels = 32; uCount = 1;

for iChannel = 1:nChannels
    if ismember(iChannel,chansWithSpikes)%if template was found in this channel
        chanClusters = unique(spikeClusters(spikeChannels==iChannel)); %infer the channel clusters from spike channels
        clusterCount = 0;
        for iCluster = chanClusters'
            if strcmp(strtrim(clusterStatus(clusterID==iCluster,:)),'good') %if cluster not noise, make unit
                clusterCount = clusterCount+1;
                u(uCount).channel = iChannel;
                
                uSpikeTimes = spikeTimes(spikeClusters==iCluster);
                u(uCount).nSpikes = numel(uSpikeTimes);
                
%                 u(uCount).spikeTimes = cell(nTrials,1);
%                 u(uCount).spikesLogical = false(nTrials,nTimePoints);
%                 for iTrial = 1:nTrials %periON is the alignment event I use
%                     iPeriOn = periOn_30k(iTrial);
%                     iStart = iPeriOn - 30*msBeforePeriOn;
%                     iEnd = iPeriOn + 30*msAfterPeriOn;
%                     
%                     iSpikes = double(uSpikeTimes(uSpikeTimes>iStart & uSpikeTimes<iEnd));
%                     if ~isempty(iSpikes)
%                         iSpikes = round((iSpikes-iPeriOn)/30); %align and convert to ms
%                         u(uCount).spikeTimes{iTrial} = unique(iSpikes); %sometimes there will be duplicates because of short trials
%                         u(uCount).spikesLogical(iTrial,ismember(spikesT,iSpikes)) = true; % no need to treat duplicates here
%                     end%if spikes in trial
%                 end%for nTrials
                
                %keep a list also (unit, channel, unit in channel)
                unitNum(uCount,1) = uCount;
                unitNum(uCount,2) = iChannel;
                unitNum(uCount,3) = clusterCount;
                
                uCount = uCount+1;
            end% if cluster not noise
        end%for channel clusters
    else
        disp(['no template for ch' num2str(iChannel)]);
    end%if channel has units
end%for nChannels

nUnits = uCount-1;

%also make nUnits x nTrials x nTimePoints for Taka
spikes3D = false(nUnits,nTrials,nTimePoints);
for iUnit = 1:nUnits, spikes3D(iUnit,:,:) = u.spikesLogical;
end

%save
disp([session ' ' array ' sorted ' num2str(nUnits) ' units']);
% saving file
save([dataDirServer monkey session 'x' array '_units.mat'],'u','nUnits','spikesT','unitNum','spikes3D'); %save to SERVER

end%function

