clear all;
dateStart = datetime(2019,4,7, 'Format', 'MM_dd_y');
dateToday = datetime('now', 'Format', 'MM_dd_y');
dateArray = dateStart:dateToday;
monkey = 'Gilligan';
binSize=.01; % bin size in seconds
alignSegs = {'GoSignal', 'StartReach', 'StartGrasp'};
windows = {[-.15,0], [-.075,.075], [0,.15]};
sigma = 5; % smoothing window in ms
graspConds = {'Extra Small Sphere', 'Large Sphere'};
singleUnits = 1;

Kernel = (-3*sigma:3*sigma);
BinSize = length(Kernel);
Kernel = (-BinSize/2:BinSize/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

pVal = .05;
allUnits = 0;
cd(['\\univ.pitt.edu\sni\Gharbawie\Lab\',monkey,'\All Data\']);
for d = 1:length(dateArray)
    currName = [monkey,'_', datestr(dateArray(d), 'mm_dd_yyyy')];
    if(exist([currName, '\Physiology\Results'], 'dir'))
        sessionDir = dir(['\\pitt\sni\gharbawie\Lab\Gilligan\All Data\', monkey,'_',datestr(dateArray(d), 'mm_dd_yyyy'),'\Physiology\Results']);
        
        sessionDir = sessionDir(~[sessionDir.isdir]);
        [~,sortInd] = natsort({sessionDir.name});
        sessionDir = sessionDir(sortInd);
        fid  = fopen(['\\pitt\sni\gharbawie\Lab\',monkey,'\All Data\',monkey,'_'...
            ,datestr(dateArray(d), 'mm_dd_yyyy'),'\Physiology\',monkey,'_',datestr(dateArray(d), 'mm_dd_yyyy'),'_notesVProbe.txt'],'r');
        text = textscan(fid,'%s','Delimiter','','endofline','');
        text = text{1}{1};
        fid  = fclose(fid);
        tk = regexp(text,'Penetration #[\s\.=]+(\d+)','tokens');
        PN = str2double(tk{:});
        reachSig = nan(length(sessionDir),5);
        graspSig = nan(length(sessionDir),5);
        reachExSig = nan(length(sessionDir),5);
        graspExSig = nan(length(sessionDir),5);
        taskModulated = nan(length(sessionDir),5);
        assign = cell(length(sessionDir),5);
        totalUnits = 0;
        for f = 1:length(sessionDir)
            load([sessionDir(f).folder, '\', sessionDir(f).name]);
            units = size(sortedSpikeData.SpikeTimes,1);
            if(isfield(sortedSpikeData, 'label') || any(contains(who('-file',[sessionDir(f).folder, '\', sessionDir(f).name]), 'label')))
                if(~isfield(sortedSpikeData, 'label'))
                    labels = label;
                else
                    labels = sortedSpikeData.label;
                end
            end
            if(singleUnits)
                numUnits = find(arrayfun(@(a) strcmp(a,'s'), labels));
            else
                numUnits = 1:length(labels);
            end
            for u = 1:length(numUnits)
                condTrials = find(cellfun(@(a) contains(a,graspConds),sortedSpikeData.ArduinoData(:,1)));
                locTrials = find(cellfun(@(a) strcmp(a,'Close')||strcmp(a,'Far'),sortedSpikeData.ArduinoData(:,8)));
                currTrialInds = intersect(condTrials,  locTrials);
                currSpikes = sortedSpikeData.SpikeTimes(numUnits(u),currTrialInds);
                currSegs = sortedSpikeData.SegTimes(numUnits(u),currTrialInds);
                badTrials = cellfun(@(a) length(a)<50, currSpikes);
                if(sum(badTrials)<length(currTrialInds)/3)
                    totalUnits = totalUnits + 1;
                    for a = 1:length(alignSegs)
                        alignInd = find(strcmp(sortedSpikeData.ConditionSegments{1}, alignSegs(a)));
                        alignTimes = cellfun(@(a) getAlignedTimes(a,alignInd), currSegs, 'UniformOutput', false);
                        alignedSpikes = cellfun(@minus, currSpikes, alignTimes, 'UniformOutput', false);
                        secondsBeforeEvent = windows{a}(1);
                        secondsAfterEvent = windows{a}(2);
                        trialHists = cellfun(@(a) histcounts(a,...
                            secondsBeforeEvent-length(Kernel)/200:binSize:...
                            secondsAfterEvent+length(Kernel)/200)./binSize, alignedSpikes,'UniformOutput', false);
                        trialHistsSmooth = cellfun(@(a) conv(a,Kernel'), trialHists, 'UniformOutput', false);
                        trialHists = cellfun(@(a) a(length(Kernel):end-length(Kernel)+1), trialHistsSmooth,'UniformOutput', false);
                        trialHists = cellfun(@sum, trialHists);
                        %                 trialHists = cellfun(@(a) sum(a>=secondsBeforeEvent & a <=secondsAfterEvent), ...
                        %                     alignedSpikes);
                        totalSpikes(:,a) = [trialHists];
                    end
                    [~,reachSig(f,u)] = ttest2(totalSpikes(:,1), totalSpikes(:,2));
                    [~,graspSig(f,u)] = ttest2(totalSpikes(:,1), totalSpikes(:,3));
                    if(mean(totalSpikes(:,1)) < mean(totalSpikes(:,2)) || mean(totalSpikes(:,1)) < mean(totalSpikes(:,3)))
                        [~,reachExSig(f,u)] = ttest2(totalSpikes(:,2), totalSpikes(:,3), 'tail', 'right');
                        [~,graspExSig(f,u)] = ttest2(totalSpikes(:,3), totalSpikes(:,2), 'tail', 'right');
                    else
                        [~,reachExSig(f,u)] = ttest2(totalSpikes(:,2), totalSpikes(:,3), 'tail', 'left');
                        [~,graspExSig(f,u)] = ttest2(totalSpikes(:,3), totalSpikes(:,2), 'tail', 'left');
                    end                    
                    clear totalSpikes
                end
            end
        end
        correctedP = pVal/nchoosek(length(alignSegs),2);
        taskModulated = min(reachSig,graspSig);
        taskModulated(taskModulated>correctedP) = NaN;
        reachExSig(reachExSig>correctedP) = NaN;
        graspExSig(graspExSig>correctedP) = NaN;
        taskModulated(~isnan(reachExSig) | ~isnan(graspExSig)) = NaN;
        for f = 1:length(sessionDir)
            load([sessionDir(f).folder, '\', sessionDir(f).name]);
            units = size(sortedSpikeData.SpikeTimes,1);
            for u = 1:units
                if(~isnan(reachExSig(f,u)))
                    assign{f,u} = 'r';
                elseif(~isnan(graspExSig(f,u)))
                    assign{f,u} = 'g';
                elseif(~isnan(taskModulated(f,u)))
                    assign{f,u} = 't';
                else
                    assign{f,u} = 'n';
                end
            end
        end
        [num,txt,raw]=xlsread(['\\130.49.229.252\gharbawie\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx']);
        [~, siteCol] = find(contains(txt, 'Site #'));
        [row, ~ ] = find(num(:,siteCol)==PN);
        for a = 1:3
            [~, col] = find(contains(txt, [alignSegs{a}, ' Percentage']));
            if(contains(alignSegs{a}, 'Reach'))
                perc = 100*sum(sum(strcmp(assign, 'r')))/totalUnits;
            elseif(contains(alignSegs{a}, 'Grasp'))
                perc = 100*sum(sum(strcmp(assign, 'g')))/totalUnits;
            elseif(contains(alignSegs{a}, 'Go'))
                perc = 100*sum(sum(strcmp(assign, 't')))/totalUnits;
                [~, col] = find(contains(txt, ['Task', ' Percentage']));
            end
            col = char('A' + col-1);
            status = 0;
            if(isempty(row) || isempty(col)||row<1||col<7)
                disp('error');
            end
            while(~status)
            status = xlswrite(['\\130.49.229.252\gharbawie\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx'],perc,'Sheet1',[col, num2str(row+1)]);
            end
        end
        allUnits = allUnits+ totalUnits;
        disp(['Reach:' num2str(100*sum(sum(strcmp(assign, 'r')))/totalUnits), '; Grasp:',...
            num2str(100*sum(sum(strcmp(assign, 'g')))/totalUnits), '; Task:', num2str(100*sum(sum(strcmp(assign, 't')))/totalUnits)]);
        disp(['Reach:' num2str(sum(sum(strcmp(assign, 'r')))), '; Grasp:',...
            num2str(sum(sum(strcmp(assign, 'g')))), '; Task:', num2str(sum(sum(strcmp(assign, 't')))), '; Total:' num2str(totalUnits)]);
    end
end
function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end