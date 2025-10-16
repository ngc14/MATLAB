conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
taskWindow =repmat({{[0, 0]}},1,length(conditions));
params = PhysRecording(string(conditions),.01,.15,-6,15,containers.Map(conditions,...
    {"StartReach","StartReach","StartReach","GoSignal"}));
phaseNames = categorical(["Baseline", "Go", "Reach", "Hold","Reward"],'Ordinal',true);
phaseAlignmentPoints = {["GoSignal","GoSignal","StartReach","StartHold","StartReward"],...
    ["GoSignal","GoSignal","StartReach","StartHold","StartReward"],...
    ["GoSignal","GoSignal","StartReach","StartHold","StartReward"],...
    ["GoSignal","GoSignal","StartReplaceHold","StartReward"]};
phaseWinSz = .2;
phaseWindows = repmat({{[-phaseWinSz, 0],[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz, 0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions)-1);
phaseWindows(end+1) = {{[-phaseWinSz, 0],[0, phaseWinSz], [-phaseWinSz,0],...
    [-phaseWinSz*(3/4),phaseWinSz*(1/4)]}};
pVal=0.05;
MIN_BLOCKS_FOR_UNIT = 13;
savePath = "S:\Lab\ngc14\Working\PMd\Task_Units\";
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","PMd");
clear rawSpikes
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
    tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(mean(t(:,1)+1,'omitnan'),params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t) mean(cell2mat(cellfun(@(a,n) max(1,mean(...
    a(:,n(~isnan(n)):n(~isnan(n))+(3/params.binSize),:),[2,3],'omitnan')),...
    p,t,'UniformOutput',false)),2,'omitnan'),num2cell([siteTrialPSTHS{:}],2),num2cell([allCondSegs{:}],2),"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) cellfun(@(p)p./repmat(nb,1,1,size(p,3)),...
    cp,'UniformOutput',false),num2cell([siteTrialPSTHS{:}],2),normBaseline,'Uniformoutput', false);
normPSTH = vertcat(normPSTH{:});
%normPSTH = horzcat(siteTrialPSTHS{:});
%%
trialInfo = cellfun(@(c) cellfun(@(t) t(strcmp(t(:,1),c),:),siteTrialInfo,'UniformOutput',false)',conditions,'UniformOutput',false);
siteTrialSegs = cellfun(@(c) cellfun(@(n) NaN(size(n,1),length(maxSegL)), c, 'UniformOutput',false), trialInfo,'UniformOutput',false);
trialInfo = cellfun(@(c) vertcat(c{:}),num2cell(horzcat(trialInfo{:}),2),'UniformOutput',false);
for j = 1:size(siteTrialSegs,2)
    for i = 1:length(siteTrialInfo)
        siteTrialSegs{j}{i}(:,condSegMappedInds{j}) = siteSegs{j}{i}{1};
    end
end
siteTrialSegs=cellfun(@(v) vertcat(v{:}),num2cell(cat(2,siteTrialSegs{:}),2),'UniformOutput',false);
allPSTHS = cellfun(@(r) cell2mat(reshape(r,1,1,[])),num2cell(normPSTH,2),'UniformOutput',false);
for s = 1:length(trialInfo)
    badTrials = all(isnan(siteTrialSegs{s}),2);
    siteTrialSegs{s} = siteTrialSegs{s}(~badTrials,:);
    trialInfo{s} = trialInfo{s}(~badTrials,:);
    allPSTHS{s} = allPSTHS{s}(:,:,~badTrials);
end
infoTable = table();
infoTable.SiteNo = "SiteNo" + table2array(siteDateMap(:,'Site'));
infoTable.TaskModulated = cellfun(@(c) any(horzcat(c{:}),2), num2cell(cat(2,tUnit{:}),2), 'UniformOutput',false);
infoTable.Condition = cellfun(@(t) string(t(:,1)), trialInfo,'UniformOutput',false);
infoTable.SegTimes = siteTrialSegs;
infoTable.PSTHS = cellfun(@(r,i)(i./i).*r,allPSTHS,infoTable.TaskModulated,'UniformOutput',false);
%%
condPSTHS = cell(1,length(conditions)-1);
segConds = cell(1,length(conditions)-1);
for s = 1:height(infoTable)
    currSite = infoTable(s,:);
    for c = 1:length(conditions)
        condInds = strcmp(currSite.Condition{:},conditions(c));
        condTimes = mean(currSite.SegTimes{1}(condInds,:),1,'omitnan');
        avgPSTHS = mean(currSite.PSTHS{1} .* (currSite.TaskModulated{:}./currSite.TaskModulated{:}),3,'omitnan');
        phaseInds = cellfun(@(pa,pw) findBins(condTimes(strcmp(maxSegL,pa))+pw,params.bins), phaseAlignmentPoints{c},phaseWindows{c},'UniformOutput',false);
        phaseAvgs = cellfun(@(pi) mean(avgPSTHS(:,pi(1):pi(end)),2,'omitnan'), phaseInds, 'UniformOutput',false);
        if(strcmp(conditions(c),"Rest"))
            phaseAvgs = [phaseAvgs(1:2), mean(cell2mat(phaseAvgs),2,'omitnan'), phaseAvgs(3:end)];
        end
        infoTable{s,string(phaseNames)+"_"+params.condAbbrev(conditions{c})} = phaseAvgs;
        condPSTHS{c}{s} = avgPSTHS;
        segConds{c}{s} = condTimes;
    end
    figure();
    currPhases = [{cell2mat(infoTable{s,6:10})};{cell2mat(infoTable{s,11:15})};{cell2mat(infoTable{s,16:20})}];
    condLabels = arrayfun(@(r) repmat(r,length(infoTable{s,'TaskModulated'}{1}),length(phaseNames)),...
        1:length(conditions)-1, 'UniformOutput',false)';
    phaseLabels = reshape(repmat({repmat(1:length(phaseNames),length(infoTable{s,'TaskModulated'}{1}),1)},length(conditions)-1,1),[],1);
    bx=boxchart(reshape(cell2mat(condLabels)',1,[]),reshape(cell2mat(currPhases)',1,[]),'GroupByColor',...
         reshape(cell2mat(phaseLabels)',1,[]),'Notch', 'on','MarkerStyle','none');
    plotJointPSTHS(params.bins,{cell2mat(cellfun(@(a) a{s},condPSTHS(1:length(conditions)-1),'UniformOutput',false)')},...
        {cell2mat(cellfun(@(a) repmat(a{s},length(currSite.TaskModulated{1}),1),segConds(1:length(conditions)-1),'UniformOutput',false)')},...
        cell2mat(cellfun(@(a) repmat(string(a),length(currSite.TaskModulated{1}),1),params.condAbbrev.values(conditions(1:end-1)),'UniformOutput',false)'),...
        repmat(currSite.TaskModulated{1},length(conditions)-1,1),[],{[-1 3]},[0 3],cell2struct(num2cell( distinguishable_colors(length(conditions)-1,'r'),2),string(params.condAbbrev.values(conditions(1:end-1)))));
end
condPhases = [cell2mat([infoTable{:,6:10}]),cell2mat([infoTable{:,11:15}]),cell2mat([infoTable{:,16:20}])];
condPhases = condPhases(cell2mat(infoTable.TaskModulated),:);
condLabels = cell2mat(arrayfun(@(c) repmat(c,size(condPhases,1),length(phaseNames)), 1:length(conditions)-1,'UniformOutput', false));
phaseLabels = repmat(1:length(phaseNames),size(condPhases,1),length(conditions)-1);
figure();
boxchart(reshape(condLabels,1,[]),reshape(condPhases,1,[]),'GroupByColor',reshape(phaseLabels,1,[]),'Notch','on','MarkerStyle','none');
ylim([0 3]);
allSegs = cellfun(@(u,s) repmat(s,size(u,1),1), vertcat(condPSTHS{1:length(conditions)-1}),vertcat(segConds{1:length(conditions)-1}), 'UniformOutput',false);
psthLabs = arrayfun(@(c) cellfun(@(s) repmat(c,size(s,1),1),allSegs(1,:), 'UniformOutput',false), string(params.condAbbrev.values), 'UniformOutput',false);
plotJointPSTHS(params.bins,{cell2mat(horzcat(condPSTHS{1:length(conditions)-1})')},...
    {cell2mat(reshape(allSegs',1,[])')},cell2mat(horzcat(psthLabs{1:length(conditions)-1})'),...
    cell2mat(repmat(infoTable.TaskModulated,length(conditions)-1,1)),[], {[-1 3]},[1 2],cell2struct(num2cell(...
    distinguishable_colors(length(conditions),'r'),2),string(params.condAbbrev.values)));
