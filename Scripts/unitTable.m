function tPhys = unitTable(conditions)
if(varargin==0)
    conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
end
params = PhysRecording(string(conditions),.001,.001,-1,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
plotPSTHS = false;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, ~, siteChannels, siteActiveInd,...
    siteRep,siteLocation,~,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","M1");
%%
mappedChannels = cellfun(@(ch,l) ch{2}(l(~isnan(l))), chMaps,siteChannels, 'Uniformoutput', false)';
sumSegs = cellfun(@(c) cellfun(@(n) [n{:}], c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
goSegs = cellfun(@(c) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(maxSegL,"GoSignal"))-.5,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp(:)),repmat(nb,1,size(vertcat(cp(:)),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
clear normBaseline;
%%
tPhys = [];
for c = 1:length(conditions)
    condTable = table();
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels)';
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    condTable.Unit = (1:length(mLabs))';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), 1:length(condUnitMapping),condUnitMapping','UniformOutput',false)');
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(mapSites2Units(condUnitMapping,siteRep'));
    condTable.Channel =  [mappedChannels{:}]';
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    PSTH = cellfun(@(m) num2cell(m{1},[2 3]), normPSTH{c},'UniformOutput',false);
    nanSegs = find(isnan(mean(cell2mat(sumSegs{c}(~cellfun(@(a) all(isnan(a),'all'),sumSegs{c}))),1)));
    nanSegs= nanSegs(nanSegs~=length(maxSegL));
    for a = 1:length(nanSegs)
        for n = 1:length(sumSegs{c})
            nextSeg = intersect(setdiff(1:size(sumSegs{c}{n},2),nanSegs),nanSegs(a)+1:size(sumSegs{c}{n},2));
            sumSegs{c}{n}(:,nanSegs(a)) = sumSegs{c}{n}(:,nextSeg(1));
        end
    end
    currSegMap = cell2mat(params.condSegMap.values(cellstr(conditions(c))));
    forwardTimes = cellfun(@(t) t(:,strcmp(currSegMap,"StartReplaceHold"))-t(:,strcmp(currSegMap,"GoSignal")),mapSites2Units(condUnitMapping,sumSegs{c}),'UniformOutput',false);
    goodTrialInds = cellfun(@(t,b) squeeze(sum(t,2))> b & squeeze(sum(t,2))< 200*b | isnan(b), vertcat(PSTH{:}), forwardTimes, 'UniformOutput',false);
    condTable.PSTH = cellfun(@(t,i) permute(t(:,:,i),[2 3 1]),vertcat(PSTH{:}),goodTrialInds,'UniformOutput',false);
    condTable.Segs = cellfun(@(s,i) s(i,:), mapSites2Units(condUnitMapping,sumSegs{c}),goodTrialInds, 'UniformOutput',false);
    tPhys = [tPhys;condTable];
end
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
clear PSTH condTable siteTrialPSTHS normPSTH sumSegs
%%
if(plotPSTHS)
    siteCondSegs = cell2mat(cellfun(@(c) cell2mat(cellfun(@(m) mean(m,1,'omitnan'),c,'UniformOutput',false)),sumSegs,'UniformOutput',false)');
    plotJointPSTHS(params,{cell2mat(cellfun(@(m) mean(m,2,'omitnan').*10,reshape(tPhys{:,contains(tPhys.Properties.VariableNames,"PSTH_")},[],1),'UniformOutput',false)')'},...
        {siteCondSegs},repmat(cellfun(@(s,t) s(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh),length(conditions),1),...
        true(length(conditions)*height(siteDateMap),1),[],{[min(params.bins),max(params.bins)]},[0 1],...
        cell2struct(num2cell(distinguishable_colors(length(unique(cell2mat(siteDateMap.SiteRep')))),2),unique(cell2mat(siteDateMap.SiteRep'))));
end
end