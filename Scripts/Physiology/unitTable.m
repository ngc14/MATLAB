function [tPhys, errs] = unitTable(conditions,params,unitType,domain)
if(~exist("conditions",'var'))
    conditions = params.condNames;
end
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
%% load all single units from M1 recordings that were not from the face zone
[siteDateMap, siteSegs, siteTrialPSTHS,~, siteChannels,chMaps,~,~,errs]=...
    getAllSessions(params,unitType,domain,"");
%%
% unit somatotopy and channel number
siteRep =cellfun(@(r,t) r(find(t==min(t+(double(ismember(r,"Face"))*100)),1)),siteDateMap.SiteRep,siteDateMap.Thresh,'UniformOutput', false);
siteRep(cellfun(@isempty,siteRep)) = {"Face"};
siteRep = cell2mat(siteRep);
mappedChannels = cellfun(@(ch,l) ch{2}(l(~isnan(l))), chMaps,siteChannels, 'Uniformoutput', false)';
% normalize PSTHs to baseline activity (averaged 1s before Cue)
sumSegs = cellfun(@(c) cellfun(@(n) [n{:}], c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
goSegs = cellfun(@(c) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(maxSegL,"GoSignal"))-3,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,'UniformOutput',false);
%% get unit information by condition
tPhys = [];
g = @(x,y,c)GetPointLineDistance(x,y,c(1),c(2),c(3),c(4));
for c = 1:length(conditions)
    condTable = table();
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels)';
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    condTable.Unit = (1:length(mLabs))';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), 1:length(condUnitMapping),condUnitMapping','UniformOutput',false)');
    condTable.Date = mapSites2Units(condUnitMapping,cell2mat(siteDateMap.Date));
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(mapSites2Units(condUnitMapping,siteRep));
    condTable.Channel =  [mappedChannels{:}]';
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    % transform monkey's recording coordinate plane to standardized axes
    % along the central sulcus and orthogonal to the central sulcus axis
    mInds = mLabs=="Skipper";
    unitLocation = mapSites2Units(condUnitMapping,num2cell(cell2mat(....
        arrayfun(@(x,y,m) [g(x,y,OrthogonalLines(m).RCLine), g(x,y,OrthogonalLines(m).MLLine)],...
        siteDateMap.x,siteDateMap.y,siteDateMap.Monkey,'UniformOutput',false)),2));
    condTable.XT(mInds) = cellfun(@(x) x(1), unitLocation(mInds)) - min(cellfun(@(x) x(1), unitLocation(mInds)));
    condTable.XT(~mInds) = cellfun(@(x) x(1), unitLocation(~mInds)) - min(cellfun(@(x) x(1), unitLocation(~mInds)));
    condTable.YT(mInds) = cellfun(@(x) x(end), unitLocation(mInds)) - min(cellfun(@(x) x(end), unitLocation(mInds)));
    condTable.YT(~mInds) = cellfun(@(x) x(end), unitLocation(~mInds)) - min(cellfun(@(x) x(end), unitLocation(~mInds)));
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    PSTH = cellfun(@(m) num2cell(m,[2 3]), siteTrialPSTHS{c},'UniformOutput',false);
    % handle unreliable behavioral markers for segment time calculations
    nanSegs = find(sum(isnan(cell2mat(sumSegs{c}(~cellfun(@(a) all(isnan(a),'all'),sumSegs{c})))),1)>=...
        sum(cellfun(@(z) size(z,1),sumSegs{c}(~cellfun(@(a) all(isnan(a),'all'),sumSegs{c}))))/2);
    nanSegs= nanSegs(nanSegs~=length(maxSegL));
    for a = 1:length(nanSegs)
        for n = 1:length(sumSegs{c})
            nextSeg = intersect(setdiff(1:size(sumSegs{c}{n},2),nanSegs),nanSegs(a)+1:size(sumSegs{c}{n},2));
            sumSegs{c}{n}(:,nanSegs(a)) = sumSegs{c}{n}(:,nextSeg(1));
        end
    end
    currSegMap = cell2mat(params.condSegMap.values(cellstr(conditions(c))));
    forwardTimes = cellfun(@(t) t(:,strcmp(currSegMap,"StartReplaceHold"))-t(:,strcmp(currSegMap,"GoSignal")),mapSites2Units(condUnitMapping,sumSegs{c}),'UniformOutput',false);
    goodTrialInds = cellfun(@(t,b) ~isnan(b), vertcat(PSTH{:}), forwardTimes, 'UniformOutput',false);
    condTable.PSTH = cellfun(@(t,i) permute(t(:,:,i),[2 3 1]),vertcat(PSTH{:}),goodTrialInds,'UniformOutput',false);
    condTable.Segs = cellfun(@(s,i) s(i,:), mapSites2Units(condUnitMapping,sumSegs{c}),goodTrialInds, 'UniformOutput',false);
    tPhys = [tPhys;condTable];
end
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
end
function distance = GetPointLineDistance(x3,y3,x1,y1,x2,y2)
% Find the numerator for our point-to-line distance formula.
numerator = abs((x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1));

% Find the denominator for our point-to-line distance formula.
denominator = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);

% Compute the distance.
distance = numerator ./ denominator;
end