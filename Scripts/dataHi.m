varNames = ["Unit" "SiteNum" "Monkey" "Somatotopy" "Channel" "X" "Y" "Condition"...
     "TaskUnits"];
conditions = ["Extra Small Sphere"    "Large Sphere"    "Photocell"];
phaseNames = ["Baseline", "Go", "Reach", "Hold"];
phaseAlignmentPoints = {["GoSignal","GoSignal","StartReach","StartHold"],...
    ["GoSignal","GoSignal","StartReach","StartHold"],["GoSignal","GoSignal","StartReach","StartHold"]};
phaseWinSz = .2;
phaseWindows = repmat({{[-phaseWinSz, 0],[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],[-phaseWinSz, 0]}},1,length(conditions));
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},{["GoSignal","StartHold"]}});
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
params = PhysRecording(string(conditions),.001,.001,-1,2,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
MIN_BLOCKS_FOR_UNIT = 13;
savePath = "S:\Lab\ngc14\Working\DataHi\Somatotopy\";
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
close all;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","M1");
clear rawSpikes;
%%
goSegs = cellfun(@(c) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(maxSegL,"GoSignal"))-.5,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp(:)),repmat(nb,1,size(vertcat(cp(:)),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
trialFR = cellfun(@(ct,cs,ta) cellfun(@(a,b) cell2mat(cellfun(@(m,tt) ...
    squeeze(mean(m(:,max(1,tt(1)):max(1,tt(end)),:),2,'omitnan').*(1*~all(isnan(tt)))),...
    squeeze(num2cell(a,[1,2])),cellfun(@(e) [cell2mat(e),repmat([NaN,NaN],isempty(cell2mat(e)),1)],num2cell(...
    [arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(1)),'UniformOutput',false),...
    arrayfun(@(bb) find(isalmost(params.bins,bb,params.binSize/1.99),1),b{1}(:,ta(end)),'UniformOutput',false)],...
    2),'UniformOutput',false),'UniformOutput',false)'),ct,cs,'UniformOutput',false), siteTrialPSTHS, siteSegs,...
    cellfun(@(c,t) arrayfun(@(e) find(strcmp(c,e)),t{1}), params.condSegMap.values(taskAlign.keys),taskAlign.values,...
    'UniformOutput',false),'UniformOutput',false);
trialCondInfo = arrayfun(@(c) cellfun(@(s) s(strcmp(s(:,1),c),:), siteTrialInfo, 'UniformOutput',false)',conditions,'UniformOutput',false);
%%
goodFR = cellfun(@(c) cellfun(@(s) s>2 & s<200,c,'UniformOutput',false),trialFR,'UniformOutput',false);
goodFR{1}(end-11:end) = cellfun(@(b) false(size(b)),goodFR{2}(end-11:end),'UniformOutput',false);
normPSTH{1}(end-11:end) = cellfun(@(b) {NaN(size(b{1}))},normPSTH{2}(end-11:end),'UniformOutput',false);
siteTrialPSTHS{1}(end-11:end) = cellfun(@(b) NaN(size(b)),siteTrialPSTHS{2}(end-11:end),'UniformOutput',false);
trialCondInfo{1}(end-11:end) = cellfun(@(b) num2cell(NaN(size(b))),trialCondInfo{2}(end-11:end),'UniformOutput',false);
goodUnits = cellfun(@(tn) cell2mat(cellfun(@(s)sum(s,2), tn,'UniformOutput',false)),num2cell(cat(2,goodFR{:}),2),'UniformOutput',false);
normPSTH = normPSTH;     %cellfun(@(c) cellfun(@(cp) cell2mat(cp),c,'UniformOutput',false),siteTrialPSTHS,'UniformOutput',false);
mappedChannels = cellfun(@(ch,l) ch{2}(l(~isnan(l))), chMaps,siteChannels, 'Uniformoutput', false)';
avgSeg = cellfun(@(ct) cellfun(@(ca) cellfun(@(t) mean(t,1,'omitnan'), ca, 'UniformOutput',false),...
    ct, 'UniformOutput',false),siteSegs, 'UniformOutput',false);
sumSegs = cellfun(@(c) cellfun(@(n) [n{:}], c, 'UniformOutput',false), siteSegs,'UniformOutput',false);
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[0, 0]}},1,length(conditions)),siteSegs,siteTrialPSTHS,false,true);
[~,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,avgSeg,normPSTH,false,false);
[~,tUnit] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,0.05),tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
avgPhase = cellfun(@(c) cellfun(@(a) median(cell2mat(reshape(cellfun(@cell2mat,a,'UniformOutput',false),1,1,[])),3,'omitnan'),...
    c, 'UniformOutput', false), avgPhase, 'UniformOutput',false);
taskUnits = cellfun(@(a,b) cell2mat(a) & repmat(sum(b,2)>MIN_BLOCKS_FOR_UNIT*size(b,2),1,size(b,2)),num2cell(cat(2,tUnit{:}),2),goodUnits,'Uniformoutput',false);
allPSTHS = cellfun(@(c,t) cellfun(@(r,n,i) permute(permute((any(i,2)./any(i,2)).*r{1},...
    [3 2 1]).*~isnan(cellfun(@str2double,n(:,end-1))),[3 2 1]),c,t,taskUnits,'UniformOutput',false),normPSTH,trialCondInfo,'UniformOutput',false);
%%
tPhys = [];
for c = 1:length(conditions)
    condTable = table();
    tUnits = cell2mat(cellfun(@(a) a(:,c), taskUnits, 'UniformOutput',false));
    tUnits(isnan(tUnits)) = 0;
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels)';    
    allReps =  mapSites2Units(condUnitMapping,siteRep');
    mLabs = mapSites2Units(condUnitMapping,siteDateMap.Monkey);
    mInds= contains(mLabs,"Skipper");
    unitLocation = mapSites2Units(condUnitMapping,num2cell([siteDateMap.x,siteDateMap.y],2));
    condTable.Unit = [1:length(mLabs)]';
    condTable.SiteNum = cell2mat(arrayfun(@(s,c) repmat(s,c,1), [1:length(condUnitMapping)]',condUnitMapping,'UniformOutput',false));
    condTable.Monkey = categorical(mLabs);
    condTable.Somatotopy = categorical(allReps);
    condTable.Channel =  [mappedChannels{:}]';
    condTable.X = mapSites2Units(condUnitMapping,siteDateMap.x);
    condTable.Y = mapSites2Units(condUnitMapping,siteDateMap.y);
    condTable.Condition = categorical(repmat({params.condAbbrev(conditions{c})},length(mLabs),1));
    condTable.TaskUnits =  logical(tUnits);
    condTable.PSTH = num2cell(cell2mat(cellfun(@(m) mean(m{1},3,'omitnan'), normPSTH{c},'UniformOutput',false)),2);
    nanSegs = find(isnan(mean(cell2mat(sumSegs{c}),1,'omitnan'))); 
    nanSegs= nanSegs(nanSegs~=length(maxSegL));
    for a = 1:length(nanSegs)
        for n = 1:length(sumSegs{c})
            nextSeg = intersect(setdiff(1:size(sumSegs{c}{n},2),nanSegs),nanSegs(a)+1:size(sumSegs{c}{n},2));
            sumSegs{c}{n}(:,nanSegs(a)) = sumSegs{c}{n}(:,nextSeg(1));
        end
    end
    condTable.Segs = mapSites2Units(condUnitMapping,sumSegs{c});
    tPhys = [tPhys;condTable];
end
plotNames = arrayfun(@(p) arrayfun(@(c) p+"_"+c{1}(1), conditions, 'UniformOutput', true), phaseNames, 'UniformOutput', false);
plotNames = [plotNames{:}];
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
%%
siteCondSegs = cell2mat(cellfun(@(c) cell2mat(cellfun(@(m) mean(m,1,'omitnan'),c,'UniformOutput',false)),sumSegs,'UniformOutput',false)');
plotJointPSTHS(params,{cell2mat(cellfun(@(m) mean(m,3,'omitnan'),cellfun(@(n) n{1}, vertcat(normPSTH{:}), 'UniformOutput',false),'UniformOutput',false))},{siteCondSegs},...
repmat(cellfun(@(s,t) s(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh),length(conditions),1),...
true(length(conditions)*height(siteDateMap),1),[],{[min(params.bins),max(params.bins)]},[0 10],...
cell2struct(num2cell(distinguishable_colors(length(unique(cell2mat(siteDateMap.SiteRep')))),2),unique(cell2mat(siteDateMap.SiteRep'))));
%%
model = "Arm";
colors = [[1 0 0]; [1 .75 0]; [0 .25 1]];
%colors = [[.75 .75 .75]; [.3 .3 .3]; [.5 .35 0]];
saveFig = true;
somatotopicReps = categorical(model); %unique(tPhys.Somatotopy)
somatotopicReps(somatotopicReps==categorical("Trunk")) = [];
allSegs= arrayfun(@(s) tPhys{tPhys.Somatotopy==s,contains(tPhys.Properties.VariableNames,"Segs_")}, somatotopicReps, 'UniformOutput',false);
allSegs= cellfun(@(c) mean(cell2mat(c),1,'omitnan'),num2cell([allSegs{:,:}],1), 'UniformOutput',false)';
taskPSTHD= arrayfun(@(a) tPhys{tPhys.Somatotopy==a,contains(tPhys.Properties.VariableNames,"PSTH_")},somatotopicReps,'UniformOutput',false);
taskPSTHD= cellfun(@(a) vertcat(a,repmat({NaN(size(a{1}))},max(cellfun(@length,num2cell([taskPSTHD{:,:}],1)))-length(a),1)),...
    num2cell([taskPSTHD{:,:}],1)', 'UniformOutput',false);

dHiStruct = struct('data',cellfun(@(c) cell2mat(cellfun(@(t) t(:,1:end),c,'UniformOutput',false)),taskPSTHD,'UniformOutput',false),...
    'traj','traj','epochStarts',cellfun(@(n) [1,n([2,3,6])], cellfun(@(a) findBins(a,params.bins),allSegs,'UniformOutput',false),'UniformOutput',false),...
    'condition',cellstr(conditions'),'epochColors',cellfun(@(c) cell2mat(cellfun(@(m) min(1,max(0,c-repmat(.15.*(m-1),size(c,1),1))),...
    num2cell(1:length(phaseNames)),'UniformOutput',false)'),num2cell(colors,2),'UniformOutput',false));
DataHigh(dHiStruct,'DimReduce');
%%
num_dims=4;
all_h = findall(groot,'Type','Figure');
guiFigs = all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h));
handles = guihandles(guiFigs);
data = guidata(guiFigs);
D = data.D;
Dtraj = D(ismember({D.type}, 'traj'));
conds = unique({Dtraj.condition});
figure(); tax=tiledlayout(1,num_dims);
ylimT = [min(arrayfun(@(m) min(m.data,[],'all'),Dtraj)),max(arrayfun(@(m) max(m.data,[],'all'),Dtraj))]-...
    [0,min(arrayfun(@(s) min(mean(s.data(1:num_dims,1:10),2,'omitnan')),Dtraj))];
for icond = 1:length(conds)
    for idim = 1:num_dims 
        for itrial = find(ismember({Dtraj.condition}, conds{icond}))
            epochs = [Dtraj(itrial).epochStarts size(Dtraj(itrial).data,2)];
            nexttile(idim); hold on; ylim(ylimT);
            for iepoch = 1:length(epochs)-1
                indices = epochs(iepoch):(epochs(iepoch+1));
                plot(indices,Dtraj(itrial).data(idim,indices),'Color',colors(icond,:),'LineWidth',2);
                if(iepoch>1)
                    if(iepoch==length(epochs)-1)
                        line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",':','Color',colors(icond,:)./1.5,'LineWidth',2);
                    else
                        if(icond==1)
                            line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                        end
                    end
                end
            end
        end
    end
end
if(saveFig)
    saveFigures(gcf,savePath,model+"_DimXConds",[]);
end
figure(); tax=tiledlayout(1,length(conds));
for icond = 1:length(conds)
    currColor = rgb2hsv(colors(icond,:));%([1 0 0]);
    currColor(2) = 1;
    currColor(end) = .4;
    for idim = 1:num_dims 
        if(idim>num_dims/2)
            currColor(end) = 1;
            currColor(2) = max(.25,currColor(2) - (.25*(idim-(num_dims/2))));            
        else
            currColor(end) = min(1,currColor(end) + (.4*(idim-1)));
        end
        for itrial = find(ismember({Dtraj.condition}, conds{icond}))
            epochs = [Dtraj(itrial).epochStarts size(Dtraj(itrial).data,2)];
            nexttile(icond); hold on; title(string(conds(icond))); ylim(ylimT);
            for iepoch = 1:length(epochs)-1
                indices = epochs(iepoch):(epochs(iepoch+1));
                plot(indices,Dtraj(itrial).data(idim,indices)-(mean(Dtraj(itrial).data(idim,1:10),'omitnan')),...
                    'Color',hsv2rgb(currColor),'LineWidth',2);
                if(iepoch>1)
                    line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                end
            end
        end
    end
end
if(saveFig)
    saveFigures(gcf,savePath,model+"_CondXDims",[]);
end