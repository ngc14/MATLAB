conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
phaseNames = ["Baseline", "Go", "Reach", "Hold"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},{["GoSignal","StartHold"]}});
params = PhysRecording(string(conditions),.001,.001,-1,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
MIN_BLOCKS_FOR_UNIT = 20;
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
close all;
%%
[siteDateMap, siteSegs, siteTrialPSTHS, ~, siteChannels, siteActiveInd,...
    siteRep,siteLocation,~,monkeys,vMask,conditions,chMaps,siteTrialInfo] = getAllSessions(params,"Single","M1");
trialCondInfo = arrayfun(@(c) cellfun(@(s) s(strcmp(s(:,1),c),:), siteTrialInfo, 'UniformOutput',false)',conditions,'UniformOutput',false);
%%
siteTrialPSTHS{1}(end-11:end) = cellfun(@(b) NaN(size(b)),siteTrialPSTHS{2}(end-11:end),'UniformOutput',false);
trialCondInfo{1}(end-11:end) = cellfun(@(b) num2cell(NaN(size(b))),trialCondInfo{2}(end-11:end),'UniformOutput',false);
siteSegs{1}(end-11:end) = cellfun(@(s) {NaN(size(s,3),length(maxSegL))}, siteTrialPSTHS{1}(end-11:end), 'UniformOutput',false);
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
    currSegs = mapSites2Units(condUnitMapping,sumSegs{c});
    forwardTimes = cellfun(@(t) t(:,strcmp(currSegMap,"StartReplaceHold"))-t(:,strcmp(currSegMap,"GoSignal")),currSegs,'UniformOutput',false);
    goodTrialInds = cellfun(@(t,b) squeeze(sum(t,2))> b & squeeze(sum(t,2))< 200*b | isnan(b), vertcat(PSTH{:}), forwardTimes, 'UniformOutput',false);
    condTable.PSTH = cellfun(@(t,i) permute(t(:,:,i),[2 3 1]),vertcat(PSTH{:}),goodTrialInds,'UniformOutput',false);
    condTable.Segs = cellfun(@(s,i) s(i,:), currSegs,goodTrialInds, 'UniformOutput',false);
    tPhys = [tPhys;condTable];
end
tPhys = unstack(tPhys,condTable.Properties.VariableNames(find(strcmp(condTable.Properties.VariableNames,"Condition"))+1:end),"Condition");
clear PSTH condTable
%%
siteCondSegs = cell2mat(cellfun(@(c) cell2mat(cellfun(@(m) mean(m,1,'omitnan'),c,'UniformOutput',false)),sumSegs,'UniformOutput',false)');
plotJointPSTHS(params,{cell2mat(cellfun(@(m) mean(m,2,'omitnan').*10,reshape(tPhys{:,contains(tPhys.Properties.VariableNames,"PSTH_")},[],1),'UniformOutput',false)')'},...
    {siteCondSegs},repmat(cellfun(@(s,t) s(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh),length(conditions),1),...
    true(length(conditions)*height(siteDateMap),1),[],{[min(params.bins),max(params.bins)]},[0 1],...
    cell2struct(num2cell(distinguishable_colors(length(unique(cell2mat(siteDateMap.SiteRep')))),2),unique(cell2mat(siteDateMap.SiteRep'))));
% if(contains(model,unique(siteRep)))
%     colors = [[1 0 0]; [1 .75 0]; [0 .25 1]]; dimCond = params.condAbbrev.values;
% else
%     colors = [[.75 .75 .75]; [.3 .3 .3]]; dimCond = ["Arm", "Hand"];
% end
%%
saveFig = true;
sampleTrials = 30;
model = "Reach";
type= 'Spike';
savePath = saveDir+type+"\";
phases = cellfun(@(c,t) arrayfun(@(e) find(strcmp(c,e)),t),params.condSegMap.values(conditions),repmat({"StartReach"},1,length(conditions)),'UniformOutput',false);
dimCond = reshape(["Arm", "Hand"]+ "_" +params.condAbbrev.values',1,[]);
colors = [[1 0 0]; [1 .75 0]; [0 .25 1]; ...
    vertcat(cell2mat(cellfun(@(a) hsv2rgb((rgb2hsv(a).*[1 0 0])+[0 .85 .75]), {[1 0 0]; [1 .75 0]; [0 .25 1]}, 'UniformOutput',false)))];
if(~exist(savePath,'dir')), mkdir(savePath); end
if(strcmp(type,'Traj'))
    allSegs= arrayfun(@(s) tPhys{tPhys.Somatotopy==extractBefore(s,"_"),contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))}, dimCond, 'UniformOutput',false);
    allSegs= cellfun(@(c) mean(cell2mat(c),1,'omitnan'),allSegs, 'UniformOutput',false);
    taskPSTHD= arrayfun(@(a) tPhys{tPhys.Somatotopy==extractBefore(a,"_"),contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
    taskPSTHD= cellfun(@(a) vertcat(a,repmat({NaN(size(a{1}))},max(cellfun(@length,taskPSTHD))-length(a),1)),...
        (vertcat(taskPSTHD)), 'UniformOutput',false);
    dHiStruct = struct('data',cellfun(@(c) cell2mat(cellfun(@(t)  mean(t(:,1:2000,:),3,'omitnan'),c,'UniformOutput',false)),taskPSTHD,'UniformOutput',false)',...
        'traj', type,'epochStarts',cellfun(@(n) [1,n([2,3,6])], cellfun(@(a) findBins(a,params.bins),allSegs,'UniformOutput',false),'UniformOutput',false)',...
        'condition',cellstr(dimCond'),'epochColors',cellfun(@(c) cell2mat(cellfun(@(m) min(1,max(0,c-repmat(.15.*(m-1),size(c,1),1))),...
        num2cell(1:length(phaseNames)),'UniformOutput',false)'),num2cell(colors,2),'UniformOutput',false));
else
    trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) repmat(m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end))),1,(all(isnan(tt))*range(tw))+1),...,
        num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
        'UniformOutput',false)',ct,cs,'UniformOutput',false),num2cell(tPhys{:,contains(tPhys.Properties.VariableNames,"PSTH_")},1),...
        num2cell(tPhys{:,contains(tPhys.Properties.VariableNames,"Segs_")},1),phases,repmat({[-100 100]},1, length(phases)),'UniformOutput',false);
    trialFRMat = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
    numTrials = cellfun(@(s) size(s,2), [trialFR{:}]);
    currD = cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(t) t(:,randi(size(t,2),1,sampleTrials)),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),...
        arrayfun(@(t) trialFRMat(all(numTrials>=MIN_BLOCKS_FOR_UNIT,2) & contains(string(tPhys.Somatotopy),extractBefore(t,"_")),...
        contains(params.condAbbrev.values,extractAfter(t,"_"))),dimCond,'UniformOutput',false),'UniformOutput',false);
    numUnits = unique(cellfun(@(m) size(m{1},1),currD),"stable");
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(currD,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    currD = cellfun(@(u,i) cellfun(@(n) n(i,:),u,'UniformOutput',false), currD, [unitInds{:}],'UniformOutput',false);
    dHiStruct = struct('data',vertcat(currD{:}),'condition',cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),dimCond,'UniformOutput',false)')),...
        'epochStarts',1,'epochColors',num2cell(cell2mat(cellfun(@(r) repmat(r,size(currD{1},1),1),num2cell(colors,2),'UniformOutput',false)),2));
end
DataHigh(dHiStruct,'DimReduce');
save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
%% spike bin analysis for 10 sessions (LDA extraction)
% cumulative variance explained / ratio of variance explained.
num_dims=4;
all_h = findall(groot,'Type','Figure');
guiFigs = all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h));
handles = guihandles(guiFigs);
data = guidata(guiFigs);
D = data.D;
plotType = unique(string({D.type}));
Ddata = D(ismember({D.type}, plotType));
conds = unique({Ddata.condition});
figure(); tax=tiledlayout(2,num_dims/2);
ylimT = [0 .5];[min(arrayfun(@(m) min(m.data,[],'all'),Ddata)),max(arrayfun(@(m) max(m.data,[],'all'),Ddata))]-[0,min(arrayfun(@(s) min(mean(s.data(1:num_dims,1:10),2,'omitnan')),Ddata))];
for icond = 1:length(conds)
    for idim = 1:num_dims 
        for itrial = find(ismember({Ddata.condition}, conds{icond}))
            if(strcmp(plotType,'traj'))
                epochs = [Ddata(itrial).epochStarts size(Ddata(itrial).data,2)];
            else
                epochs = [find(ismember({Ddata.condition}, conds)),NaN];
            end
            nexttile(idim); hold on; ylim(ylimT);
            for iepoch = 1:length(epochs)-1
                if(strcmp(plotType,'traj'))
                    indices = epochs(iepoch):(epochs(iepoch+1));
                    plot(indices,Ddata(itrial).data(idim,indices),'Color',Ddata(icond).epochColors(1,:),'LineWidth',2);
                if(iepoch>1)
                    if(iepoch==length(epochs)-1)
                        line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",':','Color',Ddata(icond).epochColors(1,:)./1.5,'LineWidth',2);
                    else
                        if(icond==1)
                            line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                        end
                    end
                end
                else
                    [bins centers] = hist(Ddata(iepoch).data(idim,:));
                    bins = bins ./ sum(bins);
                    bar(centers, bins, 'FaceColor',Ddata(iepoch).epochColors);
                end
            end
        end
    end
end
if(saveFig)
    saveFigures(gcf,savePath,model+"_DimXConds",[]);
end
if(strcmp(plotType,'traj'))
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
            for itrial = find(ismember({Ddata.condition}, conds{icond}))
                epochs = [Ddata(itrial).epochStarts size(Ddata(itrial).data,2)];
                nexttile(icond); hold on; title(string(conds(icond))); ylim(ylimT);
                for iepoch = 1:length(epochs)-1
                    indices = epochs(iepoch):(epochs(iepoch+1));
                    plot(indices,Ddata(itrial).data(idim,indices),...-(mean(Ddata(itrial).data(idim,1:10),'omitnan')),...
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
end