conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
phaseNames = ["Baseline", "Go", "Reach", "Hold"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},{["GoSignal","StartHold"]}});
params = PhysRecording(string(conditions),.001,.001,-1,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
allSegs = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegs));
maxSegL = allSegs{maxSegL};
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
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
sTrials = 30;
model = "All";
dimCond = reshape(["Hand"]+"_"+params.condAbbrev.values',1,[]);
type = 'Traj';
phases = {"StartReach","StartHold"};
phaseWindows = {[-100 100], [-200 0]};
savePath = saveDir+type+"\"; % +extractBefore(model,"_")+"\"+extractAfter(model,"_")+"\";
colors =containers.Map(dimCond,{[1 0 0];[1 .85 0];[0 0 1];[1 0 .85];[.85 .5 0];[0 1 1]});
if(~exist(savePath,'dir')), mkdir(savePath); end
if(strcmp(type,'Traj'))
    tPhysTable = tPhys(tPhys.Monkey=="Gilligan" & tPhys.Somatotopy=="Hand",:);
    allSegs= arrayfun(@(s) tPhysTable{tPhysTable.Somatotopy==extractBefore(s,"_"),contains(tPhysTable.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))}, dimCond, 'UniformOutput',false);
    allSegs= cellfun(@(c) mean(cell2mat(c),1,'omitnan'),allSegs, 'UniformOutput',false);
    taskPSTHD= arrayfun(@(a) tPhysTable{tPhysTable.Somatotopy==extractBefore(a,"_"),contains(tPhysTable.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
    taskPSTHD= cellfun(@(a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(d,sTrials),...
        a(cellfun(@(s)size(s,2)>=sTrials,a)),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(taskPSTHD), 'UniformOutput',false);
    numUnits = arrayfun(@(s) min(cellfun(@(m) size(m{1},1),taskPSTHD(contains(dimCond,s)))), unique(arrayfun(@(t) extractBefore(t,"_"),dimCond)));
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(dimCond,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    taskPSTHD = cellfun(@(u,i) cellfun(@(n) n(i,:),u,'UniformOutput',false), taskPSTHD, [unitInds{:}],'UniformOutput',false);
    dHiStruct = struct('data',vertcat(taskPSTHD{:}),'epochStarts',reshape(cellfun(@(n) [1,n([2,3,6])], ...
        repmat(cellfun(@(a) findBins(a,params.bins),allSegs,'UniformOutput',false),sTrials,1),'UniformOutput',false),[],1),...
        'condition',repmat(cellstr(dimCond'),sTrials,1),'epochColors',cellfun(@hsv2rgb,cellfun(@(l) ...
        flipud([linspace(l(1),l(1),4);linspace(1,.5,4);linspace(.85,1,4)]'),cellfun(@rgb2hsv,....
        repmat(colors.values',sTrials,1),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
else
    tPhysTable = tPhys(tPhys.Monkey=="Gilligan" & tPhys.Somatotopy=="Hand",:);
    for p = 1:length(phases)
        phaseConds = cellfun(@(t) find(strcmp(phases{p},t)), params.condSegMap.values(conditions),'UniformOutput',false);
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end))),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct,cs,'UniformOutput',false),num2cell(tPhysTable{:,contains(tPhysTable.Properties.VariableNames,"PSTH_")},1),...
            num2cell(tPhysTable{:,contains(tPhysTable.Properties.VariableNames,"Segs_")},1),phaseConds,repmat({phaseWindows{p}},1,length(phaseConds)),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
    end
    currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) ...
        downsampleTrials(r,sTrials),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),arrayfun(@(t) ...
        m(n(:,contains(params.condAbbrev.values,extractAfter(t,"_")))>=sTrials & contains(string(tPhysTable.Somatotopy),...
        extractBefore(t,"_")),contains(params.condAbbrev.values,extractAfter(t,"_"))),dimCond,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    numUnits = arrayfun(@(s) min(cellfun(@(m) size(m{1},1) ,[currD{end}(contains(dimCond,s))])), unique(arrayfun(@(t) extractBefore(t,"_"),dimCond)));
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(dimCond,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    currD = cellfun(@(u,i) cellfun(@(n) n(i,:),u,'UniformOutput',false), [currD{:}], repmat([unitInds{:}],1,length(phases)),'UniformOutput',false);
    dHiStruct = struct('data',vertcat(currD{:}),'condition',cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
        cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)')),...
        'epochStarts',1,'epochColors',{[0 0 0]});
    for i = 1:length(dHiStruct)
        dHiStruct(i).epochColors = cell2mat(colors.values(params.condAbbrev.values(conditions(...
            cellfun(@(c) contains(dHiStruct(i).condition,c),colors.keys)))));
    end
end
DataHigh(dHiStruct,'DimReduce');
save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
%%
splitGroup = "Condition";
num_dims=4;
all_h = findall(groot,'Type','Figure');
D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));%handles = guihandles(guiFigs);
plotType = unique(string({D.D.type}));
Ddata = D.D(ismember({D.D.type}, plotType));
switch(splitGroup)
    case "Somatotopy"
    condInds = cellfun(@(c) contains(c,'Arm'), {Ddata.condition});
    for u = 1:length(condInds)
        if(condInds(u)==0),Ddata(u).epochColors = [.8 .8 .8];
        else,Ddata(u).epochColors = [.2 .2 .2];end
    end
    case  "Condition"
    condInds = cellfun(@(c) contains(c,'S-'), {Ddata.condition});
    condInds = condInds + cellfun(@(c) contains(c,'_E'), {Ddata.condition});
    for u = 1:length(condInds)
        if(condInds(u)==0),Ddata(u).epochColors = [0 0 1];
        elseif(condInds(u)==1),Ddata(u).epochColors = [1 .85 0];
        else,Ddata(u).epochColors = [1 0 0];end
    end
    case "Phase"
    condInds = cellfun(@(c) contains(c,'Reach'), {Ddata.condition});
    for u = 1:length(condInds)
        if(condInds(u)==0),Ddata(u).epochColors = [1 0 1];
        else,Ddata(u).epochColors = [0 1 1];end
    end
    case "Monkey"
    condInds = cellfun(@(c) contains(c,'Gilligan'), {Ddata.condition});
    for u = 1:length(condInds)
        if(condInds(u)==0),Ddata(u).epochColors = [.8 .4 0];
        else,Ddata(u).epochColors = [0 .5 0];end
    end
end
conds = unique({Ddata.condition});
figure(); tax=tiledlayout(max(1,num_dims/2),2);
ylimT = [0 .5];%[min(arrayfun(@(m) min(m.data,[],'all'),Ddata)),max(arrayfun(@(m) max(m.data,[],'all'),Ddata))]-[0,min(arrayfun(@(s) min(mean(s.data(1:num_dims,1:10),2,'omitnan')),Ddata))];
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
    saveFigures(gcf,savePath,model+"_"+splitGroup,[]);
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

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end