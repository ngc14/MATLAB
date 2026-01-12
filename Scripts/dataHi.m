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
clear PSTH condTable siteTrialPSTHS normPSTH sumSegs
%%
siteCondSegs = cell2mat(cellfun(@(c) cell2mat(cellfun(@(m) mean(m,1,'omitnan'),c,'UniformOutput',false)),sumSegs,'UniformOutput',false)');
plotJointPSTHS(params,{cell2mat(cellfun(@(m) mean(m,2,'omitnan').*10,reshape(tPhys{:,contains(tPhys.Properties.VariableNames,"PSTH_")},[],1),'UniformOutput',false)')'},...
    {siteCondSegs},repmat(cellfun(@(s,t) s(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh),length(conditions),1),...
    true(length(conditions)*height(siteDateMap),1),[],{[min(params.bins),max(params.bins)]},[0 1],...
    cell2struct(num2cell(distinguishable_colors(length(unique(cell2mat(siteDateMap.SiteRep')))),2),unique(cell2mat(siteDateMap.SiteRep'))));
%%
model = "GilliganSkipper_ArmHand";
type = 'Spike';
num_dims=4;
sTrials = 30;
plotTrials = 0;
saveFig = true;
phases = {"StartReach","StartHold"};
phaseWindows = {[-100 100], [-200 0]};
dimCond = reshape(regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')+"_"+params.condAbbrev.values',1,[]);
%"-"+regexp(s,'(?<=Start)\w*','match'),phases,'UniformOutput',false)),1,[]); 
if(~plotTrials)
    type = type+"_Avg";
end
savePath = saveDir+type+"\"+extractBefore(model,"_")+"\"+extractAfter(model,"_")+"\";
if(~exist(savePath,'dir')), mkdir(savePath); end
splitGroup = "Somatotopy";
colors = {};
switch(splitGroup)
    case "Somatotopy"
        condInds = arrayfun(@(c) contains(c,'Arm'), dimCond);
        for u = 1:length(condInds)
            if(condInds(u)==0),colors{u} = [.8 .8 .8];
            else,colors{u} = [.2 .2 .2];end
        end
    case  "Condition"
        condInds = arrayfun(@(c) contains(c,'S-'), dimCond);
        condInds = condInds + arrayfun(@(c) contains(c,'_E'), dimCond);
        for u = 1:length(condInds)
            if(condInds(u)==0),colors{u} = [0 0 1];
            elseif(condInds(u)==1),colors{u} = [1 .85 0];
            else,colors{u} = [1 0 0];end
        end
    case "Phase"
        condInds = arrayfun(@(c) contains(c,'Reach'), dimConds);
        for u = 1:length(condInds)
            if(condInds(u)==0),colors{u} = [1 0 1];
            else,colors{u} = [0 1 1];end
        end
    case "Monkey"
        condInds = arrayfun(@(c) contains(c,'Gilligan'), dimConds);
        for u = 1:length(condInds)
            if(condInds(u)==0),colors{u} = [.8 .4 0];
            else,colors{u} = [0 .5 0];end
        end
end
colors =containers.Map(dimCond,colors);
tPhysTable = tPhys(contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]),:);
if(contains(type,'Traj'))
    allSegs= arrayfun(@(s) tPhysTable{tPhysTable.Somatotopy==extractBefore(s,"_"),contains(tPhysTable.Properties.VariableNames,"Segs_"+extractBefore(s,"_"))}, dimCond, 'UniformOutput',false);
    allSegs= cellfun(@(c) mean(cell2mat(c),1,'omitnan'),allSegs, 'UniformOutput',false);
    taskPSTHD= arrayfun(@(a) tPhysTable{tPhysTable.Somatotopy==extractBefore(a,"_"),contains(tPhysTable.Properties.VariableNames,"PSTH_"+extractBefore(s,"_"))},dimCond,'UniformOutput',false);
    if(~plotTrials)
        taskPSTHD = cellfun(@(n) {cell2mat(cellfun(@(m) mean(m,2,'omitnan')',n,'UniformOutput',false))}, taskPSTHD, 'UniformOutput',false);
    else
        taskPSTHD= cellfun(@(a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(d,sTrials),...
            a(cellfun(@(s)size(s,2)>=sTrials,a)),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(taskPSTHD), 'UniformOutput',false);
    end
    numUnits = arrayfun(@(s) min(cellfun(@(m) size(m{1},1),taskPSTHD(contains(dimCond,s)))), unique(arrayfun(@(t) extractBefore(t,"_"),dimCond)));
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(dimCond,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    taskPSTHD = cellfun(@(u,i) cellfun(@(n) n(i,1:2000),u,'UniformOutput',false), taskPSTHD, [unitInds{:}],'UniformOutput',false);
    cls = cellfun(@(r) repmat({r},max(plotTrials*sTrials,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),4);...
        linspace(1,.25,4);linspace(.85,1,4)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
    dHiStruct = struct('data',vertcat(taskPSTHD{:}),'epochStarts',reshape(cellfun(@(n) [1,n([2,3,6])], ...
        repmat(cellfun(@(a) findBins(a,params.bins),allSegs,'UniformOutput',false),max(plotTrials*sTrials,1),1),'UniformOutput',false),[],1),...
        'condition',cellstr(cell2mat(cellfun(@(d) repmat(string(d),max(plotTrials*sTrials,1),1),cellstr(dimCond),'UniformOutput',false)')),'epochColors',...
        vertcat(cls{:}));
else
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
        dHiStruct(i).epochColors = cell2mat(colors.values({extractBefore(dHiStruct(i).condition,"-")}));
    end
end
DataHigh(dHiStruct,'DimReduce');
save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
%%
all_h = findall(groot,'Type','Figure');
D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));%handles = guihandles(guiFigs);
D = D.D(ismember({D.D.type}, plotType));
plotType = unique(string({D.type}));
if(length(D)<sTrials)
    plotTrials = 0;
end
conds = unique({D.condition});
figure(); tax=tiledlayout(max(1,num_dims/2),2);
ylimT = [min(arrayfun(@(m) min(m.data,[],'all'),D)),max(arrayfun(@(m) max(m.data,[],'all'),D))]-[0,min(arrayfun(@(s) min(mean(s.data(1:num_dims,1:10),2,'omitnan')),D))];
for icond = 1:length(conds)
    for idim = 1:num_dims 
            if(strcmp(plotType,'traj'))
                itrial = find(ismember({D.condition}, conds{icond}));
                epochs = [mean(vertcat(D(itrial).epochStarts),1,'omitnan'),size(D(1).data,2)];
            else
                epochs = [find(ismember({D.condition}, conds)),NaN];
            end
            nexttile(idim); hold on; title(idim); ylim(ylimT);
            dTrial = cat(3,(D(itrial).data));
            for iepoch = 1:length(epochs)-1
                if(strcmp(plotType,'traj'))
                    indices = epochs(iepoch):(epochs(iepoch+1));
                    cellfun(@(p) plot(indices,p,'Color', cell2mat(colors.values(cellstr(dimCond(icond)))),'LineWidth',2),...
                        num2cell(squeeze(mean(dTrial(idim,indices,:),max(~plotTrials*length(size(dTrial))+1,plotTrials),'omitnan')),1+~plotTrials));
                if(iepoch>1)
                    if(iepoch==length(epochs)-1)
                        line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",':','Color',cell2mat(colors.values(cellstr(dimCond(icond))))./1.5,'LineWidth',2);
                    else
                        if(icond==1)
                            line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                        end
                    end
                end
                else
                    [bins centers] = hist(D(iepoch).data(idim,:));
                    bins = bins ./ sum(bins);
                    bar(centers, bins, 'FaceColor',cell2mat(colors.values(cellstr(dimCond(icond)))));
                end
            end
    end
    if(icond==length(conds))
        l = arrayfun(@(d) plot(NaN(length(dimCond),1),'LineWidth',3,'Color',cell2mat(colors.values(cellstr(d)))),dimCond);
        legend(l,dimCond,'Autoupdate','off','FontSize',14,'Orientation','horizontal');
    end
end
if(saveFig)
    saveFigures(gcf,savePath,model+"_"+splitGroup,[]);
end
if(strcmp(plotType,'traj'))
    figure(); tax=tiledlayout(1,length(conds)); 
    for icond = 1:length(conds)
        legendColors= {};
        currColor = rgb2hsv(cell2mat(colors.values(cellstr(dimCond(icond)))));%([1 0 0]);
        currColor(2) = 1;
        currColor(end) = .4;
        for idim = 1:num_dims
            if(idim<=num_dims/2)
                currColor(end) = min(1,currColor(end) + (.4*(idim-1)));
            else
                currColor(end) = 1;
                currColor(2) = max(.25,currColor(2) - (.25*(idim-(num_dims/2))));
            end
            itrial = find(ismember({D.condition}, conds{icond}));
            dTrial = cat(3,D(itrial).data);
                epochs = [mean(vertcat(D(itrial).epochStarts),1,'omitnan'),size(D(1).data,2)];
                nexttile(icond); hold on; title(string(conds(icond))); ylim(ylimT);
                for iepoch = 1:length(epochs)-1
                    indices = epochs(iepoch):(epochs(iepoch+1));
                    cellfun(@(p) plot(indices,p,'Color', hsv2rgb(currColor),'LineWidth',2),...
                        num2cell(squeeze(mean(dTrial(idim,indices,:),max(~plotTrials*length(size(dTrial))+1,plotTrials),'omitnan')),1+~plotTrials));
                    if(iepoch>1)
                        line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                    end
                end
            legendColors{idim} = hsv2rgb(currColor);
            if(idim==num_dims)
                l = cellfun(@(d) plot(NaN(length(num_dims),1),'Color',d),legendColors);
                legend(l,string(1:num_dims),'Autoupdate','off','Location','northeast');
            end
        end
    end
    if(saveFig)
        saveFigures(gcf,savePath,model+"_CondByDim",[]);
    end
end

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end