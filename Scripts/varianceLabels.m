phaseWindowSz = 0.2;
phaseNames = ["GoSignal","StartReach","StartHold","StartWithdraw"];
eventAlign = {"StartReach"};
margNames = {'Condition-Invariant','Condition','Noise'};
savePath = "S:\Lab\ngc14\Working\Revisions\Demixed\";

params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5,...
    containers.Map(["Extra Small Sphere","Large Sphere", "Photocell"],repmat(eventAlign,1,3)));
%%
taskAlign = containers.Map(params.condNames,repmat({{["GoSignal" "StartHold"]}},1,length(params.condNames)));
phaseAlign = containers.Map(params.condNames,cellfun(@(c) num2cell(cell2mat(c)),...
    repmat({{phaseNames}},1,length(params.condNames)),'UniformOutput',false));
phaseWin = repmat({{[0, phaseWindowSz],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)],[-phaseWindowSz*(5/4), -phaseWindowSz*(1/4)],[-phaseWindowSz*(3/4),phaseWindowSz*(1/4)]}},1,length(params.condNames));
phaseWin{strcmp(params.condNames,"Photocell")}{contains(phaseNames,"Hold")} = [-phaseWindowSz/2 0];
[siteDateMap,siteSegs,siteTrialPSTHS,~,siteChannels,chMaps,~,~] = ...
    getAllSessions(params,"Single","M1","");
simpRep =  cellfun(@(r,t) r(find(t==min(t),1)),siteDateMap.SiteRep,siteDateMap.Thresh,'UniformOutput', true)';
mappedChannels =  cell2mat(cellfun(@(ch,l) ch{end}(l(~isnan(l)))', chMaps,siteChannels, 'Uniformoutput', false)');
unitSomatotopy = cellstr(mapSites2Units(cellfun(@length, siteChannels), simpRep));
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[phaseWindowSz, 0]}},1,length(params.condNames)),siteSegs,siteTrialPSTHS,false,true);
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:length(params.condNames)),taskFR(1:length(params.condNames)),'UniformOutput',false);
taskUnits = cellfun(@cell2mat, taskUnits,'UniformOutput',false);
tUnits = any([taskUnits{:}],2);
%%
goSegs = cellfun(@(c,p) cellfun(@(a) cell2mat(cellfun(@(t) findBins(t(:,strcmp(p,"GoSignal"))-4,...
    params.bins),a,'UniformOutput',false)),c,'UniformOutput',false),siteSegs,params.condSegMap.values,'UniformOutput',false);
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(s) ...
    permute(mean(a(:,max(1,s):max(1,s)+(3.5/params.binSize),:),[2],'omitnan'),[1 3 2]),...
    num2cell(n),'UniformOutput',false),[1,1,length(n)])),3,'omitnan'))],p,t,'UniformOutput',false),siteTrialPSTHS,goSegs,"UniformOutput",false);
%%
normPSTH = cellfun(@(s,b)cellfun(@(t,n) permute(permute(sqrt(t),[1 3 2])./sqrt(n),[1 3 2]), s, b, 'Uniformoutput', false), siteTrialPSTHS, normBaseline, 'UniformOutput',false);
%normPSTH = cellfun(@(s) cellfun(@(t) zscore(t,0,2), s, 'UniformOutput',false), siteTrialPSTHS, 'UniformOutput',false);
GT = cellfun(@(g) mean(cat(3,g{:}),3,'omitnan'),num2cell([normPSTH{:}],2),'UniformOutput',false);
GT = num2cell(vertcat(GT{:}),2);
normPSTH = cellfun(@(c) cellfun(@(s) resize(s,[size(s,1),size(s,2),75],'FillValue',NaN),c,'UniformOutput',false), normPSTH,'UniformOutput',false);
normPSTH = cellfun(@(c) num2cell(cat(1,c{:}),[2 3]), normPSTH, 'UniformOutput',false);
LO = cellfun(@(c) cellfun(@(g,a) mean(g-a,3,'omitnan'), c, GT, 'UniformOutput',false), normPSTH,'UniformOutput',false);
NO = cellfun(@(c,l) cellfun(@(t,g,n) t-(g+n), c,l,GT, 'UniformOutput',false), normPSTH, LO, 'UniformOutput',false);
%%
figure
nt = tiledlayout(2,3);
plotSegs = mean(cell2mat([siteSegs{1}{:}]'),1,'omitnan');
plotSegs = plotSegs(contains(params.condSegMap(params.condNames(1)),phaseNames));
for m = 1:length(margNames)
    nexttile(); hold on; title(margNames{m});
    if(m==1)
        currM = {GT};
    elseif(m==2)
        currM = LO;
    elseif(m==3)
        currM = NO;
    end
    armM = cell2mat(reshape(cellfun(@(m) vertcat(m{unitSomatotopy=="Arm" & tUnits}),currM,'UniformOutput',false),[ones(1,ndims(currM{1})),length(currM)])).^2;
    plotMeanArm = sum(armM,find(~ismember(1:ndims(armM),2)),'omitnan')./sum(~all(isnan(armM),2),'all');
    plotVarArm = std(armM,0,find(~ismember(1:ndims(armM),2)),'omitnan')./sqrt(sum(~all(isnan(armM),2),'all'));
    ss=shadedErrorBar(params.bins,plotMeanArm,plotVarArm,'lineProps',{'Color',[0 .8 0],'LineWidth',2},'patchSaturation',0.2);
    handM = cell2mat(reshape(cellfun(@(m) vertcat(m{unitSomatotopy=="Hand" & tUnits}),currM,'UniformOutput',false),[ones(1,ndims(currM{1})),length(currM)])).^2;
    plotMeanHand = sum(handM,find(~ismember(1:ndims(handM),2)),'omitnan')./sum(~all(isnan(handM),2),'all');
    plotVarHand = std(handM,0,find(~ismember(1:ndims(handM),2)),'omitnan')./sqrt(sum(~all(isnan(handM),2),'all'));
    sh=shadedErrorBar(params.bins,plotMeanHand,plotVarHand,'lineProps',{'Color',[1 0 1],'LineWidth',2},'patchSaturation',0.2);
    xlim([-.5 2]);% ylim([0 2]);
    arrayfun(@(x) plot([x,x],get(gca,'YLim'),'k--','LineWidth',1),plotSegs);
    legend([ss.mainLine,sh.mainLine],["Arm","Hand"])
end
for c = 1:length(params.condNames)
    nexttile(); hold on; title(params.condNames(c));
    armM = cell2mat(cellfun(@(m) vertcat(m{unitSomatotopy=="Arm" & tUnits}),LO(c),'UniformOutput',false)).^2;
    plotMeanArm = sum(armM,find(~ismember(1:ndims(armM),2)),'omitnan')./sum(~all(isnan(armM),2),'all');
    plotVarArm = std(armM,0,find(~ismember(1:ndims(armM),2)),'omitnan')./sqrt(sum(~all(isnan(armM),2),'all'));
    ss=shadedErrorBar(params.bins,plotMeanArm,plotVarArm,'lineProps',{'Color',[0 .8 0],'LineWidth',2},'patchSaturation',0.2);
    handM = cell2mat(cellfun(@(m) vertcat(m{unitSomatotopy=="Hand" & tUnits}),LO(c),'UniformOutput',false)).^2;
    plotMeanHand = sum(handM,find(~ismember(1:ndims(handM),2)),'omitnan')./sum(~all(isnan(handM),2),'all');
    plotVarHand = std(handM,0,find(~ismember(1:ndims(handM),2)),'omitnan')./sqrt(sum(~all(isnan(handM),2),'all'));
    sh=shadedErrorBar(params.bins,plotMeanHand,plotVarHand,'lineProps',{'Color',[1 0 1],'LineWidth',2},'patchSaturation',0.2);
    xlim([-.5 2]); %ylim([0 1]);
    arrayfun(@(x) plot([x,x],get(gca,'YLim'),'k--','LineWidth',1),plotSegs);
    legend([ss.mainLine,sh.mainLine],["Arm","Hand"])

end
