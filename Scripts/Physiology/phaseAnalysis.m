savePath = "S:\Lab\ngc14\Working\M1\";
conditions = ["Extra Small Sphere", "Large Sphere", "Photocell"];
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment);
alignLimits = {[-.5, .15], [-.20, .20],[-.15,.5]};
pVal = 0.05;
unitNames = ["Reach", "Grasp", "Both"];
phaseNames = ["Go", "Reach", "Grasp","Withdraw"];
phaseWinSz = 0.2;
taskWindow = repmat({{[phaseWinSz, 0]}},1,length(conditions));
phaseWindows = repmat({{[0, phaseWinSz],[-phaseWinSz*(3/4),phaseWinSz*(1/4)],...
    [-phaseWinSz, 0],[-phaseWinSz*(3/4),phaseWinSz*(1/4)]}},1,length(conditions));
taskAlignmentPoints = {{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]}};
phaseAlignmentPoints = {["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"],...
    ["GoSignal","StartReach","StartHold","StartWithdraw"]};
%%
[siteDateMap, siteSegs, siteTrialPSTHS, rawSpikes, siteChannels, siteActiveInd,...
    siteRep,siteLocation,siteMasks,monkeys,vMask,conditions] = getAllSessions(params,"Single","M1");
taskAlign = containers.Map(conditions,taskAlignmentPoints);
condPhaseAlign = containers.Map(conditions,cellfun(@num2cell,phaseAlignmentPoints,'UniformOutput',false));
allCondSegs = cellfun(@(c) cellfun(@(a) cellfun(@(t) findBins(t(1)-3,params.bins),a),...
    c,'UniformOutput',false),siteSegs,'UniformOutput',false);
avgSeg = cellfun(@(ct) cellfun(@(ca) cellfun(@(t) mean(t,1,'omitnan'), ca, 'UniformOutput',false),...
    ct, 'UniformOutput',false),siteSegs, 'UniformOutput',false);
clear rawSpikes
%%
[taskBaseline,taskFR] = calculatePhases(params,taskAlign,taskWindow,siteSegs,siteTrialPSTHS,false,true);
[phaseBaseline,phaseFR] = calculatePhases(params,condPhaseAlign,phaseWindows,siteSegs,siteTrialPSTHS,false,true);
%unit PSTH
normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) ...
    permute(mean(c(:,s:s+(1/params.binSize),:),[2],'omitnan'),[1 3 2]),a(~isnan(n)),...
    num2cell(n(~isnan(n))),'UniformOutput',false),[1,1,sum(~isnan(n))])),3,'omitnan'));NaN(all(isnan(n)).*size(a{1},1),1)],p,t,...
    'UniformOutput',false),siteTrialPSTHS,allCondSegs,"UniformOutput",false);
normPSTH = cellfun(@(cp,nb) num2cell(cellfun(@(p,b)permute(permute(p,[1 3 2])./repmat(b,1,1,size(p,2)),[1 3 2]),...
    vertcat(cp{:}),repmat(nb,1,size(vertcat(cp{:}),2)),'UniformOutput',false),2),siteTrialPSTHS,normBaseline,'Uniformoutput', false);
%%
[avgBaseline,avgPhase] =  calculatePhases(params,condPhaseAlign,phaseWindows,avgSeg,normPSTH,false,false);
avgBase = cellfun(@cell2mat,cellfun(@(c) cellfun(@(a) median(cell2mat(a{1}),2,'omitnan'),...
    c, 'UniformOutput', false), avgBaseline, 'UniformOutput',false),'UniformOutput',false);
avgPhase = cellfun(@(c) cellfun(@(a) median(cell2mat(reshape(cellfun(@cell2mat,a(2),'UniformOutput',false),1,1,[])),3,'omitnan'),...
    c, 'UniformOutput', false), avgPhase, 'UniformOutput',false);
%%
phaseLim = [0 4];
countLim = [1 20];
SILim = [0.3 .6];
for c = 1:length(conditions)
    %%
    close all;
    tic
    currCond = conditions{c};
    condSegs = siteSegs{c};
    condDir = savePath+string(values(params.condAbbrev,{currCond}))+"\";
    condUnitMapping = cellfun(@(si) size(si,2),siteChannels{2})';
    condRep =  mapSites2Units(condUnitMapping,siteRep)';
    unitLocation = mapSites2Units(condUnitMapping,siteLocation');
    activityInds = mapSites2Units(condUnitMapping,siteActiveInd{c}')';
    condPSTHS = num2cell(cellfun(@(m) mean(m,3,'omitnan'),vertcat(normPSTH{c}{:}),'UniformOutput',false),1);
    %% unit FR modulation per phase
    [tAUC,tUnits] = cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
        taskBaseline{c},taskFR{c},'UniformOutput', false);
    tUnits = cell2mat(tUnits)';
    tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
        s,'UniformOutput', false),phaseFR{c},'UniformOutput', false),phaseNames, 'UniformOutput',false);
    [~,bT] = cellfun(@(pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,pVal),...
        phaseBaseline{c},pc, 'UniformOutput', false),tPhase,'UniformOutput',false);
    % phase units
    [~,rgInds] = cellfun(@(r,g) ttestTrials(r,g,1,true,pVal),tPhase{strcmp(phaseNames,"Reach")},...
        tPhase{strcmp(phaseNames,"Grasp")}, 'UniformOutput',false);
    rgInds = cell2mat(rgInds);
    AUCVals = vertcat(avgPhase{c}{:});
    goAUC = AUCVals(:,strcmp(phaseNames,"Go"));
    rAUC = AUCVals(:,strcmp(phaseNames,"Reach"));
    gAUC = AUCVals(:,strcmp(phaseNames,"Grasp"));
    wAUC = AUCVals(:,strcmp(phaseNames,"Withdraw"));
    rSI = rAUC./(rAUC + gAUC);
    gSI = gAUC./(rAUC + gAUC);

    rUnits = rgInds==1 & rAUC > gAUC;
    gUnits = rgInds==1 & rAUC < gAUC;
    bothUnits = ~rUnits & ~gUnits & tUnits'==1;
    tkAUC = cell2mat(cellfun(@(tc) mean(cell2mat(tc),2,'omitnan'),tAUC,'UniformOutput', false));
    % gRUnits = gUnits & cell2mat(bT{strcmp(phaseNames,"Reach")})' ==1;
    % gUnits = gUnits & ~gRUnits;
    % segInd = cell2mat(cellfun(@(ts) find(contains(params.condSegMap(currCond), ts)), taskAlign(currCond), 'UniformOutput', false));
    % for t = 1:length(condSegs)
    %     tLength = abs(diff(fix(condSegs{t}{1}(:,segInd)./params.binSize),1,2));
    %     tkU{t} = nanmean(tLength.*tAUC{t}{1}',1);
    % end
    plotPhaseVals = cellfun(@transpose,{tkAUC,goAUC,rAUC,gAUC,wAUC,rSI,gSI},'UniformOutput',false);        
    % AUCVals = [avgBase{c}, AUCVals];
    % AUCVals = diff(AUCVals,1,2);
    for m=1:length(monkeys)+1
        mMapping = cell2mat(arrayfun(@(m,n) ones(1,m)*n,condUnitMapping,...
            1:length(condUnitMapping),'UniformOutput',false));
        if(m==length(monkeys)+1)
            monk = monkeys;
            vmMask = [];
        else
            monk = monkeys(m);
            vmMask = vMask(monk);
            vmMask = double(vmMask{1});
        end
        mInds = contains(siteDateMap.Monkey,monk);
        muInds = mapSites2Units(condUnitMapping,mInds')' ;
        mMapping(~muInds) = NaN;
        mMapping = mMapping-(nanmin(mMapping))+1;
        tuning = repmat("",1,sum(muInds));
        if(m==length(monkeys))
            tuning(bothUnits(muInds)) = "Both";
            tuning(rUnits(muInds)) = "Reach";
            tuning(gUnits(muInds)) = "Grasp";
            somatotopy = condRep(muInds);
            somatotopy = somatotopy(~strcmp(somatotopy,"Face"));
            tuning = tuning(~strcmp(somatotopy,"Face"));
            tbl = table(tuning',somatotopy');
            tbl.Properties.VariableNames = ["Tuning", "Somatotopy"];
            % tuning(gRUnits(muInds)) = "Grasp-reach";
            % [hTable,chi2,p,labels] = crosstab(tuning,somatotopy);
            % hh = heatmap(tbl,'Tuning', 'Somatotopy','Title', ...
            %      [strjoin([mon," TuningxSomatotopy (df = ",...
            %      num2str(numel(hTable)),")"],''), strjoin(["Ci^2 :",...
            %      num2str(chi2,5),", p = ", num2str(p,5)])]);
            close all;
        end
        mSiteLocation = unitLocation-min(unitLocation(muInds',:)).*ImagingParameters.px2mm;
        trialSegs =  num2cell(vertcat(condSegs{mInds}),1);
        mutInds = muInds & tUnits==1;
        currRep = condRep;
        %task units
        currRep(~muInds) = "";
       FRFig = modulatedUnitsPerRep(currRep(mutInds),AUCVals(mutInds,:)',phaseNames,strjoin([monk,"Task Units"]));
       cellfun(@(y) set(gca(y),'YLim',[-3 3]), FRFig);
       cellfun(@(f) saveFigures(f, condDir+"FR_Diff\Task\",strjoin([monk+gca(f).Title.String]),[]),FRFig);
       %%
       close all;
        if(m<=length(monkeys))
            for p = 1:length(phaseNames)
                currAUC = AUCVals(:,p);
                currInds = currAUC.*(mutInds./mutInds)';
                siteAUC = arrayfun(@(a) currInds(mMapping==a), min(mMapping):max(mMapping), 'UniformOutput',false);
                emptyInds = ~cellfun(@any,siteAUC);
                FRMapFig = mapUnitVals(vmMask,siteMasks(mInds'),cellfun(@nanmean,siteAUC),emptyInds,0,256,phaseLim);
                title(strjoin([monk+"Task",phaseNames(p)]));
                colorbar(gca,'Ticks',linspace(0,1,5),'TickLabels',num2cell(linspace(0,4,5)));
                saveFigures(FRMapFig, condDir+"Maps\FR\",monk+phaseNames(p),[])
            end 
        end
        close all;
%%
 close all
currRep(~mutInds) = "";
% figure('Units','normalized','Position',[0 0 1 1])
% for p = 2:5
%     subplot(2,2,p-1);
%     hold on
%     currAUC = max(cell2mat(plotPhaseVals(1:end-2)'),[],1);%plotPhaseVals{p};
%     x = cellfun(@sum,arrayfun(@(a) mutInds(mMapping==a), min(mMapping):max(mMapping), 'UniformOutput',false));
%     y=cellfun(@mean,arrayfun(@(a) currAUC(mMapping==a), min(mMapping):max(mMapping),'UniformOutput',false));
%     y(isnan(y)) = 0;
%     r = fitlm(x,y);
%     scatter(x,y);
%     hold on;
%     plot(r)
%     title(["Max FR (r^2="+num2str(r.Rsquared.Ordinary,'%.02f')+")"])
% end
monkeyVMap = siteMasks(mInds');
        figure();%figure('Units','normalized','Position',[0 0 1 1]);
    for t = 1:3
        if(t==1)
            tInds = muInds & tUnits==1 & rUnits'==1;
            tName = "Reach";
        elseif(t==2)
            tInds = muInds & tUnits==1 & gUnits'==1;
            tName = "Grasp";
        else
            tInds = muInds & tUnits==1 & bothUnits'==1;
            tName = "Both";
        end
        siteUnitMods = arrayfun(@(a) tInds(mMapping==a),min(mMapping):max(mMapping),'UniformOutput', false);
        siteReps = arrayfun(@(a) unique(currRep(mMapping==a & arrayfun(@strlength,currRep)>0)), min(mMapping):max(mMapping),'UniformOutput',false);
        siteLocs = cell2mat(arrayfun(@(a) unique(unitLocation(mMapping==a,:),'rows'), min(mMapping):max(mMapping),'UniformOutput',false)');
        siteReps(cellfun(@isempty,siteReps)) = {""};
        siteLocs(cellfun(@(s) strcmp(s,""),siteReps),:) = NaN(sum(cellfun(@(s) strcmp(s,""),siteReps)),2);
        siteReps = string(siteReps);
        uniqueRep = unique(currRep(arrayfun(@strlength,currRep)>0));
        div = zeros(size(siteReps));
        for d = 1:length(uniqueRep)
            div(strcmp(siteReps,uniqueRep(d))) = sum(tInds & strcmp(currRep,uniqueRep(d)));
        end
        plotSites = cellfun(@(s,d) 100*nansum(s)/nansum(tInds), siteUnitMods,num2cell(div));
        % countMapFig = mapUnitVals(vmMask,monkeyVMap,plotSites,emptyInds,true,4,[.0005 3]);
        siteLocsN = round(siteLocs./min(siteLocs(:)));
        hm = NaN(fliplr(max(siteLocsN)));
        nonEmpty = cellfun(@any,siteUnitMods);
        for s = 1:length(siteLocsN)
            if(nonEmpty(s))
                hm(siteLocsN(s,2),siteLocsN(s,1)) = plotSites(s);
            end
        end
        hm(all(isnan(hm),2),:) = [];
        hm(:,all(isnan(hm),1)) = [];
        
        subplot(2,3,t);
        scatter3(siteLocsN(:,1),siteLocsN(:,2),plotSites);%imagesc(hm);
        set(gca,'YDir','reverse');%axis image
        ylim([min(siteLocsN(:,2))-1, max(siteLocsN(:,2))]);
        xlim([min(siteLocsN(:,1))-1, max(siteLocsN(:,1))]);
        colormap([flipud(colormap('bone'))]);
        view(90,0);
        title(strjoin([strjoin(monk,'')+tName+"YZ"]));
        ylabel("Lateral to Medial");
        
        subplot(2,3,3+t);
        scatter3(siteLocsN(:,1),siteLocsN(:,2),plotSites);%imagesc(hm)
        set(gca,'YDir','reverse');% axis image
        ylim([min(siteLocsN(:,2))-1, max(siteLocsN(:,2))]);
        xlim([min(siteLocsN(:,1))-1, max(siteLocsN(:,1))]);
        colormap([flipud(colormap('bone'))])
        view(0,0);
        title(strjoin([strjoin(monk,'')+tName+"XZ"]));
        xlabel("Caudal to Rostral");
    end
    linkaxes([subplot(2,3,1),subplot(2,3,2),subplot(2,3,3)]);
    linkaxes([subplot(2,3,4),subplot(2,3,5),subplot(2,3,6)]);
    saveFigures(gcf, strjoin([condDir,"Maps\Grids\"],''),"Scatter_"+strjoin(monk,''),[]);

%%
close all;
siteUnitMods = arrayfun(@(a) mutInds(mMapping==a),...
min(mMapping):max(mMapping),'UniformOutput', false);
emptyInds = ~cellfun(@any,siteUnitMods);
newMutInds = [plotPhaseVals{end-1}].*mutInds;
countMapFig = mapUnitVals(vmMask,siteMasks(mInds'),newMutInds,emptyInds,false,255,SILim);
h = imhandles(gcf);
cb = colorbar(h(1).Parent);
cb.Ticks = [cb.Limits(1), cb.Limits(end)];
cb.TickLabels = string(SILim);
saveFigures(gcf, strjoin([condDir,"Maps\FR\"]),"Gilligan_Task"+"rSI",[])
%%
        countAUC(currRep,mutInds,plotPhaseVals,mMapping,siteMasks(mInds'),...
            vmMask,activityInds,condDir,strjoin(monk,'')+"_Task",["Task",arrayfun(@(s) strcat(s,"Phase"),phaseNames),"RSI","GSI"],...
            m<=length(monkeys),countLim,{phaseLim.*2,phaseLim,phaseLim,phaseLim,phaseLim,SILim,SILim});
        close all;
        %
        maxUnitFR = max(cell2mat(cellfun(@(c) cell2mat(cellfun(@(p) max(p,[],2), c,...
            'UniformOutput', false)), condPSTHS, 'UniformOutput',false)),[],2,'omitnan');
        unitPSTHS = cellfun(@(p) (muInds./muInds)'.*vertcat(p{:}),condPSTHS,'UniformOutput',false);
        psthJointxUnit(currRep,mutInds,mMapping,trialSegs,params.bins,unitPSTHS,...
            maxUnitFR,mSiteLocation,activityInds,alignLimits,condDir,strjoin(monk,'')+"Task",[0,6]);
        close all;
        currRep(~mutInds) = "";
        %
        for u = 1:length(unitNames)
            switch(unitNames(u))
                case "Both"
                    uInds = bothUnits;
                case "Reach"
                    uInds = rUnits;
                case "Grasp"
                    uInds = gUnits;
                otherwise
                    uInds = tUnits';
            end
            uInds = uInds'==1&muInds==1;
            %unit type count and phase modulation
            psthJointxUnit(currRep,uInds,mMapping,trialSegs,...
                 params.bins,cellfun(@(v) vertcat(v{:}), condPSTHS, 'UniformOutput',false),...
                 maxUnitFR,mSiteLocation,activityInds,alignLimits,condDir,...
                 strjoin(monk,'')+unitNames(u),[0 10])
             close all;
            %%% Functional grouping PSTHS
            uMMapping = cell2mat(arrayfun(@(m,n) ones(1,m)*n,condUnitMapping,...
                1:length(condUnitMapping),'UniformOutput',false));
            uMMapping(~muInds) = NaN;
            uMMapping = uMMapping-(nanmin(uMMapping))+1;
            sessionInds = ~isnan(uMMapping) & uInds;
            uMMapping(~sessionInds) = 0;
            [siteIndsN,siteInds] = unique(uMMapping);
            siteInds = siteInds(siteIndsN>0);
            uReps = condRep;
            uReps(~uInds) = "";
            uPSTHS = cellfun(@(p) (uInds./uInds)'.*p, unitPSTHS,'UniformOutput',false);
            % unitPSTHFig = plotJointPSTHS(params.bins, uPSTHS, trialSegs, uReps, siteInds, activityInds, alignLimits, [0 10]);
            % saveFigures(unitPSTHFig,strcat(condDir,"UnitPSTH\"),strjoin([monk,unitNames(u)],''),[]);
            uReps(uInds) = unitNames(u);
            unitPSTHAll = plotJointPSTHS(params.bins,uPSTHS,trialSegs,uReps,siteInds,[],alignLimits,[0 10]);
            saveFigures(unitPSTHAll,strcat(condDir,"UnitPSTH\"),strjoin([monk,unitNames(u),"_ALL"],''),[]);
            close all;
            countAUC(currRep,uInds,plotPhaseVals,...
                mMapping,siteMasks(mInds'),vmMask,activityInds,...
                condDir,strjoin(monk,'')+"_"+unitNames(u),["Task",arrayfun(@(s) strcat(s,"Phase"),phaseNames),"rSI","gSI"],...
                m<=length(monkeys),countLim./2,{phaseLim.*2,phaseLim,phaseLim,phaseLim,phaseLim,...
                [0.3 .9],[0.3 .9]});
        end
        disp(join([currCond, ' for ', monk, ' done.']));
        toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sessionInds = ~isnan(mMapping) & mutInds;
        siteUnitMods = mMapping;
        siteUnitMods(~sessionInds) = 0;
        [siteIndsN,siteInds] = unique(siteUnitMods);
        siteInds = siteInds(siteIndsN>0);
        phaseJointFig = plotJointPSTHS(params.bins,unitPSTHS,trialSegs,currRep,siteInds,[],...
            alignLimits,[0 8]);
        %
        % jointInds = arrayfun(@(jN) arrayfun(@(mJ) logical(cellfun(@(ss) strcmp(ss,jN),...
        %     mJ,'UniformOutput',true)),currRep(mutInds)),unique(currRep(mutInds)),'UniformOutput', false);
        % AUCVals = cell2mat({bAUC(mutInds),goAUC(mutInds),rAUC(mutInds),gAUC(mutInds)}');
        % jointVals = cellfun(@(jI) AUCVals(:,jI),jointInds, 'UniformOutput', false);
        % jointVals = cellfun(@(jI) [jI, NaN(size(jI,1),max(...
        %     cellfun(@length,jointVals))-size(jI,2))], jointVals, 'UniformOutput',false);
        % figure(m)
        % bp = boxplot2(permute(cat(3,jointVals{:}),[1 3 2]),'notch','on');
        % cm = cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)),...
        %     unique(currRep(arrayfun(@(s) ~strcmp(s,""),currRep))), 'UniformOutput', false)');
        % for ii = 1:length(jointInds)
        %     structfun(@(b) set(b(ii,:), 'markeredgecolor', cm(ii,:), ...
        %         'color', cm(ii,:), 'markerfacecolor',cm(ii,:)), bp);
        % end
        % set(bp.out, 'marker', '.')
        % xticks(1:length(pn));
        % xticklabels([pn(end), pn(1:3)]);
        % ylim([0 40])
        %
        % saveFigures(gcf, strjoin([condDir,strjoin(monk,'')+"_","AUC\"]),"Task_Box",[]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%% cond specific mappings
% reach speed and FR photocell correlations
uMapping = cellfun(@(si) size(si,2),siteChannels{1})';
sMapping = cell2mat(arrayfun(@(m,n) ones(1,m)*n,uMapping,1:length(uMapping),...
    'UniformOutput',false));
pSegs = siteSegs{strcmp(conditions, "Photocell")};
pSegs = vertcat(pSegs{:});
ptInds = matches(params.condSegMap("Photocell"), ["StartReach","StartHold"]);
[~,~,~,~,tTime] = cellfun(@(pS,uF) cellfun(@(c) regress(c', [ones(size(pS,1),1),1./range(pS(:,ptInds),2)]),num2cell(uF{1}{end},2),'UniformOutput', false),...
    pSegs(:,1), avgPhase{strcmp(conditions,"Photocell")}, 'UniformOutput', false);
tTime = cellfun(@(d) cellfun(@(us) us(1), d), tTime, 'UniformOutput',false)
%Rsq = 1 - sum((y - pred).^2)/sum((y - mean(y)).^2)
%tTime(85:88) = cellfun(@(n) NaN(length(n),1), tTime(85:88), 'UniformOutput', false);
corrsP = cell2mat(tTime)';
%%
% sphere differences in the task phase
ssInds = strcmp(conditions,"Extra Small Sphere");
lsInds = strcmp(conditions,"Large Sphere");
stInds = matches(params.condSegMap("Extra Small Sphere"),taskAlign("Extra Small Sphere"));
ltInds = matches(params.condSegMap("Large Sphere"),taskAlign("Large Sphere"));
tcSSVals = cellfun(@(sp,sc) cellfun(@(s,c) {s{end}./range(c(:,stInds),2)'}, sp,sc, 'UniformOutput', false),...
    avgPhase{ssInds},siteSegs{ssInds}, 'Uniformoutput', false);
tcLSVals = cellfun(@(sp,sc) cellfun(@(s,c) {s{end}./range(c(:,ltInds),2)'}, sp,sc, 'UniformOutput', false),...
    avgPhase{lsInds},siteSegs{lsInds}, 'Uniformoutput', false);
[spVals, spInds] = cellfun(@(s,l) ttestTrials(s,l,1,false,pVal),...
    tcSSVals,tcLSVals, 'UniformOutput', false);
spInds = [spInds{:}];
spVals = cell2mat(cellfun(@(pv) nanmean(pv,2),[spVals{:}], 'UniformOutput', false)')';

allReps =  mapSites2Units(uMapping,siteRep)';

for m=1:length(monkeys)+1
    sReps = allReps;
    mMapping = sMapping;
    if(m>length(monkeys))
        monk = monkeys;
        vmMask = [];
    else
        monk = monkeys(m);
        vmMask = vMask(monk);
        vmMask = vmMask{1};
    end
    mInds = contains(siteDateMap.Monkey,monk);
    muInds = mapSites2Units(uMapping,mInds')' ;
    muspInds = muInds & spInds==1;

    mMapping(~muInds) = NaN;
    mMapping = mMapping-(nanmin(mMapping))+1;
    sReps(~muInds) = "";

    psthJointxUnit(sReps,muspInds,mMapping,num2cell(vertcat(...
        siteSegs{ssInds}{mInds}),1),params.bins,unitPSTH{ssInds},...
        maxUnitFR, mapSites2Units(cellfun(@(si) size(si,2),siteChannels{ssInds})',...
        siteLocation'),NaN,alignLimits,strcat(savePath,"Spheres\"),...
        strjoin([monk+"ESS_Grasp"],''),[5 35]);
    psthJointxUnit(sReps,muspInds,mMapping,num2cell(vertcat(...
        siteSegs{lsInds}{mInds}),1),params.bins,unitPSTH{lsInds},...
        maxUnitFR, mapSites2Units(cellfun(@(si) size(si,2),siteChannels{lsInds})',...
        siteLocation'),NaN,alignLimits,strcat(savePath,"Spheres\"),...
        strjoin([monk,"LS_Grasp"],''),[5 35]);

    countAUC(sReps,muspInds,{spVals},mMapping,siteMasks(mInds'),vmMask,...
        siteActiveInd{2},savePath+"Spheres\",strjoin([monk],''),"Grasp",m<=length(monkeys),...
        [1 20],{[1 7]});

    mMapping = sMapping;
    mMapping(~muspInds) = 0;
    [siteIndsN,siteInds] = unique(mMapping);
    siteInds = siteInds(siteIndsN>0);
    %    peakTimeHistograms(params.bins,sReps,corrsP,pSegs(:,1),siteInds,savePath+"ReachCorrs\Histograms\",strjoin(monk,''))
    %     FRFig = modulatedUnitsPerRep(allReps(muInds),corrsP(muInds),strjoin(["AUC:",...
    %         monk,"Corrs","Task phase"],' '),strjoin([monk,"Corrs"],''));
    %     saveFigures(FRFig, strjoin([savePath,"ReachCorrs\"],''),...
    %        strjoin([monk,"Corrs"],''),[]);
    if (m<=length(monkeys))
        FRMapFig = mapUnitVals(vmMask,siteMasks(mInds'),...
            cellfun(@nanmedian,tTime(mInds)),~cellfun(@any,tTime(mInds)),...
            false,4,[0 .1]);
        title(strjoin([monk,"Corrs Task Phase"],' '));
    end
end
close all;
clear muspInds spVals mMapping currReps mInds muInds mRange  ssMapping
%%
