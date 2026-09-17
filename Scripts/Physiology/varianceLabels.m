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
normPSTH = cellfun(@(s,b)cellfun(@(t,n) permute(permute(sqrt(t),[1 3 2])-sqrt(n),[1 3 2]), s, b, 'Uniformoutput', false), siteTrialPSTHS, normBaseline, 'UniformOutput',false);
%normPSTH = cellfun(@(s) cellfun(@(t) zscore(t,0,2), s, 'UniformOutput',false), siteTrialPSTHS, 'UniformOutput',false);
GT = cellfun(@(g) mean(cat(3,g{:}),3,'omitnan'),num2cell([normPSTH{:}],2),'UniformOutput',false);
GT = num2cell(vertcat(GT{:}),2);
normPSTH = cellfun(@(c) cellfun(@(s) resize(s,[size(s,1),size(s,2),75],'FillValue',NaN),c,'UniformOutput',false), normPSTH,'UniformOutput',false);
normPSTH = cellfun(@(c) num2cell(cat(1,c{:}),[2 3]), normPSTH, 'UniformOutput',false);
LO = cellfun(@(c) cellfun(@(g,a) mean(g-a,3,'omitnan'), c, GT, 'UniformOutput',false), normPSTH,'UniformOutput',false);
NO = cellfun(@(c,l) cellfun(@(t,g,n) t-(g+n), c,l,GT, 'UniformOutput',false), normPSTH, LO, 'UniformOutput',false);
%%
plotReps = ["Hand","Arm","Hand Arm"];
plotCombs = nchoosek( plotReps,2);
winRed = [-.5, 1.5];
fi = figure();
nt = tiledlayout(3,4,'TileIndexing','columnmajor');
timePeriod = winRed(1):params.binSize:winRed(end);
for nc = 1:4
    if(nc==1)
        projT = cellfun(@(p) p(:,findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)),GT,'UniformOutput',false);
        plotName = margNames{nc};
        st = mean(cell2mat(vertcat(siteSegs{:})),1,'omitnan');
    else
        projT = cellfun(@(p) p(:,findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)),LO{nc-1},'UniformOutput',false);
        plotName = "Condition-" + params.condAbbrev(params.condNames(nc-1));
        st = mean(cell2mat(siteSegs{nc-1}),1,'omitnan');
    end
    st = findBins(st([2,3,6,7]),params.bins(findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)));
    for i = 1:3
        prevRep = [];
        nexttile();
        hold on;
        if(nc==1)
            ylabel(i)
        end
        if(i==1)
            title(plotName);
        end
        for s = length(plotReps):-1:1
            sp =  strsplit(plotReps(s)," ");
            if(isscalar(sp))
                if(strcmp(sp,"Arm"))
                    colors = [0 .85 .4];
                else
                    colors = [1 0 .8];
                end
            else
                colors = [.7 .7 .7];
            end
            somaIndex = contains(string(unitSomatotopy),sp) & ~cellfun(@(p) all(isnan(p)), projT);
            [~,splitProj] = pca(cell2mat(cellfun(@(m) m-mean(cell2mat(projT(somaIndex)),1),projT(somaIndex),...
                'UniformOutput',false))','Economy',false,'Centered','on','Algorithm','svd');
            [~,maxI] = max(abs(splitProj),[],1);
            if(isempty(prevRep))
                orientPlot = splitProj(:,i).*sign(splitProj(maxI(i),i))';
                prevRep = orientPlot;
            else
                orientPlot = {splitProj(:,i), -splitProj(:,i)};
                [~,minInd] = min(cellfun(@(s) sum((s-prevRep).^2),orientPlot));
                orientPlot = orientPlot{minInd};
            end
            plot(timePeriod,orientPlot,'Color',colors,'LineStyle','-','LineWidth',2.5+(1*isscalar(sp)));
            if(s==1)
                arrayfun(@(t) plot([timePeriod(t),timePeriod(t)],[-30 40], 'k--'), st(st<length(timePeriod)));
            end
        end
    end
end
for f = length(plotReps):-1:1
fi = figure();
currComb = plotCombs(f,:);
nt = tiledlayout(3,4,'TileIndexing','columnmajor');
for nc = 1:4
    if(nc==1)
        projT = cellfun(@(p) p(:,findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)),GT,'UniformOutput',false);
        plotName = margNames{nc};
        st = mean(cell2mat(vertcat(siteSegs{:})),1,'omitnan');
    else
        projT = cellfun(@(p) p(:,findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)),LO{nc-1},'UniformOutput',false);
        plotName = "Condition-" + params.condAbbrev(params.condNames(nc-1));
        st = mean(cell2mat(siteSegs{nc-1}),1,'omitnan');
    end
    st = findBins(st(1,[3,6]),params.bins(findBins(winRed(1),params.bins):findBins(winRed(end),params.bins)));
    prevRep = [];
    for i = 1:3
        nexttile();
        hold on;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        if(i==1)
            title(plotName);
            view(0,90);
        elseif(i==2)
            view(0,0);
        elseif(i==3)
            view(90,0);
        end
        xlim([-20 30]); ylim([-20 30]); zlim([-20 20]);
        for s = 2:-1:1
            sp = strsplit(currComb(s), " ");
            if(isscalar(sp))
                if(strcmp(sp,"Arm"))
                    colors = [0 .85 .4];
                else
                    colors = [1 0 .8];
                end
            else
                colors = [.7 .7 .7];
            end
            somaIndex = contains(string(unitSomatotopy),sp) & ~cellfun(@(p) all(isnan(p)), projT);
            [~,splitProj] = pca(cell2mat(cellfun(@(m) m-mean(cell2mat(projT(somaIndex)),1),projT(somaIndex),...
                'UniformOutput',false))','Economy',false,'Centered','on','Algorithm','svd');
            [~,maxI] = max(abs(splitProj),[],1);
            if(isempty(prevRep))
                orientPlot = splitProj(:,1:3).*sign(splitProj(maxI(1),1))';
                prevRep = orientPlot;
            else
                orientPlot = {splitProj(:,1:3), -splitProj(:,1:3)};
                [~,minInd] = min(cellfun(@(s) sum((s-prevRep).^2,'all'),orientPlot));
                orientPlot = orientPlot{minInd};
            end
            plot3(orientPlot(:,1),orientPlot(:,2),orientPlot(:,3),'Color',colors,'LineStyle','-','LineWidth',2.5);
            colors = [colors-[.45 .45 .45]; repmat(max([0 0 0],colors-[.1 .1 .1]),length(st),1)];
            scatter3(orientPlot([1,st],1) ,orientPlot([1,st],2),orientPlot([1,st],3),[120,50,50],colors,'filled','o');
        end
    end
end
end
%%
figure
nt = tiledlayout(2,3);
plotSegs = cellfun(@(s) mean(cell2mat([s{:}]'),1,'omitnan'), siteSegs, 'UniformOutput',false);
plotSegs = cellfun(@(c,s) c(contains(cell2mat(params.condSegMap.values({s})),phaseNames)), plotSegs, num2cell(params.condNames),'UniformOutput',false);
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
    plotMeanArm = mean(armM,find(~ismember(1:ndims(armM),2)),'omitnan');
    plotVarArm = std(armM,0,find(~ismember(1:ndims(armM),2)),'omitnan')./sqrt(size(armM,1));%sqrt(sum(~all(isnan(armM),2),'all'));
    ss=shadedErrorBar(params.bins,plotMeanArm,plotVarArm,'lineProps',{'Color',[0 .8 0],'LineWidth',2},'patchSaturation',0.2);
    handM = cell2mat(reshape(cellfun(@(m) vertcat(m{unitSomatotopy=="Hand" & tUnits}),currM,'UniformOutput',false),[ones(1,ndims(currM{1})),length(currM)])).^2;
    plotMeanHand = mean(handM,find(~ismember(1:ndims(handM),2)),'omitnan');
    plotVarHand = std(handM,0,find(~ismember(1:ndims(handM),2)),'omitnan')./sqrt(size(handM,1));%sqrt(sum(~all(isnan(handM),2),'all'));
    sh=shadedErrorBar(params.bins,plotMeanHand,plotVarHand,'lineProps',{'Color',[1 0 1],'LineWidth',2},'patchSaturation',0.2);
    xlim([-.5 1]); ylim([0 5]);
    arrayfun(@(x) plot([x,x],get(gca,'YLim'),'k--','LineWidth',1),mean(cell2mat(plotSegs'),1,'omitnan'));
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
    xlim([-.5 1]); ylim([0 1.5]);
    arrayfun(@(x) plot([x,x],get(gca,'YLim'),'k--','LineWidth',1),plotSegs{c});
    legend([ss.mainLine,sh.mainLine],["Arm","Hand"])
end
