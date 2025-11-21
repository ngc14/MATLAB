conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
cLim = [0 .5];
numChannels = 32;
channelRes = 1;
condColors = [[1 0 0]; [.9 .7 0]; [0 .3 1]];
corrPhases = ["Go","Reach","Grasp"];
groupLabels = ["S","M","D"];
params = PhysRecording(conditions, .01,.15, -6, 5,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
tRange = -.5:params.binSize:1;
mainDir = "S:\Lab\";
monkey = "Gilligan";
saveDir= strcat(mainDir,"ngc14\Working\Correlations\");
saveFig = false;
%%
[siteDateMap, vMask, ~] = getMonkeyInfo(mainDir,monkey,"M1",true);
sessionDates = cellfun(@(d) datetime(d,'Format','MM_dd_yyyy'), table2cell(siteDateMap(:,'Date')));
[sessionCorrs,sessionGo,sessionReach,sessionGrasp,avgTimeCorr,sessionChannels,allFR] = ...
    loadCorrelations(sessionDates,params,mainDir,monkey,false);
%%
sessionCorrs(:,length(conditions)+1) = cellfun(@(n) mean(cat(3,n{1:length(conditions)-1}),3,'omitnan'), num2cell(sessionCorrs,2), 'UniformOutput',false);
concatMatrix = cell(1,length(conditions));
noiseCorrsConds = cell(1,length(conditions)-1);
resMap = 1:channelRes:numChannels;
chNames = arrayfun(@num2str,resMap,'UniformOutput',false);
allCondFR = NaN(length(resMap),length(sessionCorrs),length(conditions)-1,length(corrPhases));
for c = 1:length(conditions)
    [sigAll,goAll,reachAll,graspAll] = deal(NaN(length(resMap),length(resMap),length(sessionCorrs)));
    for n = 1:length(sessionCorrs)
        offDiagMat = (~diag(ones(length(sessionChannels{n}),1))./~diag(ones(length(sessionChannels{n}),1)));
        if(~isempty(offDiagMat))
            sessionChannel =  discretize(sessionChannels{n},[resMap,numChannels]);
            if(c<length(conditions))
                sessGo = abs(sessionGo{n}(:,:,c)).*offDiagMat;
                sessReach = abs(sessionReach{n}(:,:,c)).*offDiagMat;
                sessGrasp = abs(sessionGrasp{n}(:,:,c)).*offDiagMat;
            end
            sessCorrs = abs(sessionCorrs{n,c}).*offDiagMat;
            for h = 1:length(sessionChannel)
                xMapped=sessionChannel==sessionChannel(h);
                if(c<length(conditions))
                    allCondFR(sessionChannel(h),n,c,:) = mean(cell2mat(cellfun(@(h) h(xMapped,:),...
                        num2cell(allFR{n}{c}(~any(isnan(allFR{n}{c}),2),:),1),'UniformOutput',false)),1,'omitnan');
                end
                for i = 1:length(sessionChannels{n})
                    yMapped=sessionChannel==sessionChannel(i);
                    if(c<length(conditions))
                        goAll(sessionChannel(h),sessionChannel(i),n) = mean(sessGo(xMapped,yMapped),'all','omitnan');
                        reachAll(sessionChannel(h),sessionChannel(i),n) = mean(sessReach(xMapped,yMapped),'all','omitnan');
                        graspAll(sessionChannel(h),sessionChannel(i),n) = mean(sessGrasp(xMapped,yMapped),'all','omitnan');
                    end
                    sigAll(sessionChannel(h),sessionChannel(i),n) =  mean(sessCorrs(xMapped,yMapped),'all','omitnan');
                end
            end
        end
    end
    concatMatrix{c} = sigAll;
    if(c<length(conditions))
        noiseCorrsConds{c} = [{goAll}, {reachAll}, {graspAll}];
    end
end
numPairs = sum(~isnan(sigAll),3);
maxPairs = ceil(mean(numPairs(:),'omitnan')+std(numPairs(:),'omitnan'));
maxFR = mean(allCondFR(:),'omitnan')+std(allCondFR(:),'omitnan');
avgConcat = mean(cat(4,concatMatrix{1:length(conditions)-1}),4,'omitnan');
goodSess = find(~cellfun(@isempty,avgTimeCorr));
[maxVals,maxInds]= cellfun(@(c) max(cat(3,c{:}),[],4), avgTimeCorr(goodSess), 'UniformOutput', false);
maxVals = cell2mat(cellfun(@(m) reshape(m(repmat(tril(true(size(m,[1 2])),-1),1,1,size(m,3))),[],3), maxVals,'UniformOutput',false));
maxInds = cell2mat(cellfun(@(m) reshape(m(repmat(tril(true(size(m,[1 2])),-1),1,1,size(m,3))),[],3), maxInds,'UniformOutput',false));
noiseCorrMatrix = cellfun(@(c) cellfun(@(m) mean(m,3,'omitnan'),c,'UniformOutput',false), noiseCorrsConds, 'UniformOutput',false);
%%
figure(); colormap([0 0 0;colormap('jet')]); tiledlayout(3,length(conditions)+2);
plotHeatMatrix(cellfun(@(m) mean(m,3,'omitnan'),[concatMatrix,{avgConcat}],'UniformOutput',false),[conditions,"Average"],chNames,cLim);
nexttile;hold on; axis ij; axis tight;title("Average Marginals");
imagesc(mean(mean(avgConcat,3,'omitnan').*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim);colorbar;

avgNoiseCorr = cellfun(@(c) mean(cat(3,c{:}),3,'omitnan'),num2cell(cat(1,noiseCorrMatrix{:}),1),'UniformOutput',false);
plotHeatMatrix(avgNoiseCorr,corrPhases,chNames,cLim./2);
for m = 1:length(avgNoiseCorr)
    nexttile; hold on; axis tight; axis ij;title(corrPhases(m)+" marginals");
    imagesc(mean(avgNoiseCorr{m}.*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
    xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim./2);
end
colorbar;nexttile([1,1]);hold on; axis image; axis ij;title("NumPairs");
imagesc(numPairs);
yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));xticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxPairs)];cb.TickLabels = num2cell([1,round(maxPairs)]);clim([0 round(maxPairs)]);

nexttile([1,1]);hold on; axis ij; axis tight;title("FR");
imagesc(mean(allCondFR(:,:,1:3),[2,3],'omitnan'));
xticks([]);xlim([.5 1]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));
cb=colorbar(gca,'southoutside');cb.Ticks=[0,(maxFR)];cb.TickLabels = num2cell([1,(maxFR)]);clim(cLim);

chSteps = 8;
ind = find(~isnan(avgConcat));
[r,c,v] = ind2sub(size(avgConcat),ind);
[~,e] = discretize(1:length(chNames),chSteps);
pairConds = cell(1,length(e));
for i = 1:length(ind)
    pairConds{find(~(abs(r(i)-c(i))>e),1)}(end+1) = avgConcat(r(i),c(i),v(i));
end
ax =nexttile([1,4]);hold on; axis ij; axis tight;title("Average Signal Corrs by Distance");
boxchart(cell2mat(cellfun(@(p,s) repmat(p,length(s),1), num2cell(1:length(pairConds)),pairConds,'UniformOutput',false)'),...
    cell2mat(pairConds),'Notch','on');
xticks(1:length(e));xlabel("Distance"); ylim([0 1]); set(ax,'YDir','normal');
xticklabels(["Same channel",arrayfun(@(s,a) strcat(num2str(s.*100),"-",num2str(a.*100)),e(1:end-1),[e(2:end)])]);
if(saveFig)
    saveFigures(gcf,saveDir,"All",[]);
end
%%
figure(); tiledlayout(10,10,'TileSpacing','tight');
for n = 1:length(goodSess)
    nexttile; hold on; axis tight;
    sessTimeCorr = squeeze(num2cell(cat(4,avgTimeCorr{goodSess(n)}),[1 2 4]));
    allTimeCorrs = cellfun(@(m) m.*(repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))./...
        repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))),[avgTimeCorr{goodSess(n)}], 'UniformOutput',false);
    cellfun(@(cm,cl) shadedErrorBar(tRange,squeeze(mean(cm,[1,2,3],'omitnan'))',squeeze(std(cm,0,[1,2,3],'omitnan'))','lineProps',{'Color',cl,'LineWidth',1.5}), ...
        squeeze(num2cell(cat(3,allTimeCorrs{:}),[1 2 4])),num2cell(condColors,2));
    plot(tRange,zeros(1,length(tRange)),'k');
    plot([0 0],[-1 1],'--k','LineWidth',1);
    title(string(sessionDates(goodSess(n))));
    xticks(tRange([1,end]));
    ylim(cLim);
    xticklabels(arrayfun(@(s) num2str(s,'%.1f'),tRange([1,end]),'UniformOutput',false));
end
if(saveFig)
    saveFigures(gcf,saveDir+"TimeNoiseCorr\","All",[]);
end
%%
uniqueCombs=bsxfun(@minus, nchoosek(1:3+2-1, 2), 0:2-1);
uniqueCombs = 1;
groupColors = distinguishable_colors(length(corrPhases));
corrGroups = cell(length(conditions)-1,length(uniqueCombs));
sp = length(corrPhases)+1;
figure(); hold on;
title("Conditions by Phase");
for c = 1:length(conditions)-1
    for u = 1:length(uniqueCombs)
        fullMat = cellfun(@(cs) cs(:,:,:),noiseCorrsConds{c},'UniformOutput',false); %1+(11*(g2-1)):min(numChannels,11*g2)
        corrGroups{c,u} = cell2mat(cellfun(@(f)f(:),fullMat,'UniformOutput',false));
    end
    f=swarmchart(repmat(linspace(-1,1,length(corrPhases))+sp*c,size(corrGroups{c},1),1),...
        corrGroups{c},[],groupColors,'.','XJitter','Density');
    arrayfun(@(m,x) plot([x-.5 x+.5], [m m], 'k','LineWidth',2), mean(corrGroups{c},1,'omitnan'), linspace(-1,1,length(corrPhases))+sp*c);
    if(c==1)
        leg = f;
    end
    f = f(1).Parent;
end
legend(leg,corrPhases);
xlim([sp-(sp/2) (sp*length(uniqueCombs)*length(conditions)-1)-1]);
xticks(sp.*[1:length(conditions)-1]);
xticklabels(conditions);
ylim(cLim./1);
if(saveFig)
    saveFigures(gcf,saveDir,"Noise Distributions",[]);
end

uniqueCombs = [1 2 3];
groupNames = cellfun(@(h) string(strcat(h{:})), num2cell(groupLabels(uniqueCombs)',2));
figure(); tiledlayout(1,length(conditions)-1);
for c = 1:length(conditions)-1
    for u = 1:length(uniqueCombs)
        g1=uniqueCombs(u);
        %g2=uniqueCombs(u,2);
        fullMat = cellfun(@(cs) cs(1+(11*(g1-1)):min(numChannels,11*g1),:,:),noiseCorrsConds{c},'UniformOutput',false); %1+(11*(g2-1)):min(numChannels,11*g2)
        corrGroups{c,u} = cell2mat(cellfun(@(f)f(:),fullMat,'UniformOutput',false));
    end
    nexttile;hold on; title(conditions(c)+" Noise Correlations");
    for a = 1:length(uniqueCombs)
        f=swarmchart(repmat(linspace(-1,1,length(corrPhases))+sp*a,size(corrGroups{c,a},1),1),...
        corrGroups{c,a},[],groupColors,'.','XJitter','Density');
        arrayfun(@(m,x) plot([x-.5 x+.5], [m m], 'k','LineWidth',2), mean(corrGroups{c,a},1,'omitnan'), ...
               linspace(-1,1,length(corrPhases))+sp*a);
        if(a==1 && c==1)
            leg = f;
        end
    end
    legend(leg,corrPhases);
    xlim([sp-(sp/2) (sp/2)+(sp*length(uniqueCombs))]);
    xticks(sp:sp:sp*length(uniqueCombs))
    xticklabels(groupNames);
    ylim(cLim./1);
end
if(saveFig)
    saveFigures(gcf,saveDir+"Compartments\","Phases_Layers",[]);
end

groupColors = distinguishable_colors(length(uniqueCombs),groupColors);
sp = length(uniqueCombs)+2;
figure(); tiledlayout(1,length(conditions)-1);
for c = 1:length(conditions)-1
    nexttile;hold on; title(conditions(c));
    for a=1:length(corrPhases)
        f=violin(cell2mat(cellfun(@(m) [m(:,a);NaN(max(cellfun(@(s) size(s,1),corrGroups(c,:)))-size(m,1),1)], corrGroups(c,:),'UniformOutput',false)),...
            'x',linspace(-1,1,length(uniqueCombs))+sp*a,'medc',[],'facecolor',groupColors,'plotlegend',0,'ax',gca);
        f=swarmchart(repmat(linspace(-1,1,length(uniqueCombs))+sp*a,size(corrGroups{c,a},1),1),...
        corrGroups{c,a},[],groupColors,'.','XJitter','Density');
        arrayfun(@(m,x) plot([x-.5 x+.5], [m m], 'k','LineWidth',2), mean(corrGroups{c,a},1,'omitnan'), ...
               linspace(-1,1,length(corrPhases))+sp*a);
        if(c==1 && a==1)
            leg = f;
        end
    end
    legend(leg,groupNames)
    xlim([sp-(sp/2) (sp/2)+(sp*length(corrPhases))]);
    xticks(sp:sp:sp*length(corrPhases))
    xticklabels(corrPhases);
    ylim(cLim./1);
end
if(saveFig)
    saveFigures(gcf,saveDir+"Compartments\","Layers_Phases",[]);
end

%%
conditions(5) = "Average";
concatMatrix{5} = avgConcat;
noiseCorrConds{4} = mean(cell2mat(cat(4,noiseCorrsConds{1:3})),4,'omitnan');
figure(); tiledlayout(1,length(conditions));
for c = 1:length(conditions)
    t=nexttile;
    hold on; title(conditions(c));
    condSig = concatMatrix{c};
    cg = cell(1,length(uniqueCombs));
    for u = 1:length(uniqueCombs)
        g1=uniqueCombs(u);
        %g2=uniqueCombs(u,2);
        fullMat = condSig(1+(11*(g1-1)):min(size(condSig,1),11*g1),:,:);%1+(11*(g2-1)):min(size(condSig,1),11*g2)
        cg{u} = fullMat(:);
    end
    condGroups{c} = cell2mat(cellfun(@(l) [l;NaN(max(cellfun(@length,cg))-length(l),size(l,2))],cg,'UniformOutput',false));
    violin(cell2mat(cellfun(@(l) [l;NaN(max(cellfun(@length,cg))-length(l),size(l,2))],cg,'UniformOutput',false)),...
        'xlabel',{groupNames},'medc',[],'facecolor',distinguishable_colors(length(uniqueCombs)),'plotlegend',0,'ax',t); ylim([0 1]);
end
if(saveFig)
    saveFigures(gcf,saveDir+"Compartments\","Signal",[]);
end
%%
figure(); tiledlayout(1,length(conditions)-1);
for c = 1:length(conditions)-1
    nexttile; hold on;
    if(c==length(conditions)-1)
        title("Average"); scatter(tRange(round(mean(maxInds,2,'omitnan'))),mean(maxVals,2,'omitnan'));
    else
        title(conditions(c)); scatter(tRange(maxInds(:,c)),maxVals(:,c))
    end
    ylim(cLim);
end
if(saveFig)
    saveFigures(gcf,saveDir+"Compartments\","Scatter",[]);
end
%%
maxGroup = max(cellfun(@(s) max(sum(~isnan(s),1)), condGroups));
condGroups = cellfun(@(m) cell2mat(cellfun(@(l) [l(~isnan(l));NaN(maxGroup-sum(~isnan(l)),1)], ...
    num2cell(m,1),'Uniformoutput',false)),condGroups, 'UniformOutput',false);
condTables = cellfun(@(a) array2table(a,'VariableNames',groupNames),condGroups,'UniformOutput',false);
allTables = cellfun(@(a,n) renamevars(a,a.Properties.VariableNames,cellfun(@(s) strcat(s,"_",n),a.Properties.VariableNames)), ...
    condTables, conditions, 'UniformOutput',false);
allTables = stack(cat(2,allTables{:}),arrayfun(@(g) (1+(length(groupNames)*(g-1))):(length(groupNames)*g),1:length(conditions),'UniformOutput',false),...
    'NewDataVariableName',conditions,'IndexVariableName','Group');
allTables.Group = categorical(groupNames(allTables.Group));
allTables = allTables(~any(isnan(allTables{:,2:end}),2),:);

allCondPhases=cell(1,length(conditions)-1);
maxPhases = max(cellfun(@(n) max(sum(~isnan(n),1)), corrGroups),[],'all');
for c =1:length(conditions)-2
    phaseCorrs = cellfun(@(l) [l;NaN(max(cellfun(@(s) size(s,1),corrGroups(c,:)))-size(l,1),size(l,2))],corrGroups(c,:),'UniformOutput',false);
    phaseGroups = cellfun(@(m) cell2mat(cellfun(@(l) [l(~isnan(l));NaN(max(...
        cellfun(@(s) max(sum(~isnan(s),1)), phaseCorrs))-sum(~isnan(l)),1)],...
        num2cell(m,1),'UniformOutput',false)), phaseCorrs, 'UniformOutput',false);
    phaseTables = cellfun(@(a) array2table(a,'VariableNames',corrPhases), phaseGroups, 'UniformOutput',false)';
    phaseTables = cellfun(@(a,n) renamevars(a,a.Properties.VariableNames, cellfun(@(s) strcat(s,"_",n),a.Properties.VariableNames)),phaseTables,groupNames, 'UniformOutput', false);
    allPhases = stack(cat(2,phaseTables{:}),arrayfun(@(g) (1+(length(corrPhases)*(g-1))):(length(corrPhases)*g),...
        1:length(phaseTables),'UniformOutput',false),'IndexVariable','Phase');
    allPhases.Phase = categorical(corrPhases(allPhases.Phase)');
    allPhases = allPhases(~any(isnan(allPhases{:,2:end}),2),:);
    allPhases = renamevars(allPhases,allPhases.Properties.VariableNames,...
        string(arrayfun(@(g) strcat(g,"_",params.condAbbrev(conditions(c))),["Phase";groupNames],'UniformOutput',false)));
    allCondPhases{c} = [allPhases];...cell2table([num2cell(categorical(repmat("",maxPhases-height(allPhases),1))),...
        %num2cell(NaN(maxPhases-height(allPhases),length(groupNames)))],'VariableNames',allPhases.Properties.VariableNames)];
end

siteMasks = repmat({},1,height(siteDateMap));
mm = MotorMapping(35);
mRefMask = vMask{1};
[verticies, vCells] = voronoin(fliplr([siteDateMap{:,'x'},siteDateMap{:,'y'}; ...
    [0 size(mRefMask,2); size(mRefMask,1) 0; 0 0;size(mRefMask,1) size(mRefMask,2)]]));
for i = 1:height(siteDateMap)
    currSite = [siteDateMap{i,'x'},siteDateMap{i,'y'}];
    tempCircle = zeros(size(mRefMask)+2*mm.tileBuffer);
    tempCircle((currSite(2)-mm.siteRadius+mm.tileBuffer):(currSite(2)+...
        mm.siteRadius+mm.tileBuffer),(currSite(1)-mm.siteRadius+mm.tileBuffer):...
        (currSite(1)+mm.siteRadius+mm.tileBuffer))= mm.poolCircle;
    tempCircle = tempCircle(mm.tileBuffer:end-(mm.tileBuffer+1),...
        mm.tileBuffer:end-(mm.tileBuffer+1));
    siteMasks{i} = tempCircle & poly2mask(verticies(vCells{i},2),...
        verticies(vCells{i},1),size(tempCircle,1),size(tempCircle,2));
end
clear tempCircle verticies vCells poolCircle; 
vM = vMask{1};
vM(:,round(size(vM,2)/2):end) = [];
figure(); tiledlayout(1,length(conditions));
for c = 1:length(conditions)
    nexttile; hold on; title(conditions(c));
    vals = squeeze(mean(concatMatrix{c},[1 2],'omitnan'));
    im = mapUnitVals(vM,siteMasks,vals,[],false,255,cLim/.75); %setdiff(1:height(siteDateMap),goodSess)
    imshow(im);
    if(c<length(conditions))
        colorbar off;
    end
end
if(saveFig)
    saveFigures(gcf,saveDir+"Maps\","SignalCorr",[]);
end
im = {}; in = 1;
for c = 1:length(conditions)-2
    for p = 1:length(corrPhases)
        vals = squeeze(mean(noiseCorrsConds{c}{p},[1 2],'omitnan'));
        im{in} = mapUnitVals(vM,siteMasks,vals,[],false,255,cLim./1.5);
        in = in+1;
    end
end
figure();
imshow(imtile(im,[0 0 0; colormap('parula')],'GridSize',[3 3],'ThumbnailSize',[768 383]));
title(strjoin(corrPhases,'\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'));hold on;
ylabel(strjoin(conditions(1:end-2),'\n\n\n\n\n\n\n\n\n\n\n\n'),'Rotation',0);
cb=colorbar;
cb.Ticks = linspace(0,1,5);
cb.TickLabels = arrayfun(@(s) num2str(s,'%0.2f'),linspace(cLim(1)/1.5,cLim(end)/1.5,5),'UniformOutput',false);
if(saveFig)
    saveFigures(gcf,saveDir+"Maps\","Phases",[]);
end

S = (squeeze(cellfun(@(b) [b(~isnan(b));NaN(291-sum(~isnan(b),'all'),1)], num2cell(allCondFR(1:11,:,:,:),[1 2]), 'UniformOutput', false)));
M = (squeeze(cellfun(@(b) [b(~isnan(b));NaN(291-sum(~isnan(b),'all'),1)], num2cell(allCondFR(12:23,:,:,:),[1 2]), 'UniformOutput', false)));
D = (squeeze(cellfun(@(b) [b(~isnan(b));NaN(291-sum(~isnan(b),'all'),1)], num2cell(allCondFR(24:32,:,:,:),[1 2]), 'UniformOutput', false)));
S = array2table(cell2mat(S(:)'),'VariableNames',cell2mat(arrayfun(@(s) cell2mat(cellfun(@(a) strcat(s,a), corrPhases, 'UniformOutput',false)), conditions(1:3),'UniformOutput',false)));
M = array2table(cell2mat(M(:)'),'VariableNames',cell2mat(arrayfun(@(s) cell2mat(cellfun(@(a) strcat(s,a), corrPhases, 'UniformOutput',false)), conditions(1:3),'UniformOutput',false)));
D = array2table(cell2mat(D(:)'),'VariableNames',cell2mat(arrayfun(@(s) cell2mat(cellfun(@(a) strcat(s,a), corrPhases, 'UniformOutput',false)), conditions(1:3),'UniformOutput',false)));
