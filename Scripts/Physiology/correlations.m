conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
evalWindow = [-1 2];
winSz = [0 .2];
pVal=0.05;
cLim = [0 .5];
numChannels = 32;
channelRes = 1;
condColors = [[1 0 0]; [.9 .7 0]; [0 .3 1]];
corrPhases = ["Go","Reach","Grasp"];
groupLabels = ["S","M","D"];
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
tRange = -.5:params.binSize:1;
[maxSegSz,maxSegL]= max(cellfun(@length,params.condSegMap.values));
maxSegL= cell2mat(params.condSegMap.values({params.condNames(maxSegL)}));
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
goInd = cellfun(@(as) find(contains(as,"Go")),allSegs,'UniformOutput',false);
reachInd = cellfun(@(as) find(contains(as,"Reach")),allSegs,'UniformOutput',false);
graspInd = cellfun(@(as) find(contains(as,"Hold"),1),allSegs,'UniformOutput',false);
mainDir = "S:\Lab\";
monkey = "Gilligan";
saveFig = false;
%%
saveDir= strcat(mainDir,"ngc14\Working\Correlations\");
[monkeyTable, ~, ~] = getMonkeyInfo(mainDir,monkey,"M1",true);
sessionDates = cellfun(@(d) datetime(d,'Format','MM_dd_yyyy'), table2cell(monkeyTable(:,'Date')));
[sessionCorrs,sessionGo,sessionReach,sessionGrasp,avgTimeCorr,sessionChannels,allFR] = ...
    loadCorrelations(sessionDates,cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false),saveFig);
%%
sessionCorrs(:,end+1) = cellfun(@(n) mean(cat(3,n{1:length(conditions)-1}),3,'omitnan'), num2cell(sessionCorrs,2), 'UniformOutput',false);
concatMatrix = cell(1,length(conditions));
resMap = 1:channelRes:numChannels;
chNames = arrayfun(@num2str,resMap,'UniformOutput',false);
for c = 1:length(conditions)
    [sigAll,goAll,reachAll,graspAll] = deal(NaN(length(resMap),length(resMap),length(sessionCorrs)));
    allCondFR = NaN(length(resMap),length(sessionCorrs));
    for n = 1:length(sessionCorrs)
        offDiagMat = (~diag(ones(length(sessionChannels{n}),1))./~diag(ones(length(sessionChannels{n}),1)));
        sessGo = abs(sessionGo{n}).*offDiagMat;
        sessReach = abs(sessionReach{n}).*offDiagMat;
        sessGrasp = abs(sessionGrasp{n}).*offDiagMat;
        sessCorrs = abs(sessionCorrs{n,c}).*offDiagMat;
        sessionChannel =  discretize(sessionChannels{n},[resMap,numChannels]);
        for h = 1:length(sessionChannel)
            xMapped=sessionChannel==sessionChannel(h);
            allCondFR(sessionChannel(h),n) = mean(cell2mat(cellfun(@(v) mean(v(xMapped,:),'all','omitnan'), allFR{n}, 'UniformOutput',false)),'omitnan');
            for i = 1:length(sessionChannels{n})
                yMapped=sessionChannel==sessionChannel(i);
                sigAll(sessionChannel(h),sessionChannel(i),n) =  mean(sessCorrs(xMapped,yMapped),'all','omitnan');
                goAll(sessionChannel(h),sessionChannel(i),n) = mean(sessGo(xMapped,yMapped),'all','omitnan');
                reachAll(sessionChannel(h),sessionChannel(i),n) = mean(sessReach(xMapped,yMapped),'all','omitnan');
                graspAll(sessionChannel(h),sessionChannel(i),n) = mean(sessGrasp(xMapped,yMapped),'all','omitnan');
            end
        end
    end
    concatMatrix{c} = sigAll;
end
numPairs = sum(~isnan(sigAll),3);
maxPairs = ceil(mean(numPairs(:),'omitnan')+std(numPairs(:),'omitnan'));
maxFR = mean(allCondFR(:),'omitnan')+std(allCondFR(:),'omitnan');
allSigCorrs = cat(4,concatMatrix{1:length(conditions)-1});
avgConcat = mean(allSigCorrs,4,'omitnan');
goodSess = find(~cellfun(@isempty,avgTimeCorr));
noiseCorrMatrix = cellfun(@(m) mean(m,3,'omitnan'), {goAll,reachAll,graspAll}, 'UniformOutput',false);
[maxVals,maxInds]= cellfun(@(c) max(cat(3,c{:}),[],4), avgTimeCorr(goodSess), 'UniformOutput', false);
maxVals = cell2mat(cellfun(@(m) reshape(m(repmat(tril(true(size(m,[1 2])),-1),1,1,size(m,3))),[],3), maxVals,'UniformOutput',false));
maxInds = cell2mat(cellfun(@(m) reshape(m(repmat(tril(true(size(m,[1 2])),-1),1,1,size(m,3))),[],3), maxInds,'UniformOutput',false));
%%
figure(); colormap([0 0 0;colormap('jet')]); tiledlayout(3,length(conditions)+2);
plotHeatMatrix(cellfun(@(m) mean(m,3,'omitnan'),[concatMatrix,{avgConcat}],'UniformOutput',false),[conditions,"Average"],chNames,cLim);
nexttile;hold on; axis ij; axis tight;title("Average Marginals");
imagesc(mean(mean(avgConcat,3,'omitnan').*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim./2);colorbar;

plotHeatMatrix(noiseCorrMatrix,corrPhases,chNames,cLim./4);
for m = 1:length(noiseCorrMatrix)
    nexttile; hold on; axis tight; axis ij;title(corrPhases(m)+" marginals");
    imagesc(mean(noiseCorrMatrix{m}.*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
    xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim./4);
end
colorbar;nexttile([1,1]);hold on; axis image; axis ij;title("NumPairs");
imagesc(numPairs);
yticks(resMap(1:2:end));yticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));xticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxPairs)];cb.TickLabels = num2cell([1,round(maxPairs)]);clim([0 round(maxPairs)]);

nexttile([1,1]);hold on; axis ij; axis tight;title("FR");
imagesc(mean(allCondFR,2,'omitnan'));
xticks([]);xlim([.5 1]);yticks(resMap(1:2:end));yticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));yticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxFR)];cb.TickLabels = num2cell([1,round(maxFR)]);clim([0 round(maxFR)]);

nexttile([1,4]);hold on; axis ij; axis tight;title("All Signal Averages");
imagesc(squeeze(mean(cat(4,concatMatrix{1:length(conditions)-1}),[2,4],'omitnan')));
yticks(resMap(1:2:end));yticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));yticklabels([]);
cb=colorbar(gca,'southoutside');clim(cLim);
if(saveFig)
    saveFigures(gcf,mainDir+"ngc14\Working\Correlations\","All",[]);
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
    saveFigures(gcf,mainDir+"ngc14\Working\Correlations\TimeNoiseCorr\","All",[]);
end
%%
conditions(5) = "Average";
allSigCorrs(:,:,:,5) = avgConcat;
figure(); tiledlayout(1,length(conditions));
uniqueCombs=bsxfun(@minus, nchoosek(1:3+2-1, 2), 0:2-1);
groupNames = cellfun(@(h) string(strcat(h{:})), num2cell(groupLabels(uniqueCombs),2));
groupColors = num2cell(distinguishable_colors(length(uniqueCombs)),2);
for c = 1:length(conditions)
    t=nexttile;
    hold on; title(conditions(c));
    corrSigs = allSigCorrs(:,:,:,c);
    corrGroups = cell(1,length(uniqueCombs));
    for u = 1:length(uniqueCombs)
        g1=uniqueCombs(u,1);
        g2=uniqueCombs(u,2);
        fullMat = corrSigs(1+(11*(g1-1)):min(size(corrSigs,1),11*g1),1+(11*(g2-1)):min(size(corrSigs,1),11*g2),:);
        corrGroups{u} = fullMat(:);
    end
    condGroups{c} = cell2mat(cellfun(@(l) [l;NaN(max(cellfun(@length,corrGroups))-length(l),size(l,2))],corrGroups,'UniformOutput',false));
    violin(cell2mat(cellfun(@(l) [l;NaN(max(cellfun(@length,corrGroups))-length(l),size(l,2))],corrGroups,'UniformOutput',false)),...
        'xlabel',{groupNames},'medc',[],'facecolor',distinguishable_colors(length(uniqueCombs)),'plotlegend',0,'ax',t); ylim(cLim);
end
if(saveFig)
    saveFigures(gcf,mainDir+"ngc14\Working\Correlations\Compartments\","Signal",[]);
end
%%
f=axes(figure());hold on; title("Phases");
corrSigs = [{goAll},{reachAll},{graspAll}];
corrGroups = cell(1,length(uniqueCombs));
for u = 1:length(uniqueCombs)
    g1=uniqueCombs(u,1);
    g2=uniqueCombs(u,2);
    fullMat = cellfun(@(cs) cs(1+(11*(g1-1)):min(size(corrSigs{1},1),11*g1),1+(11*(g2-1)):min(size(corrSigs{1},1),11*g2),:),corrSigs,'UniformOutput',false);
    corrGroups{u} = cell2mat(cellfun(@(f)f(:),fullMat,'UniformOutput',false));
end
phaseCorrs = cellfun(@(l) [l;NaN(max(cellfun(@(s) size(s,1),corrGroups))-size(l,1),size(l,2))],corrGroups,'UniformOutput',false);
sp = length(corrPhases)+1;
for a = 1:length(uniqueCombs)
    f=violin(phaseCorrs{a},'x',linspace(-1,1,length(corrPhases))+sp*a,'medc',[],'facecolor',groupColors{a},'plotlegend',0,'ax',f);
    f = f(1).Parent;
end
xlim([sp-(sp/2) (sp/2)+(sp*length(uniqueCombs))]);
xticks(sp:sp:sp*length(uniqueCombs))
xticklabels(groupNames);
ylim(cLim);
if(saveFig)
    saveFigures(gcf,mainDir+"ngc14\Working\Correlations\Compartments\","Phases",[]);
end
%%
figure(); tiledlayout(1,length(conditions)-1);
for c = 1:length(conditions)-1
    nexttile; hold on;
    if(c==length(conditions)-1)
        title("Average"); scatter(tRange(round(mean(cell2mat(maxInds),2,'omitnan'))),mean(cell2mat(maxVals),2,'omitnan'));
    else
        title(conditions(c)); scatter(tRange(maxInds{c}),maxVals{c})
    end
    ylim(cLim);
end
if(saveFig)
    saveFigures(gcf,mainDir+"ngc14\Working\Correlations\Compartments\","Scatter",[]);
end

maxGroup = max(cellfun(@(s) max(sum(~isnan(s),1)), condGroups));
condGroups = cellfun(@(m) cell2mat(cellfun(@(l) [l(~isnan(l));NaN(maxGroup-sum(~isnan(l)),1)], ...
    num2cell(m,1),'Uniformoutput',false)),condGroups, 'UniformOutput',false);
condTables = cellfun(@(a) array2table(a,'VariableNames',groupNames),condGroups,'UniformOutput',false);
allTables = cellfun(@(a,n) renamevars(a,a.Properties.VariableNames,cellfun(@(s) strcat(s,"_",n),a.Properties.VariableNames)), ...
    condTables, conditions, 'UniformOutput',false);
allTables = stack(cat(2,allTables{:}),arrayfun(@(g) (1+(length(groupNames)*(g-1))):(length(groupNames)*g),1:length(conditions)),...
    'NewDataVariableName',conditions,'IndexVariableName','Group');
allTables.Group = categorical(groupNames(allTables.Group));
allTables = allTables(~any(isnan(allTables{:,2:end}),2),:);

maxPhases = max(cellfun(@(s) max(sum(~isnan(s),1)), phaseCorrs));
phaseGroups = cellfun(@(m) cell2mat(cellfun(@(l) [l(~isnan(l));NaN(maxPhases-sum(~isnan(l)),1)],...
    num2cell(m,1),'UniformOutput',false)), phaseCorrs, 'UniformOutput',false);
phaseTables = cellfun(@(a) array2table(a,'VariableNames',corrPhases), phaseGroups, 'UniformOutput',false)';
allPhases = cellfun(@(a,n) renamevars(a,a.Properties.VariableNames, cellfun(@(s) strcat(s,"_",n),a.Properties.VariableNames)),...
    phaseTables,groupNames, 'UniformOutput', false);
allPhases = stack(cat(2,allPhases{:}), arrayfun(@(g) (1+(length(corrPhases)*(g-1))):(length(corrPhases)*g),1:length(groupNames),'UniformOutput',false), ...
    'NewDataVariableName', groupNames, 'IndexVariableName', 'Phase');
allPhases.Phase = categorical(corrPhases(allPhases.Phase)');
allPhases = allPhases(~any(isnan(allPhases{:,2:end}),2),:);
