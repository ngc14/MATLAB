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
cLim = [-1 1];
condColors = [[1 0 0]; [.9 .7 0]; [0 .3 1]];
mainDir = "S:\Lab\";
monkey = "Gilligan";
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
corrPhases = ["Go","Reach","Grasp"];
goInd = cellfun(@(as) find(contains(as,"Go")),allSegs,'UniformOutput',false);
reachInd = cellfun(@(as) find(contains(as,"Reach")),allSegs,'UniformOutput',false);
graspInd = cellfun(@(as) find(contains(as,"Hold"),1),allSegs,'UniformOutput',false);
saveFig = false;
%%
[monkeyTable, ~, ~] = getMonkeyInfo(mainDir,monkey,"M1",true);
sessionDates = cellfun(@(d) datetime(d,'Format','MM_dd_yyyy'), table2cell(monkeyTable(:,'Date'))); % randperm(height(monkeyTable),10)
sessionCorrs = cell(length(sessionDates),length(conditions));
[sessionGo,sessionReach,sessionGrasp,sessionChannels] = deal(cell(length(sessionDates),1));
[maxSegSz,maxSegL]= max(cellfun(@length,params.condSegMap.values));
maxSegL= cell2mat(params.condSegMap.values({params.condNames(maxSegL)}));
condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
hbar = parforProgress(length(sessionDates));
avgTimeCorr = cell(length(sessionDates),1);
saveDir= strcat(mainDir,"ngc14\Working\Correlations\");
parfor n = 1:length(sessionDates)
    [condPSTHS,alignedSpikes,alignedTimes,allGoodTrials,corrTrialMatrix,sigCorr,siteTrialSegs] = deal(repmat({[]},1,length(conditions)));
    [goMatrix,reachMatrix,graspMatrix] = deal([]);
    [spikes,times,~,allTrials,~,channels,~,~,chMap] = getAllTrials(mainDir+monkey+...
        "\"+"All Data\"+monkey+"_"+string(sessionDates(n))+"\Physiology\Results","Single",true);
    if(length(spikes)>=2)
        NaNGraspInd = false(size(allTrials,1),1);
        unitNames = arrayfun(@(s) arrayfun(@(u) strcat("c",num2str(s) ,"u",num2str(u)), ...
            1:sum(chMap{end}(channels)==s),'UniformOutput',false),unique(chMap{end}(channels),'stable'),'UniformOutput',false);
        unitNames = cell2mat([unitNames{:}]);
        [unitNames,sIndex] = natsort(cellstr(unitNames));
        spikes = spikes(sIndex);
        NaNGraspInd(ismember(allTrials(:,1),conditions)) = cellfun(@(s,t) sum(isnan(t))==1 & isempty(t(strcmp(params.condSegMap(s),"StartGrasp"))), ...
            allTrials(ismember(allTrials(:,1),conditions),1),times(ismember(allTrials(:,1),conditions))');
        successfulTrials = NaNGraspInd | cellfun(@(b,t) ~isnan(str2double(b)) | isempty(b),allTrials(:,end-1));
        [alignedSpikes, alignedTimes, allGoodTrials, condPSTHS] = ...
            getConditionTrialPSTHS(params,allTrials(successfulTrials,:),taskAlign,spikes,times(successfulTrials));
        baselineAlign = cellfun(@(c,u) cell2mat(cellfun(@(t,i) max(1,mean(u(:,findBins(t(1)-3,params.bins):...
            findBins(t(1)-3,params.bins)+(1/params.binSize),i),2,'omitnan')),c,num2cell(1:size(u,3)),'UniformOutput',false)),alignedTimes,condPSTHS,'UniformOutput',false);
        for b = 1:length(baselineAlign)
            condBaseline = baselineAlign{b};
            condBaseline(isnan(condBaseline)) = 1;
            baselineAlign{b} = condBaseline;
        end
        normPSTH = cellfun(@(p,b) permute(permute(p,[1 3 2])./b,[1 3 2]),condPSTHS,baselineAlign,"UniformOutput",false);
        [taskBaseline,taskFR] = calculatePhases(params,taskAlign,repmat({{[0 0]}},length(conditions),1),...
            cellfun(@(a) {{cell2mat(a')}}, alignedTimes,'Uniformoutput',false),cellfun(@(c) {c},condPSTHS,'UniformOutput',false),false,true);
        [~,taskUnits] = cellfun(@(tb,tc) cellfun(@(b,cn) ttestTrials(b,cn,1,true,pVal),...
            tb,tc,'UniformOutput',false),taskBaseline,taskFR,'UniformOutput', false);
        taskUnits = any(cell2mat([taskUnits{:}]),2);
        unitNames = unitNames(taskUnits);
        allGoodTrials = cellfun(@(a)  a(taskUnits),allGoodTrials,'UniformOutput',false);
        alignedSpikes =  cellfun(@(a) a(taskUnits), alignedSpikes, 'UniformOutput',false);
        normPSTH = cellfun(@(a) a(taskUnits,:,:),normPSTH,'UniformOutput',false);
        sessionChannels{n} = channels(taskUnits);
        allFR{n} = cellfun(@(t) t{1}{1}{1}(taskUnits,:), taskFR,'UniformOutput',false);
        numTrials = size(condPSTHS{1},3);
        for c = 1:length(conditions)
            siteTrialSegs{c} = num2cell(NaN(numTrials,maxSegSz),2);
            for i = 1:numTrials;siteTrialSegs{c}{i}(condSegMappedInds{c}) = alignedTimes{c}{i};end
        end
        sessionTimeSegs = cell2mat(cat(1,siteTrialSegs{1:length(conditions)-1}));
        moveCondsAll = cellfun(@(i) i(1):params.binSize:i(end),num2cell([min(sessionTimeSegs(:,1:find(contains(maxSegL,"Go"),1)),[],2),...
            max(sessionTimeSegs(:,1:find(contains(maxSegL,"Hold"),1)),[],2)+range(winSz)],2),'UniformOutput',false);
        [~,tRange] = max(cellfun(@length,moveCondsAll));
        tRange = -.5:params.binSize:1;%moveCondsAll{tRange};
        corrT = cell(length(tRange),1);
        for ti = 1:length(tRange)
            currWin = round(tRange(ti),2);
            currSpks = cellfun(@(as,r,at,ag) getSpikeCounts({as},{r},...
                {at},{ag},[currWin-range(winSz)/2, currWin+range(winSz)/2]),alignedSpikes(1:length(conditions)-1),...
                reachInd(1:length(conditions)-1),alignedTimes(1:length(conditions)-1),allGoodTrials(1:length(conditions)-1),'UniformOutput',false);
            condCurrMat = NaN(sum(taskUnits),sum(taskUnits),length(currSpks));
            for u = 1:sum(taskUnits)
                for i = 1:sum(taskUnits)
                    condCurrMat(u,i,:) = cell2mat(cellfun(@(cs) corr(cs(u,:)',cs(i,:)',Rows='pairwise'), currSpks, 'UniformOutput',false));
                end
            end
            corrT{ti} = condCurrMat;
        end
        figure();
        t=tiledlayout(sum(taskUnits),sum(taskUnits),'TileSpacing','none','Padding','tight');%ceil((sum(taskUnits).*((sum(taskUnits)+1)/2)-sum(taskUnits))/5),5);
        for u = 1:sum(taskUnits)
            nexttile(t,tilenum(t,1,u)); hold on; axis tight;
            margCorr = cellfun(@(cm) cm(u,:,:), corrT, 'UniformOutput',false);
            plot(tRange,zeros(1,length(tRange)),'k');
            cellfun(@(cm,cl) shadedErrorBar(tRange,mean(cm,2,'omitnan')',std(cm,0,2,'omitnan')','lineProps',{'Color',cl,'LineWidth',1.5}), ...
              squeeze(num2cell(cat(1,margCorr{:}),[1 2])),num2cell(condColors,2));
            plot([0 0],[-1 1],'--k','LineWidth',1);
            title(strcat([unitNames{u}," Marginal"]));                
            ylim([-1 1]);
            xticks(tRange(1:20:end));
            xticklabels([]);
            yticklabels([]);
            for i = 1:u-1
                tCorr = cellfun(@(cm) squeeze(cm(u,i,:)), corrT, 'UniformOutput',false);
                nexttile(tilenum(t,u,i)); hold on; axis tight;
                colororder(condColors)
                plot(tRange,zeros(1,length(tRange)),'k');
                plot(tRange,cell2mat(tCorr')','LineWidth',1.5);
                distGo = sessionTimeSegs(:,find(contains(maxSegL,"Go"),1));
                distHold= sessionTimeSegs(:,find(contains(maxSegL,"Hold"),1));
                patch(cell2mat(arrayfun(@(r) repmat(r,1,2),[mean(distGo,'omitnan')-std(distGo,0,'omitnan'),...
                    mean(distGo,'omitnan')+std(distGo,0,'omitnan')],'Uniformoutput',false)),[-1 1 1 -1],'k','FaceAlpha',.15,'EdgeColor','none');
                patch(cell2mat(arrayfun(@(r) repmat(r,1,2),[mean(distHold,'omitnan')-std(distHold,0,'omitnan'),...
                    mean(distHold,'omitnan')+std(distHold,0,'omitnan')],'Uniformoutput',false)),[-1 1 1 -1],'k','FaceAlpha',.15,'EdgeColor','none');
                plot([0 0],[-1 1],'--k','LineWidth',1);
                title(strcat([unitNames{u},',',unitNames{i}]));
                xticks(tRange(1:20:end));
                ylim([-1 1]);
                if(i~=1)
                    yticklabels([]);
                end
                if(u~=sum(taskUnits))
                    xticklabels([]);
                else
                    xticklabels(arrayfun(@(s) num2str(s,'%.1f'),tRange(1:20:end),'UniformOutput',false));
                end
            end
        end
        nexttile(t,tilenum(t,sum(taskUnits),sum(taskUnits))); hold on; axis tight;
        allTimeCorrs = cellfun(@(m) m.*(repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))./repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))),corrT, 'UniformOutput',false);
        plot(tRange,zeros(1,length(tRange)),'k');
        avgTimeCorr{n} =  squeeze(num2cell(cat(4,allTimeCorrs{:}),[1 2 4]));
        cellfun(@(cm,cl) shadedErrorBar(tRange,squeeze(mean(cm,[1,2,3],'omitnan'))',squeeze(std(cm,0,[1,2,3],'omitnan'))','lineProps',{'Color',cl,'LineWidth',1.5}), ...
            squeeze(num2cell(cat(4,allTimeCorrs{:}),[1 2 4])),num2cell(condColors,2));
        plot([0 0],[-1 1],'--k','LineWidth',1);
        title("Average");
        ylim([-1 1]);
        xticks(tRange(1:20:end));
        yticklabels([]);
        distGo = sessionTimeSegs(:,find(contains(maxSegL,"Go"),1));
        distHold= sessionTimeSegs(:,find(contains(maxSegL,"Hold"),1));
        patch(cell2mat(arrayfun(@(r) repmat(r,1,2),[mean(distGo,'omitnan')-std(distGo,0,'omitnan'),...
            mean(distGo,'omitnan')+std(distGo,0,'omitnan')],'Uniformoutput',false)),[-1 1 1 -1],'k','FaceAlpha',.15,'EdgeColor','none');
        patch(cell2mat(arrayfun(@(r) repmat(r,1,2),[mean(distHold,'omitnan')-std(distHold,0,'omitnan'),...
            mean(distHold,'omitnan')+std(distHold,0,'omitnan')],'Uniformoutput',false)),[-1 1 1 -1],'k','FaceAlpha',.15,'EdgeColor','none');
        plot([0 0],[-1 1],'--k','LineWidth',1);
        xticklabels(arrayfun(@(s) num2str(s,'%.1f'),tRange(1:20:end),'UniformOutput',false));

        if(saveFig)
            saveFigures(gcf,saveDir+"TimeNoiseCorr\",string(sessionDates(n)),[]);
        end
        %%
        if(0)%saveFig)
            figure();
            unitLabs = cell2mat(arrayfun(@(a) repmat(string(a),numTrials,1),unitNames,'UniformOutput',false)');
            for c = 1:length(conditions)-1
                plotJointPSTHS(params,{reshape(permute(normPSTH{c},[3 1 2]),[],length(params.bins),1)},...
                    {repmat(mean(cell2mat(cat(3,siteTrialSegs{1:length(conditions)-1})),3,'omitnan'),sum(taskUnits),1)},unitLabs,...
                    [any(diff(char(unitLabs)),2);1],[],{evalWindow},[0 5],cell2struct(num2cell(repmat(condColors(c,:),sum(taskUnits),1),2),unique(unitLabs,'stable')));
            end
            saveFigures(gcf,saveDir,"PSTHS",[]);
        end
        goSpk = getSpikeCounts(alignedSpikes(1:length(conditions)-1),goInd(1:length(conditions)-1),...
            alignedTimes(1:length(conditions)-1),allGoodTrials(1:length(conditions)-1),winSz);
        reachSpk = getSpikeCounts(alignedSpikes(1:length(conditions)-1),reachInd(1:length(conditions)-1),...
            alignedTimes(1:length(conditions)-1),allGoodTrials(1:length(conditions)-1),winSz);
        graspSpk = getSpikeCounts(alignedSpikes(1:length(conditions)-1),graspInd(1:length(conditions)-1),...
            alignedTimes(1:length(conditions)-1),allGoodTrials(1:length(conditions)-1),winSz);
        for u = 1:sum(taskUnits)
            xVals = cellfun(@(m) min(m(u,:)):max(m(u,:)),{goSpk,reachSpk,graspSpk}, 'UniformOutput',false);
            for i = 1:sum(taskUnits)
                avgTimes = cellfun(@(t) mean(cell2mat(t'),1,'omitnan'),alignedTimes, 'Uniformoutput',false);
                cm = cellfun(@(cp,t,n,s) (corr(mean(cp(u,findBins(t(find(ismember(s,n{1}),1)),params.bins):...
                    findBins(t(find(ismember(s,n{1}),1,'last')),params.bins),:),3,'omitnan')',mean(cp(i,...
                    findBins(t(find(ismember(s,n{1}),1)),params.bins):findBins(t(find(ismember(s,n{1}),1,'last')),params.bins),:),...
                    3,'omitnan')')),normPSTH,avgTimes,taskAlign.values,allSegs,'UniformOutput',false);
                for c = 1:length(conditions)
                    sigCorr{c}(u,i) = cm{c};
                end
                goMatrix(u,i) = (corr(goSpk(u,:)',goSpk(i,:)',Rows='pairwise'));
                reachMatrix(u,i) = (corr(reachSpk(u,:)',reachSpk(i,:)',Rows='pairwise'));
                graspMatrix(u,i) = (corr(graspSpk(u,:)',graspSpk(i,:)',Rows='pairwise'));
                if(0)%saveFig)
                    if(i==1);figure();end
                    subplot(ceil(sum(taskUnits)/5),5,i);hold on;
                    sc=cellfun(@(s,l,m) scatter(s(u,:)', s(i,:)','Marker',m,'MarkerEdgeColor',l,'LineWidth',.5,'AlphaData',1,'SizeData',20),...
                        {goSpk,reachSpk,graspSpk},{'g','m','k'},{'x','+','*'});
                    cellfun(@(p,s,l) plot(p,mean(s(i,isalmost(s(u,:),min(s(u,:)),1)),'omitnan')+(p.*s(u,i)),'Color',l,'LineWidth',2.5),...
                        xVals,{graspSpk,reachSpk,goSpk},{[.1 .1 .1], [.75 0 .75], [0 .6 0]});
                    cellfun(@(s) set(s,'XJitter','rand','YJitter','rand','XJitterWidth',.85,'YJitterWidth',.85),sc);
                    title(strcat(unitNames(i),": ",num2str(min(sum(~isnan(goSpk(u,:))),sum(~isnan(goSpk(i,:))))),...
                        " trials (",cell2mat(compose('%.2f; ',[goMatrix(u,i),reachMatrix(u,i),graspMatrix(u,i)])),")"));
                    if(i==u);legend([sc{:}],{'Go','Reach','Grasp'});end
                    if(i==sum(taskUnits));saveFigures(gcf,saveDir+"Unit_Correlations\",unitNames{u}+"_Correlations",[]);end
                end
                % ctm = {};
                % for t = 1:numTrials
                %     ctm{t} = cellfun(@(cp) corr(cp(u,binInds(1):binInds(end),t)',...
                %         cp(i,binInds(1):binInds(end),t)'),condPSTHS,'UniformOutput',true);
                % end
                % ctm = cell2mat(ctm');
                % corrTrialMatrix{c}(u,i,:) = ctm(:,c);
            end
        end
        sessionCorrs(n,:) = sigCorr;
        sessionGo{n} = goMatrix;
        sessionReach{n} = reachMatrix;
        sessionGrasp{n} = graspMatrix;
        if(saveFig)
            avgCorrs = mean(cat(3,sigCorr{1:length(conditions)-1}),3,'omitnan');
            figure();tiledlayout(2,length(conditions)+2);colormap('jet')
            plotHeatMatrix([sigCorr,avgCorrs],[conditions,"Average"],unitNames,cLim./2);
            nexttile; hold on; axis tight; axis ij;title("Average marginals");
            imagesc(mean(avgCorrs.*(~logical(diag(true(length(unitNames),1)))./~logical(diag(true(length(unitNames),1)))),2,'omitnan'));
            xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim(cLim);
            %
            noiseCorrMatrix = {goMatrix,reachMatrix,graspMatrix};
            plotHeatMatrix(noiseCorrMatrix,corrPhases,unitNames,cLim./2);
            for m = 1:length(noiseCorrMatrix)
                nexttile; hold on; axis tight; axis ij;colormap('jet');title(corrPhases(m)+" marginals");
                imagesc(mean(noiseCorrMatrix{m}.*(~logical(diag(true(length(unitNames),1)))./~logical(diag(true(length(unitNames),1)))),2,'omitnan'));
                xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim(cLim./2);
            end
            colorbar;
            saveFigures(gcf,saveDir+"Correlations\",string(sessionDates(n)),[]);
        end
    end
    send(hbar, n);
end
%%
close all;
concatMatrix = cell(1,length(conditions));
numChannels = 32;
channelRes = 1;
resMap = 1:channelRes:numChannels;
chNames = arrayfun(@num2str,resMap,'UniformOutput',false);
for c = 1:length(conditions)
    [sigAll,goAll,reachAll,graspAll] = deal(NaN(length(resMap),length(resMap),length(sessionCorrs)));
    allCondFR = NaN(length(resMap),length(sessionCorrs));
    for n = 1:length(sessionCorrs)
        sessCorrs = sessionCorrs{n,c};
        sessGo = sessionGo{n};
        sessReach = sessionReach{n};
        sessGrasp = sessionGrasp{n};
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

figure();
colormap([0 0 0;colormap('jet')]);
tiledlayout(3,length(conditions)+2);
avgConcat = mean(cat(4,concatMatrix{1:length(conditions)-1}),4,'omitnan');
plotHeatMatrix(cellfun(@(m) mean(m,3,'omitnan'),[concatMatrix,{avgConcat}],'UniformOutput',false),[conditions,"Average"],chNames,cLim);

nexttile;hold on; axis ij; axis tight;title("Average Marginals");
imagesc(mean(mean(avgConcat,3,'omitnan').*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim./2);colorbar;

noiseCorrMatrix = cellfun(@(m) mean(m,3,'omitnan'), {goAll,reachAll,graspAll}, 'UniformOutput',false);
plotHeatMatrix(noiseCorrMatrix,corrPhases,chNames,cLim./4);

for m = 1:length(noiseCorrMatrix)
    nexttile; hold on; axis tight; axis ij;title(corrPhases(m)+" marginals");
    imagesc(mean(noiseCorrMatrix{m}.*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
    xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim(cLim./4);
end
colorbar;

nexttile([1,1]);hold on; axis image; axis ij;title("NumPairs");
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

saveFigures(gcf,mainDir+"ngc14\Working\Correlations\","All",[]);
%%
figure();     tiledlayout(10,10,'TileSpacing','none')
goodSess = find(~cellfun(@isempty,avgTimeCorr));
tRange = -.5:params.binSize:1;
for n = 1:length(goodSess)
    nexttile; hold on; axis tight;
    sessTimeCorr = squeeze(num2cell(cat(4,avgTimeCorr{goodSess(n)}),[1 2 4]));
    allTimeCorrs = cellfun(@(m) m.*(repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))./...
        repmat(~diag(ones(1,size(m,1))),1,1,size(m,3))),[avgTimeCorr{goodSess(n)}], 'UniformOutput',false);
    cellfun(@(cm,cl) shadedErrorBar(tRange,squeeze(mean(cm,[1,2,3],'omitnan'))',squeeze(std(cm,0,[1,2,3],'omitnan'))','lineProps',{'Color',cl,'LineWidth',1.5}), ...
        squeeze(num2cell(cat(3,allTimeCorrs{:}),[1 2 4])),num2cell(condColors,2));
    plot(tRange,zeros(1,length(tRange)),'k');
    plot([0 0],[-1 1],'--k','LineWidth',1);
    title("Average");
    ylim([-1 1]);
    xticks(tRange([1,end]));
    yticklabels([]);
    plot([0 0],[-1 1],'--k','LineWidth',1);
    xticklabels(arrayfun(@(s) num2str(s,'%.1f'),tRange(1:20:end),'UniformOutput',false));
end
    saveFigures(gcf,mainDir,sessionDates(n),[]);

function spCounts= getSpikeCounts(spks,segInd,times,goodTrials,wSz)
spCounts = cellfun(@(c,i,t,g) cellfun(@(a,gi) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(i)+wSz(end) & s>at(i)+wSz(1)),...
    a,t,num2cell(gi)),c,g,'UniformOutput',false),spks,segInd,times,goodTrials,'UniformOutput',false);
spCounts = reshape(permute(cell2mat(permute(cat(3,spCounts{:}),[2 1 3])),[1 3 2]),size(goodTrials{1},2),[]);
end