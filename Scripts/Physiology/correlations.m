conditions = ["Extra Small Sphere", "Large Sphere", "Photocell","Rest"];
taskAlign = containers.Map(conditions,{{["GoSignal" "StartHold"]},{["GoSignal","StartHold"]},...
    {["GoSignal","StartHold"]},{["GoSignal","StartReplaceHold"]}});
binSize = .01;
smoothKernel = .15;
secBeforeAlignment = -6;
secAfterAlignment = 5;
evalWindow = [-1 2];
winSz = .2;
pVal=0.05;
mainDir = "S:\Lab\";
monkey = "Gilligan";
params = PhysRecording(conditions,binSize,smoothKernel,secBeforeAlignment,secAfterAlignment,...
    containers.Map(conditions,{"StartReach","StartReach","StartReach","GoSignal"}));
allSegs = cellfun(@(e) e(2:end-1), params.condSegMap.values,'UniformOutput',false);
corrPhases = ["Go","Reach","Grasp"];
goInd = cellfun(@(as) find(contains(as,"Go")),allSegs,'UniformOutput',false);
reachInd = cellfun(@(as) find(contains(as,"Reach")),allSegs,'UniformOutput',false);
graspInd = cellfun(@(as) find(contains(as,"Hold"),1),allSegs,'UniformOutput',false);
%%
[monkeyTable, ~, ~] = getMonkeyInfo(mainDir,monkey,"M1",true);
sessionDates = cellfun(@(d) datetime(d,'Format','MM_dd_yyyy'), table2cell(monkeyTable(:,'Date'))); % randperm(height(monkeyTable),10)
sessionCorrs = cell(length(sessionDates),length(conditions));
[sessionGo,sessionReach,sessionGrasp,sessionChannels] = deal(cell(length(sessionDates),1));
hbar = parforProgress(length(sessionDates));
for n = 1:length(sessionDates)
    dirPath = mainDir+monkey+"\"+"All Data\"+monkey+"_"+string(sessionDates(n))+"\Physiology\Results";
    saveDir = strcat(mainDir,"ngc14\Working\Correlations\",string(sessionDates(n)),"\");
    if(0)%exist(saveDir,'dir'))
        saveFig=false;
    else
        saveFig=true;
    end
    [condPSTHS,alignedSpikes,alignedTimes,allGoodTrials,corrTrialMatrix,sigCorr,siteTrialSegs] = deal(repmat({[]},1,length(conditions)));
    [goMatrix,reachMatrix,graspMatrix] = deal([]);
    [spikes,times,~,allTrials,~,channels,~,~,chMap] = getAllTrials(dirPath,"Single",true);
    NaNGraspInd = false(size(allTrials,1),1);
    if(length(spikes)>=2)
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
        %%
        if(0)%saveFig)
            figure();
            condColors = {[1 0 0]; [.9 .7 0]; [0 .3 1]};
            numTrials = size(condPSTHS{1},3);
            unitLabs = cell2mat(arrayfun(@(a) repmat(string(a),numTrials,1),unitNames,'UniformOutput',false)');
            [maxSegSz,maxSegL ]= max(cellfun(@length,params.condSegMap.values));
            maxSegL= params.condSegMap.values({params.condNames(maxSegL)});
            condSegMappedInds = cellfun(@(f) find(contains(maxSegL,f)), allSegs, 'UniformOutput', false);
            for c = 1:length(conditions)-1
                siteTrialSegs = num2cell(NaN(numTrials,maxSegSz),2);
                for i = 1:numTrials; siteTrialSegs{i}(condSegMappedInds{c}) = alignedTimes{c}{i};end
                plotJointPSTHS(params.bins,{reshape(permute(normPSTH{c},[3 1 2]),[],length(params.bins),1)},...
                    {repmat(cell2mat(vertcat(siteTrialSegs{:})),sum(taskUnits),1)},unitLabs,[any(diff(char(unitLabs)),2);1],...
                    [],{evalWindow},[0 5],cell2struct(num2cell(repmat(condColors{c},sum(taskUnits),1),2),unique(unitLabs,'stable')));
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
                cm = cellfun(@(cp,t,n,s) abs(corr(mean(cp(u,findBins(t(find(ismember(s,n{1}),1)),params.bins):...
                    findBins(t(find(ismember(s,n{1}),1,'last')),params.bins),:),3,'omitnan')',mean(cp(i,...
                    findBins(t(find(ismember(s,n{1}),1)),params.bins):findBins(t(find(ismember(s,n{1}),1,'last')),params.bins),:),...
                    3,'omitnan')')),normPSTH,avgTimes,taskAlign.values,allSegs,'UniformOutput',false);
                for c = 1:length(conditions)
                    sigCorr{c}(u,i) = cm{c};
                end
                goMatrix(u,i) = abs(corr(goSpk(u,:)',goSpk(i,:)',Rows='pairwise'));
                reachMatrix(u,i) = abs(corr(reachSpk(u,:)',reachSpk(i,:)',Rows='pairwise'));
                graspMatrix(u,i) = abs(corr(graspSpk(u,:)',graspSpk(i,:)',Rows='pairwise'));
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
            tiledlayout(2,length(conditions)+2);
            plotHeatMatrix([sigCorr,avgCorrs],[conditions,"Average"],unitNames,[0 .5]);
            nexttile; hold on; axis tight; axis ij;title("Average marginals");
            imagesc(mean(avgCorrs.*(~logical(diag(true(length(unitNames),1)))./~logical(diag(true(length(unitNames),1)))),2,'omitnan'));
            xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);colorbar;clim([0 1]);
            %
            noiseCorrMatrix = {goMatrix,reachMatrix,graspMatrix};
            plotHeatMatrix(noiseCorrMatrix,corrPhases,unitNames,[0, .5]);
            for m = 1:length(noiseCorrMatrix)
                nexttile; hold on; axis tight; axis ij;title(corrPhases(m)+" marginals");
                imagesc(mean(noiseCorrMatrix{m}.*(~logical(diag(true(length(unitNames),1)))./~logical(diag(true(length(unitNames),1)))),2,'omitnan'));
                xlim([0.5 1]);xticklabels([]);yticks(1:length(unitNames));yticklabels(unitNames);clim([0 .5]);
            end
            colorbar;
            saveFigures(gcf,saveDir,"Correlations_ABS",[]);
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
    concatMatrix{c} =mean(sigAll,3,'omitnan');
end
numPairs = sum(~isnan(sigAll),3);
maxPairs = max(numPairs,[],'all');
maxFR = mean(allCondFR,'all','omitnan')+std(allCondFR(:),'omitnan');

figure();
colormap([0 0 0;colormap('jet')]);
tiledlayout(3,length(conditions)+2);
avgConcat = mean(cat(3,concatMatrix{1:length(conditions)-1}),3,'omitnan');
plotHeatMatrix([concatMatrix,avgConcat],[conditions,"Average"],chNames,[0 .75]);

nexttile;hold on; axis ij; axis tight;title("Average Marginals");
imagesc(mean(avgConcat.*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim([0 .75]);colorbar;

noiseCorrMatrix = cellfun(@(m) mean(m,3,'omitnan'), {goAll,reachAll,graspAll}, 'UniformOutput',false);
plotHeatMatrix(noiseCorrMatrix,corrPhases,chNames,[0 .25]);

for m = 1:length(noiseCorrMatrix)
    nexttile; hold on; axis tight; axis ij;title(corrPhases(m)+" marginals");
    imagesc(mean(noiseCorrMatrix{m}.*(~logical(diag(true(length(chNames),1)))./~logical(diag(true(length(chNames),1)))),2,'omitnan'));
    xlim([0.5 1]);xticklabels([]);yticks(1:2:length(chNames));yticklabels(chNames(1:2:end));clim([0 .25]);
end
colorbar;

nexttile([1,1]);hold on; axis image; axis ij;title("NumPairs");
imagesc(numPairs);
xticks(resMap(1:2:end));xticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));yticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxPairs)];cb.TickLabels = num2cell([1,round(maxPairs)]);clim([0 round(maxPairs)]);

nexttile([1,1]);hold on; axis ij; axis tight;title("FR");
imagesc(mean(allCondFR,2,'omitnan'));
xticks([]);xlim([.5 1]);yticks(resMap(1:2:end));yticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));yticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxFR)];cb.TickLabels = num2cell([1,round(maxFR)]);clim([0 round(maxFR)]);

nexttile([1,4]);hold on; axis ij; axis tight;title("All FR");
imagesc(allCondFR);
yticks(resMap(1:2:end));yticklabels(arrayfun(@num2str,resMap(1:2:end),'UniformOutput',false));yticklabels([]);
cb=colorbar(gca,'southoutside');cb.Ticks=[0,round(maxFR)];cb.TickLabels = num2cell([1,round(maxFR)]);clim([0 round(maxFR)]);

saveFigures(gcf,mainDir+"ngc14\Working\Correlations\","All_ABS",[]);

function spCounts= getSpikeCounts(spks,segInd,times,goodTrials,wSz)
spCounts = cellfun(@(c,i,t,g) cellfun(@(a,gi) cellfun(@(s,at,gt) (gt./gt).*sum(s<at(i)+wSz & s>at(i)),...
    a,t,num2cell(gi)),c,g,'UniformOutput',false),spks,segInd,times,goodTrials,'UniformOutput',false);
spCounts = reshape(permute(cell2mat(permute(cat(3,spCounts{:}),[2 1 3])),[1 3 2]),size(goodTrials{1},2),[]);
end