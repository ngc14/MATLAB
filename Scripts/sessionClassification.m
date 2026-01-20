num_dims = 4;
cnds = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(cnds),.001,.001,-1,3,containers.Map(cnds,{"StartReach","StartReach","StartReach"}));
phases = ["StartReach","StartHold"];
phaseWindows = {[-100 100], [-200 0]};
tPhys = unitTable(cnds,params);
MIN_NUM_TRIALS = 20;
plotSessionDiscriminants = false;
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
%%
for p = 1:length(phases)
    phaseConds = cellfun(@(t) find(strcmp(phases(p),t)), params.condSegMap.values(params.condSegMap.keys),'UniformOutput',false);
    trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end))),...,
        num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
        'UniformOutput',false)',ct,cs,'UniformOutput',false),num2cell(tPhys{:,contains(tPhys.Properties.VariableNames,"PSTH_")},1),...
        num2cell(tPhys{:,contains(tPhys.Properties.VariableNames,"Segs_")},1),phaseConds,repmat({phaseWindows{p}},1,length(phaseConds)),'UniformOutput',false);
    trialFRMat{p} = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
end
clear trialFR;
%%
close all;
[newD,projMat,lat] = deal({}); condC = NaN(max(tPhys.SiteNum),2,100);
goodUnits = all(cell2mat(cellfun(@(p) all(cellfun(@(s) size(s,2),p)>=MIN_NUM_TRIALS,2),trialFRMat,'UniformOutput',false)),2);
conditions = cell2mat(arrayfun(@(r) string(params.condAbbrev.values)+ "-" + r, phases,'UniformOutput',false));
condColors = containers.Map(params.condAbbrev.values,{[1 0 0],[1 .8 0 ],[0 0 1]});
condColorKeys = condColors.keys;
mks = {'square','pentagram'};
mksz = [70,120];
for i = 1:length(unique([tPhys.SiteNum]))
    if(any(goodUnits(tPhys.SiteNum==i)))
        sTrials = min(cellfun(@(n) min(n(tPhys.SiteNum==i & goodUnits,:),[],'all'),cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false)));
        currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) ...
            downsampleTrials(r,sTrials),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),...
            cellfun(@(t) m(tPhys.SiteNum==i & contains(params.condAbbrev.values,t) & goodUnits),...
            params.condAbbrev.values,'UniformOutput',false),'UniformOutput',false),trialFRMat,...
            cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
        currD = [cellfun(@(t) cell2mat(cellfun(@(s) sum(s,2),t,'UniformOutput',false)')',[currD{:}]', 'UniformOutput',false)];
        [pri,sci,eigi] = pca(vertcat(currD{:}));
        newD{i} = permute(reshape(sci(:,1:min(num_dims,size(sci,2))),sTrials,length(currD),min(num_dims,size(sci,2))),[3,2,1]);
        projMat{i} = pri;
        lat{i} = cumsum(eigi) ./ sum(eigi);
    else
        newD{i} = {};
        projMat{i} = {};
        lat{i} = {};
    end
    projD = newD{i}; [gm,dist]= deal({});
    if(size(projD,1)>(1+plotSessionDiscriminants))
        for icond = 1:length(conditions)
            for idim = 1:min(num_dims,size(projD,1))
                dist{icond,idim} = squeeze(projD(idim,icond,:));
            end
            gm{icond} = cell2mat(dist(icond,:));
        end
        parfor n=1:100
            cvIdx=cvpartition(size(projD,3),'Holdout',0.2);
            trIdx = training(cvIdx);
            tsIdx = test(cvIdx);
            trData = cellfun(@(m) m(trIdx,:),gm,'UniformOutput',false);
            tsData = cellfun(@(m) m(tsIdx,:),gm,'UniformOutput',false);
            if(plotSessionDiscriminants)
                figure(); hold on;
                h1 = cellfun(@(s,n) scatter3(s(:,1),s(:,2),s(:,3),mksz(strcmp(extractAfter(n,"-"),phases)),'MarkerFaceColor',cell2mat(condColors.values(cellstr(...
                    extractBefore(n,"-")))),'MarkerEdgeColor','none','Marker',mks(strcmp(extractAfter(n,"-"),phases))),trData,conditions);
            end
            for s = 1:2
                if(s==1)
                    nPairs = nchoosek(1:length(condColorKeys),2);
                    condInds = cellfun(@(d) contains(conditions,d), condColorKeys, 'Uniformoutput',false);
                else
                    nPairs = nchoosek(1:length(phases),2);
                    condInds = cellfun(@(d) contains(conditions,d), phases, 'Uniformoutput',false);
                end
                gData = {};gtData = {};
                for p = 1:length(condInds)
                    gData{p} = vertcat(trData{condInds{p}});
                    gtData{p} = vertcat(tsData{condInds{p}});
                end
                MdLinear=fitcdiscr(vertcat(gData{:}),cell2mat(cellfun(@(c) repmat(string(c),1,size(gData{1},1)),...
                    cellfun(@(o) strcat(conditions{o}),condInds,'UniformOutput',false),'UniformOutput',false)));
                if(plotSessionDiscriminants)
                    for p = 1:size(nPairs,1)
                        K=MdLinear.Coeffs(nPairs(p,1),nPairs(p,2)).Const;
                        L=MdLinear.Coeffs(nPairs(p,1),nPairs(p,2)).Linear;
                        f=@(x1,x2,x3) K+L(1)*x1+L(2)*x2+L(3)*x3;
                        fs=fimplicit3(f);
                        fs.FaceColor=[.7,0.8,0.8].*[s==1,s==2,s==2];
                        fs.EdgeColor='none';
                        fs.FaceAlpha=0.15;
                    end
                    saveFigures(gcf,saveDir+"\Site_Classifiers\Sessions\","Site_"+num2str(i),[]);
                end
                CM = confusionmat(cellstr(cell2mat(cellfun(@(c) repmat(string(c),1,size(gtData{1},1)),cellfun(@(o) ...
                    strcat(conditions{o}),condInds,'UniformOutput',false),'UniformOutput',false))),MdLinear.predict(cell2mat(gtData')));
                condC(i,s,n) = sum(CM.*(diag(ones(size(CM,1),1),0)==1),'all')./sum(CM,'all');
            end
        end
    end
end
err = std(condC,0,3,'omitnan');
condC = mean(condC,3,'omitnan');
nonEmptySessions = ~cellfun(@isempty,projMat);
weights = projMat(nonEmptySessions);
[~,ui,~] = unique(tPhys.SiteNum);
siteSomatotopy = tPhys.Somatotopy(ui);
siteSomatotopy = siteSomatotopy(nonEmptySessions);
monkeys = tPhys.Monkey(ui);
monkeys = monkeys(nonEmptySessions);
lat = lat(nonEmptySessions);
reps = unique(siteSomatotopy);
somatotopyColors = containers.Map(string(reps),{[.75 .3 .75],[1 .85 0],[0 0 1]});
reps = reps(reps~="Trunk");
figure(); tiledlayout();
for d =1:num_dims
    nexttile;hold on;title(['Factor ',num2str(d)]);
    for k = 1:length(weights)
        if(siteSomatotopy(k)~="Trunk")
            sessionW = cell2mat(weights(k)); w=[];
            for n = 1:length(sessionW)
                w(n,1) = abs(sessionW(n,min(size(sessionW,2),d)));
            end
            plot(1:length(sessionW),sort(w,'descend'),'color',[cell2mat(somatotopyColors.values({string(siteSomatotopy(k))})) 0.65],'LineWidth',0.25);
        end
        ylim([0 1]);
        xlim([1 max(cellfun(@length,weights))]);
        if(d==1)
            ylabel("Abs value of projection weights");
            xlabel('Neuron number');
        end
        if(d==1)
            l1 = arrayfun(@(p) plot(NaN,NaN,'Color',cell2mat(somatotopyColors.values({char(p)}))),reps);
            legend(l1,string(reps),'AutoUpdate','off');
        end
    end
    allLines = get(gca,'Children');
    groupInds = cellfun(@(u) strcmp(cellfun(@num2str,{allLines.Color},'UniformOutput',false),u),unique(cellfun(@num2str,{allLines.Color},'UniformOutput',false)),'UniformOutput',false);
    cellfun(@(a) plot(mean(cell2mat(cellfun(@(r) resize(r,[1,max(cellfun(@length,{allLines.YData}))],'FillValue',NaN),...
        {allLines(a).YData}','UniformOutput',false)),1,'omitnan'),'Color',max([0 0 0],allLines(find(a,1)).Color-[.2 .4 .2]),'LineWidth',2.5,'LineStyle','-.'),groupInds,'UniformOutput',false);
end
nexttile();hold on;title("Variance Explained")
varEx = cell2mat(cellfun(@(r) resize(r',max(cellfun(@length,lat)),FillValue=NaN),lat(siteSomatotopy~="Trunk"),'UniformOutput',false)');
cellfun(@(v,c) plot(v,'LineWidth',1,'Color',[cell2mat(somatotopyColors.values({c})),.5]),num2cell(varEx,2),cellstr(string(siteSomatotopy(siteSomatotopy~="Trunk"))));
allLines = get(gca,'Children');
groupInds = cellfun(@(u) strcmp(cellfun(@num2str,{allLines.Color},'UniformOutput',false),u),unique(cellfun(@num2str,{allLines.Color},'UniformOutput',false)),'UniformOutput',false);
cellfun(@(a) plot(mean(cell2mat({allLines(a).YData}'),1,'omitnan'),'Color',...
    max([0 0 0],allLines(find(a,1)).Color-[.2 .4 .2]),'LineWidth',2.5,'LineStyle','-.'),groupInds,'UniformOutput',false);
l1 = cellfun(@(p) plot(NaN,NaN,'Color',p),somatotopyColors.values);
legend(l1,string(reps),'AutoUpdate','off','Location','southeast');
plot(1:size(varEx,2),repmat(0.9,1,size(varEx,2)),'k:','LineWidth',2);
ylim([0 1]);
xlim([.5 max(cellfun(@length,lat))]);
ylabel("%");
xlabel("Latent Factor");
nexttile(); hold on;title("Discriminant Analysis");
colororder([[0 .5 0];[0 0 .75]]);
e=errorbar(condC(all(~isnan(condC),2),:),err(all(~isnan(condC),2),:)./2,'LineWidth',1,'LineStyle',':');
arrayfun(@(ee) set(ee,'CapSize',0),e);
ylabel("Accuracy");
xlabel("Session");
ylim([0 1]);
yyaxis('right');
ylabel("# of units");
plot(cellfun(@length,projMat(~cellfun(@isempty,projMat))),'k-' );
ax = gca();
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'k';
ylim([0 35]);
xlim([1 max(sum(~isnan(condC),1))]);
legend(["Conditions","Phases","# of units"],'Location','southwest');
nexttile(); hold on;
for s = 1:length(reps)
    plotAcc{s} = condC(siteSomatotopy==reps(s),:);
end
[~,sig]= arrayfun(@(c) cellfun(@(a) ttest2(a(:,1),a(:,2)),{cell2mat(cellfun(@(r) resize(r,[length(siteSomatotopy),1],'FillValue',NaN),...
    arrayfun(@(r) condC(siteSomatotopy==r,c),reps,'UniformOutput',false),'UniformOutput',false)')}),1:size(condC,2));
arrayfun(@(a) text(((1+length(sig))*(find(sig==a)-1))+1.5,.15,"p= "+num2str(a,3),'HorizontalAlignment','center'),sig);
boxplotGroup(plotAcc,'groupLabelType','both','primaryLabels',string(reps),'secondaryLabels',["Conditions","Phases"],'Notch','on');
ylim([0 1]);
ylabel("Accuracy");
title("Accuracy by Somatotopy");
nexttile(); hold on; plotAcc = {};
for s = 1:size(condC,2)
    currGroup = arrayfun(@(a) condC(siteSomatotopy==a,s),reps,'UniformOutput',false);
    plotAcc{s} = cell2mat(cellfun(@(r) resize(r,max(cellfun(@length,currGroup)),'FillValue',NaN),currGroup, 'UniformOutput',false)');
end
boxplotGroup(plotAcc,'groupLabelType','both','primaryLabels',["Conditions","Phases"],'SecondaryLabels',string(reps),'Notch','on');
ylim([0 1]);
ylabel("Accuracy")
title("Accuracy by Classifier");
saveFigures(gcf,saveDir+"\Site_Classifiers\","Summary",[]);
%mus = cellfun(@(m) mean(m,1,'omitnan'),trData,'UniformOutput',false);
%cv =  vertcat(trData{:})-cell2mat(cellfun(@(m,g) repmat(m,size(g,1),1),mus,trData,'UniformOutput',false)');
%LDA=makecdiscr(cell2mat(mus'),cv'*(cv/(size(cv,1)-length(gm))),'ClassNames',conditions);
%bins = histcounts(squeeze(projD(idim,icond,:)),[xbins,xbins(end)+1],'Normalization','percentage');
%pdf(fitdist(squeeze(projD(idim,icond,:)),'Normal'),xbins)
%all_h = findall(groot,'Type','Figure');
%D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));
function arr = downsampleTrials(r,sTrials)
if(sTrials>0)
    sz = size(r,2)-mod(size(r,2),2);
    trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
    nPairs = round(sz/sTrials)*(sTrials-trials);
    arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
        mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
else
    arr = NaN(size(r,1),1);
end
end