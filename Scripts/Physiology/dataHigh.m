model = "GilliganSkipper_ArmHand";
type = 'State';
saveDir = "S:\Lab\ngc14\Working\DataHi\";
saveFig = true;
num_dims=5;
sTrials = 20;
plotTrials = 1;
timeBins = [-.5, 1];
conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
savePath = saveDir+type+"\Full_Population\";
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
dimCond = reshape(extractAfter(model,"_")+"_"+params.condAbbrev.values',1,[]);
colors =containers.Map(dimCond,{[.7 0 0],[1 .65 0 ],[0 0 .75]}');%regexp(model,'[A-Z]+[^A-Z]+','match')
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
tPhys = unitTable(conditions,params);
%%
tableInds = contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]);
somaTable = tPhys{tableInds,"Somatotopy"};
allLocations = tPhys{tableInds,["XT","YT"]};
allSegs= arrayfun(@(s) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))},dimCond,'UniformOutput',false);%
tablePSTHD= arrayfun(@(a) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
numUnits = all(cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), tablePSTHD,'UniformOutput',false))>=sTrials,2);%
tablePSTHD = cellfun(@(c) c(numUnits), tablePSTHD, 'UniformOutput',false);
unitInds = randperm(sum(numUnits));%repmat({}',1,size(currD,2)/(length()));
somaTable = somaTable(numUnits);
allLocations = allLocations(numUnits,:);
allSegs = cellfun(@(s) s(numUnits), allSegs, 'UniformOutput', false);
somaReps = unique(somaTable);
clear tPhys tableInds numUnits;
%% smooth data and remove non-modulated units
binWidth = 10;smoothWin = 150;
%avgTrace = cellfun(@(c) mean(cell2mat(cellfun(@(u)mean(u,2,'omitnan'),c,'UniformOutput',false)'),2,'omitnan'),tablePSTHD,'UniformOutput',false);%zeros(size(params.bins))',taskPSTHD,'UniformOutput',false);%
if(strcmpi(type,"state"))
    phases = ["StartReach","StartHold"]; winSz = 200;
    phaseWin = repmat({{[-winSz*(3/4),winSz*(1/4)],[-winSz*(5/4), -winSz*(1/4)]}},1,length(conditions));
    phaseWin{end}{2} = [-winSz/2 0];
    ms_bins = binWidth;
    for p = 1:length(phases)
        phaseConds = cellfun(@(t) find(strcmp(phases{p},t)), params.condSegMap.values(params.condSegMap.keys),'UniformOutput',false);
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) sum(m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end))))/range(tw),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct(unitInds),cs(unitInds),'UniformOutput',false),tablePSTHD,allSegs,phaseConds,cellfun(@(n) n{p},phaseWin,'UniformOutput',false),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
    end
    currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) ...
        downsampleTrials(r,sTrials),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),arrayfun(@(t) ...
        m(all(n>=sTrials,2) & ismember(string(somaTable),string(somaReps)),contains(conditions,t)),conditions,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    smoothedData =  [currD{:}];
    taskPSTHD = vertcat(smoothedData{:})';
    mv = any(cell2mat(cellfun(@(a) ~all(cell2mat(reshape([a{:}],1,size(a{1},1),[]))==0,[2,3]), currD, 'UniformOutput',false)),2);
    condPhase = cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
        cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)'));
else
    ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
    if(~plotTrials)
        taskPSTHD = cellfun(@(n,a) {cell2mat(cellfun(@(m) mean(max(0,m-a),2,'omitnan')',n,'UniformOutput',false))},tablePSTHD,avgTrace,'UniformOutput',false);
    else
        taskPSTHD= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
    end
    taskPSTHD =  cellfun(@(a)cellfun(@(u)max(0,u(unitInds,:)),a,'UniformOutput',false), taskPSTHD,'UniformOutput',false);
    if(plotTrials)
        mv = sum(cell2mat(cellfun(@(m)mean(cell2mat(m),2,'omitnan').*1000>1,num2cell([taskPSTHD{:}],2),'UniformOutput',false)'),2)>sTrials/2;
    else
        mv = mean(cell2mat(cellfun(@(n)cell2mat(n),taskPSTHD,'UniformOutput',false)),2,'omitnan').*1000>1;
    end
    trialLength = floor(size(taskPSTHD{1}{1}, 2) / binWidth);
    for n = 1:length(taskPSTHD)
        smoothedData{n} = repmat({NaN(sum(mv),trialLength)},max(1,sTrials*plotTrials),1);
        for s = 1:length(smoothedData{n})
            for t = 1:trialLength
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth *t;
                smoothedData{n}{s}(:,t) = sum(taskPSTHD{n}{s}(mv,iStart:iEnd),2);%normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
            end
        end
    end
    smoothedData = cellfun(@(c) cellfun(@(s) (conv2(resize(s,[size(s,1),size(s,2)-1+floor(smoothWin/binWidth)],'Pattern','edge','side','both'),...
        transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),smoothedData,'UniformOutput',false);
    taskPSTHD = cellfun(@(s) cell2mat(cellfun(@(t) t(:,unique(round(ms_bins./binWidth))),s,'UniformOutput',false)'),smoothedData, 'UniformOutput',false);
end
somaLabs = somaTable(mv);
location = allLocations(mv,:);
segInds = cellfun(@(s) fix(s(:,~all(isnan(s),1))),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),...
    cellfun(@(aa)cell2mat(cellfun(@(a) findBins(mean(a,1,'omitnan'),params.bins),aa(unitInds),'UniformOutput',false)),...
    allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
%% plot single dimension PCA
cls = [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0;];
cpSteps = fix(270/size(cls,1)); cp=NaN(270-cpSteps,3);
for l = 1:size(cls,1)-1
    cp([1:cpSteps]+(cpSteps*(l-1)),:) = cell2mat(arrayfun(@(a) linspace(cls(l,a),cls(l+1,a),cpSteps),1:3,'UniformOutput',false)')';
end
[loadings, scores, eig,~,exp] = pca(cell2mat(taskPSTHD)','Economy',false,'Centered','off','NumComponents',num_dims,'Algorithm','eig');%%cell2mat(smoothedData')
%somaProj = cellfun(@(s)arrayfun(@(c) scores((1+(c-1)*size(scores,1)/length(dimCond)):c*size(scores,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
somaProj = cellfun(@(s)cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end)...
    .*loadings(somaLabs==s,n),1,'omitnan'),1:num_dims,'UniformOutput',false)').*binWidth,c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
for c = 1:length(dimCond)
    for s =1:length(somaReps)
        weightedPSTHS = cell2mat(somaProj{s}{c});
        avgDim{s,c} = squeeze(mean(permute(reshape(weightedPSTHS,size(weightedPSTHS,1)/...
            (max(1,plotTrials*sTrials)),[],size(weightedPSTHS,2)),[3 2 1]),2,'omitnan'));
    end
end
statePlot = contains(type,'state','IgnoreCase',1);
if(statePlot)
    pf = length(phases);
    allSomaProj = somaProj;
else
    pf = 1;
end
%%
for pI = 1:pf
    if(statePlot)
        somaProj = cellfun(@(a) a((length(conditions)*(pI-1))+[1:length(conditions)]), allSomaProj, 'UniformOutput',false);
    end
    figure(); tl1=tiledlayout(6,4);
    nexttile(tl1,[1,1]); hold on; title("Variance Explained");plot(cumsum(exp),'LineWidth',2);ylim([0 100]);xlim([0 sum(mv)]);
    plot([0 sum(mv)],[90 90],'r--','LineWidth',1);xlim([0,200]);
    bx=nexttile(tl1,[1 1]);

    boxplotGroup(bx,arrayfun(@(s) loadings(somaLabs==s,:),somaReps,'UniformOutput',false)','PrimaryLabels',...
        repmat({''},num_dims*length(somaReps),1),'Symbol','','SecondaryLabels',arrayfun(@num2str,1:num_dims,'Uniformoutput',false),'Notch','on');
    hold on; title("Loadings Distributions");ylim([-.1 .1]); xlabel("Factor");
    ax = nexttile(tl1,[1,2]); hold on;
    if(contains(type,"Traj"))
        title("PSTHS");
        psthGrouped = cellfun(@(s)cellfun(@(c)cellfun(@(t) mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end),1,'omitnan').*binWidth,...
            c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'Uniformoutput',false);
        cellfun(@(s,l) cellfun(@(m,p) plot(cell2mat(m)','LineStyle',l,'LineWidth',1.5,'Color',p),s,colors.values),psthGrouped,repmat({'-'},size(psthGrouped,1),size(psthGrouped,2)));
        condSegs = cellfun(@(i) resize(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):max(params.bins)),[1,length(epochSegs)],'FillValue',NaN),segInds,'UniformOutput',false);
        cellfun(@(x,s) plot(ax,[x;x], repmat(get(ax,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
            cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs,'UniformOutput',false)),2)',colors.values(cellstr(dimCond)));
        condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
        arrayfun(@(x) plot(ax,[x,x], get(ax,'YLim'),'Color','k','LineStyle',':','LineWidth',1.5),condSegs(1:end-2));
        set(ax,'XTickLabels',linspace(timeBins(1),max(params.bins),length(get(ax,'XTick'))));
        set(ax,'XLim',[0 trialLength-min(round(ms_bins./binWidth))]);
        set(ax,'XTick',0:50:trialLength-min(round(ms_bins./binWidth)));
    else
        title("Spikes")
        xVals = cell2mat(arrayfun(@(a) repmat(a,[1,ceil(length(loadings)/10)]),1:10,'UniformOutput',false));
        cellfun(@(p,n,cl) scatter(xVals(1:size(p{1},1))+(n*10),mean(cell2mat(p').*binWidth,2,'omitnan'),10,cl,'filled','o'),...
            smoothedData(length(conditions)*(pI-1)+[1:length(conditions)]),num2cell(1:length(conditions)),colors.values);
    end

    [locRC,sortR] = sort(location(:,1).*ImagingParameters.px2mm);
    [locML,sortM] = sort(location(:,2).*ImagingParameters.px2mm);
    minLoadings = mean(loadings)-1*std(loadings,0,1);
    maxLoadings = mean(loadings)+1*std(loadings,0,1);

    locRC = min(locRC,4);
    ct=nexttile(tl1,[num_dims 1]); hold on; title("Loadings R/C"); imagesc(cell2mat(arrayfun(@(s) loadings(sortR(somaLabs(sortR)==s),:),somaReps,'UniformOutput',false)));
    axis ij;axis tight;colormap(ct,'jet');clim([-.1 .1]);
    binGroups =arrayfun(@(s) cell2mat(arrayfun(@(f) find(fix(locRC(somaLabs==s))==f,1),1:length(unique(fix(locRC)))-1,'UniformOutput',false)')',somaReps,'UniformOutput',false);
    arrayfun(@(i,ls) line([0+.5 num_dims+.5],repmat(i,1,2),'Color','k','LineWidth',2,'LineStyle',ls),...
        [sum(somaLabs=="Arm"),binGroups{1},sum(somaLabs=="Arm")+binGroups{end}],["--",repmat(":",1,length(cell2mat(binGroups')))]);colorbar;
    for d = 1:num_dims
        bx=nexttile(tl1,tilenum(tl1,1+d,2),[1 1]); hold on; title("Binned R/C Loadings "+num2str(d));
        boxplotGroup(bx,arrayfun(@(s) cell2mat(arrayfun(@(n) resize(loadings(fix(locRC(somaLabs==s))==n,d),[sum(somaLabs==s),1],'FillValue',NaN),...
            unique(fix(locRC)),'UniformOutput',false)'),somaReps,'UniformOutput',false)',...
            'Symbol','','SecondaryLabels',arrayfun(@num2str,unique(fix(locRC)),'Uniformoutput',false),'Notch','on','PrimaryLabels',repmat({''},length(somaReps)*length(unique(fix(locRC))),1));
        ylim([-.05 .05]+(d==1*[.05 .05]));
    end
    locML = min(max(1,locML),7);
    ct=nexttile(tl1,[num_dims 1]); hold on; title("Loadings M/L"); imagesc(cell2mat(arrayfun(@(s) loadings(sortM(somaLabs(sortM)==s),:),somaReps,'UniformOutput',false)));
    axis ij;axis tight;colormap(ct,'jet');clim([-.1 .1]);
    binGroups =arrayfun(@(s) cell2mat(arrayfun(@(f) find(fix(locML(somaLabs==s))==f,1),1:length(unique(fix(locML)))-1,'UniformOutput',false)')',somaReps,'UniformOutput',false);
    arrayfun(@(i,ls) line([0+.5 num_dims+.5],repmat(i,1,2),'Color','k','LineWidth',2,'LineStyle',ls),...
        [sum(somaLabs=="Arm"),binGroups{1},sum(somaLabs=="Arm")+binGroups{end}],["--",repmat(":",1,length(cell2mat(binGroups')))]);colorbar;
    for d = 1:num_dims
        bx=nexttile(tl1,tilenum(tl1,1+d,4),[1 1]); hold on; title("Binned M/L Loadings "+num2str(d));
        boxplotGroup(bx,arrayfun(@(s) cell2mat(arrayfun(@(n) resize(loadings(fix(locML(somaLabs==s))==n,d),[sum(somaLabs==s),1],'FillValue',NaN),...
            unique(fix(locML)),'UniformOutput',false)'),somaReps,'UniformOutput',false)',...
            'Symbol','','SecondaryLabels',arrayfun(@num2str,unique(fix(locML)),'Uniformoutput',false),'Notch','on','PrimaryLabels',repmat({''},length(somaReps)*length(unique(fix(locML))),1));
        ylim([-.05 .05]+(d==1*[.05 .05]));
    end
    if(saveFig)
        saveFigures(gcf,savePath,phases(pI)+"_Variance+Loadings",[]);
    end
    tl1 = tiledlayout(4,2);
    [~,rankedWeight] = arrayfun(@(s) sort(loadings(somaLabs==s,:),1,'ascend'), somaReps, 'UniformOutput',false);
    sortR = cell2mat(cellfun(@(r,s) cell2mat(cellfun(@(c,t) t(c),num2cell(r,1),repmat({sortR(somaLabs==s)},[1,size(r,2)]),'UniformOutput',false)), rankedWeight,num2cell(somaReps),'UniformOutput',false));
    sortM = cell2mat(cellfun(@(r,s) cell2mat(cellfun(@(c,t) t(c),num2cell(r,1),repmat({sortM(somaLabs==s)},[1,size(r,2)]),'UniformOutput',false)), rankedWeight,num2cell(somaReps),'UniformOutput',false));
    rankedWeight = cell2mat(cellfun(@(r) (r-1)./(size(r,1)-1),rankedWeight,'UniformOutput',false));
    rcLoadings = rankedWeight(sortR);
    %rcLoadings = 2*((loadingsSplit-minLoadings)./(maxLoadings-minLoadings))-1;
    ct=nexttile(tl1,[4 1]); hold on; title("Ranked Loadings R/C"); imagesc(rcLoadings);axis ij;axis tight;colormap(ct,'jet');clim([0 1]);
    line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
    mlLoadings = rankedWeight(sortM);
    %mlLoadings = 2*((loadingsSplit-minLoadings)./(maxLoadings-minLoadings))-1;
    ct=nexttile(tl1,[4 1]); hold on; title("Ranked Loadings M/L"); imagesc(mlLoadings);axis ij;axis tight;colormap(ct,'jet');clim([0 1]);
    line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
    if(saveFig)
        saveFigures(gcf,savePath,phases(pI)+"_RankedLoadings",[]);
    end
    %%
    for n = 1:3
        figure();
        lc = {'-','-'};ax = {};plotSegs={};
        if(n==1)
            tileorder = 'columnmajor';
        else
            tileorder = 'rowmajor';
        end
        tl=tiledlayout(max(2,n),num_dims-(2*(n==1)),'TileIndexing',tileorder);
        for c = 1:length(dimCond)
            for s =1:length(somaReps)
                co = repmat(colors(dimCond(c)),num_dims,1);
                weightedPSTHS = (somaProj{s,c})';%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
                weightedPSTHS = permute(reshape(weightedPSTHS,size(weightedPSTHS,1)/(max(1,plotTrials*sTrials)),[],size(weightedPSTHS,2)),[3 2 1]);
                if(n==1)
                    ax{end+1}=nexttile();hold on;title(conditions(c)+"- " + string(somaReps(s)));
                    plotSegs{end+1} = c;
                    cl = rgb2hsv(colors(dimCond(c)));
                    co = hsv2rgb([linspace(cl(1),cl(1),num_dims);linspace(1,.5,num_dims);linspace(.3,1,num_dims)]');
                    %avgDim{s,c} = squeeze(mean(weightedPSTHS,2,'omitnan'));
                elseif(n==2)
                elseif(n==3)
                    lc = {'--',':'};
                    if(s==length(somaReps))
                        co = rgb2hsv(colors(dimCond(c)));
                        co = repmat(hsv2rgb(co(1), 1, .5),num_dims,1);
                    end
                    weightedPSTHS = permute(cell2mat(cellfun(@(p) cell2mat(reshape(p{c}',1,1,[])), somaProj(s), 'UniformOutput',false)),[2 3 1]);
                end
                for d = 1:size(weightedPSTHS,3)
                    if(n>=2)
                        if(n==2)
                            ax{end+1}=nexttile(((s-1)*tl.GridSize(end))+d);
                            plotSegs{end+1} = 1:length(dimCond);
                            titleName = "Dim " + num2str(d) + " - " + string(somaReps(s));
                        else
                            ax{end+1}=nexttile((c-1)*tl.GridSize(end)+d);
                            plotSegs{end+1} = c;
                            titleName = "Dim "+ num2str(d);
                        end
                        hold on;title(titleName);
                    end
                    if(~statePlot)
                        p=plot(squeeze(weightedPSTHS(:,:,d)),'LineWidth',.5,'LineStyle',lc{s});
                        arrayfun(@(pc) set(pc,'Color',[co(d,:),.35*s]),p);
                    else
                        [bins centers] = hist(squeeze(weightedPSTHS(:,:,d)),linspace(...
                            min(weightedPSTHS(:,:,d),[],'all'),...%mean(weightedPSTHS(:,:,d),'all')-std(weightedPSTHS(:,:,d),0,'all')
                            max(weightedPSTHS(:,:,d),[],'all'),10));
                        bins = bins ./ sum(bins);
                        bar(centers, bins, 'FaceColor',co(d,:));
                    end
                end
            end
        end
        if(~statePlot)
            condSegs = cellfun(@(i) resize(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):max(params.bins)),[1,length(epochSegs)],'FillValue',NaN),segInds,'UniformOutput',false);
            cellfun(@(a,p) cellfun(@(x,s) plot(a,[x;x], repmat(get(a,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
                cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(dimCond(p)))),ax,plotSegs);
            condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
            cellfun(@(a) arrayfun(@(x) plot(a,[x,x], get(a,'YLim'),'k--','LineWidth', 1.5),condSegs(1:end-2)),ax);
            cellfun(@(a) set(a,'XLim',[0 trialLength-min(round(ms_bins./binWidth))]),ax);
            cellfun(@(a) set(a,'XTick',0:50:trialLength-min(round(ms_bins./binWidth))),ax);
            cellfun(@(a) set(a,'XTickLabels',linspace(timeBins(1),max(params.bins),length(get(a,'XTick')))),ax);
        end
        if(saveFig)
            rootSave = "WeightedUnits";
            if(n==1)
                fileNameSave = rootSave+"_DIM";
            elseif(n==2)
                fileNameSave = rootSave+"_COND";
            else
                fileNameSave=rootSave+"_SOMA";
            end
            saveFigures(gcf,savePath,fileNameSave+"_Sm"+num2str(smoothWin),[]);
            save(savePath+"AvgDimTraj_Sm"+num2str(smoothWin),'avgDim');
        end
    end
end
%% DATAHIGH Dim Reduce
cls = cellfun(@(r) repmat({r},max(plotTrials*sTrials,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),5);...
    linspace(1,.25,5);linspace(.85,1,5)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
dHiStruct = struct('data',cellfun(@(t) t,vertcat(taskPSTHD{:}),'UniformOutput',false),'epochStarts',...
    cellfun(@(r) repmat([1,r],max(1,sTrials*plotTrials),1), segInds,'UniformOutput',false),'condition',cellstr(cell2mat(...
    cellfun(@(d) repmat(string(d),max(plotTrials*sTrials,1),1),cellstr(dimCond),'UniformOutput',false)')),'epochColors',cellfun(@cell2mat,cls,'UniformOutput',false));
DataHigh(dHiStruct,'DimReduce');
if(saveFig)
    save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
end
%% Get DATAHIGH data
all_h = findall(groot,'Type','Figure');handles = guihandles(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));
D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));D= D.D;
% Plot DATAHIGH single dimension
conds = unique({D.condition});
figure(); tax=tiledlayout(1,max(1,num_dims/2)*2);
ylimT = [min(arrayfun(@(m) min(m.data,[],'all'),D)),max(arrayfun(@(m) max(m.data,[],'all'),D))];
for icond = 1:length(conds)
    for idim = 1:num_dims
        itrial = find(contains({D.condition}, extractAfter(conds{icond},"_")));
        epochs = [round(mean(vertcat(D(itrial).epochStarts),1,'omitnan')),size(D(1).data,2)];
        nexttile(idim); hold on; title(idim); ylim(ylimT);
        dTrial = cat(3,(D(icond).data));
        for iepoch = 1:length(epochs)-1
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
            if(icond==length(conds) && idim==1)
                l = arrayfun(@(d) plot(NaN(length(dimCond),1),'LineWidth',3,'Color',cell2mat(colors.values(cellstr(d)))),dimCond);
                legend(l,dimCond,'Autoupdate','off','FontSize',14,'Orientation','vertical');
            end
        end
    end
end
%%
figure(); tax=tiledlayout(1,length(conds));
for icond = 1:length(conds)
    legendColors= {};
    currColor = rgb2hsv(cell2mat(colors.values(cellstr(dimCond(icond)))));currColor(2) = 1;currColor(end) = .4;
    for idim = 1:num_dims
        if(idim<=num_dims/2)
            currColor(end) = min(1,currColor(end) + (.4*(idim-1)));
        else
            currColor(2) = max(.25,currColor(2) - (.25*(idim-(num_dims/2))));currColor(end) = 1;
        end
        itrial = find(ismember({D.condition}, conds{icond}));
        dTrial = cat(3,D(itrial).data);
        epochs = [mean(vertcat(D(itrial).epochStarts),1,'omitnan'),size(D(1).data,2)];
        nexttile(icond); hold on; title(string(conds(icond))); ylim(ylimT);
        for iepoch = 1:length(epochs)-1
            indices = epochs(iepoch):(epochs(iepoch+1));
            cellfun(@(p) plot(indices,p,'Color', hsv2rgb(currColor),'LineWidth',2),...
                num2cell(squeeze(mean(dTrial(idim,indices,:),max(~plotTrials*length(size(dTrial))+1,plotTrials),'omitnan')),1+~plotTrials));
            if(iepoch>1);line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');end
        end
        legendColors{idim} = hsv2rgb(currColor);
        if(idim==num_dims)
            l = cellfun(@(d) plot(NaN(length(num_dims),1),'Color',d),legendColors);
            legend(l,string(1:num_dims),'Autoupdate','off','Location','northeast');
        end
    end
end

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end