conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
tPhys = unitTable(conditions,params);
%%
model = "GilliganSkipper_ArmHand";
type = 'State';
saveDir = "S:\Lab\ngc14\Working\DataHi\";
saveFig = false;
num_dims=5;
sTrials = 20;
plotTrials = 1;
timeBins = [-.5, 1];
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
dimCond = reshape(extractAfter(model,"_")+"_"+params.condAbbrev.values',1,[]);
colors =containers.Map(dimCond,{[.7 0 0],[1 .65 0 ],[0 0 .75]}');%regexp(model,'[A-Z]+[^A-Z]+','match')
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
savePath = saveDir+type+"\Full_Population\";
tableInds = contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]);
somaTable = tPhys{tableInds,"Somatotopy"};
allLocations = tPhys{tableInds,["X","Y"]};
allSegs= arrayfun(@(s) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))},dimCond,'UniformOutput',false);%
tablePSTHD= arrayfun(@(a) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
numUnits = all(cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), tablePSTHD,'UniformOutput',false))>=sTrials,2);%
tablePSTHD = cellfun(@(c) c(numUnits), tablePSTHD, 'UniformOutput',false);
unitInds = randperm(sum(numUnits));%repmat({}',1,size(currD,2)/(length()));
somaTable = somaTable(numUnits);
allLocations = allLocations(numUnits,:);
allSegs = cellfun(@(s) s(numUnits), allSegs, 'UniformOutput', false);
somaReps = unique(somaTable);
clear tPhys;
%% smooth data and remove non-modulated units
binWidth = 10;smoothWin = 150;
ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
%avgTrace = cellfun(@(c) mean(cell2mat(cellfun(@(u)mean(u,2,'omitnan'),c,'UniformOutput',false)'),2,'omitnan'),tablePSTHD,'UniformOutput',false);%zeros(size(params.bins))',taskPSTHD,'UniformOutput',false);%
if(strcmpi(type,"state"))
    phases = ["StartReach","StartHold"]; winSz = 200;
    phaseWin = repmat({{[-winSz*(3/4),winSz*(1/4)],[-winSz*(5/4), -winSz*(1/4)]}},1,length(conditions));
    phaseWin{end}{2} = [-winSz/2 0];
    for p = 1:length(phases)
        phaseConds = cellfun(@(t) find(strcmp(phases{p},t)), params.condSegMap.values(params.condSegMap.keys),'UniformOutput',false);
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) sum(m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end)))),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct(unitInds),cs(unitInds),'UniformOutput',false),tablePSTHD,allSegs,phaseConds,cellfun(@(n) n{p},phaseWin,'UniformOutput',false),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
    end
    currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) ...
        downsampleTrials(r,sTrials),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),arrayfun(@(t) ...
        m(all(n>=sTrials,2) & ismember(string(somaTable),string(somaReps)),contains(conditions,t)),conditions,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    currD =  [currD{:}];
    taskPSTHD = vertcat(currD{:});
    condTask = cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
        cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)'));
else
    if(~plotTrials)
        taskPSTHD = cellfun(@(n,a) {cell2mat(cellfun(@(m) mean(max(0,m-a),2,'omitnan')',n(numUnits),'UniformOutput',false))},tablePSTHD,avgTrace,'UniformOutput',false);
    else
        taskPSTHD= cellfun(@(v,a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
            v(numUnits),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD),avgTrace, 'UniformOutput',false);
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
    clear smoothedData;
end
somaLabs = somaTable(mv);
location = allLocations(mv,:);
[loadings, scores, eig,~,exp] = pca(cell2mat(taskPSTHD)','Economy',false,'Centered','off','NumComponents',num_dims,'Algorithm','eig');%%cell2mat(smoothedData')
%% plot single dimension PCA
%somaProj = cellfun(@(s)arrayfun(@(c) scores((1+(c-1)*size(scores,1)/length(dimCond)):c*size(scores,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
somaProj = cellfun(@(s)cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end)...
    .*loadings(somaLabs==s,n),1,'omitnan'),1:num_dims,'UniformOutput',false)').*binWidth,c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
psthGrouped = cellfun(@(s)cellfun(@(c)cellfun(@(t) mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end),1,'omitnan').*binWidth,...
    c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'Uniformoutput',false);
segInds = cellfun(@(s) fix(s(:,~all(isnan(s),1))),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),...
    cellfun(@(aa)cell2mat(cellfun(@(a) findBins(mean(a,1,'omitnan'),params.bins),aa(unitInds),'UniformOutput',false)),...
    allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
figure(); tl1=tiledlayout(5,4);
nexttile(tl1,[1,2]); hold on; title("Variance Explained");plot(cumsum(exp),'LineWidth',2);ylim([0 100]);xlim([0 sum(mv)]);
plot([0 sum(mv)],[90 90],'r--','LineWidth',1);
bx=nexttile(tl1,[1 1]); 
boxplotGroup(bx,arrayfun(@(s) loadings(sortR(somaLabs(sortR)==s),:),somaReps,'UniformOutput',false)','PrimaryLabels',...
    repmat({''},num_dims*length(somaReps),1),'Symbol','','SecondaryLabels',arrayfun(@num2str,1:num_dims,'Uniformoutput',false),'Notch','on');
hold on; title("Loadings Distributions");ylim([-.1 .1]); xlabel("Factor");
ax = nexttile(tl1,[1,1]); hold on; title("PSTHS");
cellfun(@(s,l) cellfun(@(m,p) plot(cell2mat(m)','LineStyle',l,'LineWidth',1.5,'Color',p), s,colors.values), psthGrouped,repmat({'-'},size(psthGrouped,1),size(psthGrouped,2)));
set(ax,'XTickLabels',linspace(timeBins(1),max(params.bins),length(get(ax,'XTick'))));
condSegs = cellfun(@(i) resize(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):max(params.bins)),[1,length(epochSegs)],'FillValue',NaN),segInds,'UniformOutput',false);
cellfun(@(x,s) plot(ax,[x;x], repmat(get(ax,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
    cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs,'UniformOutput',false)),2)',colors.values(cellstr(dimCond)));
condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
arrayfun(@(x) plot(ax,[x,x], get(ax,'YLim'),'Color','k','LineStyle',':','LineWidth',1.5),condSegs(1:end-2));
set(ax,'XLim',[0 trialLength-min(round(ms_bins./binWidth))]);
set(ax,'XTick',0:50:trialLength-min(round(ms_bins./binWidth)));

[~,sortR] = sort(location(:,1));
[~,sortM] = sort(location(:,2));
minLoadings = mean(loadingsSplit)-1*std(loadingsSplit,0,1);
maxLoadings = mean(loadingsSplit)+1*std(loadingsSplit,0,1);

loadingsSplit = cell2mat(arrayfun(@(s) loadings(sortR(somaLabs(sortR)==s),:),somaReps,'UniformOutput',false));
nexttile(tl1,[4 1]); hold on; title("Loadings R/C"); imagesc(loadingsSplit);axis ij;axis tight;colormap('jet');clim([-.1 .1]);
line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','k','LineWidth',1.5,'LineStyle','--');colorbar;
loadingsSplit = cell2mat(arrayfun(@(s) loadings(sortM(somaLabs(sortM)==s),:),somaReps,'UniformOutput',false));
nexttile(tl1,[4 1]); hold on; title("Loadings M/L"); imagesc(loadingsSplit);axis ij;axis tight;colormap('jet');clim([-.1 .1]);
line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','k','LineWidth',1.5,'LineStyle','--');colorbar;

loadingsSplit = cell2mat(arrayfun(@(s) loadings(sortR(somaLabs(sortR)==s),:),somaReps,'UniformOutput',false));
rcLoadings = 2*((loadingsSplit-minLoadings)./(maxLoadings-minLoadings))-1;
nexttile(tl1,[4 1]); hold on; title("Normalized Loadings R/C"); imagesc(rcLoadings);axis ij;axis tight;colormap('jet');clim([-1 1]);
line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
loadingsSplit = cell2mat(arrayfun(@(s) loadings(sortM(somaLabs(sortM)==s),:),somaReps,'UniformOutput',false));
mlLoadings = 2*((loadingsSplit-minLoadings)./(maxLoadings-minLoadings))-1;
nexttile(tl1,[4 1]); hold on; title("Normalized Loadings M/L"); imagesc(mlLoadings);axis ij;axis tight;colormap('jet');clim([-1 1]);
line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
if(saveFig)
    saveFigures(gcf,savePath,"Variance+Loadings",[]);
end
avgDim = {};
for n = 1:3
    figure();
    lc = {'-','-'};ax = {};plotSegs={};
    tl=tiledlayout(n,num_dims+(n==1),'TileIndexing','rowmajor');
    for c = 1:length(dimCond)
        for s =1:length(somaReps)
            co = repmat(colors(dimCond(c)),num_dims,1);
            weightedPSTHS = cell2mat(somaProj{s}{c});%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
            weightedPSTHS = permute(reshape(weightedPSTHS,size(weightedPSTHS,1)/(max(1,plotTrials*sTrials)),[],size(weightedPSTHS,2)),[3 2 1]);
            if(n==1)
                ax{end+1}=nexttile();hold on;title(conditions(c)+"- " + string(somaReps(s)));
                plotSegs{end+1} = c;
                cl = rgb2hsv(colors(dimCond(c)));
                co = hsv2rgb([linspace(cl(1),cl(1),num_dims);linspace(1,.5,num_dims);linspace(.3,1,num_dims)]');
                avgDim{s,c} = squeeze(mean(weightedPSTHS,2,'omitnan'));
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
                p=plot(squeeze(weightedPSTHS(:,:,d)),'LineWidth',.5,'LineStyle',lc{s});
                arrayfun(@(pc) set(pc,'Color',[co(d,:),.35*s]),p);
            end
        end
    end
    condSegs = cellfun(@(i) resize(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):max(params.bins)),[1,length(epochSegs)],'FillValue',NaN),segInds,'UniformOutput',false);
    cellfun(@(a,p) cellfun(@(x,s) plot(a,[x;x], repmat(get(a,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
        cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(dimCond(p)))),ax,plotSegs);
    condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
    cellfun(@(a) arrayfun(@(x) plot(a,[x,x], get(a,'YLim'),'k--','LineWidth', 1.5),condSegs(1:end-2)),ax);
    cellfun(@(a) set(a,'XLim',[0 trialLength-min(round(ms_bins./binWidth))]),ax);
    cellfun(@(a) set(a,'XTick',0:50:trialLength-min(round(ms_bins./binWidth))),ax);
    cellfun(@(a) set(a,'XTickLabels',linspace(timeBins(1),max(params.bins),length(get(a,'XTick')))),ax);
    if(saveFig)
        rootSave = "WeightedPSTHS";
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