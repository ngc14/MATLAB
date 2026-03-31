conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
tPhys = unitTable(conditions,params);
%%
model = "GilliganSkipper_ArmHand";
type = 'Traj';
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
saveFig = true;
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
allSegs= arrayfun(@(s) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))},dimCond,'UniformOutput',false);%
tablePSTHD= arrayfun(@(a) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
numUnits = all(cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), tablePSTHD,'UniformOutput',false))>=sTrials,2);%
unitInds = randperm(sum(numUnits));%repmat({}',1,size(currD,2)/(length()));
segInds = cellfun(@(s) fix(mean(s(:,~all(isnan(s),1)),1,'omitnan')),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),...
    cellfun(@(aa)cell2mat(cellfun(@(a) findBins(mean(a,1,'omitnan'),params.bins),aa(unitInds),'UniformOutput',false)),...
    allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
somaTable = somaTable(unitInds);
somaReps = unique(somaTable);
clear tPhys allSegs;
%% smooth data and remove non-modulated units
binWidth = 10;smoothWin = 150;
ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
avgTrace = cellfun(@(c) mean(cell2mat(cellfun(@(u)mean(u,2,'omitnan'),c,'UniformOutput',false)'),2,'omitnan'),tablePSTHD,'UniformOutput',false);%zeros(size(params.bins))',taskPSTHD,'UniformOutput',false);%
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
somaLabs = somaTable(mv);
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
[loadings, scores, eig,~,exp] = pca(cell2mat(cellfun(@(s) cell2mat(cellfun(@(t) t(:,unique(round(ms_bins./binWidth))),s,'UniformOutput',false)'),...
    smoothedData, 'UniformOutput',false))','Economy',false,'Centered','on','NumComponents',num_dims,'Algorithm','eig');%%cell2mat(smoothedData')
clear taskPSTHD;
%% plot single dimension PCA
%somaProj = cellfun(@(s)arrayfun(@(c)s((1+(c-1)*size(s,1)/length(dimCond)):c*size(s,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),...
% cellfun(@(s,e) scores(repmat(somaLabs,length(dimCond),1)==s,:),num2cell(somaReps),'UniformOutput',false),'UniformOutput',false);
%somaProj = cellfun(@(s)arrayfun(@(c) scores((1+(c-1)*size(scores,1)/length(dimCond)):c*size(scores,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
somaProj = cellfun(@(s)cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end).*loadings(somaLabs==s,n),1,'omitnan'),...
    1:num_dims,'UniformOutput',false)').*1000,c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
psthGrouped = cellfun(@(s)cellfun(@(c)cellfun(@(t) mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end),1,'omitnan').*binWidth,...
    c,'UniformOutput',false),smoothedData,'UniformOutput',false),num2cell(somaReps),'Uniformoutput',false);
bc = num2cell([1 0 0;1 .7 0; 0 0 1],2);
lc = {'-',':'};
ax = {};
for n = 1:3
    figure();
    ax = {};plotSegs={};
    tl=tiledlayout(n+1,num_dims+1,'TileIndexing','rowmajor');
    nexttile([1,1]); hold on; title("Variance Explained");plot(cumsum(exp),'LineWidth',2);ylim([0 100]);xlim([0 sum(mv)]);
    plot([0 sum(mv)],[90 90],'r--','LineWidth',1);
    loadingsSplit = arrayfun(@(s) loadings(somaLabs==s,:),somaReps,'UniformOutput',false);
    armCount = sum(somaLabs=="Arm");
    nexttile([1 2]); hold on; title("Loadings"); imagesc(cell2mat(loadingsSplit)); axis ij; axis tight; colormap('turbo'); clim([-.1 .1]);
    line([0+.5 num_dims+.5],repmat(armCount,1,2),'Color','k','LineWidth',2,'LineStyle',':');
    colorbar;
    bx=nexttile([1 2]); hold on; title("Loadings Distributions");
    boxplotGroup(bx,loadingsSplit','PrimaryLabels',repmat({''},num_dims*length(somaReps),1),'Symbol','',...
        'SecondaryLabels',arrayfun(@(s) "Factor "+num2str(s),1:num_dims),'Notch','on'); ylim([-.1 .1]);
    ax{end+1} = nexttile(tl,[1,1]); hold on; title("PSTHS");
    plotSegs{end+1} = 1:length(dimCond);
    cellfun(@(s,l) cellfun(@(m,p) plot(cell2mat(m)','LineStyle',l,'LineWidth',1,'Color',p), s,colors.values), psthGrouped,lc');
    for c = 1:length(dimCond)
        for s =1:length(somaReps)
            co = repmat(colors(dimCond(c)),num_dims,1);
            weightedPSTHS = cell2mat(somaProj{s}{c})';%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
            weightedPSTHS = permute(reshape(weightedPSTHS',size(weightedPSTHS,2)/(max(1,plotTrials*sTrials)),[],size(weightedPSTHS,1)),[3 2 1]);
            if(n==1)
                ax{end+1}=nexttile();hold on;title("Cond " + num2str(c));
                plotSegs{end+1} = c;
                cl = rgb2hsv(colors(dimCond(c)));
                co = hsv2rgb([linspace(cl(1),cl(1),num_dims);linspace(1,.5,num_dims);linspace(.3,1,num_dims)]');
            elseif(n==2)
            elseif(n==3)
                weightedPSTHS = permute(cell2mat(cellfun(@(p) cell2mat(reshape(p{c}',1,1,[])), somaProj(s), 'UniformOutput',false)),[2 3 1]);
            end
            for d = 1:size(weightedPSTHS,3)
                if(n>=2)
                    if(n==2)
                        ax{end+1}=nexttile(((s)*tl.GridSize(end))+d);
                        plotSegs{end+1} = 1:length(dimCond);
                    else
                        ax{end+1}=nexttile((c)*tl.GridSize(end)+d);
                        plotSegs{end+1} = c;
                    end
                    hold on;title("Dim " + num2str(d));
                end
                p=plot(squeeze(weightedPSTHS(:,:,d)),'LineWidth',.5+(s-1),'LineStyle',lc{s});
                arrayfun(@(pc) set(pc,'Color',[co(d,:),min(1,.15*(s+1))]),p);
            end
        end
    end
    condSegs = cellfun(@(i) resize(findBins(params.bins(i),timeBins(1):1/(1000/binWidth):max(params.bins)),[1,length(epochSegs)],'FillValue',NaN),segInds,'UniformOutput',false);
    cellfun(@(a,p) cellfun(@(x,s) plot(a,[x;x], repmat(get(a,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',2),num2cell(cell2mat(...
        cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(dimCond(p)))),ax,plotSegs);
    condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
    cellfun(@(a) arrayfun(@(x) plot(a,[x,x], get(a,'YLim'),[char('k'-(4*double(x==max(0)))),'--']),condSegs(1:end-2)),ax);
    cellfun(@(a) set(a,'XLim',[0 250]),ax);
    cellfun(@(a) set(a,'XTick',0:50:250),ax);
    cellfun(@(a) set(a,'XTickLabels',linspace(timeBins(1),timeBins(end),length(get(a,'XTick')))),ax);
    if(saveFig)
        rootSave = "NOMEAN_WeightedPSTHS";
        if(n==1)
            fileNameSave = rootSave+"_DIM";
        elseif(n==2)
            fileNameSave = rootSave+"_COND";
        else
            fileNameSave=rootSave+"_SOMA";
        end
        saveFigures(gcf,savePath+"Full_Signal\",fileNameSave+"_Sm"+num2str(smoothWin),[]);
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
if(saveFig)
    saveFigures(gcf,savePath,model,[]);
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
if(saveFig)
    saveFigures(gcf,savePath,model+"_CondByDim",[]);
end

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end