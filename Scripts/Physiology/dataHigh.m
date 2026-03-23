conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
tPhys = unitTable(conditions,params);
%%
model = "GilliganSkipper_ArmHand";
type = 'Traj';
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
saveFig = true;
num_dims=4;
sTrials = 20;
plotTrials = 0;
timeBins = [-.5, 2.5];
splitGroup = "Condition";
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
savePath = saveDir+type+"\PCA_Time\";
dimCond = reshape(extractAfter(model,"_")+"_"+params.condAbbrev.values',1,[]);%regexp(model,'[A-Z]+[^A-Z]+','match')
colors = {[.7 0 0],[1 .65 0 ],[0 0 .75],[1 0 .3],[1 1 0],[0 .6 1]}';
colors =containers.Map(dimCond,colors(1:length(dimCond)));
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
tableInds = contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]);
somaTable = tPhys{tableInds,"Somatotopy"};
allSegs= arrayfun(@(s) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))},dimCond,'UniformOutput',false);%
taskPSTHD= arrayfun(@(a) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
%avgTrace = mean(cell2mat(vertcat(taskPSTHD{:})'),2,'omitnan');
avgTrace = zeros(length(params.bins),1);
numUnits = all(cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), taskPSTHD,'UniformOutput',false))>=sTrials,2);
somaTable = somaTable(numUnits);
if(~plotTrials)
    taskPSTHD = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-avgTrace),2,'omitnan')',n(numUnits),'UniformOutput',false))}, taskPSTHD, 'UniformOutput',false);
else
    taskPSTHD= cellfun(@(a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-avgTrace),sTrials),...
        a(numUnits),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(taskPSTHD), 'UniformOutput',false);
end
unitInds = repmat({randperm(sum(numUnits))},1,length(taskPSTHD));%repmat({}',1,size(currD,2)/(length()));
ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
taskPSTHD =  cellfun(@(a,i,b)cellfun(@(u)max(0,u(i,ms_bins)),a,'UniformOutput',false), ...
    taskPSTHD, unitInds,allSegs,'UniformOutput',false);%cellfun(@(a,b) max(0,mean(b(max(1,findBins(mean(a(:,2))-5,params.bins)):max(1/params.binSize,findBins(mean(a(:,2))-4,params.bins))))),s(i),num2cell(u(i,:),2))),...
segInds = cellfun(@(s) fix(mean(s,1,'omitnan')),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),...
    cellfun(@(aa,i)cell2mat(cellfun(@(a) findBins(mean(a,1,'omitnan'),params.bins(ms_bins)),aa(i),'UniformOutput',false)),...
    allSegs,unitInds,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
cls = cellfun(@(r) repmat({r},max(plotTrials*sTrials,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),5);...
    linspace(1,.25,5);linspace(.85,1,5)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
%% smooth data and remove non-modulated units
timePCA =  1; binWidth = 10;smoothWin = 150;
trialLength = floor(size(taskPSTHD{1}{1}, 2) / binWidth);
mv = mean(cell2mat(cellfun(@(n) cell2mat(n),taskPSTHD,'UniformOutput',false)),2,'omitnan').*1000>1;
if(length(taskPSTHD)<sTrials)
    plotTrials = 0;
end
for n = 1:length(taskPSTHD)
    smoothedData{n} = NaN(sum(mv),trialLength);
    for t = 1:trialLength
        iStart = binWidth * (t-1) + 1;
        iEnd   = binWidth *t;
        smoothedData{n}(:,t) = sum(taskPSTHD{n}{1}(mv,iStart:iEnd),2);
    end
end
smoothedData = cellfun(@(s) conv2(resize(s,[size(s,1),size(s,2)-1+floor(smoothWin/binWidth)],'Pattern','edge','side','both'),...
    transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid'),smoothedData,'UniformOutput',false);
somaLabs = somaTable(mv);
somaReps = unique(somaLabs);
if(timePCA)
    pcaMatrix = cell2mat(smoothedData');
else
    pcaMatrix = cell2mat(smoothedData)';
end
[loadings, scores, eig] = pca(pcaMatrix,'Economy',false,'Centered','on','NumComponents',num_dims);
%smoothedData=cellfun(@(s)s(:,8+1:end-8),smoothedData,'UniformOutput',false);normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
%% plot single dimension PCA
if(timePCA)
    somaProj = cellfun(@(s)arrayfun(@(c)s((1+(c-1)*size(s,1)/length(dimCond)):c*size(s,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),...
        cellfun(@(s,e) scores(repmat(somaLabs,length(dimCond),1)==s,:),num2cell(somaReps),'UniformOutput',false),'UniformOutput',false);
else
    %somaProj = cellfun(@(s)arrayfun(@(c) scores((1+(c-1)*size(scores,1)/length(dimCond)):c*size(scores,1)/length(dimCond),:)',1:length(dimCond),'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
    somaProj = cellfun(@(s)arrayfun(@(c) cell2mat(arrayfun(@(n)mean(pcaMatrix((1+(c-1)*size(pcaMatrix,1)/length(dimCond)):c*size(pcaMatrix,1)/length(dimCond),somaLabs==s)'.*...
        loadings(somaLabs==s,n),1,'omitnan'),1:num_dims,'UniformOutput',false)').*1000,1:length(dimCond),'UniformOutput',false),num2cell(somaReps),'UniformOutput',false);
end
bc = num2cell([1 0 0;1 .7 0; 0 0 1],2);
lc = flipud(num2cell([0 0 0; .7 .7 .7;],2));
figure(); tiledlayout(2,1+num_dims);
for n = 0:num_dims+1
    if(n>num_dims);nexttile([1,num_dims+1]); hold on;
    else;nexttile(); hold on;
    end
    for i = 1:length(dimCond)
        for s =1:length(somaReps)
            if(s==1)
                if(n>0 && n<=num_dims);title("Dim " + num2str(n));
                elseif(n>0);title("Weighted PSTHS");
                else;title("PSTHS");
                end
                ls = '-';
            else;ls = ':';
            end
            if(n==0)
                weightedPSTHS = mean(smoothedData{i}(somaLabs==somaReps(s),:)',2,'omitnan');
                plot(weightedPSTHS,'LineWidth',s,'Color',cell2mat(colors.values(cellstr(dimCond(i)))),'LineStyle',ls);
            elseif(n==num_dims+1)
                weightedPSTHS = boxchart(reshape(repmat((i-1)*length(somaReps)+(s/2+0:10:num_dims*10),size(somaProj{s}{i},2),1),1,[]),...
                    reshape(somaProj{s}{i},1,[]),'WhiskerLineStyle','-','Notch','on','BoxWidth',.5,'BoxFaceColor',bc{i},'BoxEdgeColor',lc{s},'MarkerStyle','none');
            else
                condSomaInd = ismember(1:length(pcaMatrix),((1+(i-1)*size(somaLabs,1)):i*size(somaLabs,1)));
                weightedPSTHS = mean(cell2mat(cellfun(@(s) s(n,:), somaProj{s}(i),'Uniformoutput',false)),1+timePCA,'omitnan');%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
                if(timePCA)
                    weightedPSTHS = loadings(:,n).*weightedPSTHS;
                end
                plot(weightedPSTHS,'LineWidth',s,'Color',cell2mat(colors.values(cellstr(dimCond(i)))),'LineStyle',ls);
            end
        end
    end
    if(n<num_dims+1)
        condSegs = findBins(params.bins(ms_bins(cell2mat(segInds))),timeBins(1):1/(1000/binWidth):timeBins(end));
        cellfun(@(x,s) plot([x;x], repmat(get(gca,'YLim'),size(x,2),1)','LineStyle','--','Color',s),num2cell(condSegs(:,contains(epochSegs,"Withdraw")),2)',colors.values)
        condSegs = round(mean(condSegs(1,:),1,'omitnan'));
        arrayfun(@(x) plot([x,x], get(gca,'YLim'),[char('k'-(4*double(x==max(condSegs)))),'--']),condSegs(1:end-1));
    end
end
legend(reshape(cell2mat(arrayfun(@(a) a+"_"+(params.condAbbrev.values), string(somaReps), 'UniformOutput', false)),1,[]),...
    'Autoupdate','off','FontSize',14,'Orientation','horizontal');
xticks(3:10:40);
xticklabels(1:4);ylim([-2 2]);
if(saveFig)
    if(timePCA)
        fileNameSave="FactorScoreAvgs_Sm"+num2str(smoothWin);
    else
        fileNameSave="WeightedPSTHSAvgs_Sm"+num2str(smoothWin);
    end
    saveFigures(gcf,savePath,fileNameSave,[]);
end
%% DATAHIGH Dim Reduce
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
    saveFigures(gcf,savePath,model+"_"+splitGroup,[]);
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