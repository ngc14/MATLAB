conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
tPhys = unitTable(conditions,params);
%%
model = "GilliganSkipper_ArmHand";
type = 'Traj';
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
saveFig = true;
num_dims=4;
sTrials = 30;
plotTrials = 1;
timeBins = [-.5, 2.5];
splitGroup = "Condition";
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
savePath = saveDir+type+"\Full_Population\";
dimCond = reshape(extractAfter(model,"_")+"_"+params.condAbbrev.values',1,[]);%regexp(model,'[A-Z]+[^A-Z]+','match')
colors = {[.7 0 0],[1 .65 0 ],[0 0 .75],[1 0 .3],[1 1 0],[0 .6 1]}';
if(~exist(savePath,'dir')), mkdir(savePath); end
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
avgTrace = zeros(size(avgTrace,1),size(avgTrace,2));
numUnits = cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), taskPSTHD,'UniformOutput',false));
if(~plotTrials)
    taskPSTHD = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-avgTrace),2,'omitnan')',n,'UniformOutput',false))}, taskPSTHD, 'UniformOutput',false);
else
    taskPSTHD= cellfun(@(a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-avgTrace),sTrials),...
        a(all(numUnits>=sTrials,2)),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(taskPSTHD), 'UniformOutput',false);
end
unitInds = repmat({randperm(sum(all(numUnits>=sTrials,2)))},1,length(taskPSTHD));%repmat({}',1,size(currD,2)/(length(numUnits)));
taskPSTHD =  cellfun(@(a,i,b)cellfun(@(u)max(0,u(i,findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins))),a,'UniformOutput',false), ...
    taskPSTHD, unitInds,allSegs,'UniformOutput',false);%cellfun(@(a,b) max(0,mean(b(max(1,findBins(mean(a(:,2))-5,params.bins)):max(1/params.binSize,findBins(mean(a(:,2))-4,params.bins))))),s(i),num2cell(u(i,:),2))),...
cls = cellfun(@(r) repmat({r},max(plotTrials*sTrials,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),5);...
    linspace(1,.25,5);linspace(.85,1,5)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
segInds = cellfun(@(s) fix(mean(s,1,'omitnan')),cellfun(@(n) [ones(size(n,1),1),n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),...
    epochSegs(1:end-1))),min(n(:,strcmp(maxSegL,epochSegs(end))),length(timeBins(1):params.binSize:timeBins(end))-1)],...
    cellfun(@(aa,i)cell2mat(cellfun(@(a) findBins(mean(a,1,'omitnan'),params.bins(findBins(timeBins(1),params.bins):...
    findBins(timeBins(end),params.bins))),aa(i),'UniformOutput',false)),allSegs,unitInds,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
segInds = cellfun(@(r) repmat(r,sTrials*plotTrials,1), segInds,'UniformOutput',false);
dHiStruct = struct('data',vertcat(taskPSTHD{:}),'epochStarts',vertcat(segInds{:}),'condition',cellstr(cell2mat(...
    cellfun(@(d) repmat(string(d),max(plotTrials*sTrials,1),1),cellstr(dimCond),'UniformOutput',false)')),'epochColors',vertcat(cls{:}));
%%
binWidth = 10;smoothWin = 100;
trialLength = floor(size(dHiStruct(1).data, 2) / binWidth);
mv = mean([dHiStruct.data],2) * 1000;
somaLabs = somaTable(find(all(numUnits>=sTrials,2)).*double(mv>0));
for n = 1:length(dHiStruct)
    smoothedData{n} = NaN(sum(mv>0),trialLength);
    for t = 1:trialLength
        iStart = binWidth * (t-1) + 1;
        iEnd   = binWidth *t;
        smoothedData{n}(:,t) = sum(dHiStruct(n).data(mv>0,iStart:iEnd),2);
    end
end
smoothedData = cellfun(@(s) sqrt(resize(conv2(s,transpose(gausswin(ceil(smoothWin/binWidth))./sum(gausswin(ceil(smoothWin/binWidth)))),'valid')...
    ./(binWidth/1000),size(s),'Pattern','edge','Side','both')),smoothedData,'UniformOutput',false);
%smoothedData=cellfun(@(s)s(:,8+1:end-8),smoothedData,'UniformOutput',false);normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
[loadings, scores, eig] = pca([cell2mat(smoothedData)]','Economy',false,'Centered','off');
index = 0;D = dHiStruct;
for i=1:length(smoothedData)
    D(i).data= scores(index + (1:size(smoothedData{i},2)),:)';
    D(i).epochStarts = max(1,fix(D(i).epochStarts/binWidth));
    index = index + size(smoothedData{i},2);
end
%%
DataHigh(dHiStruct,'DimReduce');
if(saveFig)
    save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
end
% all_h = findall(groot,'Type','Figure');handles = guihandles(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));D= D.D;
% D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));
%%
somaReps = unique(somaLabs);
if(length(D)<sTrials)
    plotTrials = 0;
end
somaProj = cellfun(@(d)cellfun(@(s,e) e(somaLabs==s,:),num2cell(unique(somaLabs)),repmat({loadings*d},length(somaReps),1),'UniformOutput',false),...
    smoothedData,'UniformOutput',false);somaProj = reshape([somaProj{:}]',[],1);
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
                weightedPSTHS = smoothedData{i}(somaLabs==somaReps(s),:)';
            elseif(n==num_dims+1)
                weightedPSTHS = somaProj{i+(length(dimCond)*(s-1))}';
            else
                weightedPSTHS = loadings(:,n).*[smoothedData{i}];
                weightedPSTHS = weightedPSTHS(somaLabs==somaReps(s),:)';
            end
            plot(mean(weightedPSTHS,2,'omitnan'),'LineWidth',s,'Color',cell2mat(colors.values(cellstr(dimCond(i)))),'LineStyle',ls);
        end
    end
end
%l = arrayfun(@(d) arrayfun(@(l) plot(NaN(1,1),'LineStyle',l,'Color',cell2mat(colors.values(cellstr(d))))',["-",":"]),dimCond,'UniformOutput',false);
legend(reshape(cell2mat(arrayfun(@(a) a+"_"+(params.condAbbrev.values), string(somaReps), 'UniformOutput', false)),1,[]),...
    'Autoupdate','off','FontSize',14,'Orientation','horizontal');
if(saveFig)
    saveFigures(gcf,savePath,"Projections_"+num2str(num_dims)+"_Bsz"+num2str(binWidth)+"_Sm"+num2str(smoothWin),[]);
end
%%
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
    currColor = rgb2hsv(cell2mat(colors.values(cellstr(dimCond(icond)))));%([1 0 0]);
    currColor(2) = 1;currColor(end) = .4;
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
            if(iepoch>1)
                line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
            end
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
% %%
% switch(splitGroup)
%     case "Somatotopy"
%         condInds = arrayfun(@(c) contains(c,'Arm'), dimCond);
%         for u = 1:length(condInds)
%             if(condInds(u)==0),colors{u} = [.8 .8 .8];
%             else,colors{u} = [.2 .2 .2];end
%         end
%     case  "Condition"
%         condInds = arrayfun(@(c) contains(c,'S'), dimCond);
%         condInds = condInds + arrayfun(@(c) contains(c,'_E'), dimCond);
%         for u = 1:length(condInds)
%             if(condInds(u)==0),colors{u} = [0 0 1];
%             elseif(condInds(u)==1),colors{u} = [1 .85 0];
%             else,colors{u} = [1 0 0];end
%         end
%     case "Phase"
%         condInds = arrayfun(@(c) contains(c,'Reach'), dimConds);
%         for u = 1:length(condInds)
%             if(condInds(u)==0),colors{u} = [1 0 1];
%             else,colors{u} = [0 1 1];end
%         end
%     case "Monkey"
%         condInds = arrayfun(@(c) contains(c,'Gilligan'), dimConds);
%         for u = 1:length(condInds)
%             if(condInds(u)==0),colors{u} = [.8 .4 0];
%             else,colors{u} = [0 .5 0];end
%         end
% end
function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end