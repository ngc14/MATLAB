model = "GilliganSkipper_ArmHand";
type = 'State';
num_dims=5;
sTrials = 20;
plotTrials = 1;
timeBins = [-.5, 1];
conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
params = PhysRecording(string(conditions),.001,.001,-6,3,containers.Map(conditions,repmat({"StartReach"},1,length(conditions))));
dimCond = reshape(extractAfter(model,"_")+"_"+params.condAbbrev.values',1,[]);
colors =containers.Map(dimCond,{[.7 0 0],[1 .65 0 ],[0 0 .75] }');% regexp(model,'[A-Z]+[^A-Z]+','match')
allSegsL = params.condSegMap.values;
[~,maxSegL]= max(cellfun(@length,allSegsL));
maxSegL = allSegsL{maxSegL};
tPhys = unitTable(conditions,params);
%%
tableInds = contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]);
somaTable = tPhys{tableInds,"Somatotopy"};
allChannels = tPhys{tableInds,"Channel"};
allLocations = tPhys{tableInds,["XT","YT"]};
allSegs= arrayfun(@(s) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))},dimCond,'UniformOutput',false);%
tablePSTHD= arrayfun(@(a) tPhys{tableInds,contains(tPhys.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
numUnits = all(cell2mat(cellfun(@(a) cellfun(@(s) size(s,2),a), tablePSTHD,'UniformOutput',false))>=sTrials,2);%
tablePSTHD = cellfun(@(c) c(numUnits), tablePSTHD, 'UniformOutput',false);
allSegs = cellfun(@(s) s(numUnits), allSegs, 'UniformOutput', false);
somaTable = somaTable(numUnits);
allChannels = allChannels(numUnits);
allLocations = allLocations(numUnits,:);
unitInds = randperm(sum(numUnits));%repmat({}',1,size(currD,2)/(length()));
somaTable = somaTable(unitInds);
allChannels = allChannels(unitInds);
allLocations = allLocations(unitInds,:);
clear tPhys tableInds numUnits;
%% smooth data and remove non-modulated units
binWidth = 10;smoothWin = 150;
%avgTrace = cellfun(@(c) mean(cell2mat(cellfun(@(u)mean(u,2,'omitnan'),c,'UniformOutput',false)'),2,'omitnan'),tablePSTHD,'UniformOutput',false);%zeros(size(params.bins))',taskPSTHD,'UniformOutput',false);%
mv = (somaTable=="Hand" | somaTable=="Arm");
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
        m(all(n>=sTrials,2),contains(conditions,t)),conditions,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    smoothedPSTH =  [currD{:}];
    mv = mv & any(cell2mat(cellfun(@(a) ~all(cell2mat(reshape([a{:}],1,size(a{1},1),[]))==0,[2,3]), currD, 'UniformOutput',false)),2);
    taskPSTHD = vertcat(smoothedPSTH{:})';
    condPhase = cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
        cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)'));
else
    ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
    if(~plotTrials)
        msBinPSTH = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-0),2,'omitnan')',n(unitInds),'UniformOutput',false))},tablePSTHD,'UniformOutput',false);
    else
        msBinPSTH= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
    end
    msBinPSTH =  cellfun(@(a)cellfun(@(u)max(0,u(unitInds,:)),a,'UniformOutput',false), msBinPSTH,'UniformOutput',false);
    if(plotTrials)
        mv = mv & sum(cell2mat(cellfun(@(m)mean(cell2mat(m),2,'omitnan').*1000>1,num2cell([msBinPSTH{:}],2),'UniformOutput',false)'),2)>sTrials/2;
    else
        mv = mv & mean(cell2mat(cellfun(@(n)cell2mat(n),msBinPSTH,'UniformOutput',false)),2,'omitnan').*1000>1;
    end
    msBinPSTH = cellfun(@(c) cellfun(@(s) s(mv,:), c, 'UniformOutput',false), msBinPSTH, 'UniformOutput',false);
    trialLength = floor(size(msBinPSTH{1}{1}, 2) / binWidth);
    for n = 1:length(msBinPSTH)
        smoothedPSTH{n} = repmat({NaN(sum(mv),trialLength)},max(1,sTrials*plotTrials),1);
        for s = 1:length(smoothedPSTH{n})
            for t = 1:trialLength
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth *t;
                smoothedPSTH{n}{s}(:,t) = sum(msBinPSTH{n}{s}(:,iStart:iEnd),2);%normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
            end
        end
    end
    smoothedPSTH = cellfun(@(c) cellfun(@(s) (conv2(resize(s,[sum(mv),trialLength-1+floor(smoothWin/binWidth)],'Pattern','edge','side','both'),...
        transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),smoothedPSTH,'UniformOutput',false);
end
somaLabs = somaTable(mv);
channels = allChannels(mv);
location = allLocations(mv,:);
somaReps = unique(somaLabs);
segInds = cellfun(@(s) fix(s(mv,~all(isnan(s),1))),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),cellfun(@(aa)cell2mat(cellfun(@(a) ...
    findBins(mean(a,1,'omitnan'),params.bins),aa(unitInds),'UniformOutput',false)),allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
%% plot single dimension PCA
if(1)
    taskPSTHD = cellfun(@(s) cell2mat(cellfun(@(t) t(:,unique(round(ms_bins./binWidth))),s,'UniformOutput',false)'),smoothedPSTH, 'UniformOutput',false);
    [loadings, scores, eig,~,exp] = pca(cell2mat(taskPSTHD)','Economy',false,'Centered','off','Algorithm','eig');
    scores = arrayfun(@(c)reshape(scores([1:size(taskPSTHD{1},2)]+(size(taskPSTHD{1},2)*c),:),...
        [size(taskPSTHD{1},2)/max(1,plotTrials*sTrials),max(1,plotTrials*sTrials),size(scores,2)]),0:length(conditions)-1,'UniformOutput',false);
    somaProj = arrayfun(@(s) cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)mean(t(somaLabs==s,unique(round(ms_bins./binWidth)))...
        .*loadings(somaLabs==s,n),1,'omitnan'),1:num_dims,'UniformOutput',false)').*binWidth,c,'UniformOutput',false),smoothedPSTH,'UniformOutput',false),somaReps,'UniformOutput',false);
else
    taskPSTHD = cellfun(@(s) cell2mat(cellfun(@(t) t(:,unique(round(ms_bins./binWidth))'),s,'UniformOutput',false)),smoothedPSTH, 'UniformOutput',false);
    [loadings, scores, eig,~,exp] = pca(cell2mat(taskPSTHD'),'Economy',false,'Centered','on','Algorithm','eig');
    scores = arrayfun(@(c)reshape(scores([1:size(taskPSTHD{1},1)]+(size(taskPSTHD{1},1)*c),:),...
        [size(taskPSTHD{1},1)/max(1,plotTrials*sTrials),max(1,plotTrials*sTrials),size(scores,2)]),0:length(conditions)-1,'UniformOutput',false);
    somaProj = arrayfun(@(s) cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)loadings(:,n)*squeeze(mean(t(somaLabs==s,:,n),1,'omitnan')),1:num_dims,...
        'UniformOutput',false))',num2cell(c,[1 3])', 'UniformOutput',false),scores,'UniformOutput',false), somaReps,'UniformOutput',false);
    loadings = cellfun(@(m) cell2mat(cellfun(@(l) squeeze(mean(m(:,:,l),2,'omitnan')),...*l,
        num2cell((1:num_dims),1),'UniformOutput',false)),scores, 'UniformOutput',false);
    loadings = mean(cat(3,loadings{:}),3,'omitnan');
end
segVals = cellfun(@(i) findBins(params.bins(i),timeBins(1):1/(1000/binWidth):timeBins(end)),segInds,'UniformOutput',false);
for p = 0:1
    savePath = "S:\Lab\ngc14\Working\DataHigh\Non-Centered\"+type+"\";
    if(~plotTrials)
        savePath = savePath{1}(1:end-1)+"_AVG\";
    end
    if(~contains(type,"traj"))
        savePath = savePath+phases(p+1)+"\";
    end
    plotProj(loadings,exp,eig,segVals,cellfun(@(s) s([1:length(conditions)]+(length(conditions)*p)),somaProj,'UniformOutput',false),...
        somaLabs,location,num_dims,conditions,timeBins,true,savePath);
end
    %%
[epochGroups,dataGrouped,condNames]= deal(cell(length(somaReps),2,length(conditions)));
for s = 1:length(somaReps)
    for l=0:1
        currInds = somaLabs==somaReps(s) & (channels<=(16*(l+1)) & channels>(16*l));
        dataGrouped(s,l+1,:) = cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n) mean(t(currInds,min(round(ms_bins./binWidth)):end).*...
            loadings(currInds,n),1,'omitnan'),1:dims,'UniformOutput',false)').*binWidth,c,'UniformOutput',false),smoothedPSTH,'UniformOutput',false);
        epochGroups(s,l+1,:) = cellfun(@(s) s(currInds,:), segInds,'UniformOutput',false);
        channelSpan = [min([1:16]+(16*l)), max([1:16]+(16*l))];
        condNames(s,l+1,:) = arrayfun(@(r) repmat(cellstr(r),max(1,sTrials*plotTrials),1), string(somaReps(s))+"_"+num2str(channelSpan)+"_"+conditions,'UniformOutput',false);
    end
end
dataGrouped = reshape(dataGrouped,1,[]);epochGroups = reshape(epochGroups,1,[]);condNames = reshape(condNames,1,[]);
struct('data',vertcat(dataGrouped{:}),'epochStarts',num2cell(cell2mat(epochGroups'),2),'condition',vertcat(condNames{:}));

if(contains(type,"Traj"))
    title("PSTHS");
    psthGrouped = cellfun(@(s)cellfun(@(c)cellfun(@(t) mean(t(somaLabs==s,min(round(ms_bins./binWidth)):end),1,'omitnan').*binWidth,...
        c,'UniformOutput',false),smoothedPSTH,'UniformOutput',false),num2cell(somaReps),'Uniformoutput',false);
    cellfun(@(s,l) cellfun(@(m,p) plot(cell2mat(m)','LineStyle',l,'LineWidth',1.5,'Color',p),s,colors.values),psthGrouped,repmat({'-'},size(psthGrouped,1),size(psthGrouped,2)));
    condSegs = cellfun(@(i) resize(i,[1,length(epochSegs)],'FillValue',NaN),indBins,'UniformOutput',false);
    cellfun(@(x,s) plot(ax,[x;x], repmat(get(ax,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
        cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs,'UniformOutput',false)),2)',colors.values(cellstr(dimCond)));
    condSegs = round(mean(cell2mat(condSegs),1,'omitnan'));
    arrayfun(@(x) plot(ax,[x,x], get(ax,'YLim'),'Color','k','LineStyle',':','LineWidth',1.5),condSegs(1:end-2));
    set(ax,'XTickLabels',linspace(timeBins(1),timeBins(end),length(get(ax,'XTick'))));
    set(ax,'XLim',[0 trialLength-min(round(ms_bins./binWidth))]);
    set(ax,'XTick',0:50:trialLength-min(round(ms_bins./binWidth)));
else
    title("Spikes")
    xVals = cell2mat(arrayfun(@(a) repmat(a,[1,ceil(length(loadings)/10)]),1:10,'UniformOutput',false));
    cellfun(@(p,n,cl) scatter(xVals(1:size(p{1},1))+(n*10),mean(cell2mat(p').*binWidth,2,'omitnan'),10,cl,'filled','o'),...
        smoothedPSTH,num2cell(1:length(smoothedPSTH)),repmat(colors.values,1,length(phases)));%(length(conditions)*(pI-1)+[1:length(conditions)])
end
%% DATAHIGH Dim Reduce
if(1)
    trialTaskPSTHD = cellfun(@(s) cellfun(@(t) t(:,unique(round(ms_bins./binWidth))),s,'UniformOutput',false), smoothedPSTH, 'UniformOutput',false);
    trialTaskPSTHD = cellfun(@(s) arrayfun(@(r) s(:,:,somaLabs==r),somaReps, 'UniformOutput',false),scores, 'UniformOutput',false);
    trialTaskPSTHD = cellfun(@(c) cellfun(@(t) cellfun(@squeeze,num2cell(t(:,:,randperm(size(t,3),389)),[1 3]),'UniformOutput',false),c, 'UniformOutput', false), trialTaskPSTHD, 'UniformOutput',false);
    trialTaskPSTHD = cellfun(@(c) transpose([c{:}]),trialTaskPSTHD,'UniformOutput',false)
else
    trialTaskPSTHD = cellfun(@(s) cellfun(@(t) t(:,unique(round(ms_bins./binWidth)))',s,'UniformOutput',false), smoothedPSTH, 'UniformOutput',false)';
    trialTaskPSTHD=cellfun(@(s) arrayfun(@(r) cellfun(@squeeze,num2cell(s(somaLabs==r,:,:),[1 3]),'UniformOutput',false),...
        somaReps, 'UniformOutput',false)',scores,'UniformOutput',false)';
    trialTaskPSTHD = cellfun(@transpose,[trialTaskPSTHD{:}],'UniformOutput',false);
    trialTaskPSTHD = cellfun(@(c) cellfun(@(t) t(randperm(size(t,1),389),:),c, 'UniformOutput', false), trialTaskPSTHD, 'UniformOutput',false);
end
cls = cellfun(@(r) repmat({r},max((plotTrials*2*sTrials)/4,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),4);...
    linspace(1,.25,4);linspace(.85,1,4)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
epochStarts = num2cell(cell2mat(cellfun(@(s) repmat(round(mean(min((s-min(ms_bins))./binWidth,150),1,'omitnan')),2*sTrials,1), segInds,'UniformOutput',false)),2);
dHiStruct = struct('data',vertcat(trialTaskPSTHD{:}),'epochStarts',epochStarts,'condition',...
    cellstr(cell2mat(arrayfun(@(s) cell2mat(cellfun(@(d) repmat([d+string(s)],max(plotTrials*sTrials,1),1),cellstr(conditions),'UniformOutput',false)'),somaReps,'UniformOutput',false))),...
    'epochColors',cellfun(@(r) repmat(r,4,1),num2cell(cell2mat(cellfun(@cell2mat,cls,'UniformOutput',false)),2),'UniformOutput',false));
DataHigh(dHiStruct);
%save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
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