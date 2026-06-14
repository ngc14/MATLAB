model = "GilliganSkipper_ArmHand";
type = 'Traj';
centered="on";
PCATime = false;
plotTrials = false;
num_dims=5;
sTrials = 20;
timeBins = [-.5, 1];
conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
epochSegs = ["GoSignal","StartReach","StartHold","StartWithdraw"];
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
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) sum(m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end)))),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct,cs,'UniformOutput',false),tablePSTHD,allSegs,phaseConds,cellfun(@(n) n{p},phaseWin,'UniformOutput',false),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) round(downsampleTrials(cat(2,m{~cellfun(@isempty,m)}),sTrials)), [trialFR{:}], 'UniformOutput',false);
        %trialFRMat{p} = cellfun(@(m) round(mean(cat(2,m{~cellfun(@isempty,m)}),2,'omitnan')), [trialFR{:}], 'UniformOutput',false);
    end
    fullResPSTH = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(cell2mat(cellfun(@(r)reshape(r,1,1,[]),c,'UniformOutput',false)),...
        [1 2])),arrayfun(@(t) m(all(n>=sTrials,2),contains(conditions,t)),conditions,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    mv = mv & all(cell2mat(cellfun(@(a) ~all(cell2mat(reshape([a{:}],1,size(a{1},1),[]))==0,[2,3]), fullResPSTH, 'UniformOutput',false)),2);
    fullResPSTH = cellfun(@(c) cellfun(@(t) t(mv,:),c,'UniformOutput',false),[fullResPSTH{:}],'Uniformoutput',false);
    smoothedPSTH = fullResPSTH;
    % condPhase = cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
    %     cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)'));
else
    ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
    if(~plotTrials)
        fullResPSTH = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-0),2,'omitnan')',n,'UniformOutput',false))},tablePSTHD,'UniformOutput',false);
        dSegs = cellfun(@(v) {squeeze(cell2mat(reshape(cellfun(@(d) mean(d,1,'omitnan'),...
            v,'Uniformoutput',false),1,1,[])))'},vertcat(allSegs), 'UniformOutput',false);
    else
        fullResPSTH= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
        dSegs = cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(d',sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(allSegs), 'UniformOutput',false);
    end
    fullResPSTH =  cellfun(@(a)cellfun(@(u)max(0,u),a,'UniformOutput',false), fullResPSTH,'UniformOutput',false);
    if(plotTrials)
        mv = mv & sum(cell2mat(cellfun(@(m)mean(cell2mat(m),2,'omitnan').*1000>1,num2cell([fullResPSTH{:}],2),'UniformOutput',false)'),2)>sTrials/2;
    else
        mv = mv & mean(cell2mat(cellfun(@(n)cell2mat(n),fullResPSTH,'UniformOutput',false)),2,'omitnan').*1000>1;
    end
    trialLength = floor(size(fullResPSTH{1}{1}, 2) / binWidth);
    fullResPSTH = cellfun(@(c) cellfun(@(t) t(mv,:),c,'UniformOutput',false),fullResPSTH,'Uniformoutput',false);
    smoothedPSTH = cell(1,length(fullResPSTH));
    for n = 1:length(fullResPSTH)
        smoothedPSTH{n} = repmat({NaN(sum(mv),trialLength)},max(1,sTrials*plotTrials),1);
        for s = 1:length(smoothedPSTH{n})
            for t = 1:trialLength
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth *t;
                smoothedPSTH{n}{s}(:,t) = sum(fullResPSTH{n}{s}(:,iStart:iEnd),2);%normpdf(ceil(3*smoothWin/binWidth)*binWidth:binWidth:binWidth*ceil(3*smoothWin/binWidth),0,smoothWin)
            end
        end
    end
    smoothedPSTH = cellfun(@(c) cellfun(@(s) (conv2(resize(s,[size(s,1),trialLength-1+floor(smoothWin/binWidth)],'Pattern','edge','side','both'),...
        transpose(gausswin(ceil(smoothWin/binWidth)))./sum(gausswin(ceil(smoothWin/binWidth))),'valid')),c,'UniformOutput',false),smoothedPSTH,'UniformOutput',false);
    % normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) permute(mean(c(:,findBins(s-5,params.bins(...
    %     1:binWidth:end)):findBins(s-(1+rand(1)*3),params.bins(1:binWidth:end))),2,'omitnan').*binWidth,[1 3 2]),squeeze(num2cell(a,[2])),...
    %     num2cell([ones(1,all(isnan(n(:,1))));n(~isnan(n(:,1)),1)]),'UniformOutput',false),[1,max(1,sum(~isnan(n(:,1))))])),3,'omitnan'));],p,t,'UniformOutput',false),smoothedPSTH,dSegs,"UniformOutput",false);
    %smoothedPSTH = cellfun(@(c,n) cellfun(@(s,b) s(:,:), c,n,'UniformOutput', false),smoothedPSTH, normBaseline, 'UniformOutput',false);
end
somaLabs = string(somaTable(mv));
channels = allChannels(mv);
location = allLocations(mv,:);
segInds = cellfun(@(s) fix(s(mv,~all(isnan(s),1))),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),cellfun(@(aa)cell2mat(cellfun(@(a) ...
    findBins(mean(a,1,'omitnan'),params.bins),aa,'UniformOutput',false)),allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
segVals = cellfun(@(i) findBins(params.bins(i),timeBins(1):1/(1000/binWidth):timeBins(end)),segInds,'UniformOutput',false);
%% plot single dimension PCA
saveDir = "S:\Lab\ngc14\Working\";
somaIndex = contains(somaLabs,["Arm","Hand"]);
somaReps = unique(somaLabs(somaIndex));
if(PCATime & ~strcmpi(type,"state"))
    savePre = saveDir + "PCA_Time\";
    taskPSTHD = cellfun(@(s) cellfun(@(t) t(somaIndex,unique(round(ms_bins./binWidth))'),s,'UniformOutput',false),smoothedPSTH, 'UniformOutput',false);
    pcaMat = cell2mat(cellfun(@cell2mat,taskPSTHD','UniformOutput',false));
    rowWeights = repmat(arrayfun(@(s) 1/(sum(strcmp(somaLabs,s))),somaLabs)./2,max(1,plotTrials*sTrials)*length(conditions),1)./max(1,plotTrials*sTrials)*length(conditions);
    variableWeights = ones(1,size(pcaMat,2));
else
    savePre = saveDir + "DataHigh\";
    taskPSTHD= cellfun(@(s) cellfun(@(t) t(somaIndex,unique(round(ms_bins./binWidth))),... cell2mat(cellfun(@(s) s(randperm(length(s))),num2cell(
        s,'UniformOutput',false)',smoothedPSTH, 'UniformOutput',false);
    pcaMat = cell2mat(cellfun(@cell2mat,taskPSTHD,'UniformOutput',false))';
    pcaMat = bsxfun(@minus, pcaMat, mean(pcaMat,2));%zscore(taskPSTHD);
    variableWeights = ones(1,size(pcaMat,2));%arrayfun(@(s) 1/(sum(strcmp(somaLabs(somaIndex),s))),somaLabs(somaIndex))./2;
    rowWeights = ones(1,size(pcaMat,1));
end
[loadings,scores,eig,ts,exp] = pca(pcaMat,'Economy',false,'Centered',centered,...
    'Algorithm','svd','VariableWeights',variableWeights,'Weights',rowWeights);somaProj={};
%
if(PCATime & ~strcmpi(type,"state"))
    scores= permute(reshape(scores,sum(mv),max(1,sTrials*plotTrials),length(taskPSTHD),size(loadings,1)),[1 4 2 3]);
    somaProj = cellfun(@(m) cellfun(@transpose,squeeze(num2cell(m,[1 2]))','UniformOutput',false),squeeze(num2cell(scores,[1:length(size(scores))-1]))','UniformOutput',false);
    somaProj = arrayfun(@(s)cellfun(@(c)cellfun(@transpose,permute(num2cell(cell2mat(reshape(cellfun(@(t)cell2mat(arrayfun(@(n) ...
        loadings(:,n)*mean(t(n,strcmp(somaLabs(somaIndex),s)).*binWidth,2,'omitnan')',1:num_dims,'UniformOutput',false)),c,'UniformOutput',false),1,1,[])),1),[2 3 1]),'UniformOutput',false),...
        somaProj,'UniformOutput',false),somaReps, 'UniformOutput',false);
    loadings = squeeze(mean(scores,[3,4]));
else
    mZ =mean(pcaMat)'; sZ=std(pcaMat)';%bsxfun(@rdivide,t,sZ);
    taskPSTHD = cellfun(@(c) cellfun(@(t) bsxfun(@minus,t,mZ),c,'UniformOutput',false),taskPSTHD, 'UniformOutput',false);
    condScores= arrayfun(@(c) reshape(scores([1:(max(1,sTrials*plotTrials)*size(taskPSTHD{1}{1},2))]+(c*(max(1,sTrials*plotTrials)...
        *size(taskPSTHD{1}{1},2))),:),size(taskPSTHD{1}{1},2),[],size(scores,2)),0:length(conditions)-1,'UniformOutput',false);
    somaProj = projectData(somaLabs(somaIndex),taskPSTHD,num_dims,loadings);
    first3Proj = cellfun(@(c) cellfun(@(t) transpose(t'*loadings(:,1:3)),c,'UniformOutput',false),taskPSTHD,'UniformOutput',true);
    if(strcmpi(type,'state'))
        allClusters = cellfun(@(s) cellfun(@cell2mat,s,'UniformOutput',false),somaProj,'UniformOutput',false);
        allClusters = cellfun(@(c) cell2mat(cellfun(@(s) mean(s,2,'omitnan'),c,'UniformOutput',false)),allClusters,'UniformOutput',false);
        % somatotopy = mean(abs(somaProj{1}-somaProj{2}),2);
        phase = mean([abs(allClusters{1}(:,1:3)-allClusters{1}(:,4:6))],2);%,abs(allClusters{2}(:,1:3)-allClusters{2}(:,4:6))],2);
        condition = mean([abs(allClusters{1}(:,1)-allClusters{1}(:,2)),abs(allClusters{1}(:,1)-allClusters{1}(:,3)),abs(allClusters{1}(:,2)-allClusters{1}(:,3)),abs(allClusters{1}(:,4)-allClusters{1}(:,6)),abs(allClusters{1}(:,4)-allClusters{1}(:,5)),abs(allClusters{1}(:,5)-allClusters{1}(:,6)),...
            ],2);%abs(allClusters{2}(:,1)-allClusters{2}(:,2)),abs(allClusters{2}(:,1)-allClusters{2}(:,3)),abs(allClusters{2}(:,2)-allClusters{2}(:,3)),abs(allClusters{2}(:,4)-allClusters{2}(:,6)),abs(allClusters{2}(:,4)-allClusters{2}(:,5)),abs(allClusters{2}(:,5)-allClusters{2}(:,6))],2);
    end
                %somaProj{s}{c}(n,:) =  cellfun(@(w)w'*mean(loadings(somaLabs==somaReps(s),n),1)',squeeze(num2cell(condScores{c}(:,:,n),1)),'UniformOutput',false)';
end
if(strcmpi(centered,'on'));savePre = savePre+"Centered\";else;savePre = savePre+"Non-Centered\";end
newSaveDir = savePre + type +"\";
if(~plotTrials & strcmpi(type,'traj'));newSaveDir = newSaveDir{1}(1:end-1)+"_AVG\";end
for p = 0:(double(~contains(type,"traj",'IgnoreCase',true)))
    if(~contains(type,"traj",'IgnoreCase',true));savePath = newSaveDir+phases(p+1)+"\";else; savePath = newSaveDir;end
    plotProj(loadings,exp,segVals,cellfun(@(s) s([1:length(conditions)]+(length(conditions)*p)),somaProj,'UniformOutput',false),...
        somaLabs(somaIndex),num_dims,conditions,timeBins,binWidth,false,savePath);
end
%% DATAHIGH Dim Reduce
somaColors = {[0 1 .2],[.8 0 1]}; somaReps={["Arm","Hand"]};
if(PCATime & ~strcmpi(type,"state"))
    projectUnits = arrayfun(@(s) cellfun(@(n) num2cell(n,1),squeeze(num2cell(scores(contains(somaLabs(somaIndex),somaReps(s)),:,:,:),[1 2])), 'UniformOutput',false),1:length(somaReps), 'UniformOutput',false);
    projectUnits = cellfun(@(s) cellfun(@(t) cell2mat(cellfun(@(n,l) loadings(:,l)*mean(n),t,num2cell(1:length(t)),'UniformOutput',false))',s,'UniformOutput',false),projectUnits,'UniformOutput',false);
    projectUnits{end+1} = cellfun(@(a) cell2mat(cellfun(@(n,l) loadings(:,l)*mean(n),a,num2cell(1:length(a)),'UniformOutput',false))',cellfun(@(n) num2cell(n,1),squeeze(num2cell(scores,[1 2])), 'UniformOutput',false),'UniformOutput',false);
else
    taskPSTHD= cellfun(@(s) cellfun(@(t) t(somaIndex,:),s,'UniformOutput',false)',fullResPSTH, 'UniformOutput',false);
    pcaMat = (cell2mat(cellfun(@cell2mat,taskPSTHD,'UniformOutput',false))');
    projectUnits = arrayfun(@(r)  squeeze(cellfun(@squeeze,num2cell(num2cell(reshape((pcaMat(:,contains(somaLabs(somaIndex),somaReps{r}))...*loadings(contains(somaLabs(somaIndex),somaReps(r)),:)
        )',sum(contains(somaLabs(somaIndex),somaReps{r})),size(taskPSTHD{1}{1},2),max(1,plotTrials*sTrials),[]),[1 2]),3),'UniformOutput',false))',1:length(somaReps),'UniformOutput',false);
    if(strcmpi(type,'state'))
        projectUnits = cellfun(@(s) squeeze(num2cell(reshape(cellfun(@(t) cell2mat(t'),s','UniformOutput',false),[],length(conditions),length(phases)),1)), projectUnits, 'UniformOutput',false);
        projectUnits = cellfun(@(s) reshape(cellfun(@(t) cell2mat(t'),s,'UniformOutput',false),[],1),projectUnits, 'UniformOutput',false);
    else
        projectUnits = cellfun(@(s) cellfun(@(t) cell2mat(cellfun(@(m)conv(m,transpose(gausswin(smoothWin))./sum(gausswin(smoothWin)),'same')...
            ./(1/1000),num2cell(t,2),'UniformOutput',false)),vertcat(s{:}),'Uniformoutput',false), projectUnits,'Uniformoutput',false);
        projectUnits = cellfun(@(s) cellfun(@(t) t(:,findBins(timeBins(1),params.bins):end),s, 'UniformOutput',false),projectUnits, 'UniformOutput',false);
    end
end
if(length(projectUnits)>length(somaReps))
    % somaReps(end+1) = "Arm+Hand"; cls{end+1} = colors.values';
    % projectUnits{end+1} =cellfun(@(s) {vertcat(s{:})}, squeeze(cellfun(@squeeze,num2cell(num2cell(reshape((pcaMat(:,contains(somaLabs(somaIndex),somaReps))*...
    %     loadings(contains(somaLabs(somaIndex),somaReps),:))',size(pcaMat,2),size(taskPSTHD{1}{1},2),max(1,plotTrials*sTrials),[]),[1 2]),3),'Uniformoutput',false))', 'UniformOutput',false);
end
cls = cellfun(@(r) repmat({r},max(strcmpi(type,'traj')*(plotTrials*0),1),1),somaColors','UniformOutput',false);
cls = reshape(cellfun(@(v) repmat({v},length(conditions),1),vertcat(cls{:}),'UniformOutput',false),[],1);
cTrials = arrayfun(@(d) {repmat(d,max(strcmpi(type,'traj')*plotTrials*0,1),1)},conditions,'UniformOutput',false)';
cTrials = cellfun(@(p) cellfun(@(c) cell2mat(cellfun(@(t) string(t)+"-"+join(p),c,'UniformOutput',false)),cTrials,'UniformOutput',false),somaReps,'UniformOutput',false);
epochStarts = repmat(cellfun(@(e) e(:,all(e<size(taskPSTHD{1}{1},2),1)),num2cell(cell2mat(cellfun(@(s) round(mean((s-min(ms_bins)),1,'omitnan')),...%./binWidth,1,'omitnan')),...
    segInds,'UniformOutput',false)),2),'UniformOutput',false),length(somaReps),1);
if(strcmpi(type,'traj'))
    dHiStruct = struct('data',vertcat(projectUnits{:}),'condition',cellstr(cell2mat(vertcat(cTrials{:}))),'epochStarts',epochStarts,...
        'epochColors',cellfun(@(r,e)cell2mat(arrayfun(@(n) [hsv2rgb([r(1)*(1*n~=4),min(1,.5+(.2*n))*(1*n~=4),max(1*n==4,1.25-(.25*n))])],1:length(e),"UniformOutput",false)'),...
        num2cell(rgb2hsv(cell2mat(cellfun(@cell2mat,cls,'UniformOutput',false))),2),epochStarts,'UniformOutput',false));
else
    cTrials = cellfun(@(c) cell2mat(arrayfun(@(r) c+"-"+r,phases,'UniformOutput',false)'),cTrials,'UniformOutput',false);
    cls={};cls{end+1} = colors.values';
    cls = cellfun(@(p) repmat(p,length(phases),1),cls,'UniformOutput',false);
    dHiStruct = struct('data',vertcat(projectUnits{:}),'condition',cellstr(vertcat(cTrials{:})),'epochColors',cellfun(@rgb2hsv,vertcat(cls{:}),'UniformOutput',false));
    [dHiStruct.type]=deal('state');
    [dHiStruct.epochStarts]=deal(1);
    for a = 1:length(dHiStruct)
        dHiStruct(a).epochColors = hsv2rgb(max(min(dHiStruct(a).epochColors,1),0));%+([0 .25 -.25]*2*(contains(dHiStruct(a).condition,"Reach")-.5)),1),0));
    end
end
% DataHigh(dHiStruct); 
% tStruct = struct("A",cellfun(@transpose,{dHiStruct.data},'UniformOutput',false),'condition',{dHiStruct.condition},'epochStarts',{dHiStruct.epochStarts},'epochColors',{dHiStruct.epochColors});
% [Q_out, T_out] = tangleAnalysis(tStruct, params.binSize,'softenNorm',5 ,'timeStep',20,'withinConditionsOnly',true,'numPCs',20);tangle_visualize(T_out);
% formattedScatter(Q_Arm, Q_Hand, {'Q_{Arm}','Q_{Hand}'});
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
function sp = projectData(somaLabs,rawData,num_dims,loadings)
somaReps = unique(somaLabs);
sp = repmat({cell(1,length(rawData))},length(somaReps),1);
for s = 1 :length(somaReps)
    for c = 1:length(rawData)
        for n = 1:num_dims
            sp{s}{c}(n,:) =  cellfun(@(t) transpose(t(strcmp(somaLabs,somaReps(s)),:)'*loadings(strcmp(somaLabs,somaReps(s)),n)),rawData{c},'UniformOutput',false);
        end
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