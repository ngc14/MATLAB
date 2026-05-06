model = "GilliganSkipper_ArmHand";
type = 'Traj';
centered="on";
PCATime = false;
num_dims=5;
sTrials = 20;
plotTrials = 0;
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
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) sum(m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end)))),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct(unitInds),cs(unitInds),'UniformOutput',false),tablePSTHD,allSegs,phaseConds,cellfun(@(n) n{p},phaseWin,'UniformOutput',false),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) round(downsampleTrials(cat(2,m{~cellfun(@isempty,m)}),sTrials)), [trialFR{:}], 'UniformOutput',false);
        %trialFRMat{p} = cellfun(@(m) round(mean(cat(2,m{~cellfun(@isempty,m)}),2,'omitnan')), [trialFR{:}], 'UniformOutput',false);
    end
    currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(cell2mat(cellfun(@(r)reshape(r,1,1,[]),c,'UniformOutput',false)),...
        [1 2])),arrayfun(@(t) m(all(n>=sTrials,2),contains(conditions,t)),conditions,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    mv = mv & all(cell2mat(cellfun(@(a) ~all(cell2mat(reshape([a{:}],1,size(a{1},1),[]))==0,[2,3]), currD, 'UniformOutput',false)),2);
    smoothedPSTH = cellfun(@(cp) cellfun(@(t) t(mv,:),cp,'Uniformoutput',false),[currD{:}],'Uniformoutput',false);
    % condPhase = cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
    %     cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)'));
else
    ms_bins = findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins);
    if(~plotTrials)
        msBinPSTH = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-0),2,'omitnan')',n(unitInds),'UniformOutput',false))},tablePSTHD,'UniformOutput',false);
        dSegs = cellfun(@(v) {squeeze(cell2mat(reshape(cellfun(@(d) mean(d,1,'omitnan'),...
            v,'Uniformoutput',false),1,1,[])))'},vertcat(allSegs), 'UniformOutput',false);
    else
        msBinPSTH= cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-0),sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(tablePSTHD), 'UniformOutput',false);
        dSegs = cellfun(@(v) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(d',sTrials),...
            v,'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(allSegs), 'UniformOutput',false);
    end
    msBinPSTH =  cellfun(@(a)cellfun(@(u)max(0,u(unitInds,:)),a,'UniformOutput',false), msBinPSTH,'UniformOutput',false);
    if(plotTrials)
        mv = mv & sum(cell2mat(cellfun(@(m)mean(cell2mat(m),2,'omitnan').*1000>1,num2cell([msBinPSTH{:}],2),'UniformOutput',false)'),2)>sTrials/2;
    else
        mv = mv & mean(cell2mat(cellfun(@(n)cell2mat(n),msBinPSTH,'UniformOutput',false)),2,'omitnan').*1000>1;
    end
    msBinPSTH = cellfun(@(c) cellfun(@(s) s(mv,:), c, 'UniformOutput',false), msBinPSTH, 'UniformOutput',false);
    trialLength = floor(size(msBinPSTH{1}{1}, 2) / binWidth);
    smoothedPSTH = cell(1,length(msBinPSTH));
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
    normBaseline = cellfun(@(p,t)cellfun(@(a,n) [max(1,median(cell2mat(reshape(cellfun(@(c,s) ...
        permute(mean(c(:,findBins(s-5,params.bins(1:binWidth:end)):findBins(s-(1+rand(1)*3),params.bins(1:binWidth:end))),2,'omitnan').*binWidth,[1 3 2]),...
        squeeze(num2cell(a,[2])),num2cell([ones(1,all(isnan(n(mv,1))));n(~isnan(n(mv,1)),1)]),'UniformOutput',false),...
        [1,max(1,sum(~isnan(n(mv,1))))])),3,'omitnan'));],p,t,'UniformOutput',false),smoothedPSTH,dSegs,"UniformOutput",false);
    %smoothedPSTH = cellfun(@(c,n) cellfun(@(s,b) s./b', c,n,'UniformOutput', false),smoothedPSTH, normBaseline, 'UniformOutput',false);
end
somaLabs = string(somaTable(mv));
channels = allChannels(mv);
location = allLocations(mv,:);
somaReps = unique(somaLabs);
segInds = cellfun(@(s) fix(s(mv,~all(isnan(s),1))),cellfun(@(n) n(:,arrayfun(@(c)find(strcmp(maxSegL,c)),epochSegs)),cellfun(@(aa)cell2mat(cellfun(@(a) ...
    findBins(mean(a,1,'omitnan'),params.bins),aa(unitInds),'UniformOutput',false)),allSegs,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false)';
segVals = cellfun(@(i) findBins(params.bins(i),timeBins(1):1/(1000/binWidth):timeBins(end)),segInds,'UniformOutput',false);
%% plot single dimension PCA
saveDir = "S:\Lab\ngc14\Working\";
if(PCATime)
    savePre = saveDir + "PCA_Time\";
    taskPSTHD = cellfun(@(s) cellfun(@(t) t(:,unique(round(ms_bins./binWidth))'),s,'UniformOutput',false),smoothedPSTH, 'UniformOutput',false);
    pcaMat = cell2mat(cellfun(@cell2mat,taskPSTHD','UniformOutput',false));
else
    savePre = saveDir + "DataHigh\";
    taskPSTHD= cellfun(@(s) cellfun(@(t) t(:,unique(round(ms_bins./binWidth))),s,'UniformOutput',false)',smoothedPSTH, 'UniformOutput',false);
    pcaMat = zscore(cell2mat(cellfun(@cell2mat,taskPSTHD,'UniformOutput',false))');
    variableWeights = arrayfun(@(s) 1/(sum(strcmp(somaLabs,s))),somaLabs);
    rowWeights = ones(1,size(pcaMat,1));
end
[loadings,scores,eig,ts,exp] = pca(pcaMat,'Economy',false,'Centered',centered,'Algorithm','eig','VariableWeights',variableWeights,'Weights',rowWeights);somaProj={};
if(PCATime)
    scores = permute(reshape(scores,sum(mv),max(1,sTrials*plotTrials),length(taskPSTHD),size(loadings,1)),[1 4 2 3]);
    somaProj = cellfun(@(m) cellfun(@transpose,squeeze(num2cell(m,[1 2]))','UniformOutput',false),squeeze(num2cell(scores,[1:length(size(scores))-1]))','UniformOutput',false);
    somaProj = arrayfun(@(s)cellfun(@(c)cellfun(@transpose,permute(num2cell(cell2mat(reshape(cellfun(@(t)cell2mat(arrayfun(@(n) ...
        loadings(:,n).*mean(t(n,strcmp(somaLabs,s)).*binWidth,2,'omitnan')',1:num_dims,'UniformOutput',false)),c,'UniformOutput',false),1,1,[])),1),[2 3 1]),'UniformOutput',false),...
        somaProj,'UniformOutput',false),somaReps, 'UniformOutput',false);
    loadings = squeeze(mean(scores,[3,4]));
else
    mZ =mean(pcaMat)'; sZ=std(pcaMat)';
    taskPSTHD = cellfun(@(c) cellfun(@(t) bsxfun(@rdivide,bsxfun(@minus,t,mZ),sZ),c,'UniformOutput',false),taskPSTHD, 'UniformOutput',false);
    condScores= arrayfun(@(c) reshape(scores([1:(max(1,sTrials*plotTrials)*size(taskPSTHD{1}{1},2))]+(c*(max(1,sTrials*plotTrials)...
        *size(taskPSTHD{1}{1},2))),:),size(taskPSTHD{1}{1},2),[],size(scores,2)),0:length(conditions)-1,'UniformOutput',false);
    somaProj = projectData(somaLabs,taskPSTHD,num_dims,loadings);
    for s = 1:length(somaReps)
        for c = 1:length(conditions)
            for n = 1:num_dims
                %somaProj{s}{c}(n,:) =  cellfun(@(w)w'*mean(loadings(somaLabs==somaReps(s),n),1)',squeeze(num2cell(condScores{c}(:,:,n),1)),'UniformOutput',false)';
            end
        end
    end 
end
if(strcmpi(centered,'on'));savePre = savePre+"Centered";else;savePre = savePre+"Non-Centered";end
if(~PCATime);savePre=savePre+"\Z-Score";end
newSaveDir = savePre +"\"+type+"\";
if(~plotTrials & strcmpi(type,'traj'));newSaveDir = newSaveDir{1}(1:end-1)+"_AVG\";end
savePath = newSaveDir;
for p = 0:(double(~contains(type,"traj",'IgnoreCase',true)))
    if(~contains(type,"traj",'IgnoreCase',true));savePath = newSaveDir+phases(p+1)+"\";end
%    plotProj(loadings,exp,segVals,cellfun(@(s) s([1:length(conditions)]+(length(conditions)*p)),somaProj,'UniformOutput',false),...
%        somaLabs,num_dims,conditions,timeBins,binWidth,false,savePath);
end
%% DATAHIGH Dim Reduce
somaColors = {[0 1 .2],[.8 0 1]};%
if(PCATime)
    projectUnits = arrayfun(@(r) cellfun(@(t) cellfun(@(s) s(strcmp(somaLabs,r),:)*loadings,t,'UniformOutput',false),taskPSTHD, 'UniformOutput',false),somaReps, 'UniformOutput',false);
    projectUnits = cellfun(@(a) cellfun(@(d) d(:,:)',vertcat(a{:}),'Uniformoutput',false), projectUnits, 'UniformOutput',false);
else
    projectUnits = arrayfun(@(r)  squeeze(cellfun(@squeeze,num2cell(num2cell(reshape((pcaMat(:,contains(somaLabs,somaReps(r)))*...
        loadings(contains(somaLabs,somaReps(r)),:))',size(pcaMat,2),size(taskPSTHD{1}{1},2),max(1,plotTrials*sTrials),[]),[1 2]),3),'UniformOutput',false))',1:2,'UniformOutput',false);
    if(strcmpi(type,'state'))
        projectUnits = cellfun(@(s) squeeze(num2cell(reshape(cellfun(@(t) cell2mat(t'),s','UniformOutput',false),[],length(conditions),length(phases)),1)), projectUnits, 'UniformOutput',false);
        projectUnits = cellfun(@(s) reshape(cellfun(@(t) cell2mat(t'),s,'UniformOutput',false),[],1),projectUnits, 'UniformOutput',false);
    else
        projectUnits = cellfun(@(s) vertcat(s{:}), projectUnits,'Uniformoutput',false);
        projectUnits = cellfun(@(s) cellfun(@(t,n) t(:,1:end-((n==length(conditions)*1)*2.5*binWidth)),s,num2cell(1:length(s))', 'UniformOutput',false), projectUnits, 'UniformOutput',false);
        projectUnits{end+1} =cellfun(@(s) vertcat(s{:}), squeeze(cellfun(@squeeze,num2cell(num2cell(reshape((pcaMat(:,contains(somaLabs,somaReps))*...
            loadings(contains(somaLabs,somaReps),:))',size(pcaMat,2),size(taskPSTHD{1}{1},2),max(1,plotTrials*sTrials),[]),[1 2]),3),'Uniformoutput',false)), 'UniformOutput',false);
    end
end
cls = cellfun(@(r) repmat({r},max(strcmpi(type,'traj')*(plotTrials*sTrials),1),1),somaColors','UniformOutput',false);
cls = reshape(cellfun(@(v) repmat({v},length(conditions),1),vertcat(cls{:}),'UniformOutput',false),[],1);
cls{end+1} = colors.values';
cTrials = arrayfun(@(d) {repmat(d,max(strcmpi(type,'traj')*plotTrials*sTrials,1),1)},conditions,'UniformOutput',false)';
cTrials = arrayfun(@(p) cellfun(@(c) cell2mat(cellfun(@(t) string(t)+"-"+string(p),c,'UniformOutput',false)),cTrials,'UniformOutput',false),somaReps,'UniformOutput',false);
epochStarts = repmat(cellfun(@(e) e(:,all(e<size(taskPSTHD{1}{1},2),1)),num2cell(cell2mat(cellfun(@(s) repmat(round(mean((s-min(ms_bins))./binWidth,1,'omitnan')),...
    max(1,strcmpi(type,'traj')*plotTrials*sTrials),1),segInds,'UniformOutput',false)),2),'UniformOutput',false),sum(cellfun(@length,cls))/length(conditions),1);
%%
if(strcmpi(type,'traj'))
    dHiStruct = struct('data',vertcat(projectUnits{:}),'condition',cellstr(cell2mat(vertcat(cTrials{:}))),'epochStarts',epochStarts,...
        'epochColors',cellfun(@(r,e)cell2mat(arrayfun(@(n) [hsv2rgb([r(1)*(1*n~=4),min(1,.5+(.2*n))*(1*n~=4),max(1*n==4,1.25-(.25*n))])],1:length(e),"UniformOutput",false)'),...
        num2cell(rgb2hsv(cell2mat(cellfun(@cell2mat,cls,'UniformOutput',false))),2),epochStarts,'UniformOutput',false));
else
    cTrials = arrayfun(@(r)cellfun(@(c) c+"-"+r,cTrials,'UniformOutput',false),phases,'UniformOutput',false);
    dHiStruct = struct('data',vertcat(projectUnits{:}),'condition',cellstr(vertcat(cTrials{:})));
    [dHiStruct.type]=deal('state');
    [dHiStruct.epochStarts]=deal(1);
    [dHiStruct(cellfun(@(s) contains(s,'Arm'),{dHiStruct.condition})).epochColors] = deal(rgb2hsb(somaColors{1}));
    [dHiStruct(cellfun(@(s) contains(s,'Hand'),{dHiStruct.condition})).epochColors] = deal(rgb2hsb(somaColors{end}));
    for a = 1:length(dHiStruct)
        dHiStruct(a).epochColors = hsb2rgb(max(min(dHiStruct(a).epochColors+([0 .5 .0]*2*(contains(dHiStruct(a).condition,"Reach")-.5)),1),0));
    end
end
DataHigh(dHiStruct);
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