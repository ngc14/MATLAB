conditions = ["Extra Small Sphere","Large Sphere","Photocell"];
params = PhysRecording(string(conditions),.001,.001,-1,3,containers.Map(conditions,{"StartReach","StartReach","StartReach"}));
tPhys = unitTable(conditions,params);
%%
model = "GilliganSkipper_ArmHand";
type = 'Traj';
saveDir = "S:\Lab\ngc14\Working\DataHi\Combined\";
saveFig = true;
num_dims=4;
sTrials = 30;
plotTrials = 0;
timeBins = [-.5, 2.5];
phases = {"StartReach","StartHold"};
phaseWindows = {[-100 100], [-200 0]};
dimCond = reshape(regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')+"_"+params.condAbbrev.values',1,[]);%"-"+regexp(s,'(?<=Start)\w*','match'),phases,'UniformOutput',false)),1,[]);
splitGroup = "Condition";
if(~plotTrials & strcmp(type,"Traj"))
    type = type+"_Avg";
end
savePath = saveDir+type+"\";
if(~exist(savePath,'dir')), mkdir(savePath); end
colors = containers.Map(dimCond,{[.7 0 0],[1 .65 0 ],[0 0 .75],[1 0 .3],[1 1 0],[0 .6 1]}');
if(~contains(type,'Traj'))
    switch(splitGroup)
        case "Somatotopy"
            condInds = arrayfun(@(c) contains(c,'Arm'), dimCond);
            for u = 1:length(condInds)
                if(condInds(u)==0),colors{u} = [.8 .8 .8];
                else,colors{u} = [.2 .2 .2];end
            end
        case  "Condition"
            condInds = arrayfun(@(c) contains(c,'S'), dimCond);
            condInds = condInds + arrayfun(@(c) contains(c,'_E'), dimCond);
            for u = 1:length(condInds)
                if(condInds(u)==0),colors{u} = [0 0 1];
                elseif(condInds(u)==1),colors{u} = [1 .85 0];
                else,colors{u} = [1 0 0];end
            end
        case "Phase"
            condInds = arrayfun(@(c) contains(c,'Reach'), dimConds);
            for u = 1:length(condInds)
                if(condInds(u)==0),colors{u} = [1 0 1];
                else,colors{u} = [0 1 1];end
            end
        case "Monkey"
            condInds = arrayfun(@(c) contains(c,'Gilligan'), dimConds);
            for u = 1:length(condInds)
                if(condInds(u)==0),colors{u} = [.8 .4 0];
                else,colors{u} = [0 .5 0];end
            end
    end
    colors =containers.Map(dimCond,colors);
end
tPhysTable = tPhys(contains(string(tPhys.Monkey),[regexp(extractBefore(model,"_"),'[A-Z]+[^A-Z]+','match')]) & contains(string(tPhys.Somatotopy),...
    [regexp(extractAfter(model,"_"),'[A-Z]+[^A-Z]+','match')]),:);
if(contains(type,'Traj'))
    allSegs= arrayfun(@(s) tPhysTable{contains(string(tPhysTable.Somatotopy),(extractBefore(string(s),"_"))),contains(tPhysTable.Properties.VariableNames,"Segs_"+extractAfter(s,"_"))}, dimCond, 'UniformOutput',false);
    allSegs= cellfun(@(c) mean(cell2mat(c),1,'omitnan'),allSegs, 'UniformOutput',false);
    taskPSTHD= arrayfun(@(a) tPhysTable{contains(string(tPhysTable.Somatotopy),(extractBefore(string(a),"_"))),contains(tPhysTable.Properties.VariableNames,"PSTH_"+extractAfter(a,"_"))},dimCond,'UniformOutput',false);
    avgTrace = mean(cell2mat(vertcat(taskPSTHD{:})'),2,'omitnan');
    avgTrace = zeros(size(avgTrace,1),size(avgTrace,2));
    if(~plotTrials)
        taskPSTHD = cellfun(@(n) {cell2mat(cellfun(@(m) mean(max(0,m-avgTrace),2,'omitnan')',n,'UniformOutput',false))}, taskPSTHD, 'UniformOutput',false);
    else
        taskPSTHD= cellfun(@(a) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(d) downsampleTrials(max(0,d-avgTrace),sTrials),...
            a(cellfun(@(s)size(s,2)>=sTrials,a)),'Uniformoutput',false),1,1,[])),[3 1 2]),[1,2])),vertcat(taskPSTHD), 'UniformOutput',false);
    end
    numUnits = arrayfun(@(s) min(cellfun(@(m) size(m{1},1),taskPSTHD(contains(dimCond,s)))), unique(arrayfun(@(t) extractBefore(t,"_"),dimCond)));
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(dimCond,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    taskPSTHD = cellfun(@(u,i) cellfun(@(n) n(i,findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins)),u,'UniformOutput',false),...
        taskPSTHD, [unitInds{:}],'UniformOutput',false);
    cls = cellfun(@(r) repmat({r},max(plotTrials*sTrials,1),1),cellfun(@hsv2rgb,cellfun(@(l) flipud([linspace(l(1),l(1),5);...
        linspace(1,.25,5);linspace(.85,1,5)]'),cellfun(@rgb2hsv,colors.values','UniformOutput',false),'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
    dHiStruct = struct('data',vertcat(taskPSTHD{:}),'epochStarts',reshape(cellfun(@(n) [1,n([2,3,6]),min(n(7),length(timeBins(1):params.binSize:timeBins(end))-1)], ...
        repmat(cellfun(@(a) findBins(a,params.bins(findBins(timeBins(1),params.bins):findBins(timeBins(end),params.bins))),allSegs,'UniformOutput',false),...
        max(plotTrials*sTrials,1),1),'UniformOutput',false),[],1),'condition',cellstr(cell2mat(cellfun(@(d) repmat(string(d),max(plotTrials*sTrials,1),1),cellstr(dimCond),'UniformOutput',false)')),...
        'epochColors',vertcat(cls{:}));
else
    for p = 1:length(phases)
        phaseConds = cellfun(@(t) find(strcmp(phases{p},t)), params.condSegMap.values(params.condSegMap.keys),'UniformOutput',false);
        trialFR = cellfun(@(ct,cs,ta,tw) cellfun(@(a,b) cellfun(@(m,tt) m(max(1,tt+tw(1)):max(range(tw)+1,tt+tw(end))),...,
            num2cell(a,1)',arrayfun(@(bb) [find(isalmost(params.bins,bb,params.binSize/1.99),1),NaN(isnan(bb),1)],b(:,ta),'UniformOutput',false),...
            'UniformOutput',false)',ct,cs,'UniformOutput',false),num2cell(tPhysTable{:,contains(tPhysTable.Properties.VariableNames,"PSTH_")},1),...
            num2cell(tPhysTable{:,contains(tPhysTable.Properties.VariableNames,"Segs_")},1),phaseConds,repmat({phaseWindows{p}},1,length(phaseConds)),'UniformOutput',false);
        trialFRMat{p} = cellfun(@(m) cat(2,m{~cellfun(@isempty,m)}), [trialFR{:}], 'UniformOutput',false);
    end
    currD = cellfun(@(m,n)cellfun(@(c)squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) ...
        downsampleTrials(r,sTrials),c,'UniformOutput',false),1,1,[])),[3 1 2]),[1 2])),arrayfun(@(t) ...
        m(n(:,contains(params.condAbbrev.values,extractAfter(t,"_")))>=sTrials & contains(string(tPhysTable.Somatotopy),...
        extractBefore(t,"_")),contains(params.condAbbrev.values,extractAfter(t,"_"))),dimCond,'UniformOutput',false),'UniformOutput',false),...
        trialFRMat,cellfun(@(p) cellfun(@(s) size(s,2),p), trialFRMat,'UniformOutput',false),'UniformOutput',false);
    numUnits = arrayfun(@(s) min(cellfun(@(m) size(m{1},1) ,[currD{end}(contains(dimCond,s))])), unique(arrayfun(@(t) extractBefore(t,"_"),dimCond)));
    unitInds = arrayfun(@(u) repmat({randi(u,min(numUnits),1)},1,size(dimCond,2)/length(numUnits)), numUnits,'UniformOutput',false);%repmat({}',1,size(currD,2)/(length(numUnits)));
    currD = cellfun(@(u,i) cellfun(@(n) n(i,:),u,'UniformOutput',false), [currD{:}], repmat([unitInds{:}],1,length(phases)),'UniformOutput',false);
    dHiStruct = struct('data',vertcat(currD{:}),'condition',cellstr(cell2mat(cellfun(@(r) repmat(string(r),size(currD{1},1),1),...
        cell2mat(cellfun(@(s) dimCond+"-"+extractAfter(s,'Start'),phases,'UniformOutput',false)),'UniformOutput',false)')),...
        'epochStarts',1,'epochColors',{[0 0 0]});
    for i = 1:length(dHiStruct)
        dHiStruct(i).epochColors = cell2mat(colors.values({extractBefore(dHiStruct(i).condition,"-")}));
    end
end
DataHigh(dHiStruct,'DimReduce');
if(saveFig)
    save(savePath+"DStruct_"+model+".mat",'dHiStruct','-v7.3');
end
%%
all_h = findall(groot,'Type','Figure');
D = guidata(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));
handles = guihandles(all_h(arrayfun(@(s) strcmp(s.Name,'DataHigh'),all_h)));
D = D.D;
conds = unique({D.condition});
if(length(D)<sTrials)
    plotTrials = 0;
end
figure(); tax=tiledlayout(1,max(1,num_dims/2)*2);
ylimT = [min(arrayfun(@(m) min(m.data,[],'all'),D)),max(arrayfun(@(m) max(m.data,[],'all'),D))];
for icond = 1:length(conds)
    for idim = 1:num_dims
        if(contains(type,'Traj'))
            itrial = find(contains({D.condition}, extractAfter(conds{icond},"_")));
            epochs = [round(mean(vertcat(D(itrial).epochStarts),1,'omitnan')),size(D(1).data,2)];
        else
            epochs = [find(ismember({D.condition}, conds)),NaN];
        end
        nexttile(idim); hold on; title(idim); ylim(ylimT);
        dTrial = cat(3,(D(icond).data));
        for iepoch = 1:length(epochs)-1
            if(contains(type,'Traj'))
                indices = epochs(iepoch):(epochs(iepoch+1));
                cellfun(@(p) plot(indices,p,'Color', cell2mat(colors.values(cellstr(dimCond(icond)))),'LineWidth',2),...
                    num2cell(squeeze(mean(dTrial(idim,indices,:),max(~plotTrials*length(size(dTrial))+1,plotTrials),'omitnan')),1+~plotTrials));
                if(iepoch>1 && icond<=length(conds)/2)
                    if(iepoch==length(epochs)-1)
                        line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",':','Color',cell2mat(colors.values(cellstr(dimCond(icond))))./1.5,'LineWidth',2);
                    else
                        if(icond==1)
                            line([epochs(iepoch),epochs(iepoch)],ylimT,"LineStyle",'--','Color','k');
                        end
                    end
                end
            else
                [bins centers] = hist(D(iepoch).data(idim,:));
                bins = bins ./ sum(bins);
                bar(centers, bins, 'FaceColor',cell2mat(colors.values(cellstr(dimCond(icond)))));
            end
        end
    end
    if(icond==length(conds))
        l = arrayfun(@(d) plot(NaN(length(dimCond),1),'LineWidth',3,'Color',cell2mat(colors.values(cellstr(d)))),dimCond);
        legend(l,dimCond,'Autoupdate','off','FontSize',14,'Orientation','horizontal');
    end
end
if(saveFig)
    saveFigures(gcf,savePath,model+"_"+splitGroup,[]);
end
if(contains(type,'Traj'))
    figure(); tax=tiledlayout(1,length(conds));
    for icond = 1:length(conds)
        legendColors= {};
        currColor = rgb2hsv(cell2mat(colors.values(cellstr(dimCond(icond)))));%([1 0 0]);
        currColor(2) = 1;
        currColor(end) = .4;
        for idim = 1:num_dims
            if(idim<=num_dims/2)
                currColor(end) = min(1,currColor(end) + (.4*(idim-1)));
            else
                currColor(end) = 1;
                currColor(2) = max(.25,currColor(2) - (.25*(idim-(num_dims/2))));
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
end

function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = uint8(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+uint8(trials~=0)):sz-uint8(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end