muscles = cellstr({'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor',...
    'Digit Extensor','Digit Flexor'});
groupings = {string(muscles(1:3)),string(muscles(4:end))};
smoothWin = .15;
num_dims = 4;
saveFig = false;
plotTrials = 1;
sTrials = 50;
savePath = "S:\Lab\ngc14\Working\DataHi\EMG\";
load("C:\Users\ngc14\Desktop\EMG_MaxNorm.mat");
colors =containers.Map(Conditions,{[.7 0 0],[1 .65 0 ],[0 0 .75]}');
%%
trialLength = cellfun(@(tc) cellfun(@(ti,cs) max(cellfun(@(t) t(strcmp(cs,"StartHold"))-(t(strcmp(cs,"GoSignal"))-.5*Fs),ti))+(smoothWin*Fs),tc,ConditionSegs),sessionSegs, 'UniformOutput',false);
psthGrouped = cellfun(@(mc,tc,tl) cellfun(@(sc,ti) downsampleTrials(cell2mat(cellfun(@(m,t) resize(conv(m(:,t(1)-(.5*Fs):t(end-2)),...
    gausswin(smoothWin*Fs)./sum(gausswin(smoothWin*Fs)),'valid'),[size(m,1),max(tl)],'Pattern','edge','side','both'),...
    sc,ti,'UniformOutput',false)')',sTrials)',mc,tc,'UniformOutput',false),sessionSigs,sessionSegs,trialLength,'UniformOutput',false);
%%
dHiStruct = cellfun(@(c) squeeze(num2cell(permute(cell2mat(reshape(cellfun(@(r) r(:,1:length(-.5:(1/1000):1.5)*(Fs/1000)),c,'UniformOutput',false),...
    1,1,[])),[3 2 1]),[1 2])),num2cell([psthGrouped{:}],2),'UniformOutput',false);
[loadings, scores, eig,~,exp]=pca(cell2mat(vertcat(dHiStruct{:})')','Centered','on','NumComponents',num_dims);
somaProj = cellfun(@(s)cellfun(@(c) cellfun(@(t) cell2mat(arrayfun(@(n)t.*loadings(strcmp(muscles,s),n),...
    1:num_dims,'UniformOutput',false)'),c,'UniformOutput',false),psthGrouped,'UniformOutput',false),muscles,'UniformOutput',false);
%%
lc = {'-','-','-',':',':',':',':'};
ax = {};
for n = 1:3
    figure();
    ax = {};plotSegs={};
    tl=tiledlayout(n+1,num_dims+1,'TileIndexing','rowmajor');
    nexttile([1,1]); hold on; title("Variance Explained");plot(cumsum(exp),'LineWidth',2);ylim([0 100]);xlim([0.5 num_dims+.5]);
    plot([0 length(muscles)],[90 90],'r--','LineWidth',1);
    loadingsSplit = cellfun(@(s) loadings(ismember(string(muscles),s),:),groupings,'UniformOutput',false)';
    nexttile([1 2]); hold on; title("Loadings"); imagesc(cell2mat(loadingsSplit)); axis ij; axis tight; colormap('turbo'); clim([-.1 .1]);
    line([0+.5 num_dims+.5],repmat(3,1,2),'Color','k','LineWidth',2,'LineStyle',':');
    colorbar;
    bx=nexttile([1 2]); hold on; title("Loadings Distributions");
    boxplotGroup(bx,loadingsSplit','PrimaryLabels',repmat({''},num_dims*length(groupings),1),'Symbol','',...
        'SecondaryLabels',arrayfun(@(s) "Factor "+num2str(s),1:num_dims),'Notch','on'); ylim([-1 1]);
    ax{end+1} = nexttile(tl,[1,1]); hold on; title("PSTHS");
    plotSegs{end+1} = 1:length(Conditions);
    cellfun(@(s,l) cellfun(@(m,p) plot(mean(m,1,'omitnan')','LineStyle',l,'LineWidth',1,'Color',p), s',colors.values), psthGrouped,lc);
    for c = 1:length(Conditions)
        for s =1:length(muscles)
            co = repmat(colors(Conditions(c)),num_dims,1);
            weightedPSTHS = cell2mat(somaProj{s}{c})';%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
            weightedPSTHS = permute(reshape(weightedPSTHS',size(weightedPSTHS,2)/(max(1,plotTrials*sTrials)),[],size(weightedPSTHS,1)),[3 2 1]);
            if(n==1)
                ax{end+1}=nexttile();hold on;title("Cond " + num2str(c));
                plotSegs{end+1} = c;
                cl = rgb2hsv(colors(Conditions(c)));
                co = hsv2rgb([linspace(cl(1),cl(1),num_dims);linspace(1,.5,num_dims);linspace(.3,1,num_dims)]');
            elseif(n==2)
            elseif(n==3)
                weightedPSTHS = permute(cell2mat(cellfun(@(p) cell2mat(reshape(p{c}',1,1,[])), somaProj(s), 'UniformOutput',false)),[2 3 1]);
            end
            for d = 1:size(weightedPSTHS,3)
                if(n>=2)
                    if(n==2)
                        ax{end+1}=nexttile(((s)*tl.GridSize(end))+d);
                        plotSegs{end+1} = 1:length(Conditions);
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
        cellfun(@(pc) pc(:,contains(epochSegs,["Hold","Withdraw"])),condSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(Conditions(p)))),ax,plotSegs);
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


function arr = downsampleTrials(r,sTrials)
sz = size(r,2)-mod(size(r,2),2);
trials = double(mod(sz,sTrials)~=0)*(sTrials-mod(sz,sTrials));
nPairs = round(sz/sTrials)*(sTrials-trials);
arr=[mean(reshape(r(:,max(1,sz-trials+double(trials~=0)):sz-double(trials==0)),size(r,1),trials,[]),3,'omitnan'),...
    mean(reshape(r(:,1:nPairs),size(r,1),[],round(sz/sTrials)),3,'omitnan')];
end