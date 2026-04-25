function plotProj(loadings,exp,epochSegs,somaProj,somaLabs,location,num_dims,conditions,timeBins,saveFig,savePath)
colors =containers.Map(conditions,{[.7 0 0],[1 .65 0 ],[0 0 .75] }');% regexp(model,'[A-Z]+[^A-Z]+','match')
pg = [0 .8 .4;.4 0 .5; 0.5 0.5 0.5];
somaReps = unique(somaLabs);
conditions = colors.keys;

figure(); tl1=tiledlayout(2,1);
nexttile(tl1,[1,1]); hold on; title("Variance"); yyaxis right; plot(cumsum(exp),'LineWidth',2);ylim([0 100]);
plot([0 length(somaLabs)],[90 90],'r--','LineWidth',1);xlim([0,15]);
yyaxis left; bar(exp); ylim([0 100])

bx=nexttile(tl1,[1 1]);
bg=boxplotGroup(bx,arrayfun(@(s) loadings(somaLabs==s,1:2*num_dims),somaReps,'UniformOutput',false)','PrimaryLabels',...
    repmat({''},2*num_dims*length(somaReps),1),'Symbol','','SecondaryLabels',arrayfun(@num2str,1:2*num_dims,'Uniformoutput',false),...
    'Notch','on','colors',pg(1:length(somaReps),:),'GroupType','BetweenGroups');
hold on; title("Loadings Distributions");ylim([min(loadings(:,1:2*num_dims),[],'all'), max(loadings(:,1:2*num_dims),[],'all')]); xlabel("Factor");
if(0)
    nexttile(tl1,[1,2]); hold on;
    [locRC,sortR] = sort(location(:,1).*ImagingParameters.px2mm);
    [locML,sortM] = sort(location(:,2).*ImagingParameters.px2mm);
    locRC = min(locRC,4);
    ct=nexttile(tl1,[num_dims 1]); hold on; title("Loadings R/C"); imagesc(cell2mat(arrayfun(@(s) loadings(sortR(somaLabs(sortR)==s),:),somaReps,'UniformOutput',false)));
    axis ij;axis tight;colormap(ct,'jet');clim([-.1 .1]);
    binGroups =arrayfun(@(s) cell2mat(arrayfun(@(f) find(fix(locRC(somaLabs==s))==f,1),1:length(unique(fix(locRC)))-1,'UniformOutput',false)')',somaReps,'UniformOutput',false);
    arrayfun(@(i,ls) line([0+.5 num_dims+.5],repmat(i,1,2),'Color','k','LineWidth',2,'LineStyle',ls),...binGroups{1},repmat({":"},1,length(binGroups{1})));%...
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
    arrayfun(@(i,ls) line([0+.5 num_dims+.5],repmat(i,1,2),'Color','k','LineWidth',2,'LineStyle',ls),...binGroups{1},repmat({":"},1,length(binGroups{1})))%...
        [sum(somaLabs=="Arm"),binGroups{1},sum(somaLabs=="Arm")+binGroups{end}],["--",repmat(":",1,length(cell2mat(binGroups')))]);colorbar;
    for d = 1:num_dims
        bx=nexttile(tl1,tilenum(tl1,1+d,4),[1 1]); hold on; title("Binned M/L Loadings "+num2str(d));
        boxplotGroup(bx,arrayfun(@(s) cell2mat(arrayfun(@(n) resize(loadings(fix(locML(somaLabs==s))==n,d),[sum(somaLabs==s),1],'FillValue',NaN),...
            unique(fix(locML)),'UniformOutput',false)'),somaReps,'UniformOutput',false)',...
            'Symbol','','SecondaryLabels',arrayfun(@num2str,unique(fix(locML)),'Uniformoutput',false),'Notch','on','PrimaryLabels',repmat({''},length(somaReps)*length(unique(fix(locML))),1));
        ylim([-.05 .05]+(d==1*[.05 .05]));
    end
end
if(saveFig)
    saveFigures(gcf,savePath,"Variance+Loadings",[]);
end
if(0)
    tl1 = tiledlayout(4,2);
    [~,rankedWeight] = arrayfun(@(s) sort(loadings(somaLabs==s,:),1,'ascend'), somaReps, 'UniformOutput',false);
    sortR = cell2mat(cellfun(@(r,s) cell2mat(cellfun(@(c,t) t(c),num2cell(r,1),repmat({sortR(somaLabs==s)},[1,size(r,2)]),'UniformOutput',false)), rankedWeight,num2cell(somaReps),'UniformOutput',false));
    sortM = cell2mat(cellfun(@(r,s) cell2mat(cellfun(@(c,t) t(c),num2cell(r,1),repmat({sortM(somaLabs==s)},[1,size(r,2)]),'UniformOutput',false)), rankedWeight,num2cell(somaReps),'UniformOutput',false));
    rankedWeight = cell2mat(cellfun(@(r) (r-1)./(size(r,1)-1),rankedWeight,'UniformOutput',false));
    rcLoadings = rankedWeight(sortR);
    ct=nexttile(tl1,[4 1]); hold on; title("Ranked Loadings R/C"); imagesc(rcLoadings);axis ij;axis tight;colormap(ct,'jet');clim([0 1]);
    line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
    mlLoadings = rankedWeight(sortM);
    ct=nexttile(tl1,[4 1]); hold on; title("Ranked Loadings M/L"); imagesc(mlLoadings);axis ij;axis tight;colormap(ct,'jet');clim([0 1]);
    line([0+.5 num_dims+.5],repmat(sum(somaLabs=="Arm"),1,2),'Color','w','LineWidth',1.5,'LineStyle','--');colorbar;
    if(saveFig)
        saveFigures(gcf,savePath,"RankedLoadings",[]);
    end
end
%%
lg={};
plotTraj=size(somaProj{1}{1}{1},2)~=1;
for n = 1:3
    figure();
    lc = {'-','-'};ax = {};plotSegs={};
    if(n==1)
        tileorder = 'rowmajor';
    else
        tileorder = 'columnmajor';
    end
    tl=tiledlayout(length(conditions)*(n==1)+((n>1)*num_dims),max(2,n),'TileIndexing',tileorder,'TileSpacing','none','Padding','tight');
    for c = 1:length(conditions)
        for s =1:length(somaReps)+double(n==3)
            weightedPSTHS = cell2mat(somaProj{s}{c});%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
            weightedPSTHS = permute(reshape(weightedPSTHS,num_dims,[],size(weightedPSTHS,2)),[3 2 1]);
            if(n==1)
                ax{end+1}=nexttile();hold on;if(c==1);title(conditions(c)+"- " + string(somaReps(s)));end;
                plotSegs{end+1} = c;
                cl = rgb2hsv(colors(conditions{c}));
                co = hsv2rgb([linspace(cl(1),cl(1),num_dims);linspace(1,.5,num_dims);linspace(.3,1,num_dims)]');
                avgDim{s,c} = squeeze(mean(weightedPSTHS,2,'omitnan'));
            else
                co = repmat(colors(conditions{c}),num_dims,1);
                if(n==2)
                elseif(n==3)
                    lc = {'-','-','-'};
                    co = repmat(pg(s,:,:),num_dims,1);
                    % if(s==length(somaReps))
                    %     co = rgb2hsv(colors(conditions{c}));
                    %     co = repmat(hsv2rgb(co(1), 1, .5),num_dims,1);
                    % end
                    weightedPSTHS = permute(cell2mat(cellfun(@(p) cell2mat(reshape(p{c}',1,1,[])), somaProj(s), 'UniformOutput',false)),[2 3 1]);
                end
            end
            for d = 1:size(weightedPSTHS,3)
                if(n>=2)
                    if(n==2)
                        ax{end+1}=nexttile(((s-1)*tl.GridSize(1))+d);
                        plotSegs{end+1} = 1:length(conditions);
                        titleName = "Dim " + num2str(d) + " - " + string(somaReps(s));
                    else
                        ax{end+1}=nexttile((c-1)*tl.GridSize(1)+d);
                        plotSegs{end+1} = c;
                        titleName = "Dim "+ num2str(d);
                    end
                    hold on;if(d==1);title(titleName);end;
                    if(n==3 & s==1 & d==1 & c==1)
                        lg = cellfun(@(rg) plot([NaN,NaN],[NaN,NaN],'Color',rg),num2cell(pg,2),'UniformOutput',false);
                        legend([lg{:}],[somaReps;"All"],'AutoUpdate','off');
                    end
                end
                if(plotTraj)
                    % p=plot(mean(meanP,2,'omitnan'),'LineWidth',1.5*double(size(squeeze(weightedPSTHS(:,:,d)),2)==1)+.5,'LineStyle',lc{s});
                    % arrayfun(@(pc) set(pc,'Color',[co(d,:),.9-(.35*(double(n==3)*(s-1)))]),p);
                    meanP = squeeze(weightedPSTHS(:,:,d));
                    shadedErrorBar(1:size(meanP,1),meanP',std(meanP,0,2),'lineProps',...
                        {'Color',[co(d,:),.95-(.35*double(n==3)*(s-1))],'LineWidth',1.5*double(size(meanP,2)==1)+.5,'LineStyle',lc{s}});
                else
                    [bins centers] = hist(squeeze(weightedPSTHS(:,:,d)),linspace(...
                        min(weightedPSTHS(:,:,d),[],'all'),...
                        max(weightedPSTHS(:,:,d),[],'all'),10));
                    bins = bins ./ sum(bins);
                    bar(centers, bins, 'FaceColor',co(d,:));ylim([0 .5]);
                end
            end
        end
    end
    if(plotTraj)
        [~,tn,~]= unique(cellfun(@(a) tilenum(a), ax)); ex = ax(tn);
        cellfun(@(a) set(a,'XTick',0:50:size(weightedPSTHS,1)),ax);
        cellfun(@(a) set(a,'XTickLabels',linspace(timeBins(1),timeBins(end),length(get(a,'XTick')))),ax);
        cellfun(@(a) set(a,'XTickLabels',[]),ex(setdiff(1:prod(tl.GridSize),tilenum(tl,tl.GridSize(1),1:tl.GridSize(end)))));
        cellfun(@(a) set(a,'XLim',[0 size(weightedPSTHS,1)]),ax);
        cellfun(@(a) set(a,'YAxisLocation','right'),ex(tilenum(tl,1:tl.GridSize(1),tl.GridSize(end))));
        %cellfun(@(a,p) cellfun(@(x,s) plot(a,[x;x], repmat(get(a,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
        %    cellfun(@(pc) round(mean(pc(:,end-1:end),1,'omitnan')),epochSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(conditions(p)))),ax,plotSegs);
        condSegs = round(mean(cell2mat(epochSegs),1,'omitnan'));
        cellfun(@(a) arrayfun(@(x) plot(a,[x,x], get(a,'YLim'),'k--','LineWidth', 1.5),condSegs(1:end-1)),ax);
        cellfun(@(a,s) set(a,'DataAspectRatio',[1 (1+tl.GridSize(1)/tl.GridSize(end))*(1/(2*s{1}(1))),s{1}(3)],'PlotBoxAspectRatio',[1 2.5*s{2}(2) 1]),ax,...
            cellfun(@(a) get(a,{'DataAspectRatio','PlotBoxAspectRatio'}),ax,'UniformOutput',false));
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
        saveFigures(gcf,savePath,fileNameSave,[]);
        save(savePath+"AvgDimTraj",'avgDim');
    end
end
end