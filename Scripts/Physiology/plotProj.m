function plotProj(loadings,exp,epochSegs,somaProj,somaLabs,num_dims,conditions,timeBins,binWidth,saveFig,savePath)
colors = {[.7 0 0],[1 .65 0 ],[0 0 .75] };
colors =containers.Map(conditions,colors(1:length(conditions)));
somaReps = unique(somaLabs);
pg = [0 .8 .4;.4 0 .5];
if(length(somaReps)<length(somaProj))
    pg(end+1,:) = [0.25 0.25 0.25];
end
conditions = colors.keys;
plotTraj=size(somaProj{1}{1}{1},2)~=1;

figure(); tl1=tiledlayout(2,1);
nexttile(tl1,[1,1]); hold on; title("Variance"); yyaxis right; plot(cumsum(exp),'LineWidth',2);ylim([0 100]);
plot([0 length(somaLabs)],[90 90],'r--','LineWidth',1);xlim([0,15]);
yyaxis left; bar(exp); ylim([0 100])

bx=nexttile(tl1,[1 1]);
weightedLoadings = arrayfun(@(s) loadings(somaLabs==s,1:2*num_dims).^2/sum(somaLabs==s),somaReps,'UniformOutput',false)';
bg = bar(cell2mat(cellfun(@(r) sum(r./sum(cell2mat(weightedLoadings'))),weightedLoadings,'UniformOutput',false)')');
bg(1).FaceColor=pg(1,:); bg(2).FaceColor=pg(2,:);
%bg=boxplotGroup(bx,arrayfun(@(s) loadings(somaLabs==s,1:2*num_dims),somaReps,'UniformOutput',false)','PrimaryLabels',...
%    repmat({''},2*num_dims*length(somaReps),1),'Symbol','','SecondaryLabels',arrayfun(@num2str,1:2*num_dims,'Uniformoutput',false),'Notch','on','colors',pg(1:length(somaReps),:),'GroupType','BetweenGroups');
hold on; title("Loadings Distributions");ylim([0 1]); xlabel("Factor");
if(saveFig)
    saveFigures(gcf,savePath,"Variance+Loadings",[]);
end
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
%%c
lg={};
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
        for s =1:length(somaReps)+(n==3 & length(somaProj)>length(somaReps))
            weightedPSTHS = cell2mat(somaProj{s}{c});%(pcaMatrix.*loadings(:,n)').*(condSomaInd./condSomaInd),2,'omitnan');
            weightedPSTHS = permute(reshape(weightedPSTHS,num_dims,[],size(somaProj{s}{c},2)),[3 2 1]);
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
                    weightedPSTHS = permute(cell2mat(cellfun(@(p) squeeze(cell2mat(reshape(p{c},num_dims,1,1,[]))), somaProj(s), 'UniformOutput',false)),[3 2 1]);
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
                    hold on;if(d==1);title(titleName);end
                    if(n==3 && s==1 && d==1 && c==1)
                        lg = cellfun(@(rg) plot([NaN,NaN],[NaN,NaN],'Color',rg),num2cell(pg,2),'UniformOutput',false);
                        if(length(somaReps)<length(lg))
                            legend([lg{:}], [somaReps;"All"],'AutoUpdate','off','Location','eastoutside');
                        else
                            legend([lg{:}],somaReps,'AutoUpdate','off','Location','eastoutside');
                        end
                    end
                end
                if(plotTraj)
                    shadedErrorBar(1:size(weightedPSTHS,2),mean(weightedPSTHS(:,:,d)./binWidth,1),std(weightedPSTHS(:,:,d)./binWidth,0,1),'lineProps',...
                        {'Color',[co(d,:),.75-(0)],'LineWidth',1.5*double(size(weightedPSTHS,1)==1)+.5,'LineStyle',lc{s}});
                else
                    [bins centers] = hist(squeeze(weightedPSTHS(:,:,d)),linspace(...
                        quantile(weightedPSTHS(:,:,d),.20,'all'),quantile(weightedPSTHS(:,:,d),.80,'all'),10));
                    bins = bins ./ sum(bins);
                    allQuantiles= cellfun(@(m) cell2mat(cellfun(@(n) quantile([n{:}],[.01,.99],2),m,'Uniformoutput',false)'),somaProj,'UniformOutput',false);
                    allQuantiles = vertcat(allQuantiles{:});
                    bar(centers, bins, 'FaceColor',co(d,:));ylim([0 .5]);xlim([min(allQuantiles(:,1)),max(allQuantiles(:,end))]);
                end
            end
        end
    end
    if(plotTraj)
        [~,tn,~]= unique(cellfun(@(a) tilenum(a), ax)); ex = ax(tn);
        condSegs = round(mean(cell2mat(cellfun(@(m) mean(m(:,1:end-1),1,'omitnan'),epochSegs,'UniformOutput',false)),1,'omitnan'));
        cellfun(@(a) set(a,'XTick',0:50:size(weightedPSTHS,2)),ax);
        cellfun(@(a) set(a,'XTickLabels',linspace(timeBins(1),timeBins(end),length(get(a,'XTick')))),ax);
        cellfun(@(a) set(a,'XTickLabels',[]),ex(setdiff(1:prod(tl.GridSize),tilenum(tl,tl.GridSize(1),1:tl.GridSize(end)))));
        cellfun(@(a) set(a,'XLim',[0 size(weightedPSTHS,2)]),ax);
        cellfun(@(a) set(a,'YAxisLocation','right'),ex(tilenum(tl,1:tl.GridSize(1),tl.GridSize(end))));
        cellfun(@(a) arrayfun(@(x) plot(a,[x,x], get(a,'YLim'),'k--','LineWidth', 1.5),condSegs(1:2+(n==3))),ax);
        if(n<3)
            cellfun(@(a,p) cellfun(@(x,s) plot(a,[x;x], repmat(get(a,'YLim'),size(x,2),1)','LineStyle',':','Color',s,'LineWidth',1.5),num2cell(cell2mat(...
                cellfun(@(pc) round(mean(pc(:,end-1:end),1,'omitnan')),epochSegs(p),'UniformOutput',false)),2)',colors.values(cellstr(conditions(p)))),ax,plotSegs);
        end
        cellfun(@(a,s) set(a,'DataAspectRatio',[s{1}(1),2*s{1}(2),s{1}(3)],'PlotBoxAspectRatio',[s{2}(1)/2 s{2}(2) 1]),ax,...
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