figure(); tax=tiledlayout(max(1,num_dims/2),2);hold on;
projD = cat(3,D.data);
xbins = min(projD,[],'all'):max(projD,[],'all');
colorCondKeys = colors.keys; dist = {};
for icond = 1:length(conds)
    for idim = 1:num_dims
         nexttile(idim); hold on; title(idim);
        bins = histcounts(squeeze(projD(idim,:,icond)),[xbins,xbins(end)+1],'Normalization','percentage');
        bar(xbins, bins, 'FaceColor',cell2mat(colors.values(colorCondKeys(cellfun(@(D) contains(conds(icond),D),colors.keys)))),'FaceAlpha',.5);
        plot(xbins,pdf(fitdist(squeeze(projD(idim,:,icond))','Normal'),xbins).*100,'LineWidth',2,'Color',cell2mat(colors.values(colorCondKeys(cellfun(@(D) contains(conds(icond),D),colors.keys)))));
        dist{icond,idim} = squeeze(projD(idim,:,icond));
    end
    gm{icond} = cell2mat(dist(icond,:)')';
end
cvIdx=cvpartition(sTrials,'Holdout',0.1);
trIdx = training(cvIdx);
tsIdx = test(cvIdx);
trData = cellfun(@(m) m(trIdx,:),gm,'UniformOutput',false);
mus = cellfun(@(m) mean(m,1,'omitnan'),trData,'UniformOutput',false);
cv =  vertcat(trData{:})-cell2mat(cellfun(@(m,g) repmat(m,size(g,1),1),mus,trData,'UniformOutput',false)');
LDA=makecdiscr(cell2mat(mus'),cv'*(cv/(size(cv,1)-length(gm))),'ClassNames',conds);
CM = confusionmat(cellstr(cell2mat(arrayfun(@(r) repmat(string(r),sum(tsIdx),1),conds,'UniformOutput',false)')),...
    LDA.predict(cell2mat(cellfun(@(m) m(tsIdx,:),gm,'UniformOutput',false)')));
[c1,c2] = find(CM<sum(tsIdx) & CM~=0);
somC = arrayfun(@(s1,s2) strcmp(extractBefore(conds(s1),"_"),extractBefore(conds(s2),"_")), c1,c2)
condC = arrayfun(@(s1,s2) strcmp(extractBefore(extractAfter(conds(s1),"_"),"-"),extractBefore(extractAfter(conds(s2),"_"),"-")), c1,c2);



Dh = handles.orig_data;
origin = -mean([handles.orig_data.data],2);
[u sc lat] = pca([handles.orig_data.data]');
next_interval_start = 1;
for trial = 1:length(handles.orig_data)
    Dh(trial).data = sc(next_interval_start:next_interval_start+size(handles.orig_data(trial).data,2)-1, 1:num_dims)';
    next_interval_start = next_interval_start + size(handles.orig_data(trial).data,2);
end