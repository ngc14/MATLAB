load("C:\Users\ngc14\Desktop\EMG_MaxNorm.mat");
muscles = cellstr({'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor','Digit Extensor','Digit Flexor'});
smoothWin = .15;
num_dims = 5;
timeBins = [-.5, 1.5];
binWidth = 10;
maxLags = 50/binWidth;
somaReps = ["Arm","Hand"];
avgDim = load("S:\Lab\ngc14\Working\DataHi\Traj\Full_Population\AvgDimTraj_Sm150.mat");
savePath = "S:\Lab\ngc14\Working\DataHi\EMG\";
%%
smoothK = smoothWin*Fs;
allSegs = cellfun(@(a) cell2mat([a{:}]'), num2cell([sessionSegs{:}],2), 'UniformOutput',false);
allSessions = unique(cell2mat(cellfun(@(u) unique(string(u)),[sessionDates{:}],'UniformOutput',false)));
for m = 1:length(muscles)
    for s = 1:length(allSessions)
        allCondSig = cellfun(@(as,d) as(string(d)==allSessions(s)),sessionSigs{m}, sessionDates{m}, 'UniformOutput', false);
        allCondSegs = cellfun(@(g,d) g(string(d)==allSessions(s)),sessionSegs{m}, sessionDates{m}, 'UniformOutput', false);
        allCondSegs = cellfun(@(c) cellfun(@(t) t(1),c), allCondSegs,'UniformOutput',false);
        normBaseline = cellfun(@(c,s) cellfun(@(t,a) max(1,max(t(a:a+2*Fs))),...
            c,num2cell(s),'UniformOutput', false),allCondSig,allCondSegs,'UniformOutput',false);
        maxSession =mean(maxk(cell2mat(cellfun(@cell2mat,normBaseline,'UniformOutput',false)'),3));
        mTrials{s} =cellfun(@(c) cellfun(@(n)(n./maxSession),c,'UniformOutput', false),allCondSig,'UniformOutput',false);
        allCondSig = [];    allCondSegs = [];
    end
    normTrials(:,m) = cellfun(@(c) [c{:}], num2cell([mTrials{:}],2), 'UniformOutput',false);
    mTrials = {};
end
clear sessionDates sessionSigs 
trialLength = cellfun(@(tc) cellfun(@(ti,cs) max(cellfun(@(t) t(strcmp(cs,"StartReplaceSuccess"))-(t(strcmp(cs,"GoSignal"))+...
    timeBins(1)*Fs),ti))+smoothK,tc,ConditionSegs),sessionSegs, 'UniformOutput',false);
smoothedData = cellfun(@(mc,tc) cellfun(@(sc,ti) mean(cell2mat(cellfun(@(m,t) resize(conv(m(:,(-smoothK/2)+t(1)+(timeBins(1)*Fs):...
    t(end-2)+(.5*smoothK)),gausswin(smoothK)./sum(gausswin(smoothK)),'valid'),[size(m,1),max(cell2mat(trialLength),[],'all')],'Pattern','edge','side','trailing'),...
    sc,ti,'UniformOutput',false)'),1,'omitnan'),mc,tc,'UniformOutput',false),num2cell(normTrials,1),sessionSegs,'UniformOutput',false);
clear normTrials sessionSegs
avgEMG = [smoothedData{:}];
%%
avgEMGInt = cellfun(@(d) interp1(1:length(d),d,linspace(1,length(d),size(avgDim{1},1))),avgEMG,'UniformOutput',false)';
[r,l]=cellfun(@(s) cellfun(@(u,ma) cellfun(@(m) arrayfun(@(n) xcorr(u(:,n)',m,maxLags,'normalized'),transpose(1:num_dims), 'UniformOutput',false),...
    ma, 'UniformOutput',false),s,num2cell(avgEMGInt,1),'UniformOutput', false), num2cell(avgDim,2), 'UniformOutput',false);
maxCorrs = cellfun(@(t) cellfun(@(s)cell2mat(cellfun(@(u) cellfun(@(v) max(v(maxLags:end)), u),s,'UniformOutput',false)'),t,'UniformOutput',false),r, 'UniformOutput', false);
maxCorrs = vertcat(maxCorrs{:});
tc = tiledlayout(size(maxCorrs,1), size(maxCorrs,2));
for s = 1:size(maxCorrs,1)
    for c = 1:size(maxCorrs,2)
       % psthXCorr = cellfun(@(m) xcorr(mean(cell2mat(psthGrouped{s}{c}),1,'omitnan'),m,maxLags,'normalized'),avgEMGInt(:,c),'UniformOutput',false);psthXCorr = cellfun(@(p) max(p(maxLags:end)), psthXCorr)';
        nexttile([1,1]);imagesc([maxCorrs{s,c}]); colormap('jet'); colorbar; clim([-1 1]);
        hold on; title(string(somaReps(s)) + " - " + conditions(c));
        yticks(1:num_dims+1); yticklabels([arrayfun(@num2str,1:num_dims,'UniformOutput',false)]);xticklabels(muscles);
    end
end
saveFigures(gcf,savePath,"EMGxDim_Corr",[]);
%%
tc = tiledlayout(2,size(maxCorrs,2),'TileIndexing','columnmajor');
co = colororder({'g','b','c','r','m','y'});
co(end+1,:) = [1 .5 0];
colororder(co);
for c = 1:size(maxCorrs,2)
    nexttile([1,1]); hold on; cellfun(@(p) plot(p,'LineWidth',1.5),avgEMG(c,:));
    title(conditions(c)); legend(muscles,'AutoUpdate','off');
    eventM = mean(allSegs{c},1,'omitnan');
    eventM = eventM-(eventM(1)+(timeBins(1)*Fs));
    arrayfun(@(l) line([l,l],get(gca,'YLim'),'Color','k','LineStyle','--'), eventM);
    xlim([eventM(1),eventM(end-2)]);
    xticklabels(get(gca,'XTick')./Fs);
    emgCorr = NaN(size(avgEMG,2),size(avgEMG,2));
    for a = 1:size(emgCorr,1)
        for b = 1:size(emgCorr,2)
            emgCorr(a,b)=corr(avgEMG{c,a}',avgEMG{c,b}');
        end
    end
    nexttile([1,1]); hold on; imagesc(emgCorr); xticks(1:length(muscles)); yticks(1:length(muscles)); xticklabels(muscles);yticklabels(muscles);
    colormap('jet'); clim([0 1]); axis tight; colorbar;
end
saveFigures(gcf,savePath,"EMG_Corr",[]);
