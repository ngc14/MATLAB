load("C:\Users\ngc14\Desktop\EMG_MaxNorm.mat");
muscles = cellstr({'Deltoid','Biceps', 'Triceps', 'Wrist Extensor', 'Wrist Flexor','Digit Extensor','Digit Flexor'});
smoothWin = .15;
num_dims = 5;
timeBins = [-.5, 1.5];
maxLags = 200/binWidth;
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
    t(end-2)+(.5*smoothK)),gausswin(smoothK)./sum(gausswin(smoothK)),'valid'),[size(m,1),max(cell2mat(trialLength),[],'all')],'Pattern','edge','side','both'),...
    sc,ti,'UniformOutput',false)'),1,'omitnan'),mc,tc,'UniformOutput',false),num2cell(normTrials,1),sessionSegs,'UniformOutput',false);
clear normTrials sessionSegs
avgEMG = [smoothedData{:}];
%%
[r,l]=cellfun(@(s) cellfun(@(u,ma) cellfun(@(m) arrayfun(@(n) xcorr(u(:,n)',m,maxLags,'normalized'),transpose(1:num_dims), 'UniformOutput',false),...
    ma, 'UniformOutput',false),s,num2cell(cellfun(@(d) interp1(1:length(d),d,linspace(1,length(d),size(avgDim{1},1))),...
    avgEMG,'UniformOutput',false),2)','UniformOutput', false), num2cell(avgDim,2), 'UniformOutput',false);
maxCorrs = cellfun(@(t) cellfun(@(s)cell2mat(cellfun(@(u) cellfun(@(v) max(v(maxLags:end)), u),s,'UniformOutput',false)),t,'UniformOutput',false),r, 'UniformOutput', false);
maxCorrs = vertcat(maxCorrs{:});
tc = tiledlayout(size(maxCorrs,1), size(maxCorrs,2));
for s = 1:size(maxCorrs,1)
    for c = 1:size(maxCorrs,2)
        nexttile([1,1]);imagesc(maxCorrs{s,c}); colormap('jet'); colorbar; clim([-1 1]);
        hold on; title(string(somaReps(s)) + " - " + conditions(c));
        yticks(1:num_dims); yticklabels(arrayfun(@num2str,1:num_dims,'UniformOutput',false));xticklabels(muscles);
    end
end
saveFigures(gcf,savePath,"EMG_Corr",[]);