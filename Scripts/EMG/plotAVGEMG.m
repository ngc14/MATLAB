monkeys = ["Gilligan","Skipper"];
drive = "S:\Lab\"
alignments = [{'GoSignal'},{'StartReach'},{'StartHold'},{'StartWithdraw'}];
alignWindows = {[-.2 .2],[-.5 .15],[-.20 .20],[-.15 .5]};
phaseWindows = {[0 0.2],[-.15 .05],[-.2 0.0],[-.15 .05]};
gap = .1;
groupings = {{"Deltoid.mat","Biceps.mat","Triceps.mat"},...
    {"Wrist Extensor.mat","Wrist Flexor.mat","Digit Extensor.mat","Digit Flexor.mat"}};
muscles = string([groupings{:}]);
groupInds = cellfun(@(g) contains(muscles,string(g)), groupings, 'UniformOutput',false);
groupNames = cellfun(@(g) cellfun(@(s) string(s{1}(1:end-4)),g), groupings,'UniformOutput',false);
pVal = 0.01;
condLabels = {'E', 'L', 'P'};
saveFigs = true;
savePath = "C:\Users\ngc14\Desktop\";
if(saveFigs && ~exist(savePath,'dir'))
    mkdir(savePath);
    mkdir(savePath+"Phase_Boxplot\")
end
%%
sessionSigs =repmat({repmat({[]},length(condLabels),1)},1,length(muscles));
sessionSegs =repmat({repmat({[]},length(condLabels),1)},1,length(muscles));
sessionDates =repmat({repmat({[]},length(condLabels),1)},1,length(muscles));
tic
for m =1:length(monkeys)
    mPath = drive+monkeys(m)+"\Sorted EMG Data\Averages_Raw\";
    mFiles = arrayfun(@(l)  matfile(mPath+l), muscles,'UniformOutput', false);
    Conditions = string(mFiles{1}.Conditions);
    ConditionSegs = mFiles{1}.ConditionSegs;
    Conditions = Conditions(~strcmp(Conditions,"Rest"));
    ConditionSegs = ConditionSegs(~strcmp(Conditions,"Rest"));
    Fs = mFiles{1}.Fs;
    monkSigs = cellfun(@(m) m.rawTrials, mFiles, 'UniformOutput',false);
    monkSegs = cellfun(@(m) m.segs, mFiles, 'UniformOutput',false);
    monkDates = cellfun(@(m) m.dates, mFiles, 'UniformOutput',false);
    sessionSigs = cellfun(@(s,m) cellfun(@(sc,mc) horzcat(sc,mc), ...
        s(1:length(Conditions),1),m(1:length(Conditions),1),'UniformOutput',false), sessionSigs, monkSigs, 'UniformOutput',false);
    sessionSegs = cellfun(@(s,m) cellfun(@(sc,mc) horzcat(sc,mc), ...
        s(1:length(Conditions),1),m(1:length(Conditions),1),'UniformOutput',false),sessionSegs, monkSegs, 'UniformOutput',false);
    sessionDates = cellfun(@(s,m) cellfun(@(sc,mc) horzcat(sc,mc), ...
        s(1:length(Conditions),1),m(1:length(Conditions),1),'UniformOutput',false), sessionDates,  monkDates,'UniformOutput',false);
end
toc
allSessions = cellfun(@(u) unique(string(u)),[sessionDates{:}],'UniformOutput',false);
allSessions = unique([allSessions{:}]);
[alignWindows,phaseWindows] = cellfun(@(a,p) deal(a.*Fs,p.*Fs),alignWindows,phaseWindows,'UniformOutput',false);
clear monkSigs monkSegs monkDates
%%
bVal = {};
for m = 1:length(muscles)
    mTrials = {};
    normBaseline = {};
    mSigs = sessionSigs{m};
    mDates = sessionDates{m};
    mSegs = sessionSegs{m};
    for s = 1:length(allSessions)
        allCondSig = cellfun(@(as,d) as(string(d)==allSessions(s)),mSigs, mDates, 'UniformOutput', false);
        allCondSegs = cellfun(@(g,d) g(string(d)==allSessions(s)),mSegs, mDates, 'UniformOutput', false);
        allCondSegs = cellfun(@(c) cellfun(@(t) t(1),c), allCondSegs,'UniformOutput',false);
        % normBaseline{s} = mean(cellfun(@(c,s) max(1,mean(c(s:s+Fs),2,'omitnan')),...
        %     [allCondSig{:}],num2cell([allCondSegs{:}])));
        % mTrials{s} =cellfun(@(c) cellfun(@(n)(n./normBaseline{s}),c,'UniformOutput', false),...
        %     allCondSig,'UniformOutput',false);
        normBaseline = cellfun(@(c,s) cellfun(@(t,a) max(1,max(t(a:a+2*Fs))),...
            c,num2cell(s),'UniformOutput', false),allCondSig,allCondSegs,'UniformOutput',false);
        maxSession =mean(maxk(cell2mat(cellfun(@cell2mat,normBaseline,'UniformOutput',false)'),3));
        mTrials{s} =cellfun(@(c) cellfun(@(n)(n./maxSession),c,'UniformOutput', false),...
            allCondSig,'UniformOutput',false);
    end
    allCondSig = [];
    bVal{m} = normBaseline;
    normTrials(:,m) = cellfun(@(c) [c{:}], num2cell([mTrials{:}],2), 'UniformOutput',false);
end
allSessionSegs = horzcat(sessionSegs{:});
%%
[mTrials, mSigs, mDates, mSegs] = deal({});
for c = 1:length(Conditions)
    for a = 1:length(alignments)
        wind = alignWindows{a};
        alignInd = find(strcmp(ConditionSegs{c},alignments{a}));
        alignedSig{c,a} = cellfun(@(s,t) cell2mat(cellfun(@(ts,ss) ts(...
            ss(alignInd)+wind(1):ss(alignInd)+wind(end)),vertcat(s(:))',t,'UniformOutput',false)'),...
            normTrials(c,:),allSessionSegs(c,:),'UniformOutput',false);
        allSegs{c,a} = cellfun(@(s) mean(cell2mat(cellfun(@(t) t-t(alignInd),s,'UniformOutput',false)'),1,'omitnan'),...
            allSessionSegs(c,:),'UniformOutput',false);
    end
end
%% BOX PLOTS %%
rawActivity = {};
FStat = [];
stats = {};
cmpNames = join([nchoosek(Conditions,2)],',')';
cmpColors = dictionary(cmpNames,{[1 0 1], [0 1 1], [1 .7 0]});
for g = 1:length(groupInds)
    currGroup = cellfun(@(ms) ms(groupInds{g}), alignedSig, 'UniformOutput', false);
    for a = 1:length(alignments)
        p = fix(phaseWindows{a}-alignWindows{a});
        %each sigDiff (# of groups by # of alignments) is cond by muscles
        rawActivity{g,a} = cellfun(@(mu) cellfun(@(t) mean(t(:,p(1)+1:...
            min(arrayfun(@(a) a.*double((a>p(1))./(a>p(1))),[size(t,2),size(t,2)+p(end)-1]))),2,'omitnan'),...
            mu, 'UniformOutput',false), currGroup(:,a), 'UniformOutput',false);
        rawActivity{g,a} = cellfun(@(mu) cell2mat(cellfun(@(t) vertcat(mean(t(:,p(1)+1:...
            size(t,2)+p(end)-1),2),NaN(max(cellfun(@(s) size(s,1),mu))-size(t,1),1)),...
            mu, 'UniformOutput',false)), currGroup(:,a), 'UniformOutput',false);
    end
end
%%
for n = 1:numel(rawActivity)
    [r,c] = ind2sub(size(rawActivity),n);
    combinedConds=cellfun(@(a) cellfun(@(t) mean(vertcat(t,NaN(...
        max(cellfun(@length,a)-length(t)),size(t,2))),2,'omitnan')',a,'UniformOutput',false),...
        rawActivity(r,c),'UniformOutput',false);
    combinedConds=combinedConds{1};
    % combinedConds=cellfun(@(a)vertcat(a{:}).*cell2mat(cellfun(@(t)...
    %     repmat(sum(~isnan(t)),1,length(t)),a,'UniformOutput', false)),rawActivity(r,c),'UniformOutput',false);
    % [FStat(r,c),~,stats{r,c}] = anova1(horzcat(combinedConds{:}),cell2mat(cellfun(@(l,n)...
    %     ones(1,length(l)).*n, combinedConds',num2cell(1:length(Conditions)),'UniformOutput',false)),...
    %     'off');
    figure('Units','normalized','Position',[0 0 1 1]);
    hold on;
    title([strjoin([groupings{r}{:}],",")+"- "+alignments{c}]);
    outliers = horzcat(combinedConds{:});
    outliers = outliers(isoutlier(horzcat(combinedConds{:}),'quartiles'));
    b=boxplot(horzcat(combinedConds{:}),cell2mat(cellfun(@(n) ones(1,size(combinedConds{1},2)).*n,...
        num2cell(1:length(Conditions)),'UniformOutput',false)),'Notch','on');
    xticklabels(Conditions);
    xlim([0.5 3.5])
    yl = get(gca,'YLim');
    ylim([0 yl(end)]);
    if(saveFigs)
        exportgraphics(gcf,savePath+"Phase_Boxplot\"+strjoin([groupings{r}{:}],'_')+"_"+alignments{c}+'.png','Resolution',300,'ContentType','vector')
        exportgraphics(gcf,savePath+"Phase_Boxplot\"+strjoin([groupings{r}{:}],'_')+"_"+alignments{c}+'.eps','Resolution',300,'ContentType','vector')
    end
end
%% TIMECOURSES %%
close all;
pColors =([.7 0 0; .8 .4 0; 0 0 .7]);
fx = arrayfun(@gca,arrayfun(@(f) figure('Units','normalized','Position',[0 0 1 1]),...
    1:length(groupings)),'UniformOutput', false);
cellfun(@(f) hold(f,'on'),fx);
l={};
for c = 1:length(Conditions)
    currEventSegs = ConditionSegs{c};
    numSegs = length(currEventSegs);
    plotted = false(length(groupInds),length(currEventSegs));
    plotted(:,contains(string(currEventSegs),["StartGrasp","StartReplaceHold"])) = true;
    mSegs = cellfun(@(mc) cellfun(@(g) mean(mc(g,:),1,'omitnan'), groupInds, 'UniformOutput', false),...
        cellfun(@(am) cell2mat(am'),allSegs(c,:),'UniformOutput', false),'UniformOutput', false);
    xTicks = [];
    for align = 1:length(alignments)
        % plotted = false(1,length(currEventSegs));
        % plotted(contains(string(currEventSegs),["StartGrasp","StartReplaceHold"])) = true;
        wind = alignWindows{align};
        alignInd = find(strcmp(string(currEventSegs),alignments{align}));
        alignSession = cellfun(@(gI) cellfun(@(sm) sm(1:min(cellfun(@(s) size(s,1), alignedSig{c,align}(gI))),:),...
            alignedSig{c,align}(gI),'UniformOutput', false), groupInds, 'UniformOutput',false);
        alignSession = cellfun(@(a) mean(cat(3,a{:}),3,'omitnan'),alignSession, 'UniformOutput', false);
        %fx{c}=gca(figure(c)); hold on; legend('AutoUpdate','off');
        if(align==1)
            plotStart = 0;
        else
            plotStart =(gap*Fs) + xTicks{align-1}(end)+abs(wind(1));
        end

        % l{end+1}=cellfun(@(n,r,f) shadedErrorBar(plotStart+wind(1):plotStart+wind(2),...
        %     mean(n,1,'omitnan'),.5*std(n,1,1,'omitnan'),.../sqrt(size(n,1))
        %     'lineprops',{'Color',r,'LineWidth',2},...
        %     'ax',f),alignSession,repmat(num2cell(pColors(c,:),2),1,length(groupings)),fx);

        l{end+1}=cellfun(@(n,r,f) plot(f,plotStart+wind(1):plotStart+wind(2),...
            n(1:100,:),'Color',r,'LineWidth',1.25),alignSession,repmat(num2cell(pColors(c,:),2),1,...
            length(groupings)),fx,'UniformOutput',false);

        cellfun(@(f,n) plot(f,plotStart+wind(1):plotStart+wind(2),median(n,1,'omitnan'),'Color',pColors(c,:).*.85,'LineWidth',4), fx,alignSession);
        %plot(fx{c},plotStart+wind(1):plotStart+wind(end),mean(cell2mat(alignSession),1,'omitnan'),'k', 'LineWidth', 4,'LineStyle','-.');
        averageSegs = cellfun(@(ms) mean(ms,1,'omitnan'),mSegs{align},'UniformOutput',false);
        beforeSegs = cellfun(@(b)plotStart + cumsum(diff(b(alignInd:-1:1))), averageSegs,'UniformOutput', false);
        afterSegs = cellfun(@(a) plotStart + cumsum(diff(a(alignInd:end))), averageSegs, 'UniformOutput', false);
        avgSegs = cellfun(@(b,a) [fliplr(b) plotStart a], beforeSegs, afterSegs, 'UniformOutput', false);
        for s = 1:length(avgSegs{1})
            for g = 1:length(groupInds)
                if(avgSegs{g}(s)>=plotStart+wind(1) && avgSegs{g}(s)<=plotStart+wind(end) && ...
                        (~plotted(g,s) || avgSegs{g}(s)==plotStart))
                    plotted(g,s) = true;
                    if(avgSegs{g}(s)==plotStart)
                        plotColor = [.4 .4 .4];
                    elseif(contains(Conditions(c), 'Sphere') && s==4)
                        plotColor = [.7 0 .7];
                    else
                        plotColor = pColors(c,:);%[.7 .7 .7];%
                    end
                    plot(fx{g},[avgSegs{g}(s) avgSegs{g}(s)],[-5 500],'Color',plotColor,'LineStyle','--');
                end
            end
        end
        xTicks{align} = plotStart+[wind(1) 0 wind(end)];
    end
end
yL = [-5 200];
xlim([fx{:}],[xTicks{1}(1), xTicks{end}(end)]);
ylim([fx{:}],yL);
cellfun(@(a,s) title(a,s), fx, groupings);%num2cell(fNames)')%Conditions(1:end-1)');
set([fx{:}], 'XTick', sort([xTicks{:}]));
alignNames = cellfun(@(x,w) [num2str(x(1)/Fs),string(w),...
    num2str(x(end)/Fs)],alignWindows,alignments,'UniformOutput',false);
set([fx{:}], 'XTickLabels', [alignNames{:}]);
l=[l{1:length(alignments):end}];
l=l(1:length(groupings):end);
%cellfun(@(f) legend(f,[l.mainLine],Conditions),fx)%,cellfun(@(n) string(n(1:end-4)), fNames))
if(saveFigs)
    cellfun(@(f,n) exportgraphics(f,savePath+strjoin(n+".png",'_'),'Resolution',300,'ContentType','vector'),fx,groupings);
    cellfun(@(f,n) exportgraphics(f,savePath+strjoin(n+".eps",'_'),'Resolution',300,'ContentType','vector'),fx,groupings);
end
%cellfun(@(f,n) saveas(f,drive+"Lab\ngc14\Working\"+monkey+"\EMG\Muscles\Forelimb"+cellfun(@(cn) cn(1), regexp(n,' ','split'))+'png'),fx,Conditions(1:end-1)'