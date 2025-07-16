clear all; close all;
monkeys = ["Gilligan"];
conds = ["[ExtraSmallSphere]","[LargeSphere]", "[Photocell]","[Rest]"];
condComb = nchoosek(1:length(conds),2);
HPk = 250;
LPk = 5;
drivePath = "D:\";

condComNames = nchoosek(string(conds),2);
condComb(end+1:end+length(conds),:) = repmat(1:length(conds),[2 1])';
condComNames(end+1:end+length(conds),:) = repmat(string(conds),[2 1])';
[mcIms,eIms,oIms] = deal(cell(length(monkeys),length(conds)));
for m = 1:length(monkeys)
    monkey = monkeys(m);
    if(strcmp(monkey,"Gilligan"))
        mm = MotorMapping(42);
        dates = ["12_05_2018", "12_07_2018","12_09_2018", "12_10_2018", "12_12_2018", "12_13_2018","12_14_2018","01_07_2019"];
        runs = [0, 0, 0, 1, 0, 0, 1, 0];
        frames = [49:53];
    elseif(strcmp(monkey,"Skipper"))
        mm = MotorMapping(56);
        dates = ["10_30_2020","11_09_2020","11_23_2020", "11_24_2020", "11_25_2020", "11_27_2020", "11_30_2020","12_01_2020","12_02_2020"];
        runs = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        frames = [48:52];
    end
    %frames = 32:70;
    [datesTable, masksA, ~] = getMonkeyInfo(drivePath,monkey,["M1", "PMd"],false);
    refMask = masksA{1};
    siteMask = getVoronoiMask(datesTable,mm,refMask,["Arm","Hand"]);
    monkeySiteMask{m} = siteMask./siteMask;
    %%
    for c = 1:length(conds)
        [currTrials,sidx] = loadAllFrames(drivePath,monkey,dates,runs,LPk,HPk,conds{c},false,frames);
        startIdx = [1 cumsum([cellfun(@length,sidx(1:end-1))]+1)];
        tDates = arrayfun(@(s,e) currTrials(:,:,s:e), startIdx, cumsum([cellfun(@length, sidx)]), 'UniformOutput',false);
        eDates = cellfun(@(t) mean(t(:,:,2:2:end),3,'omitnan'),tDates,'UniformOutput',false);
        oDates = cellfun(@(t) mean(t(:,:,1:2:end),3,'omitnan'),tDates,'UniformOutput',false);
        mDates = cellfun(@(t) mean(t,3,'omitnan'),tDates,'UniformOutput',false);
    end
    currTrials = [];
    mcIms{m,c} = mDates;
    eIms{m,c} = eDates;
    oIms{m,c} = oDates;
end
%%
% for m = 1;length(monkeys)
%     sessionIms = cellfun(@(i,r) i(:,:,r{d}), mcIMs(m,:), 'UniformOutput',false);
%     normR{d}{1} = nanmean(cat(3,sessionIms{:}),3);
%     normR{d}{2} = nanstd(cat(3,sessionIms{:}),0,3);
%     if(any(normR{d}{2}==0,'all'))
%         [r,c] = find(normR{d}{2}==0);
%         normR{d}{2}(r,c) = nanmin(abs(normR{d}{2}(normR{d}{2}~=0)));
%     end
% end
pairN = nchoosek(["e1","e2","o1","o2"],2);
condPairInds = ~(sum(contains(pairN, num2str(1)),2)==2 | sum(contains(pairN, num2str(2)),2)==2);
corrCondFrames = cell(1,length(condComb));
for m = 1:length(monkeys)
    figure();
    dispstat('','init');
    for c = 1:length(conds)
        if(condComb(c,2)==condComb(c,1))
            corrPairs = {[oIms{m,condComb(c,1)}],[eIms{m,condComb(c,2)}]};
        else
            c1E = {eIms{m,condComb(c,1)}};
            c2E = {eIms{m,condComb(c,2)}};
            c1O = {oIms{m,condComb(c,1)}};
            c2O = {oIms{m,condComb(c,2)}};
            corrPairs = nchoosek([c1E, c2E, c1O, c2O],2);
            corrPairs = corrPairs(condPairInds,:);
        end
        cMat = {};
        parfor p =1:size(corrPairs,1)
            c1 = corrPairs{p,1};
            c2 = corrPairs{p,2};
            c1N = cellfun(@isempty,c1);
            c2N = cellfun(@isempty,c2);
            c1(c1N) = {NaN(size(monkeySiteMask{m}))};
            c2(c2N) = {NaN(size(monkeySiteMask{m}))};
            c1 = cellfun(@(ee) ee.*monkeySiteMask{m}, c1, 'UniformOutput',false);
            c2 = cellfun(@(ee) ee.*monkeySiteMask{m}, c2, 'UniformOutput',false);
            cMat{p} = cell2mat(cellfun(@(cs1) cellfun(@(cs2) xcorr(cs1(~isnan(cs1)),...
                cs2(~isnan(cs2)),0,'normalized'),c2,'UniformOutput',true), c1,'UniformOutput',false)');
            cMat{p}(cMat{p}==0) = NaN;
        end
        corrCondFrames{m,c} = cat(3,cMat{:});
    end
    corrCondFrames(m,:) = cellfun(@(cg) mean(cg,3,'omitnan'),corrCondFrames(m,:)', 'UniformOutput',false);
end
gil = cell2mat(cellfun(@(g) g(logical(eye(size(g,1)))), corrCondFrames(1,:), 'UniformOutput', false));
skip = cell2mat(cellfun(@(s) s(logical(eye(size(s,1)))), corrCondFrames(2,:), 'UniformOutput', false));
allDiag = cat(1,gil,skip);

%%
% gilligan = load('S:\Lab\ngc14\Working\Gilligan\XCorr\XCorr_AVG_EO.mat');
% gDates = size(gilligan.corrMatSessions{1},1);
% figure('Units','normalized', 'Position',[0 0 1 1]);
% for c  = 1 : length(condComNames)
%     subplot(2,5,c)
%     corrMatSessionIm = gilligan.corrMatSessions{c};
%     corrMatSessionIm(corrMatSessionIm<0) = 0;
%     imagesc(corrMatSessionIm)
%     colormap(turbo)
%     clim([0 1])
%
%     title([condComNames(c,1) + condComNames(c,2),'Med & Mean', ...
%         [num2str(nanmedian(corrMatSessionIm,'all'),'%0.2f'),', ',...
%         num2str(nanmean(corrMatSessionIm,'all'),'%0.2f')],...
%         ['(Diag: ',num2str(nanmedian(corrMatSessionIm.*eye(gDates)./eye(gDates),'all'),'%0.2f'),', ',...
%         num2str(nanmean(corrMatSessionIm.*eye(gDates)./eye(gDates),'all'),'%0.2f'), ')']]);
% end
% saveas(gcf,"S:\Lab\ngc14\Working\Gilligan\XCorr\Session_Heatmaps_Norm",'png');
% saveas(gcf,"S:\Lab\ngc14\Working\Gilligan\XCorr\Session_Heatmaps_Norm",'epsc');

allDiag(allDiag<0) = 0;
groups = ones(size(allDiag));
groups(:,strcmp(condComNames(:,1), condComNames(:,2))) = 2;
groups(:,sum(strcmp(condComNames,"[Rest]"),2)==1) = 3;
groups(:,sum(strcmp(condComNames,"[Rest]"),2)==2) = 4;
[kw,kwt,kws] = kruskalwallis(allDiag(groups~=4 & ~isnan(allDiag)), groups(groups~=4 & ~isnan(allDiag)),'off');
[mc,y] = multcompare(kws,"Alpha",0.001,"CriticalValueType","bonferroni","Display","off");
xff = figure('Units','normalized','Position', [0 0 1 1 ]);
subplot(2,1,1);
boxchart(groups(~isnan(allDiag)),allDiag(~isnan(allDiag)),'Notch', 'on');
xticks(1:max(groups,[],'all'));
xlim([0,max(groups,[],'all')+1]);
xticklabels(["Movement","Within", "Rest", "RestRest"]);
ylim([0 1]);

[~,rG] = find(groups==1);
mGroups = unique(rG);
groupsM = repmat([1:length(mGroups)],size(groups,1),1);
boxP = allDiag(:,mGroups);
[kwM,kwtM,kwsM] = kruskalwallis(boxP(~isnan(boxP)), groupsM(~isnan(boxP)),'off');
[mcM,yM] = multcompare(kwsM,"Alpha",0.001,"CriticalValueType","bonferroni","Display","off");
figure(xff);
subplot(2,1,2);
boxchart(groupsM(~isnan(boxP)),boxP(~isnan(boxP)),'Notch', 'on');
xticks(1:max(groupsM,[],'all'));
xlim([0,max(groupsM,[],'all')+1]);
condPNames = condComNames(mGroups,1) + condComNames(mGroups,2);
xticklabels(condPNames);
ylim([0 1]);

saveas(gcf,["S:\Lab\ngc14\Working\Both_XCorr_"+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Diag_Box_With_Even_Odd"],'png');
saveas(gcf,["S:\Lab\ngc14\Working\Both_XCorr_"+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Diag_Box_With_Even_Odd"],'eps');