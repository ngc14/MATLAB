monkey = "Skipper";
conds = ["[ExtraSmallSphere]","[LargeSphere]", "[Photocell]","[Rest]"];
HPk = 250;
LPk = 5;
drivePath = "X:\Lab";
saveFold = "S:\Lab\ngc14\Working\"+monkey+"\XCorr\";
if(strcmp(monkey,"Gilligan"))
    mm = MotorMapping(42);
    dates = ["12_05_2018", "12_07_2018","12_09_2018", "12_10_2018", "12_12_2018", "12_13_2018","12_14_2018","01_07_2019"];
    runs = [0, 0, 0, 1, 0, 0, 1, 0];
    frames = {49:53};
elseif(strcmp(monkey,"Skipper"))
    mm = MotorMapping(56);
    dates = ["10_30_2020","11_09_2020","11_23_2020", "11_24_2020", "11_25_2020", "11_27_2020", "11_30_2020","12_01_2020","12_02_2020"];
    runs = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    frames = {48:52};
end
frames = num2cell(1:70);
[datesTable, masks, ~] = getMonkeyInfo(drivePath+"\",monkey,["M1", "PMd"],false);
monkeySiteMask = getVoronoiMask(datesTable,mm,masks{1},["Arm","Hand"]);
monkeySiteMask = (monkeySiteMask & masks{1}) ./ monkeySiteMask;

condComb = nchoosek(1:length(conds),2);
condComb(3:4,:) = flipud(condComb(3:4,:));
condComNames = conds(condComb);
% condComb(end+1:end+length(conds),:) = repmat(1:length(conds),[2 1])';
%%
for cc = 1:length(conds)
    condFrameIntensity = cell(1,length(frames));
    for f = 1:length(frames)
        tic
        currTrials= loadAllTrials(drivePath,monkey,dates,runs,LPk,HPk,conds(cc),frames{f});
        condFrameIntensity{f} = squeeze(mean(cat(3,currTrials{:}),3,'omitnan'));
        disp(toc);
    end
    forelimbIntensity{cc}= condFrameIntensity;
end
corrComps = NaN(length(condComb),length(frames));
%%
for cm = 1:length(condComb)
    for f = 1 :length(frames)
        if(condComb(cm,1)~=condComb(cm,2))
            currCondIms1 = forelimbIntensity{condComb(cm,1)};
            currCondIms2 = forelimbIntensity{condComb(cm,2)};
            corrComps(cm,f) = xcorr(currCondIms1{f}(~isnan(monkeySiteMask)),currCondIms2{f}(~isnan(monkeySiteMask)),0,'normalized');
        else
            % currCondIms1 = evens{condComb(cm,1)};
            % currCondIms2 = odds{condComb(cm,2)};
        end
    end
end
trialCorrs = cell(length(dates),length(condComb));
%% correlations
for c =1:length(condComb)
    currCondIms1 = forelimbIntensity{condComb(c,1)};
    currCondIms2 = forelimbIntensity{condComb(c,2)};
    ci = datetime('now');
    disp(strjoin(["Correlating ", condComNames(c,:), " trials"],''));
    for d = 1:length(dates)
        [~,cd1] = evalc('cell2mat(gather(currCondIms1(1,d)))');
        [~,cd2] = evalc('cell2mat(gather(currCondIms2(1,d)))');
        corrMat = NaN(size(cd1,3),size(cd2,3));
        parfor x = 1:size(cd1,3)
            corrX = NaN(1,size(cd2,3));
            currCondImsX = cd1(:,:,x);
            currCondImsX = currCondImsX(~isnan(currCondImsX));
            for y = 1:size(cd2,3)
                currCondImsY = cd2(:,:,y);
                currCondImsY = currCondImsY(~isnan(currCondImsY));
                corrX(1,y) = xcorr(currCondImsX, currCondImsY,0,'normalized');
            end
            corrMat(x,:) = corrX;
        end
        trialCorrs{d,c} = corrMat;
    end
    disp(strjoin(["Correlations done. (", string(duration(datetime('now')-ci)), ")"], ''));
end
save(saveFold+"XCorr_"+num2str(frames(1))+"_"+num2str(frames(end))+".mat",...
    'trialCorrs','-mat');
%% heatmaps
figure('Units','normalized','Position',[0 0 1 1]);
plotCorrMat = trialCorrs;
nanDates = cellfun(@(n) all(n==0,'all')|all(isnan(n),'all'), plotCorrMat);
plInd = 1;
nBox = cell(size(plotCorrMat,1),3);
for c = 1:length(condComb)
    if(c<7)
        condB = cell2mat(cellfun(@(dd) dd(:), plotCorrMat(~nanDates(:,c),c),'UniformOutput',false));
        wCond1 = all(condComb == condComb(c,1),2);
        wCond2 = all(condComb == condComb(c,2),2);
        cond1 = cell2mat(cellfun(@(d1,c1) d1(tril(true(size(d1)),-1)),plotCorrMat(~nanDates(:,c),wCond1), 'UniformOutput',false));
        cond2 = cell2mat(cellfun(@(d2,c2) d2(tril(true(size(d2)),-1)),plotCorrMat(~nanDates(:,c),wCond2), 'UniformOutput', false));
        catCond = [NaN(length(condB)-(length(cond1)+length(cond2)),1);cond1;cond2];
        if(c<4)
            for d = 1:size(plotCorrMat,1)
                if(nanDates(d,c))
                    nBox{d,c} = NaN;
                else
                    nBox{d,c} = plotCorrMat{d,c};
                end
            end
        end
    else
        condB = cell2mat(cellfun(@(dd) dd(tril(true(size(dd)),-1)), plotCorrMat(~nanDates(:,c),c), 'UniformOutput',false));
    end
    condB(condB<0) = NaN;
    for d = 1:length(dates)
        subplot(length(condComb),length(dates), plInd);
        plInd = plInd +1;
        currMat = plotCorrMat{d,c};
        currMat(currMat<=0) = 0;
        imagesc(currMat);
        colormap(jet)
        clim([0 1])
        if(d==1)
            xPos = get(gca,'Position');
            text(xPos(1)+0.2, 1, condComNames(c,1) + condComNames(c,2),'HorizontalAlignment', 'right');
        end
        if(d==length(dates))
            text(55,10,[num2str(median(cellfun(@(n) median(n.*double(n>0),'all','omitnan'),...
                plotCorrMat(~nanDates(:,c),c)),'omitnan'),'%0.2f'),', ',num2str(mean(...
                cellfun(@(n) mean(n.*double(n>0),'all','omitnan'),plotCorrMat(~nanDates(:,c),c)),'omitnan'),...
                '%0.2f')],'HorizontalAlignment', 'left');
        end
        title(['Med&Mean: ', num2str(median(currMat,'all','omitnan'),'%0.2f'), ', ',...
            num2str(mean(currMat,'all','omitnan'),'%0.2f')]);
    end
end
saveas(gcf,saveFold+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Heatmaps",'png');
saveas(gcf,saveFold+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Heatmaps",'eps');
%% statistical testing
tsCorrs = cellfun(@(p) p(tril(true(size(p)))), plotCorrMat, 'UniformOutput', false);
groups = ones(size(tsCorrs));
groups(:,strcmp(condComNames(:,1), condComNames(:,2))) = 2;
groups(:,sum(strcmp(condComNames,"[Rest]"),2)==1) = 3;
groups(:,sum(strcmp(condComNames,"[Rest]"),2)==2) = 4;

tsCorrsV = vertcat(tsCorrs{:});
groupTrials =  cellfun(@(r,b) repmat(r,length(b),1),num2cell(groups),tsCorrs,'UniformOutput', false);
groupTrials = vertcat(groupTrials{:});
[~,rG] = find(groups==1);
mGroups = reshape(repmat(1:length(unique(rG)),size(groups,1),1),[],1);
groupTrialsM =  cellfun(@(r,b) repmat(r,length(b),1),num2cell(mGroups),tsCorrs(groups==1),'UniformOutput', false);
groupTrialsM = vertcat(cell2mat(groupTrialsM));

[kw,kwt,kws] = kruskalwallis(tsCorrsV(groupTrials~=4),groupTrials(groupTrials~=4),'off');
[mc,y] = multcompare(kws,"Alpha",0.001,"CriticalValueType","bonferroni","Display","off");
[kwM,kwtM,kwsM] = kruskalwallis(tsCorrsV(groupTrials==1), groupTrialsM,'off');
[mcM,yM] = multcompare(kwsM,"Alpha",0.001,"CriticalValueType","bonferroni","Display","off");

figure('Units','normalized','Position',[0 0 1 1]);
plotDists(1,groupTrials(groupTrials~=4),tsCorrsV(groupTrials~=4),["Movement","Within", "Rest", "RestRest"],mc);
plotDists(2,groupTrialsM,tsCorrsV(groupTrials==1),condComNames(unique(rG),1) + condComNames(unique(rG),2),mcM);

saveas(gcf,saveFold+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Violin",'png');
saveas(gcf,saveFold+ num2str(frames(1))+ "_"+num2str(frames(end))+"_Violin",'eps');

function plotDists(sI,groupLab,corrVals,xLabs,stp)
gcf();

nGroups = length(unique(groupLab));
subplot(2,2,sI);
hold on;
title('Negative')
reorderCorr = cell(1,nGroups);
for n = 1:nGroups
    reorderCorr{n} = corrVals(groupLab==n);
end
[ax,L,mX,mS,bw] = violin(reorderCorr,'xlabel',num2cell(xLabs),'mc',[],...
    'medc','k','facecolor',['m','g','b']','edgecolor',[],'plotlegend',0);
for cm = 1:size(stp,1)
    if(sign(stp(cm,3)) == sign(stp(cm,5)))
        line(stp(cm,1:2),repmat(-1+((cm-1) * .05),1,2), 'Color', [.25 .25 .25], 'LineWidth',2)
    end
end
%boxchart(groupLab,corrVals,'Notch', 'on');
xticks(1:nGroups);
xticklabels(xLabs);
ylim([-1 1.2]);

subplot(2,2,sI+2);
hold on;
title('Clipped')
clippedCorrs = cellfun(@(r) r.*double(r>0), reorderCorr, 'UniformOutput',false);
[ax,L,mX,mS,bw] = violin(clippedCorrs,'xlabel',num2cell(xLabs),'mc',[],...
    'medc','k','facecolor',['m','g','b']','edgecolor',[],'plotlegend',0);
%boxchart(groupLab,clippedCorrs,'Notch', 'on');
xticks(1:nGroups);
xticklabels(xLabs);
ylim([0 1.2]);
end