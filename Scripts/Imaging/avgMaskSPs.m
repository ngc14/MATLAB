%close all;
clear all;
monkey = 'Skipper';
indiTrials = 0;
if(strcmp(monkey, 'Gilligan'))
    dates = {'12_03_2018', '12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    dates = {'12_19_2018'};
    runs = {[1], [0], [0], [0], [1], [0], [0], [1], [0]};
    runs = {[0]};
    HPk = 550;
else
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    dates = {'12_09_2020','12_18_2020'};
    runs = {[0],[0]};
    HPk = 250;
end

dates = {'12_09_2020', '12_09_2020', '12_18_2020', '12_18_2020'};
runs = {[0], [1], [0], [1]};
stimSelect = [4, 5];                   % stim #s to do superpixel for
stims = {'LargeSphere', 'Photocell'};
SPData = {};
for s = 1:length(stimSelect)
    for d = 1:length(dates)
        if(indiTrials)
            SPDataLoad = load(['S:\Lab\', monkey, '\All Data\', monkey, '_', dates{d}, '\Imaging\run0', num2str(runs{d}), '\Results\Superpixel\LP5HP',num2str(HPk),'C0\SPData_',num2str(stimSelect(s)),'.mat']);
            if d==1
                SPData{s} = smoothdata(SPDataLoad.SPData.maskMeanAll{1},2);
            else
                SPData{s} = cat(3,SPData{s}, smoothdata(SPDataLoad.SPData.maskMeanAll{1},2));
            end
        else
            SPDataLoad = load(['S:\Lab\', monkey, '\All Data\', monkey, '_', dates{d}, '\Imaging\run0', num2str(runs{d}), '\Results\Superpixel\LP5HP',num2str(HPk),'C0\SPData.mat']);
            if(d==1)
                SPData{s} = SPDataLoad.SPData.maskMean_stim_Avg2{s};
                if(s==1)
                    SPData{length(stimSelect)+1} = SPDataLoad.SPData.maskMean_blank_Avg2{1};
                end
            elseif (d==length(dates) && strcmp(monkey, 'Gilligan'))
                SPData{s} = cat(3,SPData{s},SPDataLoad.SPData.maskMean_stim_Avg2{s}(:,1:70));
            else
                SPData{s} = cat(3,SPData{s}, SPDataLoad.SPData.maskMean_stim_Avg2{s});
                if (s==1)
                    SPData{length(stimSelect)+1} = cat(3,SPData{length(stimSelect)+1}, SPDataLoad.SPData.maskMean_blank_Avg2{1});
                end
            end
        end
        
    end
end
% SPData1_2 = load('S:\Lab\Skipper\All Data\Skipper_12_03_2020\Imaging\run00\Results\Superpixel\LP5HP250C0\SPData_2.mat');
% SPData2_2 = load('S:\Lab\Skipper\All Data\Skipper_12_04_2020\Imaging\run00\Results\Superpixel\LP5HP250C0\SPData_2.mat');
% SPData1 = cat(3,SPData1_2.SPData.maskMeanAll{1}, SPData2_2.SPData.maskMeanAll{1});
% SPData1 = smoothdata(SPData1,2);

maskFiles = dir(['S:\Lab\',monkey,'\Mapping\SPMasks\*.bmp']);
STDWeight = .5;
rmsVals = {};
for i = 1:length(maskFiles)
    lgnd{i} = ['Mask ',maskFiles(i).name(1:end-4)];
end

subplotSize = numSubplots(length(maskFiles));
colors = {'r', 'g','k'};
for s = 1:length(stimSelect)
    for m = 1:size(SPData{1},1)
        subplot(subplotSize(1), subplotSize(2),m);
        hold on;
        meanSP = smooth(mean(SPData{s}(m,:,:),3));
        plot(meanSP,colors{s},'linewidth',2);
        % plot(meanSP-(STDWeight.*smooth(std(SPData{s}(m,:,:),1,3))), 'r--', 'linewidth', 1);
        % plot(meanSP+(STDWeight.*smooth(std(SPData{s}(m,:,:),1,3))), 'r--', 'linewidth', 1);
        % rmsVals{s}(m,:) = rms(squeeze(SPData{s}(m,10:70,:)));
        rmsVals{s}(m,:) = rms(squeeze(SPData{s}(m,:,:)));
        if m==1
            legend([stims, 'Rest']);
        end
        title(['Mask ', lgnd{m}]);
        xlabel('Frame #');
        ylabel('R/dR %');
        set(gca,'fontsize',14);
        ylim([-.1 .1])
        xticks([0:10:70]);
        xticklabels([-1:1:6]);
    end
end
%%
mFig = figure();
for m = 1:size(SPData{1},1)
    rmValues = cellfun(@(a) a(m,:), rmsVals, 'UniformOutput', false);
    rmGroups = arrayfun(@(a,b) repmat(a,b,1), stimSelect, cellfun(@length, rmValues),'UniformOutput', false);
    rmValueArr = cell2mat(rmValues);
    rmGroupArr = reshape(cell2mat(rmGroups'),size(rmValueArr));
    [~,~,stats] = anova1(rmValueArr,rmGroupArr);
    hg1 = get(gca,'Children');
    figm = get(hg1,'Children');
    figure();
    allStats = multcompare(stats);
    sigInds = find(allStats(:,end)<0.05);
    comps = figure(3);
    plots = figure(4);
    statFig = figure(5);
    figure(mFig);
    aFig = subplot(subplotSize(1), subplotSize(2),m);
    copyobj(figm,aFig);
    hold on;
    for s = 1:length(sigInds)
        sigline([allStats(sigInds(s),1), allStats(sigInds(s),2)],'p<0.05',max(rmValueArr));
    end
    close([plots;statFig;comps]);
end
%%
minFig = figure();
for m = 1:size(SPData{1},1)
    %     [~,inds1] = min(squeeze(SPData{s}(m,20:end,:)));
    %     inds1 = inds1 + 20;
    minVals = cellfun(@(a) returnMinInd(a,m),SPData  ,'UniformOutput', false);
    minGroups = arrayfun(@(a,b) repmat(a,b,1), stimSelect,  cellfun(@length, minVals),'UniformOutput', false);
    minValArr = cell2mat(minVals);
    minGroupArr = reshape(cell2mat(minGroups'),size(minValArr));
    [~,~,stats] = anova1(minValArr,minGroupArr);
    hg1 = get(gca,'Children');
    figm = get(hg1,'Children');
    figure();
    allStats = multcompare(stats);
    %allStats = allStats([1,3],:);
    sigInds = find(allStats(:,end)<0.05);
    comps = figure(4);
    plots = figure(5);
    statFig = figure(6);
    figure(minFig);
    aFig = subplot(subplotSize(1), subplotSize(2),m);
    copyobj(figm,aFig);
    hold on;
    for s = 1:length(sigInds)
        sigline([allStats(sigInds(s),1), allStats(sigInds(s),2)],'p<0.05',max(minValArr));
    end
    close([plots;statFig;comps]);
end

function indices = returnMinInd(a,m)
[~,indices] = min(squeeze(a(m,:,:)));
end