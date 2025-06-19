close all;
clear all;
monkey = 'Gilligan';
areaNames = {'Active M1','Active PMd', 'Inactive M1', 'Inactive PMd', 'All M1', 'All PMd'};
MMNames = {'Arm', 'Hand', 'Face', 'Trunk'};
maskInd = 1;
for a = 1:length(areaNames)
    for m = 1:length(MMNames)
        maskNames{maskInd} = [areaNames{a}, ': ', MMNames{m}];
        maskInd = maskInd + 1;
    end
end
% clear maskNames;
% maskNames{1} = 'Active M1 Arm';
% maskNames{2} = 'Active M1 Hand';
% maskNames{3} = 'Active PMd Arm';
% maskNames{4} = 'Active PMd Hand';
% maskNames{5} = 'Active Arm';
% maskNames{6} = 'Active Hand';
% maskNames{7} = 'Active M1';
% maskNames{8} = 'Active PMd';
if strcmp(monkey, 'Gilligan')
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
else
dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0]};
end
stimNames = {'ExtraSmallSphere','LargeSphere', 'Photocell'};

refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'], 'bmp')>200;
consecutiveFrames = 7;
windowSmoothSize = 5;
stdThresh = 0;
areaThresh = 0.25;

colors = {'r', 'g', 'b', 'k'};
for c = 1:length(stimNames)+1
    toPlot = {};
    for d = 1:length(dates)
        if(strcmp(monkey,'Gilligan'))
            sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP550C0\SPData_Best.mat']);
        else
            sp=load(['S:\Lab\',monkey,'\All Data\',monkey,'_',dates{d},'\Imaging\run0',num2str(runs{d}),'\Results\Superpixel\LP5HP250C0\SPData_Best.mat']);
        end
        sp = sp.SPData;
        if(strcmp(monkey, 'Skipper') && c<length(stimNames)+1 && ((d==1 && c>1) || (d==2)))
            if(d==1)
                toPlot{d} = [];
            elseif(d==2)
                if(c==1)
                    toPlot{d} =[];
                else
                    toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c-1},2,'movmean',windowSmoothSize);
                end
            end
        elseif(c==length(stimNames)+1)
            toPlot{d} = smoothdata(sp.maskMean_blank_Avg2{1},2,'movmean',windowSmoothSize);
        else
            toPlot{d} = smoothdata(sp.maskMean_stim_Avg2{c},2,'movmean',windowSmoothSize);
        end
        if(strcmp(monkey, 'Gilligan') && d == length(dates))
            toPlot{d} = toPlot{d}(:,1:70);
        end
        
        allMasksMeans = vertcat(sp.maskMean_stim_Avg2{:});
        allMasksMeans = [allMasksMeans; sp.maskMean_blank_Avg2{1}];
        meanVals = nanmean(allMasksMeans(:));
        stdVals = nanstd(allMasksMeans(:));
        %toPlot{d} = (toPlot{d} - meanVals) ./ stdVals;
        
        if(strcmp(monkey, 'Skipper') && d ==length(dates))
            toPlot = toPlot(~cellfun(@isempty,toPlot));
        end
       
    end
    clusterPath = ['S:\Lab\',monkey,'\Mapping\Clustering\15px_Grid\XCorr_Elbow\Biased\Norm_No_Clip\'];
    if c==1
        maskFiles = openfig([clusterPath, 'ExtraSmallSphere_Grid.fig'], 'invisible');
    elseif c==2
        maskFiles = openfig([clusterPath, 'LargeSphere_Grid.fig'], 'invisible');
        
    elseif c==3
        maskFiles = openfig([clusterPath, 'Photocell_Grid.fig'], 'invisible');
        
    else
        maskFiles = openfig([clusterPath, 'Rest_Grid.fig'], 'invisible');
        
    end
    maskFiles = unique(maskFiles.Children(1).Children(2).CData);
    maskFiles(maskFiles==0) = [];
    colorP = colormap(flipud(jet(length(maskFiles))));
    
%     maskNames = {maskFiles.name};
%     [maskNames, sortIdx] = natsort(maskNames);
    sortIdx = 1:length(maskFiles);
    
    %toPlot = cellfun(@(a) normalize(a), toPlot, 'UniformOutput', false);
    

    allMasks = reshape(cell2mat(toPlot),[size(toPlot{1},1), size(toPlot{1},2), length(toPlot)]);
    subplotSize = numSubplots(length(maskFiles));
    allMasks = allMasks(sortIdx,:,:);
        
    figure();
    for m = 1:size(allMasks,1)
        masks =  nanmean(allMasks,3);
        errs = nanstd(allMasks,0,3);
        timecourse = masks(m,:);
        subplot(subplotSize(1), subplotSize(2),m);
        hold on;
        plot(timecourse, 'Color',colorP(m,:));
        shadedErrorBar(1:length(errs(m,:)),timecourse,errs(m,:),'lineprops',{'Color', [colorP(m,:)]});
        ylim([-.1 .1]);
        if(m==1)
            ytickformat('%.2f');
        else
            yticklabels([]);
        end
        xticks([0,10,20,30,40,50,60,70]);
        set(gca,'XTickLabel',arrayfun(@(a) sprintf('%.1f',a), [-1,0,1,2,3,4,5,6], 'UniformOutput', 0));
        xlim([0, 71]);

        %if(c==length(stimNames)+1)
        %title(maskNames{m});
        %             if(m==1)
        %                 legend([{'Extra Small Sphere'}, {'Large Sphere'}, {'Photocell'}, {'Rest'}]);
        %             end
        % end
    end
    
    [~,allMins] = min(allMasks,[],2);
    for d = 1:size(allMasks,3)
        for m = 1:size(allMasks,1)
            diffVals = find(diff(allMasks(m,:,d))<0)+1;
            areas = cumsum(abs(allMasks(m,diffVals,d)));
            areas = trapz(abs(allMasks(m,diffVals(ismember(diffVals,1:14)),d)));
            %areas = [];
            minValInd = strfind(diff(allMasks(m,:,d))<0, ones(1, consecutiveFrames));
            minValInd = [];
            if(~isempty(minValInd) || ~isempty(areas))
                % minValInd = diffVals(areas > areaThresh*areas(end));
                % allMins(m,d) = minValInd(1);
                %allMins(m,d) = areas;
            else
                allMins(m,d) = NaN;
            end
        end
    end
    %     minVals(c,:) = nanmean(squeeze(allMins)');
    %     errors(c,:) = minVals(c,:)./nanstd(squeeze(allMins)');
    %     allVals{c} = squeeze(allMins)';
end
%%
for c = 1:length(stimNames)
    figure();
    % plotOrg(1,1:2) = minVals(c,1:2);
    % plotOrg(2,1:2) = minVals(c,3:4);
    % plotErrors(1,1:2) = errors(c,1:2);
    % plotErrors(2,1:2) = errors(c,3:4);
    % hbar = bar(plotOrg);
    %
    % for p = 1:size(plotOrg,2)
    %     ctr(p,:) = bsxfun(@plus, hbar(1).XData, hbar(p).XOffset);
    % end
    % hold on
    % errorbar(ctr', plotOrg, plotErrors, '.k')
    % title(stimNames{c});
    % xticklabels({'M1', 'PMd'});
    % legend({'Arm','Hand'});
    
    %val = allVals{c};
    %boxplot(val,'Notch','on', 'labels', maskNames,'Symbol','');
    %ylim([-.2 1]);
    title(stimNames{c});
    %saveas(gcf, ['S:\Lab\',monkey,'\Mapping\SPMasks\Consecutive Frames\',stimNames{c},'_Box.fig']);
end
