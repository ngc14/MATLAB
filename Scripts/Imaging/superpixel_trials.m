monkey = 'Gilligan';                   % animal name
if(strcmp(monkey,'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
    mm = MotorMapping(42);
elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};    
    mm = MotorMapping(56);
end
saveVar = false;

conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
maskNames = ["FOV", "M1Arm","M1Hand","M1Forelimb","PMd"];
LPk = 5;
HPk = 250;
%%
[datesTable, masksR, ~] = getMonkeyInfo('S:\Lab\',monkey,["M1", "PMd"],false);
refMask = masks{1};
siteMask = getVoronoiMask(datesTable,mm,refMask,["Arm","Hand"]);
%
allMMs = {};
allMMs{end+1} = siteMask & ~refMask;
allMMs{end+1} = siteMask & ~refMask & masksR{2} & any(siteMasks(:,:,strcmp(simpRep, "Arm")),3);
allMMs{end+1} = siteMask & ~refMask & masksR{2} & any(siteMasks(:,:,strcmp(simpRep, "Hand")),3);
allMMs{end+1} = siteMask & ~refMask & masksR{2};
allMMs{end+1} = siteMask & ~refMask & masksR{3};
mks = dir(['S:\Lab\ngc14\Working\',monkey,'\*.bmp']);
maskNames = [maskNames, mks.name];
for m = 1:length(mks)
    masks = im2uint8(imread([mks(m).folder,'\',mks(m).name]))>50;
    masks = masks(:,:,1)./masks(:,:,1);
    sizeDiffX = (size(refMask,1) - size(masks,1))/2;
    sizeDiffY = (size(refMask,2) - size(masks,2))/2;
    masks = masks(1+ceil(abs(sizeDiffX)):end+ceil(sizeDiffX),:);
    masks = masks(:,1+ceil(abs(sizeDiffY)):end+ceil(sizeDiffY));
    masks(1:2,:) = NaN;
    masks(:,1:2) = NaN;
    masks(end-1:end,:) = NaN;
    masks(:,end-1:end) = NaN;
    masks = double(masks);
    masks(masks==0) = NaN;
    allMMs{end+1} = masks;
end
allMMs = cellfun(@(m) double(m./m), allMMs, 'UniformOutput',false);
clear siteMask siteMasks tempCircle masksR vCells verticies siteLocation simpRep datesTable repThreshInd masks
%%
imA = cell(1,length(conds));
for cn =1:length(conds)
    ti = datetime('now');
    trialAvgFile = ['S:\Lab\ngc14\Working\',monkey,'\',conds{cn},'\TrialAVG.mat'];
    if(exist(trialAvgFile,'file'))
        imAll = load(trialAvgFile,'frames');
        imA{cn} = imAll.frames;
    else
        imAll= loadAllFrames(monkey,dates,runs,LPk,HPk,conds{cn},false,1:70);
        imA{cn} = cell2mat(reshape(cellfun(@gather, imAll, 'UniformOutput', false),1,1,[]));
    end
    %imsSP{f} = cellfun(@(m) nanmean(allTrials.*repmat(m,[1,1,size(allTrials,3)]),'all'), allMMs, 'UniformOutput',true);%
    %imsVar{f} = cellfun(@(m) nanstd(allTrials.*repmat(m,[1,1,size(allTrials,3)]),0,'all'), allMMs, 'UniformOutput', false);
    disp(strjoin(["Loaded ", conds{cn}," trials. (", num2str(minutes(datetime('now')-ti))," min)"],''));
end

%save(['S:\Lab\ngc14\Working\', monkey,'\Delay_SPs\12_18_2018.mat'],'imA');
%%
close all;
for m = 1:length(allMMs)
figure();
plot(movmean(cell2mat(cellfun(@(cf) squeeze(nanmean(cf.*repmat(allMMs{m},[1 1 size(m,3)]),[1 2])), imA, 'UniformOutput',false)),5,1));
hold on;
ylim([-0.25 0.25])
title(maskNames(m))
legend(conds)
end
%%
savePlots(conds,allMMs,imA,maskNames,'mean',monkey);
%savePlots(conds,allMMs,imA,maskNames,'median',monkey);
function savePlots(conds,allMMs,imA,maskNames,avgType,monkey)
colorsP = distinguishable_colors(length(maskNames));
figure('Units','normalized','Position',[0 0 1 1])
smoothROIs = cell(1,length(imA));

for c = 1:length(imA)
    subplot(2,2,c)
    hold on;
    title(conds{c})
    if(strcmp(avgType,'mean'))
        smoothROIs{c} = cell2mat(cellfun(@(a) squeeze(mean(a.*imA{c},[1 2],'omitnan')),allMMs, 'UniformOutput', false));
    elseif(strcmp(avgType,'median'))
        smoothROIs{c} = cell2mat(cellfun(@(a) squeeze(median(a.*imA{c},[1 2],'omitnan')),allMMs,'UniformOutput',false));
    end
    plot(movmean(smoothROIs{c},3,1),'LineWidth',2)
    lns = get(gca,'Children');
    cellfun(@(lh,ch) set(lh,'Color',ch), num2cell(lns),num2cell(colorsP,2))
    ylim([-.263 .25])
end
legend(maskNames)
exportgraphics(gcf,['S:\Lab\ngc14\Working\',monkey,'\ROI_',avgType,'.eps'],'ContentType','vector','BackgroundColor','none');
exportgraphics(gcf,['S:\Lab\ngc14\Working\',monkey,'\ROI_',avgType,'.png'],'ContentType','vector','BackgroundColor','none');
figure('Units','normalized','Position',[0 0 1 1]);
for m = 1:length(maskNames)
    subplot(2,ceil(length(maskNames)/2),m)
    hold on;
    title(maskNames(m));
    plot(movmean(cell2mat(cellfun(@(r) r(:,m), smoothROIs, 'UniformOutput',false)),3,1), 'LineWidth', 2);
    ylim([-.25 .25])
end
legend(conds);
exportgraphics(gcf,['S:\Lab\ngc14\Working\',monkey,'\ROIMasks_',avgType,'.eps'],'ContentType','vector','BackgroundColor','none');
exportgraphics(gcf,['S:\Lab\ngc14\Working\',monkey,'\ROIMasks_',avgType,'.png'],'ContentType','vector','BackgroundColor','none');
end