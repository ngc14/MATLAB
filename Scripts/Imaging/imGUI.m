close all;
monkey = 'Skipper';
if(strcmp(monkey,'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
end
conds = {'[ExtraSmallSphere]','[LargeSphere]', '[Photocell]','[Rest]'};
angles = linspace(0, 2*pi, 1000);
radius = 20;
dColors = distinguishable_colors(20);
lStyle = {'-', '-.', ':',':' };
refMask = imread(['S:\Lab\',monkey,'\Mapping\clean_mask_filled'],'bmp')>180;
refMask = refMask(:,:,1);
[nx,ny] = ndgrid(1:size(refMask,2),1:size(refMask,1));
activityC = cellfun(@(cn) logical(imread(['S:\Lab\ngc14\Working\',monkey,'\',cn,...
    '\tTests\NaN\HP250\tTest_Nomask48  49  50  51  52_p-0.0001.bmp'])), conds, 'UniformOutput', false);

condFrames = cellfun(@(cn) load(['S:\Lab\ngc14\Working\',monkey,'\',cn,'\TrialAVG.mat'],'frames'), conds, 'UniformOutput', false);
condFrames = cellfun(@(cn) cn.frames, condFrames, 'Uniformoutput', false);
delayFrames = load("S:\Lab\ngc14\Working\"+monkey+"\Delay_SPs\Average_Trial_Frames.mat");
delayFrames = delayFrames.imF;
%%
sm = im2uint8(imread(['S:\Lab\',monkey,'\Mapping\Simplified_MM-01.png']));
sm = sm(:,1:end-(size(sm,2)-size(refMask,2)),:);
sm(repmat(sum(sm,3)==0,1,1,size(sm,3))) = 255;
sm(repmat(~refMask,1,1,size(sm,3))) = 175;
rF = figure('Units','Normalized');
hold on;
sf1 = subplot(2,1,1);
hold on;
sf2 = subplot(2,1,2);
hold on;
axF = figure('Units','Normalized');
offset =  get(axF,'InnerPosition')- get(axF,'OuterPosition');
offset = (offset(1:2)+[0 -1]);

figure(axF);
csm = sm;
colorInds = reshape(randi(size(csm,3)+1,1,size(csm,3)*length(conds))-1,[],size(csm,3));
for c = 1:length(conds)
    currColor = colorInds(c,:);
    currColor = currColor(currColor~=0);
    for l = 1:length(currColor)
        ch = csm(:,:,currColor(l));
        ch(refMask & activityC{c}) = max(0,min(125,ch(refMask & activityC{c}) + 30*sign(mod(c,2)-.5)));
        csm(:,:,currColor(l)) = ch;
    end

end
imshow(csm);

rInd = 0;
continueROI = true;
waitfor(axF,'CurrentPoint');
while continueROI
    center = get(axF, 'CurrentPoint');
    %center = center + offset;
    rInd = rInd + 1;
    center = [center(1) 1-center(end)].*size(refMask,[2,1]);

    centerx = center(2); centery = center(1);
    x = radius * cos(angles) + centerx;
    y = radius * sin(angles) + centery;

    inon = inpolygon(nx(:),ny(:),x,y);
    lin = sub2ind(size(refMask),nx(inon), ny(inon));

    ROIIm = false(size(refMask));
    ROIIm(lin) = true;
    ROIIm = ROIIm & refMask;
    
    figure(axF)
    if(mod(rInd,2))
        hold off;
        imshow(csm);
        hold on;
        cla(sf1);
        cla(sf2);
    end
    viscircles(gca(axF),center,radius,'Color',dColors(mod(rInd,20)+~mod(rInd,20),:));
    axes(sf1);
    cellfun(@(cf,ls) plot(movmean(squeeze(mean(cf.*repmat(ROIIm./ROIIm,[1 1 size(cf,3)]),...
        [1 2],'omitnan')),3),'Color',dColors(mod(rInd,20)+~mod(rInd,20),:),'LineStyle',ls), condFrames, lStyle,'UniformOutput', false);
    xticks(0:10:70);
    xticklabels(arrayfun(@num2str,-1:1:6,'UniformOutput',false));
    ylim([-0.2 0.2]);
    axes(sf2);
    cellfun(@(cs,ls) plot(movmean(squeeze(mean(cs.*repmat(ROIIm./ROIIm,[1 1 size(cs,3)]),...
        [1 2],'omitnan')),3),'Color',dColors(mod(rInd,20)+~mod(rInd,20),:),'LineStyle',ls), delayFrames, lStyle(1:length(delayFrames)),'UniformOutput', false);
    ylim([-0.3 0.3]);

    xticks(0:10:90);
    xticklabels(arrayfun(@num2str, 0:9,'UniformOutput',false));
    figure(axF);
    waitfor(axF,'CurrentPoint');
    continueROI = ~strcmp(get(axF,'CurrentCharacter'),'q');
end
