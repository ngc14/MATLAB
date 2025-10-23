close all;
clearvars;
%% user input
animal = 'Gilligan';
RAD = 42; %28 = 0.5mm, 1.0mm^2; 42 = 0.75 mm, 2.2mm^2; 56 = 1mm, 3.4mm^2;
siteSize = 8;
REMOVE_NR = true;
USE_CORNERS = true;
SAVE = 0;
maxThresh = 40;
%% load green
green = imread(['S:\Lab\', animal, '\Mapping\M1_green_edited.bmp']);
xlimits = [1 size(green,2)];
ylimits = [1 size(green,1)];
%% load and parse excel sheet
excel = readcell(['S:\Lab\',animal, '\Mapping\', animal, '_MM_Sites.xlsx']);
headings = excel(1,:);

if REMOVE_NR
    NRi = logical(~strcmpi(excel(:,6),'None'));
else
    NRi = logical(~cellfun(@(m) all(ismissing(m)), excel(:,5)));
end
NRi(1) = false;
%NRi(cell2mat(excel(2:end,1))>288) = false;
% trunkSites = cellfun(@(a) strfind(a,'Trunk'), excel, 'UniformOutput', false);
% faceSites = cellfun(@(a) strfind(a,'Face'), excel, 'UniformOutput', false);
% [trunkSites,~] = find(~cellfun(@isempty,trunkSites));
% [faceSites,~] =find(~cellfun(@isempty,faceSites));
%
% for f = 1:length(faceSites)
%     strAll = excel{faceSites(f),4};
%     strAll = deleteWords(strAll, 'Face');
%     excel{faceSites(f),4} = strAll;
% end
%
% for t = 1:length(trunkSites)
%     strAll = excel{trunkSites(t),4};
%     strAll = deleteWords(strAll, 'Trunk');
%     excel{trunkSites(t),4} = strAll;
% end
excel = excel(NRi,:);
commas = find(cellfun(@(a) strcmp(a(end-1:end), ', '), excel(:,4)));
for m = 1:length(commas)
    excel{commas(m),4} = excel{commas(m),4}(1:end-2);
end

for i = 1:length(excel)
    siteNum(i,1) = excel{i,1};
    % x y
    siteLoc(i,:) = [excel{i,2} excel{i,3}];
    % resp strength
    switch excel{i,6}
        case 'Very Strong'
            respStrength(i,1) =  6;
        case 'Strong'
            respStrength(i,1) =  5;
        case 'Medium'
            respStrength(i,1) =  4;
        case 'Weak'
            respStrength(i,1) =  3;
        case 'Very Weak'
            respStrength(i,1) =  2;
        case 'None'
            respStrength(i,1) =  1;
    end
    if respStrength(i,1) ~= 1
        % receptive field
        tempStr = [',',excel{i,4},','];
        tempStr = tempStr(tempStr~=' ');
        commas = strfind(tempStr,',');
        for n = 1:(sum(tempStr==',')-1)
            tempTempStr = tempStr(commas(n)+1:commas(n+1)-1);
            parenth = [strfind(tempTempStr,'(') strfind(tempTempStr,')')];
            if ~isempty(strfind(tempTempStr(1:parenth(1)-1),'Digit'))
                Joint{i,n} = 'Digit';
            else
                Joint{i,n} = tempTempStr(1:parenth(1)-1);
            end
            Movement{i,n} = tempTempStr(parenth(1)+1:parenth(2)-1);
        end
    else % no response
        Joint{i,1} = 'None';
        Movement{i,1} = 'None';
    end
end
Thresh = nan(size(Movement));
for i = 1:size(excel,1)
    temp = sscanf(num2str(excel{i,5}),'%f,');
    Thresh(i,[1:length(temp)]) = temp;
end
[~, sortThresh] = sort(Thresh,2, 'ascend');
Thresh(Thresh>maxThresh) = maxThresh;
[Thresh,te] = discretize(Thresh,[0:5:maxThresh]);
%% encode primary RFs
[Joint_num, Movement_num] = deal(zeros(size(Joint)));
Joint_num(strcmpi(Joint,'None')) = 1;
Joint_num(strcmpi(Joint,'Digit')) = 2;
Joint_num(strcmpi(Joint,'Wrist')) = 3;
Joint_num(strcmpi(Joint,'Elbow') | strcmpi(Joint,'Forearm')) = 4;
Joint_num(strcmpi(Joint,'Shoulder')) = 5;
Joint_num(strcmpi(Joint,'Face')) = 6;
Joint_num(strcmpi(Joint,'Trunk')) = 7;
Movement_num(strcmpi(Movement,'None')) = 1;
Movement_num(strcmpi(Movement,'Flexion')) = 2;
Movement_num(strcmpi(Movement,'Extension')) = 3;
Movement_num(strcmpi(Movement,'Pronation')) = 4;
Movement_num(strcmpi(Movement,'Supination')) = 5;
Movement_num(strcmpi(Movement,'Abduction')) = 6;
Movement_num(strcmpi(Movement,'Adduction')) = 7;
Movement_num(strcmpi(Joint,'Face')) = 8;
Movement_num(strcmpi(Joint,'Trunk')) = 9;
Movement_num(strcmpi(Movement,'???')) = 10;

Joint_num_unique = zeros(size(Joint_num));
for i = 1:size(Joint_num,1)
    temp = unique(Joint_num(i,:));
    Joint_num_unique(i,[1:length(temp(temp~=0))]) = temp(temp~=0);
end
%% create radial mask around every site
imDims = size(green);
circle = fspecial('disk',RAD);
siteMask = zeros(imDims(1),imDims(2),length(siteNum));
for n = 1:length(siteNum)
    tempCircle = zeros(imDims(1)+400,imDims(2)+400);
    tempCircle((siteLoc(n,2)-RAD+200):(siteLoc(n,2)+RAD+200),(siteLoc(n,1)-RAD+200):(siteLoc(n,1)+RAD+200)) = circle;
    tempCircle = tempCircle(200:end-201,200:end-201);
    siteMask(:,:,n) = tempCircle>0;
end
%% striped image
stripeWidth = 6;
stripeIm = ones(imDims(1),imDims(2),5);
x = 1:stripeWidth*2:size(green,1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,2) = 2;
end
x = 1:stripeWidth*3:size(green,1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,3) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,3) = 3;
end
x = 1:stripeWidth*4:size(green,1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,4) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,4) = 3;
    stripeIm((x(i)+stripeWidth*2+1):(x(i)+stripeWidth*3),:,4) = 4;
end
x = 1:stripeWidth*5:size(green,1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,4) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,5) = 3;
    stripeIm((x(i)+stripeWidth*2+1):(x(i)+stripeWidth*3),:,5) = 4;
    stripeIm((x(i)+stripeWidth*3+1):(x(i)+stripeWidth*4),:,5) = 5;
end
%% corner sites
if USE_CORNERS
    corners(1,:) = [0 size(green,2)];
    corners(2,:) = [size(green,1) 0];
    corners(3,:) = [0 0];
    corners(4,:) = [size(green,1) size(green,2)];
end
%% voronoi
if ~USE_CORNERS
    [v, c] = voronoin(siteLoc);
else
    [v, c] = voronoin([siteLoc; corners]);
end
cmap_Joint = [0 0 0; 1 0 0; 1 .5 0; 0 1 0; 0 1 1;...
    155/255 9/255 144/255; .5 .25 .1];
cmap_Joint=[0 0 0; .35 .35 .35; .35 .35 .35; .75 .75 .75; .75 .75 .75; 155/255 9/255 144/255; .5 .25 .1];
%% construct SM image
MM_image(:,:,1) = zeros(imDims(1),imDims(2));
MM_image(:,:,2) = zeros(imDims(1),imDims(2));
MM_image(:,:,3) = zeros(imDims(1),imDims(2));
[MM_image_flat, MM_respStrength_image_flat, MM_movement_image_flat, ...
    MM_thresh_image_flat] = deal(NaN(imDims(1),imDims(2)));

[MM_joint_image_flat2, MM_movement_image_flat2, MM_thresh_image_flat2] = ...
    deal(NaN(imDims(1),imDims(2),size(Joint_num,2)));

dispstat('Creating MM image... 0%','init');
if USE_CORNERS
    numSites = length(c)-4;
else
    numSites = length(c);
end
for i = 1:numSites
    reorderTemp = sortThresh(i,:);
    Joint_num(i,:) = Joint_num(i,reorderTemp);
    Movement_num(i,:) = Movement_num(i,reorderTemp);
    bestThresh = Thresh(i,reorderTemp(1));
    for j = 1:sum(Joint_num(i,:)~=0)-1
        if(Thresh(i,j+1)==0)
            Thresh(i,j+1) = Thresh(i,j);
        end
    end
    temp = length(unique(Joint_num(i,Thresh(i,:)<=bestThresh+bestThresh*.1 & Joint_num(i,:)~=0)));
    %temp = sum(Joint_num(i,:)~=0);
    for x = 1:imDims(1)
        for y = 1:imDims(2)
            if siteMask(x,y,i) == 1
                [in, on] = inpolygon(x,y,v(c{i},2),v(c{i},1));
                if ((in || on) && temp>0)
                    MM_image(x,y,1) = cmap_Joint(Joint_num(i,stripeIm(x,y,temp)),1);
                    MM_image(x,y,2) = cmap_Joint(Joint_num(i,stripeIm(x,y,temp)),2);
                    MM_image(x,y,3) = cmap_Joint(Joint_num(i,stripeIm(x,y,temp)),3);
                    MM_image_flat(x,y) = Joint_num(i,stripeIm(x,y,temp));
                    MM_movement_image_flat(x,y) = Movement_num(i,stripeIm(x,y,temp));
                    MM_respStrength_image_flat(x,y) = respStrength(i);
                    MM_thresh_image_flat(x,y) = Thresh(i,stripeIm(x,y,temp));
                    for t = 1:temp
                        MM_joint_image_flat2(x,y,reorderTemp(t)) = Joint_num(i,t);
                        MM_movement_image_flat2(x,y,reorderTemp(t)) = Movement_num(i,t);
                        MM_thresh_image_flat2(x,y,reorderTemp(t)) = Thresh(i,t);
                    end
                end
            end
        end
    end
    percDone = (i)/(length(c)-4);
    percDone = round(percDone*10000)/100;
    dispstat(['Creating MM image... ',num2str(percDone),'%']);
end
dispstat('Creating MM image... 100%');
%% plot voronoi
figure('color',[1 1 1],'position',[150 150 1000 800]);
imagesc(green);
axis image off; hold on;
imagesc(MM_image,'alphadata',~isnan(MM_image_flat));
title([animal,' - Voronoi RF Map']);
xlim(xlimits);
ylim(ylimits);
%% plot response strength
figure('color',[1 1 1],'position',[200 150 1000 800]);
cmap_respStrength = [0 0 0; jet(5)];
cmap_cutProp = [0 0 0; 0 0 0; 1 1 1; 0.5 0.5 0.5];
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_respStrength); caxis([0 6]);
for i = 1:length(siteNum)
    t = linspace(0,2*pi);
    patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
        [1 1 1],'edgecolor','none');
    patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
        cmap_respStrength(respStrength(i),:),'edgecolor','none');
end
set(h,'ticks',(1:6)-0.5,'ticklabels',{'NR','VWR','WR','MR','GR','VGR'});
title([animal,' - Movement Strength']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% plot joint
figure('color',[1 1 1],'position',[250 150 1000 800]);
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_Joint); caxis([0 size(cmap_Joint,1)]);
for i = 1:length(siteNum)
    temp = sum(Joint_num(i,:)~=0);
    for n = 1:temp
        t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
        patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
            [1 1 1],'edgecolor','none');
        patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
            cmap_Joint(Joint_num(i,n),:),'edgecolor','none');
    end
end
set(h,'ticks',(1:size(cmap_Joint,1))-0.5,'ticklabels',{'NR','Digit','Wrist','Elbow','Shoulder','Face','Trunk'});
title([animal,' - Joint']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% plot movement
figure('color',[1 1 1],'position',[250 150 1000 800]);
% cmap_Movement = [0 0.8 1; 0 0 1; 1 0 0; 0.6 0 0; 0.1 1 0; 0.05 0.6 0];
cmap_Movement = distinguishable_colors(10); cmap_Movement(5,:) = [0.6 0 0.6];
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_Movement(1:6,:)); caxis([0 6]);
for i = 1:length(siteNum)
    temp = sum(Movement_num(i,:)~=0);
    for n = 1:temp
        if Movement_num(i,n)>=2 && Movement_num(i,n)<=7
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_Movement(Movement_num(i,n)-1,:),'edgecolor','none');
        end
    end
end
set(h,'ticks',(1:size(cmap_Movement,1))-0.5,'ticklabels',{'Flexion','Extension',...
    'Pronation','Supination','Abduction','Adduction'});
title([animal,' - Movement Type']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% plot Thresh
figure('color',[1 1 1],'position',[250 150 1000 800]);
cmap_Thresh = hot(max(Thresh,[],'all'));
cmap_Thresh(end,:) = [0.5 0.5 0.5];
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_Thresh); caxis([0 size(cmap_Thresh,1)]);
for i = 1:length(siteNum)
    temp = sum(~isnan(Thresh(i,:)));
    for n = 1:temp
        t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
        patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
            cmap_Thresh(Thresh(i,n),:),'edgecolor','none');
    end
end
title([animal,' - Threshold']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% plot digit and wrist w/ movement
figure('color',[1 1 1],'position',[250 150 1400 800]);
ax1 = subplot(1,2,1);
imagesc(green); axis image off; hold on;
cmap_digit_movement = [0 51/255 255/255; 0 230/255 255/255];
cmap_wrist_movement = [0.5 0 0; 1 0 0];
h = colorbar; colormap(ax1,[cmap_digit_movement; cmap_wrist_movement]); caxis([0 4]);
for i = 1:length(siteNum)
    temp = sum(Joint_num(i,:)~=0);
    for n = 1:temp
        if Joint_num(i,n)==2 && (Movement_num(i,n)==2 || Movement_num(i,n)==3)
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+3)*cos(t) 0],siteLoc(i,2)+[(siteSize+3)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_digit_movement(Movement_num(i,n)-1,:),'edgecolor','none');
        end
        if Joint_num(i,n)==3 && (Movement_num(i,n)==2 || Movement_num(i,n)==3)
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+3)*cos(t) 0],siteLoc(i,2)+[(siteSize+3)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_wrist_movement(Movement_num(i,n)-1,:),'edgecolor','none');
        end
    end
end
set(h,'ticks',(1:4)-0.5,'ticklabels',{'D Flex','D Extend','W Flex','W Extend'});
title([animal,' - Digit/Wrist Movement']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% plot shoulder w/ movement
ax2 = subplot(1,2,2);
imagesc(green); axis image off; hold on;
cmap_shoulder_movement = [0.2 0.2 0.2; 0.7 0.7 0.7];
cmap_elbow_movement = [0 0.5 0; 0 1 0];
h = colorbar; colormap(ax2,[cmap_elbow_movement; cmap_shoulder_movement]); caxis([0 4]);
for i = 1:length(siteNum)
    temp = sum(Joint_num(i,:)~=0);
    for n = 1:temp
        if Joint_num(i,n)==5 && (Movement_num(i,n)==2 || Movement_num(i,n)==3)
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+3)*cos(t) 0],siteLoc(i,2)+[(siteSize+3)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_shoulder_movement(Movement_num(i,n)-1,:),'edgecolor','none');
        elseif Joint_num(i,n)==4 && (Movement_num(i,n)==2 || Movement_num(i,n)==3)
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+3)*cos(t) 0],siteLoc(i,2)+[(siteSize+3)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_elbow_movement(Movement_num(i,n)-1,:),'edgecolor','none');
        end
    end
end
set(h,'ticks',(1:4)-0.5,'ticklabels',{'E Flex','E Extend','S Flex','S Extend'});
title([animal,' - Elbow/Shoulder Movement']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);
%% SAVE
if SAVE
    if ~exist(['S:\Lab\',animal,'\Mapping\Motor Maps V3'],'dir')
        mkdir(['S:\Lab\',animal,'\Mapping\Motor Maps V3']);
    end
    cd(['S:\Lab\',animal,'\Mapping\Motor Maps V3']);
    if REMOVE_NR
        NRstr = '_RNR_Recorded_MIN';
    else
        NRstr = '_All';
    end

    filename = [animal,' MM Joint ',num2str(RAD),'px',NRstr,'.png'];
    imwrite(ind2rgb(MM_image_flat,cmap_Joint),filename,'alpha',double(~isnan(MM_image_flat)));
    exportgraphics(gcf,[filename(1:end-4),'.eps'],'ContentType','vector');

    %     filename = [animal,' MM respStrength ',num2str(RAD),'px',NRstr,'.png'];
    %     imwrite(ind2rgb(MM_respStrength_image_flat,cmap_respStrength),filename,'alpha',double(~isnan(MM_respStrength_image_flat)));
    %     filename = [animal,' MM Movement ',num2str(RAD),'px',NRstr,'.png'];
    %     imwrite(ind2rgb(MM_movement_image_flat,cmap_Movement),filename,'alpha',double(~isnan(MM_movement_image_flat)));
    %     filename = [animal,' MM Threshold ',num2str(RAD),'px',NRstr,'.png'];
    %     imwrite(ind2rgb(MM_thresh_image_flat,cmap_Thresh),filename,'alpha',double(~isnan(MM_thresh_image_flat)));
    %     %MM_numbered_joints = ind2rgb(MM_image_flat,cmap_Joint);
    %     %MM_numbered_joints = cellfun(@(a,b) insertText(MM_numbered_joints,...
    %     %a(2),a(1),num2str(b),'AnchorPoint','center');
    %     %imwrite(MM_numbered_joints, cmap_Joint, [animal, ' MM Numbered ', num2str(RAD), 'px', NRstr,'.png'], 'alpha', double(~isnan(MM_image_flat)));
    %     filename = [animal,' MM All ',num2str(RAD),'px',NRstr,'.mat'];
    %     save(filename,'MM_image_flat','MM_movement_image_flat','MM_thresh_image_flat','MM_respStrength_image_flat',...
    %         'cmap_Joint', 'cmap_Movement', 'cmap_Thresh', 'cmap_respStrength',...
    %         'siteLoc', 'siteNum', 'v', 'c','siteMask');
end

function deleted = deleteWords(fullStr, deleteWord)
deleted = fullStr;
startInd = regexp(fullStr, deleteWord);
endInd = startInd + regexp(fullStr(startInd:end), ', ', 'end')-1;
if(length(startInd)>1)
    deleted = [fullStr(1:startInd(1)-1), deleteWords(fullStr(startInd(2):end), deleteWord)];
else
    if(isempty(endInd))
        endInd = length(fullStr);
    end
    deleted(startInd:endInd) = [];
end
end