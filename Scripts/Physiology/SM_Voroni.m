% function [] = unbiasedMM_V3()

close all;
clearvars;

%% user input
animal = 'Skipper';

RAD = 42; %28 = 0.5mm radius, 56 = 1mm radius
siteSize = 8;
USE_CORNERS = true;
SM_Grid = false;
SAVE = true;
xlimits = [10 1300];
ylimits = [10 1070];
maxThresh = 35;
%% load green
%green = imread(['S:\Lab\', animal, '\Mapping\clean_mask_S1_filled.bmp']);
green = imread(['S:\Lab\', animal, '\Mapping\M1_green_edited.bmp']);

%% load and parse excel sheet
[excel.num, excel.txt, excel.raw] = xlsread(['S:\Lab\',...
    animal, '\Mapping\', animal, '_SM_Sites.xlsx']);
respStrength = [];

for i = 1:size(excel.num,1)
    siteNum(i,1) = excel.num(i,1);
    siteLoc(i,:) = [excel.num(i,2) excel.num(i,3)];
    
    responses = strsplit(excel.raw{i+1,5}, ',');
    responses = cellfun(@assignRS, responses);
    col = 1;
    for r = 1:length(responses)
        respStrength(i,col) = responses(r);
        col = col+1;
    end
    
    if any(respStrength(i,:)>0)
        % receptive field
        tempStr = [',',excel.raw{i+1,4},','];
        tempStr = tempStr(tempStr~=' ');
        commas = strfind(tempStr,',');
        
        for n = 1:(sum(tempStr==',')-1)
            tempTempStr = tempStr(commas(n)+1:commas(n+1)-1);
            parenth = [strfind(tempTempStr,'(') strfind(tempTempStr,')')];
            Joint{i,n} = tempTempStr(1:parenth(1)-1);
            Response{i,n} = tempTempStr(parenth(1)+1:parenth(2)-1);
        end
        
    else % no response
        Joint{i,1} = 'None';
        Response{i,1} = 'None';
    end
end

[~, sortResp] = sort(respStrength,2, 'descend');
%% encode primary RFs
Joint_num = zeros(size(Joint));
Location_num = zeros(size(Joint));

Joint_num(strcmpi(Joint,'None')) = 1;
Joint_num(strcmpi(Joint,'Digit1')) = 2;
Joint_num(strcmpi(Joint,'Digit2')) = 3;
Joint_num(strcmpi(Joint,'Digit3')) = 4;
Joint_num(strcmpi(Joint,'Digit4')) = 5;
Joint_num(strcmpi(Joint,'Digit5')) = 6;
Joint_num(strcmpi(Joint,'Pad1') | strcmpi(Joint, 'Pad2') |...
    strcmpi(Joint,'Pad3') | strcmpi(Joint, 'Thenar') | strcmpi(Joint, 'Hypothenar') ) = 7;
Joint_num(strcmpi(Joint,'Forearm')) = 8;
Joint_num(strcmpi(Joint,'Face')) = 9;

Location_num = cellfun(@assignD, Response, 'UniformOutput', false);

Joint_num_unique = zeros(size(Joint_num));
for i = 1:size(Joint_num,1)
    temp = unique(Joint_num(i,:));
    temp = temp(temp~=0);
    Joint_num_unique(i,[1:length(temp)]) = temp;
end

Location_num_unique = zeros(size(Joint_num));
for i = 1:size(Joint_num,1)
    temp = unique([Location_num{i,:}]);
    temp = temp(temp~=0);
    Location_num_unique(i,[1:length(temp)]) = temp;
end

%% create radial mask around every site
imDims = size(green);
circle = fspecial('disk',RAD);
siteMask = zeros(imDims(1),imDims(2),length(siteNum));
for n = 1:length(siteNum)
    
    tempCircle = zeros(imDims(1)+400,imDims(2)+400);
    tempCircle((siteLoc(n,2)-RAD+200):(siteLoc(n,2)+RAD+200),(siteLoc(n,1)-RAD+200):(siteLoc(n,1)+RAD+200)) = circle;
    disp('');
    
    %     if (size(tempCircle,1) > imDims(1)) || (size(tempCircle,2) > imDims(2))
    %         tempCircle = tempCircle(1:imDims(1),1:imDims(2));
    %     end
    tempCircle = tempCircle(200:end-201,200:end-201);
    
    siteMask(:,:,n) = tempCircle>0;
end

%% striped image
stripeWidth = 5;
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
    
    %     xrange = [min(siteLoc(:,1))-35 max(siteLoc(:,1))+35];
    %     yrange = [min(siteLoc(:,2))-35 max(siteLoc(:,2))+35];
    %     corners(1,:) = [xrange(1) yrange(2)];
    %     corners(2,:) = [xrange(2) yrange(1)];
    %     corners(3,:) = [xrange(1) yrange(1)];
    %     corners(4,:) = [xrange(2) yrange(2)];
end


%% voronoi
if ~USE_CORNERS
    [v, c] = voronoin(siteLoc);
else
    [v, c] = voronoin([siteLoc; corners]);
end

cmap_Joint = [.8 .8 .8; 1 0 0; 1 165/255 0; 1 1 0; 0 1 0; 0 1 1;...
    0 0 1; 139/255 0 139/255;  0 0 0];

% figure;
% for i = 1:length(c)-4
%     if all(c{i}~=1)   % If at least one of the indices is 1,
%                       % then it is an open region and we can't
%                       % patch that.
%
%         patch(v(c{i},1),v(c{i},2),cmap2(RF_Primary_num(i,1),:),'facealpha',1,'edgecolor','none');
%         hold on;
%     end
% end

%% construct SM image
SM_image(:,:,1) = zeros(imDims(1),imDims(2));
SM_image(:,:,2) = zeros(imDims(1),imDims(2));
SM_image(:,:,3) = zeros(imDims(1),imDims(2));

SM_image_flat = NaN(imDims(1),imDims(2));
SM_respStrength_image_flat = NaN(imDims(1),imDims(2));
SM_location_image_flat = NaN(imDims(1),imDims(2));

SM_joint_image_flat2 = NaN(imDims(1),imDims(2),size(Joint_num,2));
SM_location_image_flat2 = SM_joint_image_flat2;
SM_respStrength_image_flat = SM_joint_image_flat2;

dispstat('Creating MM image... 0%','init');

if USE_CORNERS
    numSites = length(c)-4;
else
    numSites = length(c);
end

for i = 1:numSites
    tempJ = sum(Joint_num_unique(i,:)~=0);
    tempL = sum(Location_num_unique(i,:)~=0);
    reorderTemp = sortResp(i,:);
    for j = 1:sum(Joint_num(i,:)~=0)-1
        if(respStrength(i,j+1)==0)
            respStrength(i,j+1) = respStrength(i,j);
        end
    end
    
    for x = 1:imDims(1)
        for y = 1:imDims(2)
            if siteMask(x,y,i) == 1
                [in, on] = inpolygon(x,y,v(c{i},2),v(c{i},1));
                
                if ((in || on))
                    
                    if(tempJ>0)
                        SM_image(x,y,1) = cmap_Joint(Joint_num(i,stripeIm(x,y,tempJ)),1);
                        SM_image(x,y,2) = cmap_Joint(Joint_num(i,stripeIm(x,y,tempJ)),2);
                        SM_image(x,y,3) = cmap_Joint(Joint_num(i,stripeIm(x,y,tempJ)),3);
                        
                        SM_image_flat(x,y) = Joint_num(i,stripeIm(x,y,tempJ));
                    end
                    if (tempL>0)
                        SM_location_image_flat(x,y) = Location_num_unique(i,stripeIm(x,y,tempL));
                    end
                    SM_respStrength_image_flat(x,y) = respStrength(i,reorderTemp(1));
                    
                    for t = 1:tempJ
                        SM_joint_image_flat2(x,y,reorderTemp(t)) = Joint_num(i,reorderTemp(t));
                    end
                    for t = 1:tempL
                        SM_location_image_flat2(x,y,reorderTemp(t)) = Location_num_unique(i,reorderTemp(t));
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

if SM_Grid
    for x = xlimits(1):pxPerMM:xlimits(2)
        plot([x x],ylimits,'w-','linewidth',1);
    end
    
    for y = ylimits(1):pxPerMM:ylimits(2)
        plot(xlimits,[y y],'w-','linewidth',1);
    end
end

imagesc(SM_image,'alphadata',~isnan(SM_image_flat));
title([animal,' - Voronoi RF Map']);
xlim(xlimits);
ylim(ylimits);


%% plot response strength
figure('color',[1 1 1],'position',[200 150 1000 800]);
cmap_respStrength = [0 0 0; colormap(flipud(hot(6)))];
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_respStrength); caxis([0 7]);

if SM_Grid
    for x = xlimits(1):pxPerMM:xlimits(2)
        plot([x x],ylimits,'w-','linewidth',1);
    end
    
    for y = ylimits(1):pxPerMM:ylimits(2)
        plot(xlimits,[y y],'w-','linewidth',1);
    end
end

for i = 1:length(siteNum)
    t = linspace(0,2*pi);
    patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
        [1 1 1],'edgecolor','none');
    patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
        cmap_respStrength(max(respStrength(i,:))+1,:),'edgecolor','none');
end
set(h,'ticks',(1:size(cmap_respStrength,1))-0.5,'ticklabels',{'NR','IR','VWR','WR','GR','VGR','ER'});
title([animal,' - Response Strength']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);



%% plot joint
figure('color',[1 1 1],'position',[250 150 1000 800]);
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_Joint); caxis([0 size(cmap_Joint,1)]);

if SM_Grid
    for x = xlimits(1):pxPerMM:xlimits(2)
        plot([x x],ylimits,'w-','linewidth',1);
    end
    
    for y = ylimits(1):pxPerMM:ylimits(2)
        plot(xlimits,[y y],'w-','linewidth',1);
    end
end

for i = 1:length(siteNum)
    temp = sum(Joint_num_unique(i,:)~=0);
    
    for n = 1:temp
        t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
        patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
            [1 1 1],'edgecolor','none');
        patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
            cmap_Joint(Joint_num(i,n),:),'edgecolor','none');
    end
end
set(h,'ticks',(1:size(cmap_Joint,1))-0.5,'ticklabels',{'NR','Digit1','Digit2','Digit3','Digit4','Digit5','Pads','Forearm','Face'});
title([animal,' - Joint']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);

%% plot location
figure('color',[1 1 1],'position',[250 150 1000 800]);
cmap_Location = [1 0 0; 1 .5 0; 0 1 0; 0 0 1;  0 1 1; 0 .2 0];
imagesc(green); axis image off; hold on;
h = colorbar; colormap(cmap_Location); caxis([0 6]);

if SM_Grid
    for x = xlimits(1):pxPerMM:xlimits(2)
        plot([x x],ylimits,'w-','linewidth',1);
    end
    
    for y = ylimits(1):pxPerMM:ylimits(2)
        plot(xlimits,[y y],'w-','linewidth',1);
    end
end

for i = 1:length(siteNum)
    temp = sum(Location_num_unique(i,:)~=0);
    
    for n = 1:temp
        if Location_num_unique(i,n)>0
            t = linspace(2*(n-1)/temp*pi,2*n/temp*pi);
            patch(siteLoc(i,1)+[(siteSize+2)*cos(t) 0],siteLoc(i,2)+[(siteSize+2)*sin(t) 0],...
                [1 1 1],'edgecolor','none');
            patch(siteLoc(i,1)+[siteSize*cos(t) 0],siteLoc(i,2)+[siteSize*sin(t) 0],...
                cmap_Location(Location_num_unique(i,n),:),'edgecolor','none');
        end
    end
end
set(h,'ticks',(1:size(cmap_Location,1))-0.5,'ticklabels',{'Distal','Medial',...
    'Proximal','Glaborous/Ventral','Dorsal','Proprioceptive'});
title([animal,' - Location Type']); set(gca,'fontsize',16);
xlim(xlimits);
ylim(ylimits);

%% SAVE
if SAVE
    if ~exist(['S:\Lab\',animal,'\Mapping\Sensory Maps V3'],'dir')
        mkdir(['S:\Lab\',animal,'\Mapping\Sensory Maps V3']);
    end
    cd(['S:\Lab\',animal,'\Mapping\Sensory Maps V3']);
    
    
    filename = [animal,' SM Joint ',num2str(RAD),'px','.png'];
    imwrite(ind2rgb(SM_image_flat,cmap_Joint),filename,'alpha',double(~isnan(SM_image_flat)));
    
    filename = [animal,' SM respStrength ',num2str(RAD),'px','.png'];
    imwrite(ind2rgb(max(SM_respStrength_image_flat,[],3),cmap_respStrength),filename,'alpha',double(~isnan(max(SM_respStrength_image_flat,[],3))));
    
    filename = [animal,' SM Location ',num2str(RAD),'px','.png'];
    imwrite(ind2rgb(SM_location_image_flat,cmap_Location),filename,'alpha',double(~isnan(SM_location_image_flat)));
    
    % mat files
    filename = [animal,' SM All ',num2str(RAD),'px','.mat'];
    save(filename,'SM_joint_image_flat2','SM_location_image_flat2','SM_respStrength_image_flat');
    
end
function value = assignRS(response)

switch strtrim(response)
    case 'ER'
        value =  6;
    case 'VGR'
        value =  5;
    case 'GR'
        value =  4;
    case 'WR'
        value =  3;
    case 'VWR'
        value =  2;
    case 'IR'
        value =  1;
    case 'NR'
        value = 0;
end
end

function value = assignD(d)
value = [];
if ~isempty(d)
    if contains(d, 'Distal')
        value(end+1) = 1;
    end
    if contains(d,'Medial')
        value(end+1) = 2;
    end
    if contains(d,'Proximal')
        value(end+1) = 3;
    end
    if contains(d,'Glaborous')
        %value(end+1) = 4;
    end
    if contains(d,'Ventral')
        %value(end+1) = 4;
    end
    if contains(d,'Dorsal')
       % value(end+1) = 5;
    end
    if contains(d,'Proprioceptive')
       % value(end+1) = 6;
    end
end
end
% end