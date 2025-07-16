close all;
clearvars;

%% user input
animal = 'Gilligan';
RAD = 42; %28 = 0.5mm, 1.0mm^2; 42 = 0.75 mm, 2.2mm^2; 56 = 1mm, 3.4mm^2;
USE_CORNERS = true;
xlimits = [10 1080];
xlimits = [10 758];
ylimits = [10 1300];
ylimits = [10 758];
PSTH_type =  'Trial';
cond = 'ESS';
columnNumber = 18;
singleOrAll = 'Single';

%% load and parse excel sheet
[excel.num, excel.txt, excel.raw] = xlsread(['S:\Lab\',...
    animal, '\Mapping\Encoding Maps\PSTHs\',singleOrAll,'\FRs\','M1_',cond,'_',PSTH_type, '_Summary.xlsx']);
headings = excel.raw(1,:);
saveName = [headings{columnNumber},PSTH_type,'_', singleOrAll,'_',cond];


NRi = logical(~isnan(excel.num(:,find(strcmp(headings, 'Total Units')))));
excel.raw = excel.raw(NRi,:);
excel.num = excel.num(NRi,:);
excel.txt = excel.txt(NRi,:);
movementsInd = find(strcmp(headings, 'Evoked Movement(s)'));
commas = find(cellfun(@(a) strcmp(a(end-1:end), ', '), excel.raw(:,movementsInd)));

for m = 1:length(commas)
    excel.raw{commas(m),movementsInd} = excel.raw{commas(m),movementsInd}(1:end-2);
end

for i = 1:length(excel.raw)
    siteNum(i,1) = excel.raw{i,find(strcmp(headings, 'Site #'))};
    siteLoc(i,:) = [excel.raw{i,find(strcmp(headings, 'x'))},...
        excel.raw{i,find(strcmp(headings, 'y'))}];
end
%% load green
green = imread(['S:\Lab\', animal, '\Mapping\M1_green_edited.bmp']);

corners(1,:) = [0 size(green,2)];
corners(2,:) = [size(green,1) 0];
corners(3,:) = [0 0];
corners(4,:) = [size(green,1) size(green,2)];
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
dispstat('Creating MM image... 0%','init');
if USE_CORNERS
    numSites = length(c)-4;
else
    numSites = length(c);
end
plotMap = NaN(imDims(1), imDims(2));
for i = 1:numSites
    temp = 1;
    for x = 1:imDims(1)
        for y = 1:imDims(2)
            if siteMask(x,y,i) == 1
                [in, on] = inpolygon(x,y,v(c{i},2),v(c{i},1));
                if (in || on)
                    plotMap(x,y) = round(10*excel.num(i,columnNumber));
                end
            end
        end
    end
    percDone = (i)/(length(c)-4);
    percDone = round(percDone*10000)/100;
    dispstat(['Creating MM image... ',num2str(percDone),'%']);
end
%%
figure('Units', 'normalized', 'Position', [0 0 1 1]);
imagesc(plotMap);
axis image;

colormap(jet);
caxis([1,10]);
cc = colorbar;
cc.Ticks = 1:10;
cc.TickLabels = cc.Ticks/10;
saveFigures(gcf,['S:\Lab\',animal, '\Mapping\Encoding Maps\Voronoi\'],saveName,[])