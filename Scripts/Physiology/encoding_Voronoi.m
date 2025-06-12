close all;
clearvars;
%% user input
animal = 'Skipper';
RAD = 42; %28 = 0.5mm, 1.0mm^2; 42 = 0.75 mm, 2.2mm^2; 56 = 1mm, 3.4mm^2;
PSTH_type =  'Unit';
cond = 'Extra Small Sphere';
countType = {'Single'};
USE_CORNERS = true;
SAVE = true;
countZero = true;

colorSteps = 6;
totalRange = [linspace(0,60,colorSteps),Inf];
unitRange = [linspace(0,50,colorSteps),Inf];
encodingRange = [linspace(0,3.3,colorSteps),Inf];
subRange = [-Inf,linspace(-1.3,1.3,colorSteps),Inf];
subIndRange = [-Inf,linspace(-2,2,colorSteps-1),Inf];

unitType = {'Go', 'Task', 'Reach', 'Grasp'};
phaseType = {'Go', 'Task', 'Reach', 'Grasp', 'Rest'};
MMfilename = 'MM_Simp_RGB-01.png';
parentPath = ['S:\Lab\', animal,'\Mapping\'];

subTypes = cellfun(@(a) {strcat(a, ' Units Grasp Phase FR Change'),...
    strcat(a,' Units Reach Phase FR Change')},[countType, unitType], 'UniformOutput', false);
encodingType = cellfun(@(b) cellfun(@(a) [b, ' Units ', a, ' Phase FR Change'],...
    phaseType,'UniformOutput', false),[countType, unitType],'UniformOutput', false);
xcorrType = cellfun(@(a) [a, 'Units XCorr'], unitType, 'UniformOutput', false);
encodingType = [encodingType{:}];
%% load green
MM = double(imread([parentPath, 'Motor Maps V3\', MMfilename]));
% M1/PMd border is yellow
M1border = bwareafilt(MM(:,:,1)>125 & MM(:,:,2)>125 & MM(:,:,3)<125, 1);
% M1 mask
[row,col] = find(M1border);
pointsX = unique(row);
pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
pointsX(end+1:end+2) = [size(MM,2), 1];
pointsY(end+1:end+2) = [1 1];
M1mask = poly2mask(pointsY,pointsX, size(MM,1),size(MM,2));
imDims = size(MM);
clear MM green;
%% load and parse excel sheet
[excel.num, excel.txt, excel.raw] = xlsread(['S:\Lab\',...
    animal, '\Mapping\Encoding Maps\PSTHs\',countType{1},'\FRs\','M1_',PSTH_type, '_', cond,'_Summary.xlsx']);
headings = excel.raw(1,:);

NRi = logical(~isnan(excel.num(:,find(strcmp(headings, 'Total Units')))));

excel.raw = excel.raw(NRi,:);
excel.num = excel.num(NRi,:);
excel.txt = excel.txt(NRi,:);

commas = find(cellfun(@(a) strcmp(a(end-1:end), ', '), excel.raw(:,4)));
for m = 1:length(commas)
    excel.raw{commas(m),4} = excel.raw{commas(m),4}(1:end-2);
end
for i = 1:length(excel.num)
    siteNum(i,1) = excel.num(i,1);
    siteLoc(i,:) = [excel.num(i,2) excel.num(i,3)];
end
excelNum = excel.num;
excelRaw = excel.raw;
clear excel;
%% create radial mask around every site
circle = fspecial('disk',RAD);
siteMask = zeros(imDims(1),imDims(2),length(siteNum));
for n = 1:length(siteNum)
    tempCircle = zeros(imDims(1)+400,imDims(2)+400);
    tempCircle((siteLoc(n,2)-RAD+200):(siteLoc(n,2)+RAD+200),(siteLoc(n,1)-RAD+200):(siteLoc(n,1)+RAD+200)) = circle;
    tempCircle = tempCircle(200:end-201,200:end-201);
    siteMask(:,:,n) = tempCircle>0;
    M1Inds(n) = sum(sum(siteMask(:,:,n) & M1mask))>(numel(circle)/5);
end
siteNum = siteNum(M1Inds);
siteLoc = siteLoc(M1Inds,:);
excelNum = excelNum(M1Inds,:);
siteMask = siteMask(:,:,M1Inds);
%% encoding assignmets
%% Unit counts
countColoring = [linspace(.35,.35,colorSteps)', ...
    linspace(.5,1,colorSteps)', linspace(1,.6,colorSteps)'];
if(countZero)
    countColoring(end+1,:) = [.35, .3, 1];
end
countColoring = hsv2rgb(countColoring);
tMaps = cell(1,length(countType));
tMaps(1:end) = {NaN(size(siteNum))};
for t = 1:length(countType)
    for p = 1:length(excelNum)
        if(strcmp(countType{t}, 'All'))
            currValueT = excelNum(p,find(strcmp(headings, 'Total Units')));
        else
            currValueT = excelNum(p,find(strcmp(headings, [countType{t},' Units'])));
        end
        for r = 1:length(totalRange)-1
            if(countZero && currValueT == 0)
                tMaps{t}(p) = length(countColoring);
            elseif(currValueT>totalRange(r) && currValueT<=totalRange(r+1))
                tMaps{t}(p) = r;
                if(t==1 & r==length(totalRange)-1)
                    disp(currValueT);
                end
            end
        end
    end
end

countColoring2 = [linspace(.6,.6,colorSteps-countZero)' ...
    linspace(.5,1,colorSteps-countZero)', linspace(1,.6,colorSteps-countZero)'];
if(countZero)
    countColoring2(end+1,:) = [.6, .3, 1];
end
countColoring2 = hsv2rgb(countColoring2);
uMaps = cell(1,length(unitType));
uMaps(1:end) = {NaN(size(siteNum))};
for u = 1:length(unitType)
    for p = 1:length(excelNum)
        currValueU = excelNum(p,find(strcmp(headings, [unitType{u}, ' Unit Counts'])));
        for r = 1:length(unitRange)-1
            if(countZero && currValueU==0)
                uMaps{u}(p) = length(countColoring2);
            elseif(currValueU>unitRange(r) && currValueU<=unitRange(r+1))
                uMaps{u}(p) = r;
            end
        end
    end
end
%% Encoding maps
encodingColoring = [.6,1,1; linspace(0,0,colorSteps-(1+countZero))',...
    linspace(.5,1,colorSteps-(1+countZero))', ...
    linspace(1,.6,colorSteps-(1+countZero))'];
if(countZero)
    encodingColoring(end+1,:) = [0, .3, 1];
end
encodingColoring = hsv2rgb(encodingColoring);
eMaps = cell(1,length(encodingType));
eMaps(1:end) = {NaN(size(siteNum))};
[~, noChangeInd] = min(abs(diff(abs(encodingRange-1))));
for e =1:length(eMaps)
    for p = 1:length(excelNum)
        currValueE = excelNum(p,find(strcmp(headings, encodingType{e})));
        for r = 1:length(encodingRange)-1
            if(countZero && currValueE>encodingRange(noChangeInd) && ...
                    currValueE <=encodingRange(noChangeInd+1))
                eMaps{e}(p) = length(encodingRange);
            elseif(currValueE>encodingRange(r) && currValueE<=encodingRange(r+1))
                eMaps{e}(p) = r;
            end
        end
    end
end
%% Subtracted encodings
encodingColoringSub = [linspace(.15,.15,colorSteps/2)', linspace(.5,1,colorSteps/2)',...
    linspace(1,.5,colorSteps/2)';linspace(.75,.75,colorSteps/2)',...
    linspace(1,.5,colorSteps/2)', linspace(.6,1,colorSteps/2)';];
encodingColoringSub(end+1,:) = [0 0 .5];
encodingColoringSub = hsv2rgb(encodingColoringSub);
sMaps = cell(1,length(subTypes));
sMaps(1:end) = {NaN(size(siteNum))};
for s = 1:length(sMaps)
    for p = 1:length(excelNum)
        currValueS = excelNum(p,find(strcmp(headings, subTypes{s}{1}))) -...
            excelNum(p,find(strcmp(headings, subTypes{s}{2})));
        for r = 1:length(subRange)-1
            if(countZero && currValueS==0)
                sMaps{s}(p) = length(encodingColoringSub);
            elseif(currValueS>subRange(r) && currValueS<=subRange(r+1))
                sMaps{s}(p) = r;
            end
        end
    end
end
%% XCorr maps
xcorrColoring = colormap(jet(colorSteps));
xMaps = cell(1,length(xcorrType));
xMaps(1:end) = {NaN(size(siteNum))};
% for x =1:length(xcorrType)
%     for p = 1:length(excelNum)
%         xMaps{x}(p) = round(10*excelNum(p,find(strcmp(headings, xcorrType{x}))));
%     end
% end
%% striped image
stripeWidth = 6;
stripeIm = ones(imDims(1),imDims(2),5);
x = 1:stripeWidth*2:imDims(1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,2) = 2;
end
x = 1:stripeWidth*3:imDims(1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,3) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,3) = 3;
end
x = 1:stripeWidth*4:imDims(1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,4) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,4) = 3;
    stripeIm((x(i)+stripeWidth*2+1):(x(i)+stripeWidth*3),:,4) = 4;
end
x = 1:stripeWidth*5:imDims(1);
for i = 1:length(x)
    stripeIm(x(i):(x(i)+stripeWidth),:,4) = 2;
    stripeIm((x(i)+stripeWidth+1):(x(i)+stripeWidth*2),:,5) = 3;
    stripeIm((x(i)+stripeWidth*2+1):(x(i)+stripeWidth*3),:,5) = 4;
    stripeIm((x(i)+stripeWidth*3+1):(x(i)+stripeWidth*4),:,5) = 5;
end
%% corner sites
if USE_CORNERS
    corners(1,:) = [0 imDims(2)];
    corners(2,:) = [imDims(1) 0];
    corners(3,:) = [0 0];
    corners(4,:) = [imDims(1) imDims(2)];
end
%% voronoi
if ~USE_CORNERS
    [v, c] = voronoin(siteLoc);
else
    [v, c] = voronoin([siteLoc; corners]);
end
%% construct images
for t = 1:size(tMaps,2)
    tMapsFlat{t} = NaN(imDims(1), imDims(2));
end
for u = 1:size(uMaps,2)
    uMapsFlat{u} = NaN(imDims(1), imDims(2));
end
for s = 1:size(eMaps,2)
    eMapsFlat{s} = NaN(imDims(1), imDims(2));
end
for e = 1:size(sMaps,2)
    sMapsFlat{e} = NaN(imDims(1), imDims(2));
end
% for x = 1:size(xMaps,2)
%     xMapsFlat{x} = NaN(imDims(1), imDims(2));
% end
dispstat('Creating MM image... 0%','init');
if USE_CORNERS
    numSites = length(c)-4;
else
    numSites = length(c);
end
for i = 1:numSites
    temp = 1;
    for x = 1:imDims(1)
        for y = 1:imDims(2)
            if siteMask(x,y,i) == 1
                [in, on] = inpolygon(x,y,v(c{i},2),v(c{i},1));
                if (in || on)
                    for t = 1:size(tMaps,2)
                        tMapsFlat{t}(x,y) = tMaps{t}(i,stripeIm(x,y,temp));
                    end
                    for u = 1:size(uMaps,2)
                        uMapsFlat{u}(x,y) = uMaps{u}(i,stripeIm(x,y,temp));
                    end
                    for s = 1:size(sMaps,2)
                        sMapsFlat{s}(x,y) = sMaps{s}(i,stripeIm(x,y,temp));
                    end
                    for e = 1:size(eMaps,2)
                        eMapsFlat{e}(x,y) = eMaps{e}(i,stripeIm(x,y,temp));
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
%% SAVE
if SAVE
    if ~exist(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\'],'dir')
        mkdir(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\']);
    end
    cd(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\']);
    NRstr = ['_', cond];
    if(length(countType)==1)
        if(strcmp(countType, 'All'))
            animal = ['MULTI_', animal];
        end
    end
    for t = 1:length(countType)
        filename = [animal, ' ', countType{t}, ' Units ', num2str(RAD),'px',NRstr,'.eps'];
        set(gcf,'visible','off','Units', 'Pixels', 'Position', [0 0 flip(imDims(1:2))]);
        %set(gcf,'visible','off') %suppress figure
        tMapsFlat{t} = tMapsFlat{t}+1;
        countColoring = [1 1 1;countColoring];
        tMapsFlat{t}(isnan(tMapsFlat{t})) = 1;
        imagesc(ind2rgb(tMapsFlat{t},countColoring));
        axis image               % resolution based on image
        axis off                 % avoid printing axis
                set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
        set(gca,'LooseInset',get(gca,'TightInset')); % removing extra white space in figure
        saveas(gcf,filename,'epsc');
    end
    for u = 1:length(unitType)
        filename = [animal, ' ', unitType{u}, ' Units ', num2str(RAD),'px',NRstr,'.png'];
        imwrite(ind2rgb(uMapsFlat{u},countColoring2),filename,'alpha',double(~isnan(uMapsFlat{u})));
        imagesc(ind2rgb(uMapsFlat{u},countColoring2));
        saveImages(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\'],filename(1:end-4),[]);
    end
    for e = 1:length(encodingType)
        filename = [animal, ' ', encodingType{e}, ' ', num2str(RAD),'px',NRstr,'.png'];
        imwrite(ind2rgb(eMapsFlat{e},encodingColoring),filename,'alpha',double(~isnan(eMapsFlat{e})));
        imagesc(ind2rgb(eMapsFlat{e},encodingColoring));
                saveImages(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\'],filename(1:end-4),[]);

    end
    for s = 1:length(subTypes)
        filename = [animal, ' ', subTypes{s}{1}, '-',subTypes{s}{2} , ' '];
        imwrite(ind2rgb(sMapsFlat{s},encodingColoringSub),...
            [filename,num2str(RAD),'px',NRstr,'.png'],'alpha',double(~isnan(sMapsFlat{s})));
        imagesc(ind2rgb(sMapsFlat{s},encodingColoringSub));
                saveImages(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\'],filename(1:end-4),[]);

        subtractedG_R = eMapsFlat{find(strcmp(encodingType,subTypes{s}{1}))}...
            - eMapsFlat{find(strcmp(encodingType,subTypes{s}{2}))};
        indSubtractedG_R = NaN(size(subtractedG_R));
        for r = 1:length(subIndRange)-1
            indSubtractedG_R(subtractedG_R>subIndRange(r) & ...
                subtractedG_R<=subIndRange(r+1)) = r;
        end
        imwrite(ind2rgb(indSubtractedG_R, encodingColoringSub),...
            [filename,'_IND ',num2str(RAD),'px',NRstr,'.png'],...
            'alpha', double(~isnan(indSubtractedG_R)));
        imagesc(ind2rgb(indSubtractedG_R, encodingColoringSub));
                saveImages(['S:\Lab\',animal,'\Mapping\Encoding Maps\Voronoi\',PSTH_type,'_PSTHs\'],filename(1:end-4),[]);

    end
    %     for x = 1:length(xcorrType)
    %         filename = [animal, ' ', xcorrType{x}, ' ', num2str(RAD),'px',NRstr,'.png'];
    %         imwrite(ind2rgb(xMapsFlat{x},xcorrColoring),filename,'alpha',double(~isnan(xMapsFlat{x})));
    %     end
    
    %     MM_numbered_joints = ind2rgb(allUnitsIM,countColoring);
    %     for i = 1:numSites
    %        MM_numbered_joints =  insertText(MM_numbered_joints,[siteLoc(i,:)], num2str(siteNum(i)), 'AnchorPoint', 'center');
    %     end
    %     imwrite(MM_numbered_joints, countColoring, [animal, ' Encoding Numbered ', num2str(RAD), 'px', NRstr,'.png'], 'alpha', double(~isnan(allUnitsIM)));
    
    %     filename = [animal,' Encoding ',num2str(RAD),'px',NRstr,'.mat'];
    %     save(filename,'allUnitsIM','singleUnitsIM','reachUnitsIM','graspUnitsIM',...
    %         'taskUnitsIM', 'goUnitsIM','reachFrIM','graspFrIM','taskFrIM','goFrIM',...
    %         'G_R_FrIM','wTFRIM','reachTFrIM','graspTFrIM','goTFrIM','G_R_TFrIM',...
    %         'indSubtractedG_R', 'indSubtractedTG_R', 'countColoring', ...
    %         'countColoring2', 'encodingColoring', 'encodingColoringSub',...
    %          'siteLoc', 'siteNum', 'v', 'c','siteMask');
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