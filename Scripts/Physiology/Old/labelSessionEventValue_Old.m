clear all;
dateStart = datetime(2019,6,7, 'Format', 'MM_dd_y');
dateToday = datetime('now', 'Format', 'MM_dd_y');
dateArray = dateStart:dateToday;
dateArray = dateStart:datetime(2019,7,10, 'Format', 'MM_dd_y');

monkey = 'Gilligan';

singleUnits = 0;
binSize=.01; % bin size in seconds
sigma = 25; % smoothing window in ms
secondsBeforeEvent = 5;
secondsAfterEvent = 2.1;
alignSegs = {'GoSignal','StartReach', 'StartGrasp', 'StartLift'};
conds = {'Large Sphere'};
cond2 = {'Extra Small Sphere'};
loc = {'Far'};
loc2= {'Close'};
histBins = -secondsBeforeEvent:binSize:secondsAfterEvent;
pVal = 0.05;
correctedP = pVal/2;
correctedP = pVal;

Kernel = (-5*sigma:3*sigma);
BinSize = length(Kernel);
Kernel = (-BinSize/2:BinSize/2);
Factor = (1.0003)/(sigma*sqrt(2*pi));
Kernel = Factor*(exp(-(0.5*((Kernel./sigma).^2))));

% armHandMask = imread(['S:\Lab\',monkey,'\Mapping\ArmHandMask.png']);
% armMask = armHandMask(:,:,1)>100 & armHandMask(:,:,1)<255;
% handMask = armHandMask(:,:,3)>100 & armHandMask(:,:,3)<255;
% ccArm = bwconncomp(armMask);
% ccHand = bwconncomp(handMask);
% sizeCCArm = cellfun(@length,ccArm.PixelIdxList);
% sizeCCHand = cellfun(@length,ccHand.PixelIdxList);
% [~,maxArm] = max(sizeCCArm);
% [~,maxHand] = max(sizeCCHand);
% armMask = zeros(size(armMask));
% handMask = zeros(size(handMask));
% armMask(ccArm.PixelIdxList{maxArm}) = 1;
% handMask(ccHand.PixelIdxList{maxHand}) = 1;

MM = double(imread(['S:\Lab\',monkey,'\Mapping\Motor Maps V3\MM_Simp_RGB.png']));
%border is yellow
M1border = MM(:,:,1)>200 & MM(:,:,2) > 200 & MM(:,:,3) < 100;
[row,col] = find(M1border);
pointsX = unique(row);
pointsY = arrayfun(@(a) max(col(row==a)),pointsX);
pointsX(end+1:end+2) = [768, 1];
pointsY(end+1:end+2) = [1 1];
M1mask = double(poly2mask(pointsY,pointsX, 768,768));
M1mask(M1mask==0) = NaN;
M1mask = repmat(M1mask,[1,1,3]);
MM = MM.*M1mask;

%arm is red
armMask = MM(:,:,1) > 100 & MM(:,:,2) < 100 & MM(:,:,3) < 100;
armMask = imfilter(conv2(armMask, ones(5)/5^2, 'same')>.5,ones(15));
%hand is green
handMask = MM(:,:,2) > 100 & MM (:,:,1) < 100 & MM(:,:,3) < 100;
handMask = imfilter(conv2(handMask, ones(5)/5^2, 'same')>.5,ones(15));

activationMask = imresize(imread(['S:\Lab\',monkey,'\Mapping\tTests\ESS_HSV.png']),[768,768]);
activationMask = activationMask(:,:,1)>200 & activationMask(:,:,2) < 100 & activationMask(:,:,3) < 100;
ccAct = bwconncomp(activationMask);
sizeCCAct = cellfun(@length,ccAct.PixelIdxList);
smallDomains = sizeCCAct<200;
pixels = cell2mat(ccAct.PixelIdxList(smallDomains)');
activationMask(pixels) = 0;
domains = regionprops(activationMask,'PixelIdxList');
domainIm = zeros(768,768);
for i = 1:length(domains)
    domainIm(domains(i).PixelIdxList) = i;
end

cd(['S:\Lab\',monkey,'\All Data\']);
armHandVal = [];
domainInds = [];
SIVals = [];
ZVals = [];
SIR = [];
SIG = [];
ZR = [];
ZG = [];
reachVals = [];
reachAvgVals = [];
graspVals = [];
graspAvgVals = [];
[num,txt,raw]=xlsread(['S:\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx']);
for d = 1:length(dateArray)
    currName = [monkey,'_', datestr(dateArray(d), 'mm_dd_yyyy')];
    if(exist([currName, '\Physiology\Results\',currName,'_1.mat'], 'file'))
        sessionDir = dir(['S:\Lab\Gilligan\All Data\', monkey,'_',datestr(dateArray(d), 'mm_dd_yyyy'),'\Physiology\Results']);
        
        sessionDir = sessionDir(~[sessionDir.isdir]);
        [~,sortInd] = natsort({sessionDir.name});
        sessionDir = sessionDir(sortInd);
        fid  = fopen(['S:\Lab\',monkey,'\All Data\',monkey,'_'...
            ,datestr(dateArray(d), 'mm_dd_yyyy'),'\Physiology\',monkey,'_',datestr(dateArray(d), 'mm_dd_yyyy'),'_notesVProbe.txt'],'r');
        txtNote = textscan(fid,'%s','Delimiter','','endofline','');
        txtNote = txtNote{1}{1};
        fid  = fclose(fid);
        tk = regexp(txtNote,'Penetration #[\s\.=]+(\d+)','tokens');
        PN = str2double(tk{:});
        FRvalues = [];
        reachGraspSigExV = [];
        reachLags = [];
        reachLagsAvg = [];
        graspLags = [];
        graspLagsAvg = [];
        FRR = [];
        FRG = [];
        reachGraspR = [];
        reachGraspG = [];
        for f = 1:length(sessionDir)
            load([sessionDir(f).folder, '\', sessionDir(f).name]);
            if(isfield(sortedSpikeData, 'label') || any(contains(who('-file',[sessionDir(f).folder, '\', sessionDir(f).name]), 'label')))
                if(~isfield(sortedSpikeData, 'label'))
                    labels = label;
                else
                    labels = sortedSpikeData.label;
                end
                if(singleUnits)
                    numUnits = find(arrayfun(@(a) strcmp(a,'s'), labels));
                else
                    numUnits = 1:length(labels);
                end
                for u = 1:length(numUnits)
                    increase = 0;
                    %currTrialInds = cellfun(@(a) contains(a,conds),sortedSpikeData.ArduinoData(:,1));
                    currTrialInds = find(cellfun(@(a) contains(a,loc),sortedSpikeData.ArduinoData(:,8)));
                    currSpikes = sortedSpikeData.SpikeTimes(numUnits(u),currTrialInds);
                    currSegs = sortedSpikeData.SegTimes(numUnits(u),currTrialInds);
                    spikesPerTrial = cellfun(@length, currSpikes);
                    badTrials = cellfun(@(a) length(a)<10 | length(a)>2000, currSpikes);
                    currSpikes = currSpikes(~badTrials);
                    currSegs = currSegs(~badTrials);
                    
                    %currTrialInds2 = cellfun(@(a) contains(a,cond2),sortedSpikeData.ArduinoData(:,1));
                    currTrialInds2 = find(cellfun(@(a) contains(a,loc2),sortedSpikeData.ArduinoData(:,8)));
                    currSpikes2 = sortedSpikeData.SpikeTimes(numUnits(u),currTrialInds2);
                    currSegs2 = sortedSpikeData.SegTimes(numUnits(u),currTrialInds2);
                    spikesPerTrial2 = cellfun(@length, currSpikes2);
                    badTrials2 = cellfun(@(a) length(a)<10 | length(a)>2000, currSpikes2);
                    currSpikes2 = currSpikes2(~badTrials2);
                    currSegs2 = currSegs2(~badTrials2);
                    if(sum(badTrials)<length(currTrialInds)/2 && sum(badTrials2)<length(currTrialInds2))
                        alignTimesGo = cellfun(@(a) getAlignedTimes(a,1), currSegs, 'UniformOutput', false);
                        alignedSpikesGo = cellfun(@minus, currSpikes, alignTimesGo, 'UniformOutput', false);
                        bins =  -secondsBeforeEvent-length(Kernel)/(2*1/binSize):binSize:...
                            secondsAfterEvent+length(Kernel)/(2*1/binSize);
                        bins = bins(length(Kernel)/2:end-(length(Kernel)/2)+1);
                        trialHists = cellfun(@(a) histcounts(a,bins)./binSize, alignedSpikesGo,'UniformOutput', false);
                        trialHistsSmooth = cellfun(@(a) conv(a,Kernel), trialHists, 'UniformOutput', false);
                        trialHistsSmooth = cellfun(@(a) a((length(Kernel)/2):end-(length(Kernel)/2)+1), trialHistsSmooth,'UniformOutput', false);
                        totalHists = reshape(cell2mat(trialHistsSmooth)',length(trialHistsSmooth{1}),[])';
                        PSTH = mean(totalHists);
                        [~, maxFR] = max(PSTH);
                        maxFRV = histBins(maxFR);
                        [~, pv] = ttest2(reshape(totalHists(:,find(isalmost(bins,-5,.001)):find(isalmost(bins,-3,.001))),[],1), ...
                            reshape(totalHists(:,find(isalmost(bins,0,.001)):find(isalmost(bins,2,.001))),[],1));
                        if(pv<correctedP)
                            alignInd = cell2mat(cellfun(@(a) find(strcmp(sortedSpikeData.ConditionSegments{1}, a)), alignSegs, 'UniformOutput', false));
                            alignTimes = cellfun(@(a) getAlignedTimes(a,alignInd), currSegs, 'UniformOutput', false);
                            alignTimes = cellfun(@(a) [a(1), a(2)-.05, a(3:end)], alignTimes, 'UniformOutput', false);
                            baseSpikes = cellfun(@(a,b) length(a(a>b(1)-5 & a<b(1)))/5,currSpikes, alignTimes, 'UniformOutput', false);
                            reachSpikes = cellfun(@(a,b,c) abs((length(a(a>b(2) & a<b(3)))/(b(3)-b(2)))-c), currSpikes, alignTimes, baseSpikes,'UniformOutput', false);
                            graspSpikes = cellfun(@(a,b,c) abs((length(a(a>b(3) & a<b(4)))/(b(4)-b(3)))-c), currSpikes, alignTimes, baseSpikes,'UniformOutput', false);
                            FRvalues(end+1) = (mean(cell2mat(reachSpikes))-mean(cell2mat(graspSpikes))) / mean(cellfun(@(a,b) a+b, reachSpikes, graspSpikes));
                            reachGraspSigExV(end+1) = (mean(cell2mat(reachSpikes)) - mean(cell2mat(graspSpikes)))/sqrt(var(cell2mat(reachSpikes)) + var(cell2mat(graspSpikes)));
                            [~,maxPeakI,~,pi] = findpeaks(PSTH);
                            [~,maxPeakD,~,pd] = findpeaks(-PSTH);
                            peaks = [pi pd];
                            maxPeaks = [maxPeakI maxPeakD];
                            [~,largestP] = sort(peaks,'descend');
                            maxPeaks = maxPeaks(largestP(1));
                            maxPeak = min(maxPeaks);
                            [~,maxPeak] = max(PSTH);
                            [~,maxPeakTrial] = max(totalHists');
                            reachLagsAvg(end+1) = bins(maxPeak) - mean(cellfun(@(a,b) a(2)+.05-b, alignTimes,alignTimesGo));
                            reachLags(end+1) = mean(bins(maxPeakTrial)-cellfun(@(a,b) a(2)+.05-b, alignTimes,alignTimesGo));
                            graspLagsAvg(end+1) = bins(maxPeak) - mean(cellfun(@(a,b) a(3)-b, alignTimes,alignTimesGo));
                            graspLags(end+1) = mean(bins(maxPeakTrial)-cellfun(@(a,b) a(3)+-b, alignTimes,alignTimesGo));
                            if(~isempty(cond2) || ~isempty(loc2))
                                alignTimes2 = cellfun(@(a) getAlignedTimes(a,alignInd), currSegs2, 'UniformOutput', false);
                                alignTimes2 = cellfun(@(a) [a(1), a(2)-.05, a(3:end)], alignTimes2, 'UniformOutput', false);
                                baseSpikes2 = cellfun(@(a,b) length(a(a>b(1)-5 & a<b(1)))/5,currSpikes2, alignTimes2, 'UniformOutput', false);
                                reachSpikes2 = cellfun(@(a,b,c) abs((length(a(a>b(2) & a<b(3)))/(b(3)-b(2)))-c), currSpikes2, alignTimes2, baseSpikes2,'UniformOutput', false);
                                graspSpikes2 = cellfun(@(a,b,c) abs((length(a(a>b(3) & a<b(4)))/(b(4)-b(3)))-c), currSpikes2, alignTimes2, baseSpikes2,'UniformOutput', false);
                                FRR(end+1) = (mean(cell2mat(reachSpikes))-mean(cell2mat(reachSpikes2))) / (mean(cell2mat(reachSpikes))+mean(cell2mat(reachSpikes2)));
                                FRG(end+1) = (mean(cell2mat(graspSpikes))-mean(cell2mat(graspSpikes2))) / (mean(cell2mat(graspSpikes))+mean(cell2mat(graspSpikes2)));
                                reachGraspR(end+1) = (mean(cell2mat(reachSpikes)) - mean(cell2mat(reachSpikes2)))/sqrt(var(cell2mat(reachSpikes)) + var(cell2mat(reachSpikes2)));
                                reachGraspG(end+1) = (mean(cell2mat(graspSpikes)) - mean(cell2mat(graspSpikes2)))/sqrt(var(cell2mat(graspSpikes)) + var(cell2mat(graspSpikes2)));
                            end
                        end
                    end
                    
                end
                
            end
        end
        %FRvalues = (FRvalues+1)./2;
        [~, siteCol] = find(contains(txt, 'Site #'));
        [row, ~ ] = find(num(:,siteCol)==PN);
        xLoc = num(row,2);
        yLoc = num(row,3);
        %             status = 0;
        %             while(status==0)
        %             [status, message] = xlswrite(['\\130.49.229.252\gharbawie\Lab\', monkey,'\Mapping\', monkey,'_MM_Sites.xlsx'],SI,'Sheet1',[col, num2str(row+1)]);
        %             pause(.5);
        %             end
        ZVals = [ZVals; reachGraspSigExV'];
        SIVals = [SIVals; FRvalues'];
        ZG = [ZG; reachGraspG'];
        ZR = [ZR; reachGraspR'];
        SIR = [SIR; FRR'];
        SIG = [SIG; FRG'];
        reachAvgVals = [reachAvgVals; reachLagsAvg'];
        reachVals = [reachVals; reachLags'];
        graspAvgVals = [graspAvgVals; graspLagsAvg'];
        graspVals = [graspVals; graspLags'];
        domainInds(end+1:end+length(FRvalues)) = domainIm(yLoc,xLoc);
        if(armMask(yLoc,xLoc)==1)
            armHandVal(end+1:end+length(FRvalues)) = 1;
        elseif(handMask(yLoc,xLoc)==1)
            armHandVal(end+1:end+length(FRvalues)) = 2;
        else
            armHandVal(end+1:end+length(FRvalues)) = 0;
        end
        disp(datestr(dateArray(d), 'mm_dd_yyyy'));
    end
end
%%
SIFig = figure;
ZFig = figure;
reachLagFig = figure;
graspLagFig = figure;

plotHists(gca(SIFig),SIVals(armHandVal==1),SIVals(armHandVal==2),-.9,.9,-.4,.4,'SIVals',1);
plotHists(gca(ZFig),ZVals(armHandVal==1),ZVals(armHandVal==2),-2,2,-1,1,'ZVals',1);
plotHists(gca(reachLagFig),reachAvgVals(armHandVal==1),reachAvgVals(armHandVal==2),-.2,.2,NaN,NaN,'Reach Lags', 1);
plotHists(gca(graspLagFig),graspAvgVals(armHandVal==1),graspAvgVals(armHandVal==2),-.2,.2,NaN,NaN,'Grasp Lags', 1);

domainFigSI = figure();
domainFigZ = figure();
for i =0:length(domains)
    set(0,'CurrentFigure',domainFigSI)
    SIs = subplot(ceil(length(domains)/3),3,i+1);
    set(0,'CurrentFigure',domainFigZ)
    Zs = subplot(ceil(length(domains)/3),3,i+1);
    if(i==0)
       imshow(domainIm, 'Parent', SIs);
       colormap(SIs,'jet');
       caxis(SIs,[-1, length(domains)]);
       title(SIs,'SI values');
       
       imshow(domainIm, 'Parent', Zs);
       colormap(Zs,'jet');
       caxis(Zs,[-1, length(domains)]);
       title(Zs,'Z values');
    else
        valsHand = domainInds==i & armHandVal==2;
        valsArm = domainInds ==i & armHandVal==1;
        plotHists(Zs,ZVals(valsArm),ZVals(valsHand),-2,2,-1,1,num2str(i));
        plotHists(SIs,SIVals(valsArm),SIVals(valsHand),-.9,.9,-.4,.4,num2str(i));

    end
end
function aligned = getAlignedTimes(a,alignInd)
if(sum(isnan(a))==0)
    aligned = a(alignInd);
else
    aligned = NaN;
end
end
function plotHists(ax, valsArm, valsHand,minVal, maxVal,minCutoff, maxCutoff,titleInfo,legendFlag)
stepVal = abs((minVal-maxVal)/9);
[bArm,~] = hist(valsArm,minVal:stepVal:maxVal);
[bHand,e] = hist(valsHand,minVal:stepVal:maxVal);
b = bar(ax,e, [bArm./sum(bArm);bHand./sum(bHand)]');
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 0 1];
YZ = ax.YLim(2);
line(ax,[minCutoff, minCutoff], [0, YZ],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
line(ax,[maxCutoff, maxCutoff], [0, YZ],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
ax.YLim(2) = YZ;
if(exist('legendFlag','var'))
    legend(ax,sprintf('Arm (n=%d)',length(valsArm)),sprintf('Hand (n=%d)',length(valsHand)),'Location', 'best');
    title(ax,titleInfo);
else
    title(ax,[titleInfo, sprintf(', Arm (n=%d), Hand (n=%d)',length(valsArm),length(valsHand))]);
end
end