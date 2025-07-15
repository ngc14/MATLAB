%WithinSessionGraphs.m
%CREATED: 5/16/2017
%LAST EDITTED: 04/26/2018

% Create graphs for an individual session that shows the percent error,
% reaction time, reach duration, lift duration, hold duration, and
% withdrawal duration for each condition. The percent error, reaction time,
% etc. are appeneded to the EXCEL spreadsheet GilliganCummulativeData.xlsx.

% NOTE ABOUT CONDITIONS: This script does not allow the user to modify the
% conditions (e.g small sphere, photocell, etc.) as this should remain
% consistent in the GilliganCummulativeData.xlsx for subsequent analysis
% with AcrossSessionsGraphs.m.

% NOTE ABOUT SEGMENTS OF THE TASK: This script only works for sessions that
% have a reaction time, reach duration, lift duration, hold duration, and
% withdrawal duration.

clear
close all
%% VARIABLES FOR USER MODIFICATION
appendToGilliganCummulativeData='no'; %'yes' or 'no'*
monkey = 'Skipper';
sessionDate = '11_25_2020';

%splitErrors = {'Stim', ''};
splitErrors = {'Extra Small Sphere','Large Sphere','Photocell','Rest'};

%mainEffect = {'Extra Small Sphere','Stim_Extra Small Sphere','Photocell', 'Stim_Photocell'};
mainEffect = {'Extra Small Sphere','Large Sphere','Photocell'};

%sigComps = {{'Extra Small Sphere','Stim_Extra Small Sphere'},{'Photocell', 'Stim_Photocell'}}
sigComps = {{[],[]}};

% Variable needed for multiple plots and proper addition to
% GilliganCummulative Data (DO NOT CHANGE)
allConditions_withRest={'Extra Small Sphere','Large Sphere','Photocell' ,'Rest'};
%allConditions_withRest={'Extra Small Sphere','Stim_Extra Small Sphere','StimBeforeGo_Extra Small Sphere','StimBeforeGo10_Extra Small Sphere'};
% Do not append the same file information twice! If a mistake is made, the
% duplicate information in appendToGilliganCummulativeData will need to be
% deleted.
alphaVal = 0.05;
%% IMPORT DATA FROM TXT FILE
% need importdatafile function
% importdatafile works regardless of if the txt file has a header row or not

%%%%%%%%%%%%% can edit this portion so that the script runs for multiple files
% cd('S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16')
% tempallfilenames=dir;
% tempallfilenames={tempallfilenames.name};
% allfilenames={tempallfilenames{:,[89:133]}};
%
% for z=1:length(allfilenames)
%     filename=allfilenames{z};
%%%%%%%%%%%%%
filename= [monkey,'_',sessionDate]; % DO NOT INCLUDE FILE EXTENSION (e.g. '.txt')
filePathDir=['S:\Lab\', monkey,'\All Data\',filename,'\Arduino\'];

% Copy Arduino file from BehaveBackup\Monkey_Training to
% % Gilligan\All Data
if(~exist(filePathDir, 'dir'))
    mkdir(filePathDir);
end
if(~exist([filePathDir,filename,'txt'],'file'))
    if(strcmp(monkey, 'Gilligan'))
        copyfile(['S:\Lab\BehaveBackup\Monkey_Training\Gilligan_Macaque_110-16\',filename,'\',filename,'.txt'],[filePathDir,filename,'.txt']);
    elseif(strcmp(monkey, 'Skipper'))
        copyfile(['S:\Lab\BehaveBackup\Monkey_Training\Skipper_Macaque_152_17\',filename,'\',filename,'.txt'],[filePathDir,filename,'.txt']);
    end
end
% Import data
addpath('S:\Lab\Jennifer\MATLAB\GilliganGraphScripts');

session=importdatafile([filePathDir,filename,'.txt']);

% Removes pauses
session=session(~isnan(session.Trial),:);
% blockInd = randi(50,1,30);
% session = session(ismember(session.Block,blockInd),:);

% Remove trials with 1 ms durations
Idx_1ms=sum((table2array(session(:,5:9))==1),2)~=0;
trialIdx=find(~Idx_1ms);
session=session(trialIdx,:);

% Data for GilliganCummulativeData
spreadsheetData{1}=filename;

numAllConditions_withRest=length(allConditions_withRest);

% Seperate Rest conditions from all other conditions
conditions=allConditions_withRest(~strcmp(allConditions_withRest,'Rest'));
numConditions=length(conditions);
restTrialIdx=strcmp(session.Condition,'Rest');
session_noRest=session(~restTrialIdx,:);
session_Rest=session(restTrialIdx,:);

%% PLOT OF ERRORS
% NOTE: The number of each specific type of error may not match the last
% line in the text file outputted from Arduino because the rest conditions
% are dealt seperately in this script. For instance, a false start on a
% rest trial is not counted as a false start here. Rather, a false start on
% a rest trial is ONLY counted as a rest error.
%h1=figure('units','normalized','outerposition',[0 0 1 1],'Name',[filename, ': plot of errors'],'NumberTitle','off');

% Rest errors
% Convert Successful Trial column from text to matrix of numbers and NaNs
tempSuccessfulTrial=cellfun(@str2num,session.SuccessfulTrial,'UniformOutput',false);
emptyIndx=cellfun(@isempty,tempSuccessfulTrial);
tempSuccessfulTrial(emptyIndx)={NaN};
successfulTrial=cell2mat(tempSuccessfulTrial);
restTrials=successfulTrial(restTrialIdx);

numRestTrials=length(restTrials);
failedRest=sum(isnan(restTrials));

barGroup = [];
% All other errors
for s = 1:length(splitErrors)
    if(isempty(splitErrors{s}))
        sessionSplit = session(~contains(session.Condition, splitErrors(~ismember(1:length(splitErrors),s))),:);
    else
        sessionSplit = session(contains(session.Condition,splitErrors{s}),:);
    end
    falseStart=sum(strcmpi(sessionSplit.SuccessfulTrial,'False Start'));
    failedReach=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Reach'));
    failedContact=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Lift') ...
        & isnan(sessionSplit.LiftDuration)) + sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Contact'));
    failedLift=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Lift')& ~isnan(sessionSplit.LiftDuration));
    failedHold=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Hold'));
    failedReplace=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Replace'));
    failedReplaceHold=sum(strcmpi(sessionSplit.SuccessfulTrial,'Failed-to-Replace-hold'));
    barGroup = [barGroup; falseStart,failedReach,failedContact,...
        failedLift,failedHold,failedReplace,failedReplaceHold];
end
totalErrors=sum(sum(barGroup));
trials=size(session,1);
errorPercent=totalErrors/trials*100;
if(any(cellfun(@isempty, splitErrors)) && length(splitErrors)>1)
    splitErrors{find(cellfun(@isempty,splitErrors))} = 'No_Stim';
end

% subplot(1,2,1)
% bar(barGroup')
if isempty(session_Rest)
    errorType={'False Start','Failed Reach','Failed Contact','Failed Lift',...
        'Failed Hold','Failed Replace','Failed Replace Hold'};
else
    errorType={'False Start','Failed Reach','Failed Contact','Failed Lift',...
        'Failed Hold','Failed Replace','Failed Replace Hold','Failed Rest'};
end
% xlabel('Error Type')
% ylabel('Number of Occurances')
% set(gca,'XTickLabel',errorType,'XTick',1:numel(errorType),'FontSize',16,...
%     'box','off','xTickLabelRotation',45)
% legend(splitErrors);
% ylim([0,max([falseStart,failedReach,failedContact,failedLift,...
%     failedHold,failedReplace,failedReplaceHold])+5])
% title('Number of Each Error Type')
% 
% text(7,14,['Total Errors: ',num2str(totalErrors)],'FontSize',16)
% text(7,11,['Number of Trials: ',num2str(trials)],'FontSize',16)
% text(7,8,['Error Percent: ',num2str(errorPercent),'%'],'FontSize',16)

% Find number of trials for each condition
correctTrialsPerCond=[];
for i=1:numAllConditions_withRest
    tempCondIdx=strcmp(session.Condition,allConditions_withRest(i));
    ct = session.Block(tempCondIdx);
    correctTrialsPerCond=[correctTrialsPerCond,mean(arrayfun(@(s) sum(ct==s), unique(ct))-1)];
end
correctTrialsPerCond'
% Data for GilliganCummulativeData
tempCorrectTrials=cellfun(@str2num,session.SuccessfulTrial,'UniformOutput',false);
correctTrials=sum(~cellfun(@isempty,tempCorrectTrials));
spreadsheetData{2}=[correctTrials,trials,correctTrialsPerCond,falseStart,...
    failedReach,failedContact,failedLift,failedHold,failedReplace,...
    failedReplaceHold,[],totalErrors,errorPercent];
clearvars restData restTrials restTrialsIdx failedRest errorType trials...
    errorPercent totalErrors tempCorrectTrials correctTrials

%% PLOT SUCCESSFUL TRIALS

for i=1:numAllConditions_withRest
    tempCondIdx=strcmp(session.Condition,allConditions_withRest(i));
    tempNumCond=sum(tempCondIdx);
    tempNumSuccessful=sum(successfulTrial(tempCondIdx)>0);
    tempPercentSuccessful=tempNumSuccessful/tempNumCond;
    
    percentSuccessfulCummulative(i)=tempPercentSuccessful*100;
end
subplot(1,2,2)
bar(percentSuccessfulCummulative(~isnan(percentSuccessfulCummulative)),0.5)
ylabel('Percent (%)')
set(gca,'XTickLabel',allConditions_withRest(~isnan(percentSuccessfulCummulative)),...
    'XTick',1:numel(allConditions_withRest(~isnan(percentSuccessfulCummulative))),...
    'FontSize', 16, 'box','off','xTickLabelRotation',45)
ylim([0,100])
title('Percent of Successful Trials')

% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},percentSuccessfulCummulative];
clearvars i tempCondIdx tempNumCond tempNumSuccessful...
    tempPercentSuccessful percentSuccessfulCummulative

%% PLOT REACTION TIME
h2=figure('units','normalized','outerposition',[0 0 1 1],'Name',[filename, ': only correct trials included'],'NumberTitle','off');

% Isolate successful trials for graphing.
successfulTrial_noRest=cellfun(@str2num,session_noRest.SuccessfulTrial,'UniformOutput',false);
successfulIdx=find(~cellfun(@isempty,successfulTrial_noRest));
session_noRest=session_noRest(successfulIdx,:);

% Maybe not neccessary. Remove rows with a NaN in the SuccessfulTrial
% column. (Just an additional safegaurd.)
%successfulTrial_noRest=cellfun(@str2num,session_noRest.SuccessfulTrial,'UniformOutput',false);
%successfulIdx=~cellfun(@isnan,successfulTrial_noRest);
%session_noRest=session_noRest(successfulIdx,:);

for i=1:numConditions
    tempCondIdx=strcmp(session_noRest.Condition,conditions(i));
    condRxtTime=session_noRest.Reaction_time(tempCondIdx);
    averageCondRxtTime=nanmean(condRxtTime);
    semCondRxtTime=nanstd(condRxtTime)/sqrt(length(condRxtTime));
    
    avgRxtTimeCummulative(i)=averageCondRxtTime;
    semRxtTimeCummulative(i)=semCondRxtTime;
    
    rxtTimeArray{i}=condRxtTime;
end
subplot(2,3,1)
bar(avgRxtTimeCummulative(~isnan(avgRxtTimeCummulative)))
hold on;
h=errorbar(avgRxtTimeCummulative(~isnan(avgRxtTimeCummulative)),...
    semRxtTimeCummulative(~isnan(avgRxtTimeCummulative)),'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
ylabel('Milliseconds(ms)')
set(gca,'XTickLabel',conditions(~isnan(avgRxtTimeCummulative)),...
    'XTick',1:numel(conditions(~isnan(avgRxtTimeCummulative))),...
    'FontSize', 10, 'box','off','xTickLabelRotation',45)
title('Average Reaction Time')


% Overlay reaction times for individual trials
hold on;
yMax=-Inf;
xLoc=1;
for i=find(~isnan(avgRxtTimeCummulative))
    plot(xLoc*ones(length(rxtTimeArray{i}),1),rxtTimeArray{i},'or')
    tempYMax=max(rxtTimeArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
    xLoc=xLoc+1;
end
% end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.6]);

allRTs = cellfun(@(a) session_noRest.Reaction_time(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
maxLength = max(cellfun(@length, allRTs));
allRTs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allRTs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allRTs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.Reaction_time(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.Reaction_time(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgRxtTimeCummulative,semRxtTimeCummulative];
clearvars averageCondRxtTime avgRxtTimeCummulative condRxtTime h i...
    reactionTimeIdx rxtTimeArray stdCondRxtTime stdRxtTimeCummulative...
    tempCondIdx tempsession

%% PLOT REACH DURATION

for i=1:numConditions
    tempCondIdx=strcmp(session_noRest.Condition,conditions(i));
    condReachDur=session_noRest.ReachDuration(tempCondIdx);
    averageCondReachDur=nanmean(condReachDur);
    semCondReachDur=nanstd(condReachDur)/sqrt(length(condReachDur));
    
    avgReachDurCummulative(i)=averageCondReachDur;
    semReachDurCummulative(i)=semCondReachDur;
    
    reachDurArray{i}=condReachDur;
end
subplot(2,3,2)
bar(avgReachDurCummulative(~isnan(avgReachDurCummulative)))
hold on;
h=errorbar(avgReachDurCummulative(~isnan(avgReachDurCummulative)),...
    semReachDurCummulative(~isnan(avgReachDurCummulative)),'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
ylabel('Milliseconds (ms)')
set(gca,'XTickLabel',conditions(~isnan(avgReachDurCummulative)),'XTick',...
    1:numel(conditions(~isnan(avgReachDurCummulative))),'FontSize',10,...
    'box','off','xTickLabelRotation',45)
title('Average Reach Duration')

% Overlay reach durations for individual trials
hold on;
yMax=-Inf;
xLoc=1;
for i=find(~isnan(avgReachDurCummulative))
    plot(xLoc*ones(length(reachDurArray{i}),1),reachDurArray{i},'or')
    tempYMax=max(reachDurArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
    xLoc=xLoc+1;
end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.2]);

allRs = cellfun(@(a) session_noRest.ReachDuration(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
allRs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allRs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allRs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.ReachDuration(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.ReachDuration(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgReachDurCummulative,semReachDurCummulative];
clearvars averageCondReachDur avgReachDurCummulative condReachDur...
    reachDurArray reachDurationIdx stdCondReachDur...
    stdReachDurCummulative tempsession h i tempCondIdx

%% PLOT GRASP DURATION
for i=1:numConditions
    tempCondIdx=strcmp(session_noRest.Condition,conditions(i));
    condGraspDuration=session_noRest.GraspDuration(tempCondIdx);
    averageCondGraspDur=nanmean(condGraspDuration);
    semCondGraspDur=nanstd(condGraspDuration)/sqrt(length(condGraspDuration));
    
    avgGraspDurCummulative(i)=averageCondGraspDur;
    semGraspDurCummulative(i)=semCondGraspDur;
    
    graspDurArray{i}=condGraspDuration;
end
subplot(2,3,3)
bar(avgGraspDurCummulative(~isnan(avgGraspDurCummulative)))
hold on;
h=errorbar(avgGraspDurCummulative(~isnan(avgGraspDurCummulative)),...
    semGraspDurCummulative(~isnan(avgGraspDurCummulative)),'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
ylabel('Milliseconds (ms)')
set(gca,'XTickLabel',conditions(~isnan(avgGraspDurCummulative)),'XTick',...
    1:numel(conditions(~isnan(avgGraspDurCummulative))),'FontSize',10,...
    'box','off','xTickLabelRotation',45)
title('Average Grasp Duration')

% Overlay reach durations for individual trials
hold on;
yMax=-Inf;
xLoc=1;
for i=find(~isnan(avgGraspDurCummulative))
    plot(xLoc*ones(length(graspDurArray{i}),1),graspDurArray{i},'or')
    tempYMax=max(graspDurArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
    xLoc=xLoc+1;
end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.2]);

allGs = cellfun(@(a) session_noRest.GraspDuration(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
allGs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allGs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allGs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.GraspDuration(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.GraspDuration(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgGraspDurCummulative,semGraspDurCummulative];
clearvars averageCondGraspDur avgGraspDurCummulative condGraspDur...
    graspDurArray graspDurationIdx stdCondGraspDur...
    stdReachDurCummulative tempsession h i tempCondIdx
%% PLOT LIFT DURATION

% Remove conditions that do not have a lift duration
tempsession=session_noRest(~strcmp(session_noRest.Condition,'Photocell'),:);
tempsession=tempsession(~strcmp(tempsession.Condition,'Empty'),:);

tempConditions=conditions(~strcmp(conditions,'Photocell'));
tempConditions=tempConditions(~strcmp(tempConditions,'Empty'));
tempNumConditions=length(tempConditions);

for i=1:tempNumConditions
    tempCondIdx=strcmp(tempsession.Condition,tempConditions(i));
    condLiftDur=tempsession.LiftDuration(tempCondIdx);
    averageCondLiftDur=nanmean(condLiftDur);
    semCondLiftDur=nanstd(condLiftDur)/sqrt(length(condLiftDur));
    
    avgLiftDurCummulative(i)=averageCondLiftDur;
    semLiftDurCummulative(i)=semCondLiftDur;
    
    liftDurArray{i}=condLiftDur;
end
subplot(2,3,4)
bar([1:length(avgLiftDurCummulative)],avgLiftDurCummulative)
hold on;
h=errorbar([1:length(avgLiftDurCummulative)],avgLiftDurCummulative,...
    semLiftDurCummulative,'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
ylabel('Milliseconds (ms)')
set(gca,'XTickLabel',tempConditions,'XTick',1:numel(tempConditions),...
    'FontSize', 10, 'box','off','xTickLabelRotation',45)
title('Average Lift Duration')

% Overlay lift durations for individual trials
hold on;
yMax=-Inf;
for i=1:tempNumConditions
    plot(i*ones(length(liftDurArray{i}),1),liftDurArray{i},'or')
    tempYMax=max(liftDurArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.2]);

allLs = cellfun(@(a) session_noRest.LiftDuration(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
allLs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allLs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allLs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.LiftDuration(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.LiftDuration(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgLiftDurCummulative,semLiftDurCummulative];
clearvars averageCondLiftDur avgLiftDurCummulative condLiftDur h i...
    liftDurArray liftDurationIdx stdCondLiftDur stdLiftDurCummulative...
    tempCondIdx tempConditions tempNumConditions tempsession

%% PLOT HOLD DURATION

for i=1:numConditions
    tempCondIdx=strcmp(session_noRest.Condition,conditions(i));
    condHoldDur=session_noRest.HoldDuration(tempCondIdx);
    averageCondHoldDur=nanmean(condHoldDur);
    semCondHoldDur=nanstd(condHoldDur)/sqrt(length(condHoldDur));
    
    avgHoldDurCummulative(i)=averageCondHoldDur;
    semHoldDurCummulative(i)=semCondHoldDur;
    
    holdDurArray{i}=condHoldDur;
end
subplot(2,3,5)
bar(avgHoldDurCummulative(~isnan(avgHoldDurCummulative)))
hold on;
h=errorbar(avgHoldDurCummulative(~isnan(avgHoldDurCummulative)),...
    semHoldDurCummulative(~isnan(avgHoldDurCummulative)),'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
xlabel('Condition')
ylabel('Milliseconds (ms)')
set(gca,'XTickLabel',conditions(~isnan(avgHoldDurCummulative)),...
    'XTick',1:numel(conditions(~isnan(avgHoldDurCummulative))),'FontSize',...
    10, 'box','off','xTickLabelRotation',45)
title('Average Hold Duration')

% Overlay hold durations for individual trials
hold on;
yMax=-Inf;
xLoc=1;
for i=find(~isnan(avgHoldDurCummulative))
    plot(xLoc*ones(length(holdDurArray{i}),1),holdDurArray{i},'or')
    tempYMax=max(holdDurArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
    xLoc=xLoc+1;
end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.2]);

allHs = cellfun(@(a) session_noRest.HoldDuration(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
allHs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allHs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allHs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.HoldDuration(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.HoldDuration(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgHoldDurCummulative,semHoldDurCummulative];
clearvars averageCondHoldDur avgHoldDurCummulative condHoldDur h...
    holdDurArray holdDurationIdx i stdCondHoldDur stdHoldDurCummulative...
    tempCondIdx tempsession

%% PLOT WITHDRAWAL DURATION

for i=1:numConditions
    tempCondIdx=strcmp(session_noRest.Condition,conditions(i));
    condWithdrawalDur=session_noRest.WithdrawalDuration(tempCondIdx);
    averageCondWithdrawalDur=nanmean(condWithdrawalDur);
    semCondWithdrawalDur=nanstd(condWithdrawalDur)/sqrt(length(condWithdrawalDur));
    
    avgWithdrawalDurCummulative(i)=averageCondWithdrawalDur;
    semWithdrawalDurCummulative(i)=semCondWithdrawalDur;
    
    withdrawalDurArray{i}=condWithdrawalDur;
end
subplot(2,3,6)
bar(avgWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)))
hold on;
h=errorbar(avgWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)),...
    semWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)),'k','lineWidth',2);
set(h,'linestyle','none') %just keep error bars
xlabel('Condition')
ylabel('Milliseconds (ms)')
set(gca,'XTickLabel',conditions(~isnan(avgWithdrawalDurCummulative)),'XTick',1:numel(conditions(~isnan(avgWithdrawalDurCummulative))),'FontSize',...
    10, 'box','off','xTickLabelRotation',45)
title('Average Withdrawal Duration')

% Overlay withdrawal durations for individual trials
hold on;
yMax=-Inf;
xLoc=1;
for i=find(~isnan(avgWithdrawalDurCummulative))
    plot(xLoc*ones(length(withdrawalDurArray{i}),1),withdrawalDurArray{i},'or')
    tempYMax=max(withdrawalDurArray{i});
    if tempYMax>yMax
        yMax=tempYMax;
    end
    xLoc=xLoc+1;
end
ylim([0,yMax])
ylim([0,(max(h.YData) + max(h.YPositiveDelta))*1.2]);

allWs = cellfun(@(a) session_noRest.WithdrawalDuration(strcmp(session_noRest.Condition,a)),mainEffect, 'UniformOutput', false);
allWs = cell2mat(cellfun(@(a) fillNans(a,maxLength), allWs, 'UniformOutput', false));
if(isempty(mainEffect) || anova1(allWs,[],'off')<alphaVal)
    if(~isempty(mainEffect))
        correctedLength = length(sigComps);
    else
        correctedLength = 1;
    end
    for s = 1:length(sigComps)
        dist1 = session_noRest.WithdrawalDuration(strcmp(session_noRest.Condition,sigComps{s}{1}));
        dist2 = session_noRest.WithdrawalDuration(strcmp(session_noRest.Condition,sigComps{s}{2}));
        if(any(dist1) && any(dist2))
            if(length(dist1)~=length(dist2))
                dist1(end:max(length(dist1),length(dist2))) = NaN;
                dist2(end:max(length(dist1),length(dist2))) = NaN;
            end
            if(anova1([dist1,dist2],[],'off')<alphaVal/correctedLength)
                sigline([find(strcmp(sigComps{s}{1},allConditions_withRest)), ...
                    find(strcmp(sigComps{s}{2},allConditions_withRest))],['p<',num2str(alphaVal)],max(h.YData));
            end
        end
    end
end
% Data for GilliganCummulativeData
spreadsheetData{2}=[spreadsheetData{2},avgWithdrawalDurCummulative,semWithdrawalDurCummulative];
clearvars averageCondWithdrawalDur avgWithdrawalDurCummulative...
    condWithdrawalDur h i semCondWithdrawalDur semWithdrawalDurCummulative...
    tempCondIdx withdrawalDurArray withdrawalDurationIdx

%% CREATE TABLE WITH  ONLY CORRECT TRIAL INFORMATION TO SAVE DAILY
% Remove failed trials.
falseTrialIdx=cellfun(@str2num,session.SuccessfulTrial,'UniformOutput',false);
correctTrialIdx=find(~cellfun(@isempty,falseTrialIdx));
successfulTrial=session(correctTrialIdx,4:10);

%allConditions_withRest={'Large Sphere','Photocell','Empty','Rest','Extra Small Sphere','Small Sphere'};

sortedArduinoData.Date=filename;

for i=1:numAllConditions_withRest
    tempCondIdx=strcmp(successfulTrial.Condition,allConditions_withRest(i));
    if sum(tempCondIdx)~=0
        tempCondition=allConditions_withRest{i};
        tempCondition=tempCondition(~isspace(tempCondition));
        
        if strcmp(allConditions_withRest(i),'Photocell')
            sortedArduinoData.(tempCondition)=successfulTrial(tempCondIdx,[1:4,6:7]);
        else
            sortedArduinoData.(tempCondition)=successfulTrial(tempCondIdx,:);
        end
    end
end



% SAVE
if(~exist([filePathDir,'Results\'],'dir'))
    mkdir([filePathDir,'Results\'])
end
% Save sortedArduinoData
save([filePathDir,'Results\',filename],'sortedArduinoData','-v7.3')
%save(['S:\Lab\', monkey,'\Sorted Arduino Data\',filename],'sortedArduinoData','-v7.3')

% Save graph of errors and successful trials.
saveas(h1,[filePathDir,'Results\ErrorsAndSuccessGraph.png']);

% Save graphs for reaction time, reach duration, lift duration,
% hold duration, and withdrawal duration for each condition.
saveas(h2,[filePathDir,'Results\SummaryGraphs_correctTrial.png']);

% Append to excel spreadsheet GilliganCummulativeData
if strcmp(appendToGilliganCummulativeData,'yes')
    cd(['S:\Lab\', monkey,'\All Data'])
    [num,txt,raw]=xlsread([monkey, 'CummulativeData'],'B3:B10000');
    row=length(txt)+3;
    range1=['B',num2str(row),':B',num2str(row)]; %for filename
    range2=['C',num2str(row),':CD',num2str(row)]; %for data
    xlswrite([monkey, 'CummulativeData.xlsx'],spreadsheetData(1),'Session Information',range1);
    xlswrite([monkey, 'CummulativeData.xlsx'],spreadsheetData{2},range2)
else
    disp('Data not added to GilliganCummulativeData');
end

clearvars -except sortedArduinoData allConditions_withRest alphaVal sigComps session_noRest
%%
if(0);
    %%
    figure();
    conds = allConditions_withRest(~strcmp(allConditions_withRest,'Rest'));
    react_to_grasp = cellfun(@(a) table2array(session_noRest(strcmp(session_noRest.Condition,a),5:7)), conds, 'UniformOutput', false);
    maxLength = max(cellfun(@length, react_to_grasp));
    react_to_grasp = cellfun(@(a) fillNans(a,maxLength), react_to_grasp, 'UniformOutput', false);
    reach_to_grasp = cellfun(@(a) sum(a(:,2:3),2), react_to_grasp, 'UniformOutput', false);
    react_to_grasp = cellfun(@(a) a(:,1), react_to_grasp, 'UniformOutput', false);
    
    %
    subplot(1,2,1);
    bar(cellfun(@nanmean,react_to_grasp));
    hold on
    e = errorbar(cellfun(@nanmean,react_to_grasp), cellfun(@nanstd,react_to_grasp)./sqrt(cellfun(@length,react_to_grasp)), 'k', 'LineWidth', 2);
    set(e, 'linestyle', 'none');
    set(gca,'XTickLabel', conds);
    set(gca,'xTickLabelRotation', 45);
    ylabel('Milliseconds (ms)');
    title('Reaction, Reach, Grasp');
    react_to_grasp = cell2mat(react_to_grasp);
    if(anova1(react_to_grasp,[],'off')<alphaVal)
        for s = 1:length(sigComps)
            sInd1 = find(strcmp(sigComps{s}(1),conds));
            sInd2 = find(strcmp(sigComps{s}(2),conds));
            if(anova1([react_to_grasp(:,sInd1),react_to_grasp(:,sInd2)],[],'off')<alphaVal/length(sigComps))
                sigline([sInd1, sInd2],['p<',num2str(alphaVal)],max(e.YData + e.YPositiveDelta));
            end
        end
    end
    
    subplot(1,2,2);
    bar(cellfun(@nanmean,reach_to_grasp));
    hold on
    e = errorbar(cellfun(@nanmean,reach_to_grasp), cellfun(@nanstd,reach_to_grasp)./sqrt(cellfun(@length,reach_to_grasp)), 'k', 'LineWidth', 2);
    set(e, 'linestyle', 'none');
    set(gca,'XTickLabel', conds);
    set(gca,'xTickLabelRotation', 45);
    ylabel('Milliseconds (ms)');
    title('Reach, Grasp');
    reach_to_grasp = cell2mat(reach_to_grasp);
    if(anova1(reach_to_grasp,[],'off')<alphaVal)
        for s = 1:length(sigComps)
            sInd1 = find(strcmp(sigComps{s}(1),conds));
            sInd2 = find(strcmp(sigComps{s}(2),conds));
            if(anova1([reach_to_grasp(:,sInd1),reach_to_grasp(:,sInd2)],[],'off')<alphaVal/length(sigComps))
                sigline([sInd1, sInd2],['p<',num2str(alphaVal)],max(e.YData + e.YPositiveDelta));
            end
        end
    end
    
    figure()
    subplot(2,1,1);
    plot(react_to_grasp);
    legend(conds);
    title('Reaction, Reach, Grasp')
    xlabel('Trial')
    ylabel('Milliseconds (ms)')
    subplot(2,1,2);
    plot(reach_to_grasp);
    title('Reach, Grasp')
    xlabel('Trial')
    ylabel('Milliseconds (ms)')
    
    figure();
    subplot(2,1,1);
    for i = 1:floor(length(react_to_grasp)/5)
        blocksRG(i,:) = mean(react_to_grasp(5*(i-1)+1:5*(i-1)+5,:));
        blocksRGS(i,:) = std(react_to_grasp(5*(i-1)+1:5*(i-1)+5,:))./sqrt(5);
    end
    hbar = bar(blocksRG);
    legend(conds, 'Orientation','horizontal','Location', 'north');
    hold on
    for k1 = 1:size(blocksRG,2)
        ctr(k1,:) = bsxfun(@plus, hbar(1).XData, hbar(k1).XOffset');
    end
    hold on
    errorbar(ctr', blocksRG, blocksRGS, '.k')
    title('Reaction, Reach, Grasp')
    xlabel('Blocks of 5')
    ylabel('Milliseconds (ms)')
    
    subplot(2,1,2);
    for i = 1:floor(length(react_to_grasp)/5)
        blocksRHG(i,:) = mean(reach_to_grasp(5*(i-1)+1:5*(i-1)+5,:));
        blocksRGHS(i,:) = std(reach_to_grasp(5*(i-1)+1:5*(i-1)+5,:))./sqrt(5);
    end
    hbar = bar(blocksRHG);
    hold on
    for k1 = 1:size(blocksRHG,2)
        ctr(k1,:) = bsxfun(@plus, hbar(1).XData, hbar(k1).XOffset');
    end
    hold on
    errorbar(ctr', blocksRHG, blocksRGHS, '.k')
    title('Reach, Grasp')
    xlabel('Blocks of 5')
    ylabel('Milliseconds (ms)')
    
    figure();
    react_to_grasp = react_to_grasp(:,1:2);
    reach_to_grasp = reach_to_grasp(:,1:2);
    for i = 1:2
        aa(1,i) = nanmean(react_to_grasp(:,i));
        aa(2,i) = nanmean(reach_to_grasp(:,i));
        aaS(1,i) = nanstd(react_to_grasp(:,i))./sqrt(length(react_to_grasp(:,i)));
        aaS(2,i) = nanstd(reach_to_grasp(:,i))./sqrt(length(reach_to_grasp(:,i)));
    end
    hbar = bar(aa);
    hold on
    for k1 = 1:2
        cc(k1,:) = bsxfun(@plus, hbar(k1).XData, hbar(k1).XOffset);
    end
    hold on
    e =errorbar(cc', aa, aaS, '.k');

    xticklabels({'Reaction Time', 'Reach-to-Grasp Time'});
    ylabel('Time (ms)')
    if(anova1(reach_to_grasp,[],'off')<alphaVal)
        for s = 1:1
            sInd1 = find(strcmp(sigComps{s}(1),conds));
            sInd2 = find(strcmp(sigComps{s}(2),conds));
            if(anova1([reach_to_grasp(:,sInd1),reach_to_grasp(:,sInd2)],[],'off')<alphaVal/length(sigComps))
                sigline([2-hbar(1).XOffset, 2+hbar(1).XOffset],['p<',num2str(alphaVal)],min([e.YData] + [e.YPositiveDelta])+19);
            end
        end
    end
    
     if(anova1(react_to_grasp,[],'off')<alphaVal)
        for s = 1:1
            sInd1 = find(strcmp(sigComps{s}(1),conds));
            sInd2 = find(strcmp(sigComps{s}(2),conds));
            if(anova1([react_to_grasp(:,sInd1),react_to_grasp(:,sInd2)],[],'off')<alphaVal/length(sigComps))
                sigline([1-hbar(1).XOffset, 1+hbar(1).XOffset],['p<',num2str(alphaVal)],max([e.YData] + [e.YPositiveDelta]));
            end
        end
    end
end


function filled = fillNans(vec,maxLength)
filled = vec;
filled(end:maxLength,:) = NaN;
end
%%%%%%%%%%%%
% clearvars -except allfilenames
% end