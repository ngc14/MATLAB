%WithinSessionGraphs.m
%CREATED: 5/17/2017
%NEED TO DOUBLE CHECK THAT THIS WORKS

% Create graphs to compare the percent error,reaction time, reach duration,
% lift duration, hold duration, and withdrawal duration across multiple
% sessions.

%% VARIABLES FOR USER MODIFICATION
% Refer to GilliganCummulativeData to determine the date range that you
% want to graph. Modify the cell range to extract data from the desired
% sessions.
monkey = 'Gilligan';
startRow = '3';
endRow = '43';
data=xlsread(['\\130.49.229.252\gharbawie\Lab\', monkey,'\Sorted Arduino Data\',monkey, 'CummulativeData.xlsx'],['C',startRow,':CD',endRow]);
[num,filenames,raw]=xlsread(['\\130.49.229.252\gharbawie\Lab\', monkey,'\Sorted Arduino Data\',monkey, 'CummulativeData.xlsx'],['B',startRow,':B',endRow]);

%% ORGANIZE DATA INTO TABLES
startIndex = 9;
% Extract the dates selected to use for x-axis labels
for i=1:length(filenames)
    tempFilename=filenames{i};
    nameIdx=strfind(tempFilename,'Gilligan');
    if nameIdx==1
        dateString=tempFilename(10:end);
    else
        dateString=tempFilename(1:10);
    end
    dateMatrix(i,1)=datenum(dateString,'mm_dd_yyyy');
    dateStringArray{i,1}=datestr(dateMatrix(i,1));
end

% trials=data(:,1)
errors=data(:,startIndex:startIndex+7);
errors=array2table(errors,'VariableNames',{'FalseStart','FailedReach',...
    'FailedGrasp','FailedLift','FailedHold','FailedReplace',...
    'FailedReplaceHold','FailedRest'});
% totalErrors=data(:,10);
% errorPercent=data(:,11);
startIndex = startIndex + 10;

percentSuccessful=data(:,startIndex:startIndex+5);
percentSuccessful=array2table(percentSuccessful,'VariableNames',...
    {'LargeSphere','Photocell','Empty','Rest','SmallCube','SmallSphere'});
startIndex = startIndex + 6;

averageReactionTime=data(:,startIndex:startIndex+4);
averageReactionTime=array2table(averageReactionTime,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
stdAverageReactionTime=data(:,startIndex+5:startIndex+9);
stdAverageReactionTime=array2table(stdAverageReactionTime,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
startIndex = startIndex + 10;

averageReachDuration=data(:,startIndex:startIndex+4);
averageReachDuration=array2table(averageReachDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
stdAverageReachDuration=data(:,startIndex+5:startIndex+9);
stdAverageReachDuration=array2table(stdAverageReachDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
startIndex = startIndex + 10;

averageGraspDuration=data(:,startIndex:startIndex+4);
averageGraspDuration=array2table(averageGraspDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
stdAverageGraspDuration=data(:,startIndex+5:startIndex+9);
stdAverageGraspDuration=array2table(stdAverageGraspDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
startIndex = startIndex + 10;

averageLiftDuration=data(:,startIndex:startIndex+2);
averageLiftDuration=array2table(averageLiftDuration,'VariableNames',...
    {'LargeSphere','ExtraSmallSphere','SmallSphere'});
stdAverageLiftDuration=data(:,startIndex+3:startIndex+5);
stdAverageLiftDuration=array2table(stdAverageLiftDuration,'VariableNames',...
    {'LargeSphere','ExtraSmallSphere','SmallSphere'});
startIndex = startIndex + 6;

averageHoldDuration=data(:,startIndex:startIndex+4);
averageHoldDuration=array2table(averageHoldDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
stdAverageHoldDuration=data(:,startIndex+5:startIndex+9);
stdAverageHoldDuration=array2table(stdAverageHoldDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
startIndex = startIndex + 10;

averageWithdrawalDuration=data(:,startIndex:startIndex+4);
averageWithdrawalDuration=array2table(averageWithdrawalDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});
stdAverageWithdrawalDuration=data(:,startIndex+5:startIndex+9);
stdAverageWithdrawalDuration=array2table(stdAverageWithdrawalDuration,'VariableNames',...
    {'LargeSphere','Photocell','Empty','ExtraSmallSphere','SmallSphere'});

%% PLOT THE DATA

% ERRORS
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:size(errors,2)
    plot(dateMatrix,table2array(errors(:,i)),'o-')
    hold on
end
legend(errors.Properties.VariableNames)
set(gca,'FontSize',16,'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
xlabel('Session')
ylabel('Number of Occurances')
title('Number of Each Error Type')
ylim([0 200])

% PERCENT SUCCESSFUL
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:size(percentSuccessful,2)
    subplot(2,4,1)
    plot(dateMatrix,table2array(percentSuccessful(:,i)),'.-','MarkerSize',15)
    hold on
end
% legend(percentSuccessful.Properties.VariableNames)
set(gca,'FontSize',10,'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
ylabel('Percent (%)')
ylim([0 100])
title('Percent of Successful Trials')

% AVERAGE REACTION TIME
for i=1:size(averageReactionTime,2)
    subplot(2,4,2)
    plot(dateMatrix,table2array(averageReactionTime(:,i)),'.-','MarkerSize',15)
%         errorbar(table2array(averageReactionTime(:,i)),...
%             table2array(stdAverageReactionTime(:,i)),'o-')
    hold on
end
legend(averageReactionTime.Properties.VariableNames)
set(gca,'FontSize',10,'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
ylabel('Milliseconds (ms)')
ylim([0 2000])
title('Average Reaction Time')

% AVERAGE REACH DURATION
for i=1:size(averageReachDuration,2)
    subplot(2,4,3)
    plot(dateMatrix,table2array(averageReachDuration(:,i)),'.-','MarkerSize',15)
    %     errorbar(table2array(averageReachDuration(:,i)),...
    %         table2array(stdAverageReachDuration(:,i)),'o-')
    hold on
end
legend(averageReachDuration.Properties.VariableNames)
set(gca,'FontSize',10,'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
ylabel('Milliseconds (ms)')
ylim([0 800])
title('Average Reach Duration')

% AVERAGE GRASP DURATION
for i=1:size(averageGraspDuration,2)
    subplot(2,4,4)
    plot(dateMatrix,table2array(averageGraspDuration(:,i)),'.-','MarkerSize',15)
    %     errorbar(table2array(averageGraspDuration(:,i)),...
    %         table2array(stdAverageGraspDuration(:,i)),'o-')
    hold on
end
legend(averageGraspDuration.Properties.VariableNames)
set(gca,'FontSize',10,'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
ylabel('Milliseconds (ms)')
ylim([0 800])
title('Average Grasp Duration')

% AVERAGE LIFT DURATION
for i=1:size(averageLiftDuration,2)
    subplot(2,4,5)
    plot(dateMatrix,table2array(averageLiftDuration(:,i)),'.-','MarkerSize',15)
    %     errorbar(table2array(averageLiftDuration(:,i)),...
    %         table2array(stdAverageLiftDuration(:,i)),'o-')
    hold on
end
legend(averageLiftDuration.Properties.VariableNames)
set(gca,'FontSize',10, 'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
xlabel('Session')
ylabel('Milliseconds (ms)')
ylim([0 200])
title('Average Lift Duration')

% AVERAGE HOLD DURATION
for i=1:size(averageHoldDuration,2)
    subplot(2,4,6)
    plot(dateMatrix,table2array(averageHoldDuration(:,i)),'.-','MarkerSize',15)
    %     errorbar(table2array(averageHoldDuration(:,i)),...
    %         table2array(stdAverageHoldDuration(:,i)),'o-')
    hold on
end
legend(averageHoldDuration.Properties.VariableNames)
set(gca,'FontSize',10, 'box','off','XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
xlabel('Session')
ylabel('Milliseconds (ms)')
ylim([0 4000])
title('Average Hold Duration')

% AVERAGE WITHDRAWAL DURATION
for i=1:size(averageWithdrawalDuration,2)
    subplot(2,4,7)
    plot(dateMatrix,table2array(averageWithdrawalDuration(:,i)),'.-','MarkerSize',15)
    %     errorbar(table2array(averageWithdrawalDuration(:,i)),...
    %         table2array(stdAverageWithdrawalDuration(:,i)),'o-')
    hold on
end
legend(averageWithdrawalDuration.Properties.VariableNames)
set(gca,'FontSize',10, 'box','off','XTickLabel',...
    1:size(averageWithdrawalDuration,1),'XTickLabel',dateStringArray,'XTick',dateMatrix)
xtickangle(90)
xlabel('Session')
ylabel('Milliseconds (ms)')
ylim([0 800])
title('Average Withdrawal Duration')