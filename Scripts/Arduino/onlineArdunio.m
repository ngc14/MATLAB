clear all;
dateToday = datestr(datetime('now', 'Format', 'MM_dd_y'), 'mm_dd_yyyy');
tempFileName = ['C:\Users\Omar Lab\Documents\tempArduino.txt'];
filePath = ['C:\Users\Omar Lab\Documents\Monkey Training\Gilligan_Macaque_110-16\Gilligan_', dateToday, '\Gilligan_', dateToday, '.txt'];
conds = {'Extra Small Sphere','Stim_Extra Small Sphere','Photocell','Stim_Photocell','Stim_Rest'};
delimiter = '\t';
formatSpec = '%*s%*s%*s%s%s%s%s%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
h2=figure('units','normalized');

%%
if(1)
    while(true)
        copyfile(filePath, tempFileName, 'f');
        %%
        % Open the text file.
        fileID = fopen(tempFileName,'r');
        
        % Read columns of data according to the format.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
        smallestLength = min(cellfun(@length, dataArray));
        dataArray = cellfun(@(a) a(1:smallestLength), dataArray, 'UniformOutput', false);
        % Close the text file.
        fclose(fileID);
        
        % Create output variable
        recordedTrial = [dataArray{1:end-1}];
        
        % remove lines that are not task related (pause, headers)
        recordedTrial = recordedTrial(ismember(recordedTrial(:,1), conds),:);
        
        % Remove failed trials.
        falseTrialIdx=cellfun(@str2num,recordedTrial(:,8),'UniformOutput',false);
        recordedTrial = recordedTrial(find(~cellfun(@isempty,falseTrialIdx)),:);
        clf;
        %% REACTION TIME
        for i=1:length(conds)
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condRxtTime=recordedTrial(tempCondIdx,2);
            condRxtTime = cell2mat(cellfun(@str2num, condRxtTime, 'UniformOutput', false));
            averageCondRxtTime=mean(condRxtTime);
            semCondRxtTime=std(condRxtTime)/sqrt(length(condRxtTime));
            
            avgRxtTimeCummulative(i)=averageCondRxtTime;
            semRxtTimeCummulative(i)=semCondRxtTime;
        end
        subplot(2,3,1)
        bar(avgRxtTimeCummulative(~isnan(avgRxtTimeCummulative)))
        hold on;
        h=errorbar(avgRxtTimeCummulative(~isnan(avgRxtTimeCummulative)),...
            semRxtTimeCummulative(~isnan(avgRxtTimeCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        ylabel('Milliseconds(ms)')
        set(gca,'XTickLabel',conds(~isnan(avgRxtTimeCummulative)),...
            'XTick',1:numel(conds(~isnan(avgRxtTimeCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Reaction Time')
        
        %% REACH TIME
        for i=1:length(conds)
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condReachDur=recordedTrial(tempCondIdx,3);
            condReachDur = cell2mat(cellfun(@str2num, condReachDur, 'UniformOutput', false));
            averageCondReachDur=mean(condReachDur);
            semCondReachDur=std(condReachDur)/sqrt(length(condReachDur));
            
            avgReachDurCummulative(i)=averageCondReachDur;
            semReachDurCummulative(i)=semCondReachDur;
        end
        subplot(2,3,2)
        bar(avgReachDurCummulative(~isnan(avgReachDurCummulative)))
        hold on;
        h=errorbar(avgReachDurCummulative(~isnan(avgReachDurCummulative)),...
            semReachDurCummulative(~isnan(avgReachDurCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        ylabel('Milliseconds (ms)')
        set(gca,'XTickLabel',conds(~isnan(avgReachDurCummulative)),...
            'XTick',1:numel(conds(~isnan(avgReachDurCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Reach Duration')
        
        %% GRASP TIME
        for i=1:length(conds)
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condGraspDuration=recordedTrial(tempCondIdx,4);
            condGraspDuration = cell2mat(cellfun(@str2num, condGraspDuration, 'UniformOutput', false));
            averageCondGraspDur=mean(condGraspDuration);
            semCondGraspDur=std(condGraspDuration)/sqrt(length(condGraspDuration));
            
            avgGraspDurCummulative(i)=averageCondGraspDur;
            semGraspDurCummulative(i)=semCondGraspDur;
        end
        subplot(2,3,3)
        bar(avgGraspDurCummulative(~isnan(avgGraspDurCummulative)))
        hold on;
        h=errorbar(avgGraspDurCummulative(~isnan(avgGraspDurCummulative)),...
            semGraspDurCummulative(~isnan(avgGraspDurCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        ylabel('Milliseconds (ms)')
        set(gca,'XTickLabel',conds(~isnan(avgGraspDurCummulative)),...
            'XTick',1:numel(conds(~isnan(avgGraspDurCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Grasp Duration')
        
        %% LIFT TIME
        for i=1:length(conds)
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condLiftDur=recordedTrial(tempCondIdx,5);
            condLiftDur = cell2mat(cellfun(@str2num, condLiftDur, 'UniformOutput', false));
            averageCondLiftDur=mean(condLiftDur);
            semCondLiftDur=std(condLiftDur)/sqrt(length(condLiftDur));
            
            avgLiftDurCummulative(i)=averageCondLiftDur;
            semLiftDurCummulative(i)=semCondLiftDur;
        end
        subplot(2,3,4)
        bar(avgLiftDurCummulative(~isnan(avgLiftDurCummulative)))
        hold on;
        h=errorbar(avgLiftDurCummulative(~isnan(avgLiftDurCummulative)),...
            semLiftDurCummulative(~isnan(semLiftDurCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        ylabel('Milliseconds (ms)')
        set(gca,'XTickLabel',conds(~isnan(avgLiftDurCummulative)),...
            'XTick',1:numel(conds(~isnan(avgLiftDurCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Lift Duration')
        
        %% HOLD TIME
        for i=1:length(conds);
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condHoldDur=recordedTrial(tempCondIdx,6);
            condHoldDur = cell2mat(cellfun(@str2num, condHoldDur, 'UniformOutput', false));
            averageCondHoldDur=mean(condHoldDur);
            semCondHoldDur=std(condHoldDur)/sqrt(length(condHoldDur));
            
            avgHoldDurCummulative(i)=averageCondHoldDur;
            semHoldDurCummulative(i)=semCondHoldDur;
        end
        subplot(2,3,5)
        bar(avgHoldDurCummulative(~isnan(avgHoldDurCummulative)))
        hold on;
        h=errorbar(avgHoldDurCummulative(~isnan(avgHoldDurCummulative)),...
            semHoldDurCummulative(~isnan(avgHoldDurCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        xlabel('Condition')
        ylabel('Milliseconds (ms)')
        set(gca,'XTickLabel',conds(~isnan(avgHoldDurCummulative)),...
            'XTick',1:numel(conds(~isnan(avgHoldDurCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Hold Duration')
        
        %% WITHDRAWAL TIME
        for i=1:length(conds)
            tempCondIdx=strcmp(recordedTrial(:,1),conds(i));
            condWithdrawalDur=recordedTrial(tempCondIdx,7);
            condWithdrawalDur = cell2mat(cellfun(@str2num, condWithdrawalDur, 'UniformOutput', false));
            averageCondWithdrawalDur=mean(condWithdrawalDur);
            semCondWithdrawalDur=std(condWithdrawalDur)/sqrt(length(condWithdrawalDur));
            
            avgWithdrawalDurCummulative(i)=averageCondWithdrawalDur;
            semWithdrawalDurCummulative(i)=semCondWithdrawalDur;
        end
        subplot(2,3,6)
        bar(avgWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)))
        hold on;
        h=errorbar(avgWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)),...
            semWithdrawalDurCummulative(~isnan(avgWithdrawalDurCummulative)),'k','lineWidth',2);
        set(h,'linestyle','none') %just keep error bars
        xlabel('Condition')
        ylabel('Milliseconds (ms)')
        set(gca,'XTickLabel',conds(~isnan(avgWithdrawalDurCummulative)),...
            'XTick',1:numel(conds(~isnan(avgWithdrawalDurCummulative))),...
            'FontSize', 10, 'box','off','xTickLabelRotation',45)
        title('Average Withdrawal Duration')
        pause(60);
    end
end