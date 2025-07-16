close all;
unitTypeName = "unitType_E";
physTableCond =  readtable('S:\Lab\ngc14\Working\Both\Full_Baseline\physTable_Cond.xlsx');
savePath = "S:\Lab\ngc14\Working\Both\Full_Baseline\Misc\unitTypes\groupedBySomatotopy\"+unitTypeName+"\Count\";
tb = groupsummary(physTableCond(physTableCond.TaskUnits~=0,:),["Somatotopy",unitTypeName,"Condition"]);
filteredTable = tb;
filteredTable = filteredTable(filteredTable.(unitTypeName)~=0 & filteredTable.Condition == 3,:);
% Create pivoted table
pivotedData = pivot(filteredTable, Rows='Somatotopy', Columns=unitTypeName, ...
    DataVariable="GroupCount", Method="sum", RowLabelPlacement="rownames");
rowNames = reordercats(categorical(pivotedData.Properties.RowNames),["Hand","Arm", "Trunk"]);
colNames = reordercats(categorical(pivotedData.Properties.VariableNames),[1 3 2]);
data = table2array(pivotedData);
data2=data(:,[1 3 2]);
%data2=data2./sum(data,2);
bh = bar(rowNames,data2,"stacked");
hold on;
set(gca,TickLabelInterpreter="none");
% Compute the height of each segment and write text to plot
text(reshape([bh.XEndPoints],1,[]), reshape(cumsum(data2,2),1,[]), ...
    reshape(compose('%.3f',[bh.YData]),1,[]), 'Color', 'k', ...
    'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
colors = cell2mat(cellfun(@(s) str2num(s),unique(cellfun(@num2str,num2cell(vertcat(bh.CData),2),'UniformOutput',false)),"uniformoutput",false));
set(bh, 'CData', colors);
l=legend(string([1 3 2]));
title("Photocell");
saveFigures(gcf,savePath,"P_Somatotopy_Type",[]);