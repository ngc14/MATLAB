function figHandle = modulatedUnitsPerRep(repLabels,unitValues,valueLabel,condNames)
xBuff = 0.5;
axColors = [.75 .45 0; 0 .4 .4;];

repColorLabels = fieldnames(MotorMapping.repColors);
repNames = repColorLabels(ismember(repColorLabels,unique(horzcat(repLabels(:)))));
jointInds = arrayfun(@(jN) arrayfun(@(mJ) logical(cellfun(@(ss) strcmp(ss,jN),...
    mJ,'UniformOutput',true)),repLabels,'UniformOutput', false),repNames,'UniformOutput', false);

spAx = [];spInd = 1;
figHandle = {};
for c = 1:length(condNames)
    for p = 1:size(unitValues,1)
        figHandle{p} = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
        %subplotInd = ((c-1)*size(unitValues,1)) + p;
        %spAx(subplotInd) = subplot(length(condNames),size(unitValues,1),subplotInd);hold on;
        spAx(spInd) = gcf;
        currUnitTypes = unitValues(p,:);
        if(any(contains(["logical", "uint32"],class(currUnitTypes))))
            set(gca,'ColorOrder',axColors);
            xlim([xBuff,(2*length(repNames))+xBuff]);
            xticks(1:2*length(repNames))
            xticklabels(repmat(repNames,1,2));
            refreshdata(spAx(spInd));
            %count
            yyaxis left;
            b = bar(1:length(repNames), cellfun(@(jI) nansum(currUnitTypes(...
                cell2mat(jI'))),jointInds),.4,'FaceColor', 'flat', 'EdgeColor',...
                axColors(1,:), 'LineWidth',2);
            b.CData = cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)), repNames, 'UniformOutput', false));

            ylabel(join(['# ', valueLabel(p)]),'Interpreter','none');
            %percent
            yyaxis right;
            b = bar(length(repNames)+1:(2*length(repNames))+xBuff, 100.*...
                [cellfun(@(jI) abs(nansum(currUnitTypes(cell2mat(jI'))))./...
                nansum(cell2mat(jI')),jointInds)],.4,'FaceColor', 'flat',...
                'EdgeColor', axColors(2,:), 'LineWidth',2);
            ylim([0 100]);
            ylabel(join(['%',valueLabel(p)]),'Interpreter','none');
            b.CData = cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)), repNames, 'UniformOutput', false));
        else
            currUnitTypes = vertcat(currUnitTypes(:));
            xlim([xBuff,length(repNames)+xBuff]);
            xticks(1:length(repNames));
            xticklabels(repNames);
            jointVals = cell2mat(cellfun(@(jI) currUnitTypes(cell2mat(jI')),jointInds, 'UniformOutput', false))';
            jointLabels = cellfun(@(s,r) repmat(string(r),1,sum(cell2mat(s))), jointInds, repNames, 'UniformOutput',false);
            b = boxplot(jointVals,horzcat(jointLabels{:}),'GroupOrder',string(repNames),'Colors',...
                cell2mat(cellfun(@(a) MotorMapping.repColors.(a), repNames,'UniformOutput',false)),...
                'Notch','on','BoxStyle','outline','Widths',.5,'PlotStyle','traditional',...
                'Whisker',1,'Symbol','');
            ob = findobj('LineStyle','--'); set(ob, 'LineStyle','-');
        end
        if(p==1)
            %y = ylabel(condNames{c})
        end
        if(c==1)
            title(valueLabel{p},'Interpreter','none')
        end
        spInd = spInd + 1;
    end
end
%linkaxes(spAx);
end