function figHandle = modulatedUnitsPerRep2(repLabels,unitValues,valueLabel,condNames)
xBuff = 0.5;
axColors = [.75 .45 0; 0 .4 .4;];

repColorLabels = fieldnames(MotorMapping.repColors);
repNames = repColorLabels(ismember(repColorLabels,unique(horzcat(repLabels(:)))));
jointInds = arrayfun(@(jN) arrayfun(@(mJ) logical(cellfun(@(ss) strcmp(ss,jN),...
    mJ,'UniformOutput',true)),repLabels,'UniformOutput', false),repNames,'UniformOutput', false);

spAx = [];spInd = 1;
figHandle = figure('Units', 'normalized', 'Position', [0 0 1 1]);hold on;
for c = 1:length(condNames)
    for p = 1:size(unitValues,1)
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
        else
            currUnitTypes = vertcat(currUnitTypes(:));
            xlim([xBuff,length(repNames)+xBuff]);
            xticks(1:length(repNames));
            xticklabels(repNames);
            jointVals = cellfun(@(jI) currUnitTypes(cell2mat(jI')),jointInds, 'UniformOutput', false);
            b = bar(cellfun(@nanmean, jointVals),'FaceColor', 'flat');
            errorbar(cellfun(@nanmean,jointVals),cellfun(@nanstd,jointVals)./...
                cellfun(@(jv) sqrt(sum(~isnan(jv))),jointVals),'k', 'linestyle', 'none');
        end
        b.CData = cell2mat(arrayfun(@(jN) MotorMapping.repColors.(char(jN)), repNames, 'UniformOutput', false));
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