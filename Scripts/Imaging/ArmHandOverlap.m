monkey = "Gilligan";
pval = '0.0001';
tTestRange = {32:70};
activationConditions = ["[ExtraSmallSphere]", "[LargeSphere]","[Photocell]","[Rest]"];
if(strcmp(monkey,"Skipper"))
    mm = MotorMapping(56);
else
    mm = MotorMapping(42);
end
[datesTable, masks, ~] = getMonkeyInfo("S:\Lab\",monkey,["M1", "PMd"],false);
siteMask = getVoronoiMask(datesTable,mm,masks{1},["Arm","Hand"]);
maskNames = ["FOV_Forelimb", "M1","PMd"];

allMMs = {};
allMMs{end+1} = siteMask;
allMMs{end+1} = masks{2} & siteMask;
allMMs{end+1} = masks{3} & siteMask;
totals = [];
%%
for c= 1:length(activationConditions)   
    for f = 1:length(tTestRange)
        ttestName = strjoin(["S:\Lab\ngc14\Working",monkey,activationConditions{c},...
        "tTests\NaN\HP250\tTest_NoMask"+num2str(tTestRange{f})+"_p-"+pval+".bmp"],'\');
        % for ff = 1:length(tTestRange)
        %     ttestName = strjoin([frameDir,strjoin(["tTest_",...
        %           num2str(tTestRange(ff)),"_p-",pVal,".png"],'')],'');
        % end
        imR = mat2gray(imread(ttestName));
        activation = (imR==1 & siteMask==1);
        for m = 1:length(allMMs)
            currMask = allMMs{m};
            currActivation = activation & currMask;
            % totals = [armSA,handSA]./([sum(sum(armArea & ~handArea & currMask)),...
            % sum(sum(handArea & ~armArea & currMask))]+sum(sum(armHandArea & currMask))/2);
            totals(m,c,f) =  sum(currActivation(:)) / sum(sum(allMMs{m}));

        end
    end
end
res = totals.*100;
forelimbTotal = totals.*(18^2./1e6);

if(0)
    b = plotBarStackGroups([forelimbTotal(:,1) + forelimbTotal(:,3), forelimbTotal(:,2)], {'M1', 'PMd'});
    res(:,c) = forelimbTotal(:,2) ./ (sum(forelimbTotal,2));%forelimbTotal(:,1) + forelimbTotal(:,3));
    title(activationConditions{c});
    ylim([0 5]);
    barNames = cellfun(@(cc) regexp(cc,'(?=)]','split'), activationConditions,'UniformOutput',false);
    barNames = cellfun(@(cc) string(cellfun(@(c2) c2(2:end), cc(1:2), 'UniformOutput',false)),...
        barNames, 'UniformOutput',false);
    barNames = cellfun(@(s) strjoin(s,'>'),[barNames, ...
        cellfun(@fliplr, barNames, 'UniformOutput', false)],'UniformOutput', false);
    saveas(gcf, ['S:\Lab\ngc14\Figures\Imaging\', activationConditions{c}(1:end-4), '_Areas_', monkey, '.eps'], 'epsc');
end
figure('Units','normalized','Position',[0 0 1 1]);
b = plotBarStackGroups(res, maskNames);
yticks(0:10:100);
ylim([0 100]);
%arm/forelimb,'
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 1 0];
b(3).FaceColor = [0 0 1];
legend(activationConditions)
saveas(gcf,strjoin(['S:\Lab\ngc14\Working\',monkey,'\Percentage_',num2str(tTestRange(1)),'_',num2str(tTestRange(end)),'_',pval,'.eps'],''),'epsc');
saveas(gcf,strjoin(['S:\Lab\ngc14\Working\',monkey,'\Percentage_',num2str(tTestRange(1)),'_',num2str(tTestRange(end)),'_',pval,'.png'],''),'png');