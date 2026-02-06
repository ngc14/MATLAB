%rawActivity =  load("C:\Users\ngc14\Desktop\RawEMG.mat");
phaseNames = ["Go", "Reach", "Grasp","Withdraw"];
%  rAc = num2cell(cellfun(@(d) cellfun(@(c) vertcat(c), d, 'UniformOutput',false), rawActivity,'UniformOutput',false),2);
% rAc = cellfun(@(cm) num2cell((horzcat(cm{:})),1),rAc,'UniformOutput',false);
%rAc = cellfun(@(s) cellfun(@(n) cell2mat(n'),s,'UniformOutput',false),rawActivity,'UniformOutput',false);
rAc = rawActivity;
tEMG = [];
varNames = ["Trial" "LimbSeg" "Condition" "Go" "Reach" "Grasp" "Withdraw"];
repNames = ["Arm","Hand"];
conditions = [ "Extra Small Sphere","Large Sphere","Photocell"];
%%
vs = {};
for c = 1:length(conditions)
    trial=1;
    for a = 1:size(rAc,1)
        AUCrep=cellfun(@(r) r{c} ,rAc,'UniformOutput',false);
        AUCVals = cell2mat(cellfun(@(cn) mean(cn,2,'omitnan'), AUCrep(a,:), 'UniformOutput',false));
        condTable = table();
        condTable.Trial = [trial:trial + length(AUCVals)-1]';
        condTable.LimbSeg = repmat(repNames(a),length(AUCVals),1);
        condTable.Condition = categorical(cellstr(repmat(conditions{c}(1),length(AUCVals),1)));
        condTable.Go = AUCVals(:,strcmp(phaseNames,"Go"));
        condTable.Reach = AUCVals(:,strcmp(phaseNames,"Reach"));
        condTable.Grasp = AUCVals(:,strcmp(phaseNames,"Grasp"));
        condTable.Withdraw = AUCVals(:,strcmp(phaseNames,"Withdraw"));
        condTable.Properties.VariableNames = varNames;
        tEMG = [tEMG;condTable];
        trial = trial + length(AUCVals);
        vs{c,a} = AUCVals;
    end
end
tEMG = unstack(tEMG,["Go","Reach","Grasp","Withdraw"],"Condition");
varNames = tEMG.Properties.VariableNames;
%%
rAc = cellfun(@(s) s(:,2:end),num2cell(rAc,2),'UniformOutput',false);
rAc = cellfun(@(r) cellfun(@(b) b, num2cell(horzcat(r{:}),2), 'UniformOutput',false), rAc, 'UniformOutput',false);
xVals = cellfun(@(v,a) repmat([1:size(v,2)]+a,size(v,1),1),vs,repmat({[-.2 -.2 -.2 -.2],[.2 .2 .2 .2]},size(vs,1),1),'UniformOutput',false);
figure(); hold on;
cl = [0 .5 0; .5 0 .5];
limbLab = cellfun(@(p) [repmat("Arm",size(p{1},1),1);repmat("Hand",size(p{2},1),1)],num2cell(vs,2),'UniformOutput',false);
s=swarmchart(reshape(cell2mat(cellfun(@(a,b) b+vertcat(a{:}),num2cell(xVals,2),num2cell(0:5:5*size(xVals,1)-1)','UniformOutput',false)),1,[]),...
    reshape(cell2mat(cellfun(@(a) vertcat(a{:}),num2cell(vs,2),'UniformOutput',false)),1,[]),[],repmat(cl(strcmp(cell2mat(limbLab),"Hand")+1,:),4,1),'filled');
arrayfun(@(p) set(p,'XJitter','rand'),s);
arrayfun(@(p) set(p,'XJitterWidth',.1),s);
ylim([0 .6]);
xa = gca();
set(xa,'XTick',2+[1:max(xa.XLim)/length(conditions):max(xa.XLim)]);
set(xa,'XTickLabel',cellfun(@(s) s(1),cellstr(conditions),'UniformOutput',false));
[meanVals,gNames] = cellfun(@(vs,lbs) groupsummary(vs,{lbs},"median"),...
    cellfun(@(a) vertcat(a{:}),num2cell(vs,2),'UniformOutput',false),limbLab,'UniformOutput',false);
plotVals = reshape(cell2mat(meanVals),[],1);
scatter(unique(cell2mat(cellfun(@(a,b) b+vertcat(a{:}),num2cell(xVals,2),num2cell([0:5:5*size(xVals,1)-1])','UniformOutput',false)),'stable'),...
    plotVals,400,'black',"_",'LineWidth',3);
%bx=violin(condPhases,'x',length(phaseNames)*(condLabels-1+repmat([-.3 -.1 .1 .3],1,length(conditions))),'facecolor',cl(phaseLabels,:));
saveFigures(gcf,"S:\Lab\ngc14\Working\PSTHS\","EMG",[]);

%%
%condXrepXphaseEMG = cellfun(@(b) b,  num2cell(rAc,1), 'UniformOutput' ,false);
%condXrepXphaseEMG = num2cell(cellfun(@(b) b, (vertcat(rAc{:})'), 'UniformOutput' ,false),2);
condXrepXphaseEMG = num2cell(cellfun(@(b) b, (horzcat(rAc{:})'), 'UniformOutput' ,false)',2);
plotGroupedBars(condXrepXphaseEMG,"S:\Lab\ngc14\Working\Save\Phase_Boxplot\"+"EMG_AUC",true);