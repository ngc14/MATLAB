%rawActivity =  load("C:\Users\ngc14\Desktop\RawEMG.mat");
phaseNames = ["Baseline","Go", "Reach", "Grasp","Withdraw"];
%  rAc = num2cell(cellfun(@(d) cellfun(@(c) vertcat(c), d, 'UniformOutput',false), rawActivity,'UniformOutput',false),2);
% rAc = cellfun(@(cm) num2cell((horzcat(cm{:})),1),rAc,'UniformOutput',false);
%rAc = cellfun(@(s) cellfun(@(n) cell2mat(n'),s,'UniformOutput',false),rawActivity,'UniformOutput',false);
rAc = rawActivity;
tEMG = [];
varNames = ["Trial" "LimbSeg" "Condition" "Go" "Reach" "Grasp" "Withdraw"];
repNames = ["Arm","Hand"];
conditions = [ "Extra Small Sphere","Large Sphere","Photocell"];
%%
for c = 1:length(conditions)
    trial=1;
    for a = 1:size(rAc,1)
        AUCrep= rAc(a,:);
        AUCVals = cell2mat(cellfun(@(cn) mean(cn{c},2,'omitnan'), AUCrep, 'UniformOutput',false));
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
    end
end
rAc = cellfun(@(s) s(:,2:end),num2cell(rAc,2),'UniformOutput',false);
rAc = cellfun(@(r) cellfun(@(b) b, num2cell(horzcat(r{:}),2), 'UniformOutput',false), rAc, 'UniformOutput',false);
%%
%condXrepXphaseEMG = cellfun(@(b) b,  num2cell(rAc,1), 'UniformOutput' ,false);
%condXrepXphaseEMG = num2cell(cellfun(@(b) b, (vertcat(rAc{:})'), 'UniformOutput' ,false),2);
condXrepXphaseEMG = num2cell(cellfun(@(b) b, (horzcat(rAc{:})'), 'UniformOutput' ,false)',2);
plotGroupedBars(condXrepXphaseEMG,"S:\Lab\ngc14\Working\Save\Phase_Boxplot\"+"EMG_AUC",true);