load("X:/Lab/ngc14/Working/Both/Baseline_FR/phaseAnalysis.mat",...
    'taskBaseline','taskFR','siteChannels','avgPhase','phaseFR','siteRep','siteDateMap','siteLocation','siteSegs','normPSTH');
monkeys = ["Gilligan", "Skipper"];
alignLimits = {[-.5, .15], [-.20, .20],[-.15,.5]};
pVal = 0.05;
params = PhysRecording(["Extra Small Sphere","Large Sphere", "Photocell"],.01,.15,-6,5);
vMask = containers.Map(monkeys,arrayfun(@(s) im2double(im2gray(imread(strcat("S:\Lab\",s,"\Mapping\clean_reference_mask.png")))), monkeys,'UniformOutput',false));
chUnitMap = cell(height(siteDateMap),1);
chUnitMap(strcmp([siteDateMap.Monkey],"Gilligan")) = {[1:2:32,2:2:32]};
chUnitMap(strcmp([siteDateMap.Monkey],"Skipper")) = {[32:-1:1]};
condUnitMapping = cellfun(@(s) size(s,2), siteChannels{2})';
mappedChannels = cellfun(@(ch,l) ch(l(~isnan(l))), chUnitMap,siteChannels{2}, 'Uniformoutput', false);
mappedChannels = [mappedChannels{:}];
unitRep =  mapSites2Units(condUnitMapping,siteRep)';
unitLocation = mapSites2Units(condUnitMapping,siteLocation');
condDir = ["E","L","P"];
phaseLim = [0 3];
unitNames = ["Reach", "Grasp", "Both"];
phaseNames = ["Go","Reach", "Grasp"];
savePath = "S:\Lab\ngc14\Working\Both\Baseline_FR\Decoding\All\ActiveInds\";
%%
close all;
[~,taskUnits] = cellfun(@(pb,pc) cellfun(@(b,p)  ttestTrials(b,p,1,true,0.01),...
    pb,pc, 'UniformOutput', false),taskBaseline(1:3),taskFR(1:3),'UniformOutput',false);
taskUnits = cell2mat(cellfun(@cell2mat, taskUnits,'UniformOutput',false));
for c = 1:3
    tPhase = arrayfun(@(pn) cellfun(@(s) cellfun(@(t) {t{strcmp(phaseNames,pn)}},...
        s,'UniformOutput', false),phaseFR{c},'UniformOutput', false),phaseNames, 'UniformOutput',false);
    [~,rgInds] = cellfun(@(r,g) ttestTrials(r,g,1,true,pVal),tPhase{strcmp(phaseNames,"Reach")},...
        tPhase{strcmp(phaseNames,"Grasp")}, 'UniformOutput',false);
    rgInds = cell2mat(rgInds);
    AUCVals = vertcat(avgPhase{c}{:});
    goAUC = AUCVals(:,strcmp(phaseNames,"Go"));
    rAUC = AUCVals(:,strcmp(phaseNames,"Reach"));
    gAUC = AUCVals(:,strcmp(phaseNames,"Grasp"));
    rSI = rAUC./(rAUC + gAUC);
    gSI = gAUC./(rAUC + gAUC);

    rUnits = rgInds==1 & rAUC > gAUC;
    gUnits = rgInds==1 & rAUC < gAUC;
    bothUnits = ~rUnits & ~gUnits & taskUnits(:,c);
    plotPhaseVals = cellfun(@transpose,{goAUC,rAUC,gAUC,rSI,gSI},'UniformOutput',false);
    for m=3:length(monkeys)+1
        currRep = mapSites2Units(condUnitMapping,siteRep)';
        mMapping = cell2mat(arrayfun(@(m,n) ones(1,m)*n,condUnitMapping,...
            1:length(condUnitMapping),'UniformOutput',false));
        if(m==length(monkeys)+1)
            monk = monkeys;
            vmMask = [];
        else
            monk = monkeys(m);
            vmMask = double(vMask(monk));
        end
        mInds = contains(siteDateMap.Monkey,monk);
        muInds = mapSites2Units(condUnitMapping,mInds');
        mSiteLocation = unitLocation-min(unitLocation(muInds',:)).*ImagingParameters.px2mm;
        trialSegs =  num2cell(vertcat(siteSegs{c}{mInds}),1);
        mutInds = muInds & taskUnits(:,c)==1;
        mMapping(~muInds) = NaN;
        mMapping = mMapping-(nanmin(mMapping))+1;
        currRep(~mutInds) = "";
        %figure('Units', 'normalized', 'Position', [0 0 1 1]);
        currPSTH = num2cell(cellfun(@(m) mean(m,3,'omitnan'),vertcat(normPSTH{c}{:}),'UniformOutput',false),1);

        for u = 1:length(unitNames)
                                figHandle = figure('Units', 'normalized', 'Position', [0 0 1 1]); hold on;

            switch(unitNames(u))
                case "Both"
                    uInds = bothUnits;
                case "Reach"
                    uInds = rUnits;
                case "Grasp"
                    uInds = gUnits;
                otherwise
                    uInds = taskUnits(:,c);
            end
            mRep = currRep;
            mRep(~uInds) = "";
            for i = 1:2
                if(i==1)
                    channelLayer = double(mappedChannels'>16);
                    lrs = "Superficial";
                else
                    channelLayer = double(mappedChannels'<=16);
                    lrs = "Deep";
                end
                umInds = uInds==1&muInds==1& mapSites2Units(cellfun(@length,siteChannels{2}), siteActiveInd{c})==(i-1);%channelLayer==1;
                uPSTHS = cellfun(@(v) vertcat(v{:}), currPSTH,'UniformOutput', false);
                uPSTHS = cellfun(@(u) u.*(umInds./umInds), uPSTHS, 'UniformOutput',false);
                % countAUC(currRep,uInds',plotPhaseVals,mMapping,siteMasks(mInds'),vmMask,[],...
                %     savePath+condDir(c)+"\",strjoin(monk,'')+"_"+unitNames(u)+lrs,[arrayfun(@(s) strcat(s,"Phase"),phaseNames),"rSI","gSI"],...
                %     m<=length(monkeys),[1 10],{phaseLim,phaseLim,phaseLim,[0.3 .9],[0.3 .9]});
                plotJointPSTHS(params.bins,uPSTHS, num2cell(vertcat(siteSegs{c}{:}),1), ...
                    mRep,cell2mat(arrayfun(@(a) find(mMapping==a,1), unique(mMapping), 'UniformOutput', false)),...
                    [], alignLimits,[1 8], structfun(@(s) round(s./((5*i)-4),2), MotorMapping.repColors,'UniformOutput',false));
            end
            saveFigures(gcf,savePath+"PSTH\"+condLabels(c)+"\",strjoin(monk,'')+"_"+unitNames(u), []);

        end
    end
    close all;
end

function [phaseSubVals, sigs] = ttestTrials(dist1,dist2,taskPhase,paired,pVal)
if(paired)
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)>=max(1,floor(length(dist1)/2));
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}./...
        max(1,d1{taskPhase})),dist1,dist2,'UniformOutput',false)',1);
else
    sigs = nansum(cell2mat(vertcat(cellfun(@(d1,d2) ttest2(...
        d1{taskPhase},d2{taskPhase},'Alpha', pVal,'Dim',2),...
        dist1,dist2,'UniformOutput',false))),2)'>=floor(length(dist1)/2);
    smallestTrialCount = min(size(dist1{1}{taskPhase},2),size(dist2{1}{taskPhase},2));
    phaseSubVals = num2cell(cellfun(@(d1,d2) abs(d2{taskPhase}(:,1:smallestTrialCount)-...
        d1{taskPhase}(:,1:smallestTrialCount)),dist1,dist2,'UniformOutput',false)',1);
end
phaseSubVals = cellfun(@(m) median(cat(3,m{:}),3), phaseSubVals, 'UniformOutput', false);
sigs = double(sigs);

end