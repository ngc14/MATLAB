directory = dir('S:\Lab\Gilligan\All Data');
allSessions = {};
conds = ["Extra Small Sphere","Large Sphere","Photocell"];
tic
for f = 3:size(directory,1)
    resultFolder = [directory(f).folder, '\',directory(f).name,'\EMG\Results_New\'];
    if(exist(resultFolder, 'dir'))
        results = dir([resultFolder,'*.mat']);
        sortedEMGFile = arrayfun(@(r) any(strcmp(string({whos(matfile([resultFolder,results(r).name])).name}),'sortedEMGData')),1:length(results));
        sessionT = cell(1,length(results));
        parfor r = 1:length(sortedEMGFile)
            matF = matfile([results(r).folder, '\', results(r).name]);
            if(sortedEMGFile(r))
                condCheck = matF.sortedEMGData;
                condCheck = condCheck.Conditions;
            else
                condCheck = matF.Conditions;
            end
            condCheck = reshape(condCheck,[],1);
            muscle = regexp(results(r).name,'_','split');
            muscle = upper(erase(string(muscle{end}(1:end-4))," "));
            if(~any(cellfun(@(c) contains(c,"Stim"),condCheck)) && all(ismember(conds,condCheck)))
                if(r==1),disp([results(r).folder, '\', results(r).name]),end
                if(sortedEMGFile(r))
                    matF = matF.sortedEMGData;
                    sampleRate = matF.SampleRate;
                    trials = matF.EMGData;
                    segTimes = matF.SegTimes;
                    arduino = matF.ArduinoData;
                    arduino = arduino(:,1);
                    validConds = ismember(arduino,condCheck(ismember(condCheck,conds)));
                    arduino = arduino(validConds);
                    segTimes = segTimes(validConds);
                    trials = trials(validConds);
                else
                    sampleRate = matF.Fs;
                    trials = matF.rawTrials;
                    validConds = ismember(condCheck,conds);
                    trials = [trials{validConds}];
                    segTimes = matF.segs;
                    segTimes = segTimes(validConds);
                    arduino = cell2mat(cellfun(@(c,r) repmat(string(c),size(r,2),1),condCheck(validConds),segTimes,'UniformOutput',false));
                    segTimes = [segTimes{:}];
                end
                nanTrials = cellfun(@(s) all(isnan(s)),segTimes);
                badTrials = getBadTrials(trials(~nanTrials), segTimes(~nanTrials), sampleRate);
                sessionCount = zeros(3,length(conds));
                for c = 1:length(conds)
                    sessionCount(1,c) = sum(strcmp(arduino(~nanTrials),conds(c)) & ~badTrials');
                    sessionCount(2,c) = sum(strcmp(arduino(~nanTrials),conds(c)));
                    sessionCount(3,c) = sum(strcmp(arduino,conds(c)));
                end
                sessionT{r} = {muscle,sessionCount};
            end
        end
        if(~all(cellfun(@isempty,sessionT)))
            sessionT = sessionT(~cellfun(@isempty,sessionT));
            allSessions{end+1} = dictionary(cell2mat(cellfun(@(s) s{1},sessionT,'UniformOutput',false)), cellfun(@(s) s{end},sessionT,'UniformOutput',false));
        end
        toc
    end
end
%%
finalT = configureDictionary("string","cell");
for f = 1:length(allSessions)
    addKeys = allSessions{f}.keys;
    for k = 1:length(addKeys)
        currKey = addKeys{k};
        if(~finalT.isKey(currKey) && ~any(cell2mat(cellfun(@(f) max([strfind(f,currKey),strfind(currKey,f(1:end-1))]),finalT.keys,'UniformOutput',false))))
            finalT = insert(finalT, string(currKey),{allSessions{f}{currKey}});
        else
            fInd = finalT.keys;
            fStr = cellfun(@(s)  max([strfind(s,currKey),strfind(currKey,s(1:end-1))]), finalT.keys, 'UniformOutput', false);
            fStr(cellfun(@isempty, fStr)) = {Inf};
            fInd = fInd(cell2mat(fStr)==min(cell2mat(fStr)));
            finalT(fInd) = {finalT{fInd} + allSessions{f}{currKey}};
        end
    end
end
finalT("TRICEPS") = {finalT{"TRICEPS"}+finalT{"TRIEPS"}};
finalT = finalT.remove("TRIEPS");
finalT("DIGITEXTENSOR") = {finalT{"DIGITSEXTENSOR"}+finalT{"DIGITEXTENSOR"}};
finalT("DIGITFLEXOR") = {finalT{"DIGITSFLEXOR"}+finalT{"DIGITFLEXOR"}};
finalT = finalT.remove("DIGITSFLEXOR");
finalT = finalT.remove("DIGITSEXTENSOR");
finalT("DIGITFLEXOR") = {finalT{"DIGITFLEXOR"}+finalT{"DISTALFLEXOR"}};
finalT = finalT.remove("DISTALFLEXOR");
finalT("D4D5EXTENSOR") = {finalT{"D4D5EXTENSOR"}+finalT{"DIGITEXTENSOR"}};
finalT = finalT.remove("DIGITEXTENSOR");
finalT("DELTOID") = {finalT{"DELTOID"}+finalT{"DELTPOD"}};
finalT = finalT.remove("DELTPOD");
finalT("WRISTFLEXOR") = {finalT{"WRISTFLEXOR"}+finalT{"WIRSTFLEXOR"}};
finalT = finalT.remove("WIRSTFLEXOR");
mean(reshape(cell2mat(cellfun(@(a) cell2mat(cellfun(@(v) v(end,:)-v(1,:), a.values, 'UniformOutput',false)),...
    allSessions, 'UniformOutput',false)'),[],1))
save(['S:Lab\Gilligan\Sorted EMG Data\Results\count'],'allSessions' , 'finalT');
%AVG from both monkeys: 3.15 trials / 5.87% of trials