PSdate = '10_23_2018';
EMGdate = '10_11_2018';
locations = {'Close', 'Far'};
muscle = 'Deltoid';
joint = 'FE';
figure();
for loc = 1:1
    EMGdata = load(['\\univ.pitt.edu\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
        EMGdate, '\EMG\Results_New\', locations{loc}, '\',muscle, '_aligned']);
    
    PSdir = dir(['\\univ.pitt.edu\sni\gharbawie\Lab\Gilligan\All Data\Gilligan_', ...
        PSdate, '\Phasespace\Results\', locations{loc}]);
    PSdir = PSdir(~[PSdir.isdir]);
    jInd = find(contains({PSdir.name}, joint));
    PSdata = load([PSdir(jInd).folder, '\', PSdir(jInd).name]);
    EMGconds = EMGdata.conds;
    PSconds = PSdata.conditions;
    conds = {EMGconds, PSconds};
    [~,condInd] = min(cellfun(@length, conds));
    conds = conds{condInd};
    for c = 1:length(conds)
        hs = subplot(3,1,c);
        labels = [];
        xTicks = [];
        hold on;
        cond = conds{c};
        title(cond);
        EMGcond = find(strcmp(EMGconds, cond));
        PScond = find(strcmp(PSconds, cond));
        yyaxis(hs,'left');
        for align = 1:sum(~cellfun(@isempty,EMGdata.alignedDataX(EMGcond,:)))
            EMGY = EMGdata.alignedDataY{EMGcond,align};
            shadedErrorBar(EMGdata.alignedDataX{EMGcond,align}./2000, ...
                mean(EMGY),std(EMGY),'lineprops', {'b-','LineWidth',2});
            labels = [labels, [alignedDataX{EMGcond,align}(1)./2000 alignedDataX{EMGcond,align}(end)./2000]];
            for s =  1:length(EMGdata.segs{EMGcond,align})
                avgSegs = EMGdata.segs{EMGcond,align}./2000;
                if(avgSegs(s)>=EMGdata.alignedDataX{EMGcond,align}(1)/2000 && avgSegs(s)<=EMGdata.alignedDataX{EMGcond,align}(end)/2000)
                    plot([avgSegs(s), avgSegs(s)],[0 max(max(EMGY))],'b--');
                end
            end
        end
        
        yyaxis(hs, 'right');
        for align = 1:sum(~cellfun(@isempty,PSdata.alignedDataX(PScond,:)))
            PSY = PSdata.alignedDataY{PScond,align}';
            shadedErrorBar(PSdata.alignedDataX{PScond,align}./480,...
                mean(PSY),std(PSY),'lineprops', {'r-','LineWidth',2});
            
            for s =  1:length(PSdata.segs{PScond,align})
                avgSegs = PSdata.segs{PScond,align}./480;
                if(avgSegs(s)>=PSdata.alignedDataX{PScond,align}(1)/480 && avgSegs(s)<=PSdata.alignedDataX{PScond,align}(end)/480)
                    plot([avgSegs(s), avgSegs(s)],[0 max(max(PSY))],'r--');
                end
            end
        end
        set(gca, 'XTick', labels);
        set(gca, 'XTickLabel', arrayfun(@num2str, (labels.*1000), 'UniformOutput', false))
    end
end