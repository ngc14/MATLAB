sessionDate = '04_26_2018';
testNum = '_1';
directory = ['\\130.49.229.252\gharbawie\Lab\Gilligan\All Data\Gilligan_',sessionDate, '\'];
load([directory,'EMG\Results\EMGinterp_',sessionDate,'.mat']);
load([directory,'Phasespace\Results\wristInterp_',sessionDate,testNum,'.mat']);
fs = 480;
conds = fieldnames(wristAnglesInterpl);
mapConds = [3 1 2];
fig1 = figure();
fig2 = figure();
fig3 = figure();
for c = 1:3 %length(conds)
    overlay = subplot(3,1,c,axes(fig1));
    corr = subplot(3,1,c,axes(fig2));
    combinedPS = [];
    combinedEMG = cell(3,1);
    interpolatedEMG = {};
    currEMG = EMGinterp.(conds{mapConds(c)});
    currPS = wristAnglesInterpl.(conds{mapConds(c)});
    hold(overlay, 'on');
    for s = 1:size(currPS,2)
        % make single NaNs an array of NaNs the same length of the given
        % segment time in order to make cell2mat functional
        tempData=currPS(:,s);
        lengths=cellfun(@length,tempData);
        tempData(find(lengths==1))={repmat(NaN,1,max(lengths))};
        % concatenate all segments into trial-by-totalTime matrix
        combinedPS = [combinedPS cell2mat(tempData)];
        % plot dashed line for segment times
        segs(s) = plot(overlay, [size(combinedPS,2)/fs, size(combinedPS,2)/fs], [-50 50],'--k');
    end
    avgPS = nanmean(combinedPS);
    wristLegend = plot(overlay,(1:size(avgPS,2))/fs, avgPS, '-m');
    yyaxis(overlay, 'left');
    ylabel(overlay, 'Flexion/Extension');
    colors = {'-r', '-g', '-b'};
    for m = 1:size(currEMG,1)
        ax3 = subplot(3,3,3*c+m-3,axes(fig3));
        for s = 1:size(currPS,2)
            if(s==1)
                % 1 second before go cue
                currSeg = currEMG{m,s}(:,end-2000+1:end);
                combinedEMG{m} = [combinedEMG{m} currSeg];
            elseif(s==size(currPS,2))
                % 1 second after replace home
                currSeg = currEMG{m,s}(:, 1:2000);
                combinedEMG{m} = [combinedEMG{m} currSeg];
            else
                combinedEMG{m} = [combinedEMG{m} currEMG{m,s}];
            end
            
        end
        combinedEMG{m} = combinedEMG{m}/max(max(combinedEMG{m}));
        % interpolate EMG time points down to kinematic time points
        ratio = size(combinedPS,2)/(size(combinedEMG{m},2)-1);
        spacing = 1:(1/ratio):size(combinedEMG{m},2)-1;
        % interpolate each trial
        for t = 1:size(combinedEMG{m},1)
            interpolatedEMG{m}(t,:) = interp1(1:size(combinedEMG{m},2),...
                movmean(combinedEMG{m}(t,:),20),spacing);
        end
        %plot(interpolatedEMG{m}');
        interpolatedEMGAvg{m} = nanmean(interpolatedEMG{m});
        yyaxis(overlay, 'right');
        overlay.YAxis(1).Color = 'm';
        overlay.YAxis(2).Color = 'k';
        muscleLegend(m) = plot(overlay,(1:size(interpolatedEMGAvg{m},2))/fs,...
            interpolatedEMGAvg{m}, colors{m});
        
        hold(corr, 'on');
        scatter(corr, avgPS,interpolatedEMGAvg{m}, colors{m}(2));
        
        % get first half of kinematic and EMG data (start to middle of hold)
        startSeg = fs*(segs(2).XData(1)-.15);
        endSeg = fs*(segs(4).XData(1)+.1);
        phaseEMG = interpolatedEMG{m}(:,round(startSeg):round(endSeg));
        phasePS = combinedPS(:,round(startSeg):round(endSeg));
        % fit data to linear model for R^2 value
        model = fitlm(reshape(phasePS,numel(phasePS),1), reshape(phaseEMG,numel(phaseEMG),1));
        R = model.Rsquared.Adjusted;
        scatter(ax3,reshape(phasePS,numel(phasePS),1),...
            reshape(phaseEMG,numel(phaseEMG),1), colors{m}(2));
        
        % formatting
        title(ax3, [conds{mapConds(c)}, ' ', EMGinterp.Muscles{m}, ', R^2 = ' num2str(R)]);
        xlabel(ax3, 'Flexion/Extension');
        ylabel(ax3, 'Voltage');
    end
    % formatting
    legend(overlay, [wristLegend;muscleLegend'] , ['Wrist'; EMGinterp.Muscles]);
    title(overlay, conds{mapConds(c)});
    xlabel(overlay, 'Time (s)');
    ylabel(overlay, 'Voltage');
    legend(corr, EMGinterp.Muscles);
    title(corr,conds{mapConds(c)});
    xlabel(corr,'Flexion/Extension');
    ylabel(corr,'Voltage');
end