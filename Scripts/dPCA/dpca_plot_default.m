function dpca_plot_default(data, time, yspan, explVar, compNum, events, signif, condNames)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number

lines = [.8 0 0; .9 .67 0; 0 0 .8; 1 0 1; 0 1 0; 0.5 0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data, 'legend')
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

    % if there is one parameter
    elseif length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(1:numOfStimuli,:);
        hold on
        
        for f = 1:numOfStimuli
            p(f) = plot([NaN NaN], [NaN NaN], 'color', colors(f,:), 'LineWidth', 2);
            %text(1.2, f, [condNames(f)])
        end
        legend(p,condNames,'Location','northeast');
        % axis([0 3 -1 1.5+numOfStimuli])
        % set(gca, 'XTick', [])
        % set(gca, 'YTick', [])
        % set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision (decision can only have two
    % values)
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(1:numOfStimuli,:);
        hold on
        
        for f = 1:numOfStimuli
            p(f) = plot([NaN NaN],[NaN NaN], 'color', colors(f,:), 'LineWidth', 2);
            %text(1.2, f, ['Stimulus ' num2str(f)])
        end
        p(end+1) = plot([NaN NaN],[NaN NaN], 'k', 'LineWidth', 2);
        p(end+1) = plot([NaN NaN],[NaN NaN], 'k:', 'LineWidth', 2);
        %text(1.2, -2, 'Decision 1')
        %text(1.2, -3, 'Decision 2')
        legend([condNames,"Arm","Hand"])
        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
        
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) yspan])
hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

if ~isempty(events)
    plot([events; events], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    colororder(lines(1:numOfStimuli,:));
    plot(time, squeeze(data(1,:,:)), 'LineWidth', 2)    

elseif ndims(data) == 4 && size(data,3)==2
    % different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    colors = lines(1:numOfStimuli,:);

    for f=1:numOfStimuli 
        plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, :)), ':', 'color', colors(f,:), 'LineWidth', 2)
    end

else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(time, data, 'LineWidth', 2)    
end
