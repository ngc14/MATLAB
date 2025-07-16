
clc
clear all
close all

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 0.02;                  % bin size in seconds
ElectrodeChannel = 10;       % ElectrodeChannel = 1 for single electrode recordings, ElectrodeChannel = 2 for second electrode OR second contact on V-probe.
SpikeChannel = 1;            % Spike channel 1 = 1a; Spike channel 2 = 1b; Spike channel 3 = 1c; Spike channel 4 = 1d;
Sigma=1;                     % Width (ms) of Gaussian for smoothing PSTH plots
bins_from_midpoint = 5;       % number of bins before (and after) midpoint for ANOVA analysis

time_before_event = -.01;    % seconds before event
time_after_event  = .01;     % seconds before event
stimElectrodes = [6:6:32];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Select a File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fname, pathname] = uigetfile('*.plx', 'Select a Plexon file','C:\Users\OMAR\Desktop\Ms Howell\MsHowell Single Unit\Right Hemisphere');
filename = strcat(pathname, fname);
results_folder = strcat(pathname, fname);
fid = fopen(filename, 'r');
if(fid == -1)
    disp('cannot open file', 'r');
    return
end
disp(strcat('file = ', filename));


%%

tic

%%
%%%%%%%%%%%%%%%%%%%%%%%%%Read File%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% event_time_stamps is a cell array for holding event_names in the first row and
% their respecitve time stamps in the second row
event_time_stamps {2,14} = [];
% read event_names
[n,event_names] = plx_chan_names(filename);

% store events names from channel 3 to 14 into first row of cell array
% event_time_stamps
for event_channel = 1:14
    event_time_stamps (1,event_channel) = cellstr (event_names (event_channel+2,:));
end

% read time stamps for event channels 3 to 14. Store time stamps for events
% in second row of cell array event_time_stamps
for event_channel = 1:14
    [n, event_ts] = plx_event_ts(filename, event_channel+2);
    event_time_stamps (2,event_channel) = {(event_ts)}';
end

% read spike time stamps for units a, b, c, d recorded on SpikeChannel (user
% specified above). Store unit letter in first row of the array
% spike_time_stamps. Store time stamps in second row of the array
% spike_time_stamps
unit_letter = {'a','b','c','d'};
spike_time_stamps {2,4} = [];

for unit_number = 1:4
    [n, spike_ts] = plx_ts(filename, 1, unit_number);      
    spike_time_stamps (1, unit_number) = unit_letter (unit_number);
    spike_time_stamps (2, unit_number) = {(spike_ts)};
end
%%


% Determine the highest number of events (e.g. 150 occasions of "grasp_sphere") 
maximum_number_event_time_stamps = max(cellfun(@length,event_time_stamps(2,:)));

% Preallocate a 3-D array named center_time_stamps  

center_time_stamps (length(spike_time_stamps{2,SpikeChannel}),maximum_number_event_time_stamps,event_channel) = zeros;

% Trigger spike_time_stamps to event_ts and store results in
% center_time_stamps
for  index_event = 1:length(event_time_stamps)
    
    %skip events that are empty
    if  isempty(event_time_stamps{2,index_event})
        index_event = index_event +1;
    else
        
        for index_event_timestamps = 1:length(event_time_stamps{2,index_event})
            event_ts = event_time_stamps{2,index_event}(index_event_timestamps);
            center_time_stamps (:, index_event_timestamps, index_event) = (spike_time_stamps{2,SpikeChannel}) - event_ts;
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%% PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psth_edges = time_before_event:bin:time_after_event;
[row, column, dimension] = size(center_time_stamps);
psth (length(psth_edges), column) = zeros;
psth_normalized (length(psth_edges), dimension) = zeros;
psth_normalized_filtered (length(psth_edges), dimension) = zeros;
psth_std (length(psth_edges), dimension) = zeros;
psth_std_filtered (length(psth_edges), dimension) = zeros;
center_time_stamps_temporary (length(center_time_stamps),1) = zeros; 
raster {column, dimension} = []; 

for i = 1: dimension
    psth ( :, :) = zeros;
    
    for j = 1:column
        if center_time_stamps(1,j,i) == 0;
            j = j-1;
            break;
        end;
        center_time_stamps_temporary = center_time_stamps(:,j,i);
        spike_times_stamps_triggered = center_time_stamps_temporary (center_time_stamps_temporary>time_before_event & center_time_stamps_temporary<time_after_event);
        
        raster{j,i} = {(spike_times_stamps_triggered)}; 
        
        psth (:,j) = histc (spike_times_stamps_triggered, (psth_edges));
        spike_times_stamps_triggered = [];
        
    end
    
    
    %%%%%%%%%%%%%%%%% Filter with Gaussian Kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kernel=(-3*Sigma:3*Sigma);
    BinSize=length(Kernel);
    Half_BW=(BinSize-1)/2;
    Kernel=(-BinSize/2:BinSize/2);
    Factor=(1.0003)/(Sigma*sqrt(2*pi));
    Kernel=Factor*(exp(-(0.5*((Kernel./Sigma).^2))));
    Kernel=Kernel';
    %%%%%%%%%%%%%%%%%% Create PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    psth_normalized (:,i) = (sum(psth,2))/j/0.02;
    psth_normalized_filtered  = convn(psth_normalized, Kernel, 'same');
    
    psth_std (:,i) = std (psth,0,2);
    psth_std_filtered = convn(psth_std, Kernel, 'same');
    
end

%%%%%%%%%%%%%%%%%%%%%%% Plot Smoothed PSTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_positive(length(psth_edges), dimension) = zeros; 
error_negative(length(psth_edges), dimension) = zeros; 

error_positive = (psth_normalized_filtered (:,:) + 3* psth_std_filtered (:,:));
error_negative = (psth_normalized_filtered (:,:) - 3* psth_std_filtered (:,:));
 
yMax = max(max(error_positive));
 
 
figure
for subplots = 1:index_event;
    subplot (index_event/7, index_event/2, subplots);
    x = time_before_event:bin:time_after_event;
    y = psth_normalized_filtered (1:end,subplots);
   
    plot (x,y,'r','LineWidth',2)
    hold on
    plot (x,error_positive(:,subplots),'r','LineWidth',0.5)
    hold on 
    plot (x,error_negative(:,subplots),'r','LineWidth',0.5)

    title (event_time_stamps{1,subplots});
    ylim ([0 yMax])
    xlim ('auto')
    axis square
    xlabel ('seconds');
    ylabel ('spikes/sec');
    
end

%  spaceplots(gcf,[0 0 0 0], [0 0]);

saveas (gcf, (strcat(pathname, fname, '_SpikeChannel_', num2str(SpikeChannel), '_PSTH', '.fig'))); 


%%%%%%%%%%%%%%%%%%%%%%% Normalize PSTH range 0 - 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psth_normalized_0_to_1(length(psth_edges), dimension) = zeros; 

min_psth = min(min(psth_normalized_filtered));
max_psth = max(max(psth_normalized_filtered));

psth_normalized_0_to_1 = (psth_normalized_filtered - min_psth)/ (max_psth - min_psth); 

psth_filename = strcat(pathname, fname, '_SpikeChannel_', num2str(SpikeChannel), '.xlsx')


%%%%%%%%%%%%%%%%%%%%%%% Write to Excel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% xlswrite(psth_filename, (psth_edges)',  'psth_filtered_normalized_0_to_1', 'A2')
% xlswrite(psth_filename, event_time_stamps(1,:),  'psth_filtered_normalized_0_to_1', 'B1')
% xlswrite(psth_filename, psth_normalized_0_to_1,  'psth_filtered_normalized_0_to_1', 'B2')
% 
% 
% xlswrite(psth_filename, (psth_edges)',  'psth_filtered_normalized', 'A2')
% xlswrite(psth_filename, event_time_stamps(1,:),  'psth_filtered_normalized', 'B1')
% xlswrite(psth_filename, psth_normalized_filtered,  'psth_filtered_normalized', 'B2')
% 
% 
% xlswrite(psth_filename, (psth_edges)',  'psth_3xStd_filtered_normalized', 'A2')
% xlswrite(psth_filename, event_time_stamps(1,:),  'psth_3xStd_filtered_normalized', 'B1')
% xlswrite(psth_filename, 3*(psth_std_filtered),  'psth_3xStd_filtered_normalized', 'B2')


%%%%%%%%%%%%%%%%%%%%%%% Plot Rasters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

for subplots = 1:index_event;
    subplot (index_event/7, index_event/2, subplots);
    x = time_before_event:0.001:time_after_event;
    
    trial = 1;
    
    for i_raster = 1:length(raster)
        
        if ~isempty(cell2mat(raster{i_raster,subplots})) && i_raster <= length(raster)
            spike_times_raster = cell2mat(raster{i_raster, subplots});
            for tick = 1:length(spike_times_raster);
                line ([spike_times_raster(tick) spike_times_raster(tick)], [trial-1 trial], 'Color',[0 0 0]);
            end
            
            trial = trial +1;
            i_raster = i_raster+1;
            
        elseif isempty(cell2mat(raster{i_raster,subplots})) && i_raster <= length(raster)
            i_raster = i_raster + 1;
            
        end
        
    end
    
    title (event_time_stamps{1,subplots});
    ylim ([0 trial])
    xlim ('auto')
    axis square
    xlabel ('seconds');
    ylabel ('Trials');
    
end
% spaceplots(gcf,[0 0 0 0], [0 0]);

saveas (gcf, (strcat(pathname, fname, '_SpikeChannel_', num2str(SpikeChannel), '_Rasters', '.fig'))); 

%%%%%%%%%%%%%%%%%%%%%%% Run 1-way ANOVA & Bonferroni %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_in =  round(length(psth_edges)/2 - bins_from_midpoint);
time_out =  round(length(psth_edges)/2 + bins_from_midpoint);

[p_value,table,stats] = anova1(psth_normalized_filtered (time_in : time_out, :));

[pos_hoc_compare, means, h] = multcompare(stats,'ctype','bonferroni');

ylabel ('Referenced to "event-names"');
xlabel ('spikes/sec');

toc

figure(h) = view(-90,90);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


