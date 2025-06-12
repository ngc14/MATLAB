function [] = ttest_nsc15_V3_multirun(varargin)

animal = 'Skipper';
hemi = 'Right';
%dates = {'12_03_2018', '12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
%runs = {[1], [0], [0], [0], [1], [0], [0], [1], [0]};
dates = {'10_21_2020', '10_30_2020', '11_07_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
runs = {[0], [0], [0], [0], [0], [0], [0], [0],[0],[0],[0]};

frames = {[51:56]};

style = [2];
stim1 = [2 4 5];
stim2 = [0];
stimNames = containers.Map([0 2 4 5],{'Rest','ExtraSmallSphere','LargeSphere', 'Photocell'});

if ~isempty(varargin)
    % dates
    if ~isempty(varargin{1})
        dates = {varargin{1}};
    end
    
    % runs
    if ~isempty(varargin{2})
        runs = {varargin{2}};
    end
    
    % stims
    if ~isempty(varargin{3})
        stimNames = varargin{3};
        stim1 = cell2mat(stimNames.keys);
        %stim1 = stim1(stim1~=0);
    end
end


LPk = [15];
HPk = [250];
clipMethod = 2;
clip = [.3];             % does not affect ttest

number_runs = length(frames)*length(LPk)*length(HPk)*length(clip)*length(stim1)*length(style)*length(stim2)*length(dates);
for d = 1:length(dates)
    number_runs = number_runs*length(runs{d});
end;

disp(['There will be ',num2str(number_runs),' total runs of ttest_nsc15_V3.m']);
disp(' ');
startTime = clock;

k = 1;
for d = 1:length(dates)
    if d==1
        stim1 = [2];
    elseif d ==2 || d==5
        stim1 = [4 5];
    else
        stim1 = [2, 4,5];
    end
    for r = 1:length(runs{d})
        for l = 1:length(LPk)
            for h = 1:length(HPk)
                for c = 1:length(clip)
                    for s = 1:length(stim1)
                        for f = 1:length(frames)
                            for s2 = 1:length(style)
                                for s3 = 1:length(stim2)
                                    if(stim1(s)~= stim2(s3))
                                        run = ['run',num2str(runs{d}(r),'%1.2d')];
                                        ttest_ngc14_V1(animal,hemi,dates{d},...
                                            run,stim1(s),style(s2),...
                                            LPk(l),HPk(h),clipMethod,...
                                            clip(c),frames{f},stim2(s3),...
                                            stimNames)
                                        disp([num2str(round(k/number_runs*10000)/100),'% done with multirun']);
                                        k = k+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
end

totalTime = round(etime(clock,startTime)/60*10)/10;
timePerRun = round(totalTime/number_runs*10)/10;
disp(['BlockView multirun took ',num2str(totalTime),' mins (',num2str(timePerRun),' min per run)']);
% end