% function [] = ttest_nsc15_V3_multirun()
clearvars;

animal = 'Gilligan';
hemi = 'Right';
dates = {'12_03_2018', '12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
runs = {[1], [0], [0], [0], [1], [0], [0], [1], [0]};
frames = {[50:55]};

style = [2];
stim1 = [2, 4, 5];
stim2 = [2, 4, 5];
stimName = containers.Map([0 2 3 4 5],{'Rest', 'ExtraSmallSphere', 'SmallSphere', 'LargeSphere', 'Photocell'});

LPk = [15];
HPk = [550];
clipMethod = 2;
clip = [.75];             % does not affect ttest

number_runs = length(frames)*length(LPk)*length(HPk)*length(clip)*length(stim1)*length(style)*length(stim2)*length(dates);
for d = 1:length(dates)
    number_runs = number_runs*length(runs{d});
end;

disp(['There will be ',num2str(number_runs),' total runs of ttest_nsc15_V3.m']);
disp(' ');
startTime = clock;

k = 1;
for d = 1:length(dates)
    for r = 1:length(runs{d})
        if runs{d}(r)<10
            run = ['run0',num2str(runs{d}(r))];
        else
            run = ['run',num2str(runs{d}(r))];
        end
        superpixels = dir(['\\pitt\sni\Gharbawie\Lab\' animal,'\All Data\',animal,'_',dates{d},...
            '\Imaging\',run,'\Results\Superpixel\LP5HP550C0']);
        superpixels = superpixels(contains({superpixels.name}, '.mat'));
        
        for l = 1:length(LPk)
            for h = 1:length(HPk)
                for c = 1:length(clip)
                    for s = 1:length(stim1)
%                         superpixels_cond = superpixels(contains({superpixels.name}, stimName(stim1(s))));
%                         frames = {};
%                         for sc = 1:length(superpixels_cond)
%                             load([superpixels_cond(sc).folder, '\', superpixels_cond(sc).name]);
%                             midpoint = round(mean([domain_stim{:}]));
%                             frames{sc} = [midpoint-1:midpoint+1];
%                             if(midpoint==70)
%                                 frames{sc} = [69:70];
%                             end
%                         end
                        for f = 1:length(frames)
                            for s2 = 1:length(style)
                                for s3 = 1:length(stim2)
                                    if(stim1(s)~= stim2(s3))
                                        ttest_ngc14_V3(animal,hemi,dates{d},run,stim1(s),...
                                            stimName(stim1(s)),style(s2),LPk(l),HPk(h),clipMethod,clip(c),frames{f},stim2(s3))
                                        
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