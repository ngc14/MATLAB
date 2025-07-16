% function [] = ttest_nsc15_V3_multirun()
clearvars;
animal = 'Gilligan';
hemi = 'right';
dates = {'12_12_2018'};
frames = {[12:14]};
run = 'run01';

style = [4];
stim1 = [3 4];
stim2 = [6 6];
stimName = containers.Map([0 2 3 4 5],{'Rest', 'ExtraSmallSphere', 'SmallSphere', 'LargeSphere', 'Photocell'});
stimName = containers.Map([0 2 3 4 6 7],{'Blank', 'D1_75','pD2', 'P1','D1_50','D1_25'});

LPk = [10];
HPk = [350,550];
clipMethod = 2;
clip = [.75];             % does not affect ttest

number_runs = length(frames)*length(LPk)*length(HPk)*length(clip)*length(stim1)*length(style)*length(stim2)*length(dates);

disp(['There will be ',num2str(number_runs),' total runs of ttest_nsc15_V3.m']);
disp(' ');
startTime = clock;
k=0;
for l = 1:length(LPk)
    for h = 1:length(HPk)
        for c = 1:length(clip)
            for s = 1:length(stim1)
                for f = 1:length(frames)
                    for s2 = 1:length(style)
                        for s3 = 1:length(stim2)
                            if(stim1(s)~=stim2(s3))
                                ttest_ngc14_V3(animal,hemi,dates{1},run,stim1(s),...
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

totalTime = round(etime(clock,startTime)/60*10)/10;
timePerRun = round(totalTime/number_runs*10)/10;
disp(['BlockView multirun took ',num2str(totalTime),' mins (',num2str(timePerRun),' min per run)']);
% end