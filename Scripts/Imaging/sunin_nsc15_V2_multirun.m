% function [] = sunin_nsc15_V2_multirun()
clearvars;

animal = 'Bilbo';
hemi = 'Left';
dates = {'07_24_2018'};
runs = {[1 2 3]};
frames = {[13:15]};

style = [4];
stim1 = [2 3];
stimName = {'Blue LED 1500ms pulse','Blue LED train'}; %{'150BiPulses_60microAmps'};

LPk = [10];
HPk = [300 350];
clipMethod = 2;
clip = [1.5];             % does not affect ttest

number_runs = length(frames)*length(LPk)*length(HPk)*length(clip)*length(stim1)*length(style);
for d = 1:length(dates)
    number_runs = number_runs*length(runs{d});
end;

disp(['There will be ',num2str(number_runs),' total runs of Sunin_nsc15_V2.m']);
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
        
        for l = 1:length(LPk)
            for h = 1:length(HPk)
                for c = 1:length(clip)
                    for f = 1:length(frames)
                        for s = 1:length(stim1)
                            for s2 = 1:length(style)

                            sunin_nsc15_V2(animal,hemi,dates{d},run,stim1(s),...
                                stimName{s},style(s2),LPk(l),HPk(h),clipMethod,clip(c),frames{f})

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