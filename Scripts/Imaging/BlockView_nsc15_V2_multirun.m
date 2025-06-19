% function [] = BlockView_nsc15_V2_multirun()
clearvars;

animal = 'Skipper';
hemi = 'Right';
dates = {'12_18_2020'};
runs = {[1]};
frames = {[50:55]};
excludeFrames = [];

style = [2];
stim1 = [0,2,4,5];
stimName = {'Rest', 'Extra Small Sphere', 'Large Sphere', 'Photocell'};
condStimMap = containers.Map([0 2 4 5],{'Rest', 'Extra Small Sphere', 'Large Sphere', 'Photocell'});

LPk = [5];
HPk = [250];
clipMethod = 6;
clip = [1];             % does not affect ttest
stim = containers.Map([0 2 4 5],{'Rest','ExtraSmallSphere','LargeSphere', 'Photocell'});


number_runs = length(frames)*length(LPk)*length(HPk)*length(clip)*length(stim1)*length(style);
for d = 1:length(dates)
    number_runs = number_runs*length(runs{d});
end;

disp(['There will be ',num2str(number_runs),' total runs of BlockView_nsc15_V2.m']);
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

                            BlockView_ngc14_V3(animal,hemi,dates{d},run,stim1(s),...
                                condStimMap(stim1(s)),style(s2),LPk(l),HPk(h),clipMethod,clip(c),frames{f},excludeFrames, [], [],stim)

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