function BlockView_ngc14_V2_multirun()
monkey = 'Skipper';
if(strcmp(monkey,'Gilligan'))
    dates = {'12_05_2018', '12_07_2018','12_09_2018', '12_10_2018', '12_12_2018', '12_13_2018','12_14_2018','01_07_2019'};
    runs = {[0], [0], [0], [1], [0], [0], [1], [0]};
elseif(strcmp(monkey,'Skipper'))
    dates = {'10_30_2020','11_09_2020','11_23_2020', '11_24_2020', '11_25_2020', '11_27_2020', '11_30_2020','12_01_2020','12_02_2020'};
    dates = {'12_18_2020', '12_18_2020'};
    runs = {[0], [0], [0], [0], [0], [0], [0], [0], [0]};
    runs = {[0], [1]};
end
LPk = 5;
HPk = 250;
allConds = [4 5];
% delete(gcp('nocreate'));
% parpool('threads');
totalRuns = length(LPk)*length(HPk)*length(allConds)*length(dates);
disp(['There will be ',num2str(totalRuns),' total runs of BlockView_ngc14_V2.m']);
startTime = datetime('now');
for d =1:length(dates)
    stim1 = allConds;
    stimName = containers.Map([0 2 4 5],{'Rest','ExtraSmallSphere','LargeSphere', 'Photocell'});
    if(strcmp(dates{d},"12_09_2020"))
        stimName = containers.Map([0 4 5],{'Rest','LargeSphere', 'Photocell'});
    end
%     if(strcmp(monkey,"Skipper"))
%         if d ==1
%             stim1 = [0 2];
%         elseif d==2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                2
%             stim1 = [0 5];
%             stimName = containers.Map([0 4 5],{'Rest', 'LargeSphere', 'Photocell'});
%         end
%     end

    for l = 1:length(LPk)
        for h = 1:length(HPk)
            for s = 1:length(stim1)
                startrun = datetime('now');
                BlockView_ngc14_V1(monkey,dates{d},['run0',num2str(runs{d})],stim1(s), stimName,LPk(l),HPk(h))
                disp([num2str(time2num(datetime('now')-startrun,'minutes')), ' min to run']);
                disp(['BlockView multirun took ',num2str(time2num(datetime('now')-startrun,'minutes')),' mins']);
            end
        end
    end
end
totalTime = time2num((datetime('now')-startTime),'hours');
timePerRun = (totalTime*60)/totalRuns;
disp(['BlockView multirun took ',num2str(totalTime),' hrs (',num2str(timePerRun),' min per run)']);
end