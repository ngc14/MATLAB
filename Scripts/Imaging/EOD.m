date = '12_18_2020';
run = 0;
stimConds = containers.Map([0 2 4 5],{'Rest',  'ExtraSmallSphere','LargeSphere', 'Photocell'});

imageAlignment_ngc14(date,run,cell2mat(stimConds.keys));
pause(60);
BlockView_ngc14_V2_multirun(date,run,stimConds);
pause(60);
ttest_nsc15_V3_multirun(date,run,stimConds);

pause(60);
date = '12_18_2020';
run = 1;
stimConds = containers.Map([0 2 4 5],{'Rest',  'ExtraSmallSphere','LargeSphere', 'Photocell'});

imageAlignment_ngc14(date,run,cell2mat(stimConds.keys));
pause(60);
BlockView_ngc14_V2_multirun(date,run,stimConds);
pause(60);
ttest_nsc15_V3_multirun(date,run,stimConds);