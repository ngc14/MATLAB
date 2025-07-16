function gsbinder

% conbine several goodstim file
% 080526: Made by Hisashi

clear all;
datadriver = 'K:\';     % Data disk name
datafolder = 'expt1\';   % Data folder name on data disk, results will be saved in 'expresult'
expname = '101231JerL2G_Atn\'; % Exp folder name (in both data folder and result folder)
runname = 'run00\';      % Run foler name (in both data folder and result folder)

gsfilename1='LF101231t0116PM_Sorted_reward_gs_wCt.txt';
gsfilename2='_shiftlog.2_ff_gs.txt';
gsfilename3='_shiftlog.3_within_gs.txt';
outputgsfilename='goodstim00_wCt_150302.txt';

NStim=9; %size(goodstim, 2)
% NBlock=94; %size(goodstim, 1)

%%
fidgoodstim1=fopen(strcat(gsfilename1), 'r'); %Hisashi 071023
goodstim1=fscanf(fidgoodstim1, '%d');
fclose(fidgoodstim1);
NBlock=size(goodstim1,1)/NStim;
if NBlock~=floor(NBlock)
    NBlock
    fprintf('\r*** Error: Check NStim first! ***\n');
end
goodstim1=reshape(goodstim1, NStim, NBlock);
goodstim1=goodstim1';

% if gsfilename2 ~=''
    fidgoodstim2=fopen(strcat(gsfilename2), 'r'); %Hisashi 071023
    goodstim2=fscanf(fidgoodstim2, '%d');
    fclose(fidgoodstim2);
    goodstim2=reshape(goodstim2, NStim, NBlock);
    goodstim2=goodstim2';
% end
% if gsfilename3 ~=''
    fidgoodstim3=fopen(strcat(gsfilename3), 'r'); %Hisashi 071023
    goodstim3=fscanf(fidgoodstim3, '%d');
    fclose(fidgoodstim3);
    goodstim3=reshape(goodstim3, NStim, NBlock);
    goodstim3=goodstim3';
% end

goodstim=goodstim1.*goodstim2.*goodstim3;

fidoutputgsfilename=fopen(strcat(outputgsfilename), 'w');



for i=1:NBlock
    for j=1:NStim
        fprintf(fidoutputgsfilename, '%d\t', goodstim(i,j));
    end
    fprintf(fidoutputgsfilename, '\r\n');
end
fclose(fidoutputgsfilename);

fprintf('\r*** The gs files were conbinded! ***\n');

return 