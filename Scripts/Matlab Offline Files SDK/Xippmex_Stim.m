% PATH TO XIPPMEX MIGHT NEED TO BE CHANGED IN LINE 74
% THE FILENAME WHERE THE DATA IS SAVED CAN BE CHANGED IN LINES 85 AND 94
% THIS SCRIPT REQUIRES THE FUNCTION importblock.m

%LAST EDITTED: 9/7/2017
monkey = 'Gilligan';
%% VARIABLES FOR USER MODIFICATION
trainFrequency = 2;                         %frequency of train period (Hz)

% Stimulation parameters 1
numberOfPulses         = [1,2];                 %number of pulses in train
frequency_Hz           = 200;               %frequency of pulse period (Hz)
aDur                   = 0.2;               %duration of anode phase, positive (ms) (leave at .1 ms for no phase)
cDur                   = 0.2;               %duration of cathode phase, negative (ms) (leave at .1 ms for no phase)
aAmp                   = 40;                 %height of anode phase's current (headstage steps - [0, 127]) (uA)
cAmp                   = 40;                %height of cathode phase's current (headstage steps - [0, 127]) (uA)
stimElectrodes = {8,16,24};
iterations = 150;
condIter = 1;
for p = 1:length(numberOfPulses)
    for f = 1:length(frequency_Hz)
        for a = 1:length(aDur)
            for c = 1:length(cDur)
                for aa = 1:length(aAmp)
                    for ca = 1:length(cAmp)
                        for e = 1:length(stimElectrodes)
                        stimParams(condIter,:) = [p,f,a,c,aa,ca,e];
                        condIter = condIter + 1;
                        end
                    end
                end
            end
        end
    end
end
%% OTHER VARIABLES AND PARAMETERS

%one indexed list of wanted electrodes - must be integers
%Port B starts at number 129
% 129 - nano2+stim electrode labeled 1
% 130 - nano2+stim electrode labeled 2
% 131 - nano2+stim electrode labeled 3
FEMap = [129:160];
FEMap = [1:32];
probeMap = [1:2:32, 2:2:32];
mappedElectrodes = FEMap(probeMap);

%for bipolar stimulation: 1 - cathodic first, 0 - anodic first
polarity               = 1;

%delay - delay for interleaving (ms)
electrodeDelay_ms      = 0;

%duration of time to fast settle after the full biphasic pulse (ms)
fastSettle = 1;

%% NIP Hardware Initialization

% Initialize xippmex
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

% Give the NIP some time to process any commands we have sent
pause(0.5)

%% GENERATE THE RANDOM ORDER THAT THE STIMULATION PATTERNS WILL BE DELIVERED

% Path to Xippmex
path(path, 'C:\Program Files (x86)\Ripple\Trellis\Tools\xippmex');
xippmex('stim', 'enable', '1')
cd(['C:\Users\ripple\Documents\Data\110-16_Gilligan_Macaque\',monkey,'_',datestr(now,'mm_dd_yyyy')]);
% Ask the user if they want to continue a new run or start a new one
choice = questdlg('Do you wish to continue a previous run or start a new run?', ...
    '', ...
    'Continue previous run','Start new run','Quit program','Quit program');

if strcmp(choice,'Start new run')
    % Create new file
    runNumber = inputdlg('Enter the new run number:',...
        'Run Number', [1 30]);
    sessionFilename=[monkey,'_',datestr(now,'mm_dd_yyyy'),'_stim_',runNumber{:},'.txt'];
    fileID = fopen(sessionFilename,'a');
    fprintf(fileID,'%s \t %s \t %s \t %s \t %s \t %s\r\n','Electrodes','Block',...
        'Number of Pulses', 'Frequency (Hz)', 'Pulse Width (Cathode,Anode ms)','Amplitude (Cathode,Anode mA)');
    block=1;
elseif strcmp(choice,'Continue previous run')
    % Open existing file
    runNumber = inputdlg('Enter the run number you wish to continue:',...
        'Run Number', [1 30]);
    sessionFilename=['XippmexStimulation_run',runNumber{:},'_',date,'.txt'];
    previousBlock=importblocks(sessionFilename);
    block=previousBlock(end,1)+1;
    fileID = fopen(sessionFilename,'a');
    fprintf(fileID,'%s\r\n','---Paused---');
end
%%
% The rest of this script is not executed if 'Quit program is chosen'
if strcmp(choice,'Start new run')||strcmp(choice,'Continue previous run')
    while(block<=iterations)
        conditionOrder = randperm(size(stimParams,1));
        for cond = 1:length(conditionOrder)
            t0 = clock;
            trialParams = stimParams(conditionOrder(cond),:);
            
            pulses = numberOfPulses(trialParams(1));
            freq = frequency_Hz(trialParams(2));
            anDur = aDur(trialParams(3));
            caDur = cDur(trialParams(4));
            anAmp = aAmp(trialParams(5));
            caAmp = cAmp(trialParams(6));
            elecs = mappedElectrodes(stimElectrodes{trialParams(7)});
            
            disp(['Electrode(s): ',num2str(stimElectrodes{trialParams(7)})])
            disp(['Block: ',num2str(block)])
            sizeE = size(elecs);
            % Calculate the train length from the number of pulses and the frequency
            trainLength_ms=1000*(pulses/freq);
            
            % Generate stimulation string
            elect_str = strcat('Elect=', sprintf('%d,', elecs));
            tl_str = strcat('TL=', sprintf('%.3f,', repmat(trainLength_ms,sizeE)));
            freq_str = strcat('Freq=', sprintf('%.0f,', repmat(freq,sizeE)));
            cDur_str = strcat('CathDur=', sprintf('%.3f,', repmat(caDur,sizeE)));
            aDur_str = strcat('AnodDur=', sprintf('%.3f,', repmat(anDur,sizeE)));
            cAmp_str = strcat('CathAmp=', sprintf('%d,', repmat(caAmp,sizeE)));
            aAmp_str = strcat('AnodAmp=', sprintf('%d,', repmat(anAmp,sizeE)));
            delay_str = strcat('TD=', sprintf('%.3f,', repmat(electrodeDelay_ms,sizeE)));
            pol_str = strcat('PL=', sprintf('%d,', repmat(polarity,sizeE)));
            fs_str = strcat('FS=', sprintf('%.4f,', repmat(fastSettle,sizeE)));
            
            stimString = strcat(elect_str, ';', tl_str, ';', freq_str, ';',...
                cDur_str, ';', aDur_str, ';', cAmp_str, ';', aAmp_str, ';', delay_str, ';',...
                pol_str, ';', fs_str, ';');
            
            % Execute stimulation
            xippmex('stim',stimString);
            
            % Write to txt file
            fprintf(fileID,[repmat('%1.0f ', 1, length(stimElectrodes{c})),'\t %1.0f \t %1.0f \t %5.2f \t %4.2f \t %4.2f \t %5.2f \t %5.2f\r\n'],...
                stimElectrodes{trialParams(7)},block,pulses,freq,caDur,anDur,caAmp,anAmp);
            
            %toc(forloop)
            t1 = clock;
            
            pause(((1000*(1/trainFrequency)) - (etime(t1,t0)+trainLength_ms))/1000);
            
        end
        block=block+1;
    end
    fclose(fileID)
end

%% %%%%%%%%%%%%%%% NOTES %%%%%%%%%%
%check electode 132 is labeled 134
