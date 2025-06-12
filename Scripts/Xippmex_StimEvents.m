clear all;
numberOfPulses         = 33;                 %number of pulses in train
frequency_Hz           = 333;               %frequency of pulse period (Hz)
aDur                   = 0.2;               %duration of anode phase, positive (ms) (leave at .1 ms for no phase)
cDur                   = 0.2;               %duration of cathode phase, negative (ms) (leave at .1 ms for no phase)
aAmp                   = 5;                 %height of anode phase's current (headstage steps - [0, 127]) (uA)
cAmp                   = 5;                %height of cathode phase's current (headstage steps - [0, 127]) (uA)
stimElectrodes = [2:2:16];


%for bipolar stimulation: 1 - cathodic first, 0 - anodic first
polarity               = 1;

%delay - delay for interleaving (ms)
electrodeDelay_ms      = 0;

%duration of time to fast settle after the full biphasic pulse (ms)
fastSettle = 1;

FEMap = [1:32];
probeMap = [1:2:16, 2:2:16, 17:32];
mappedElectrodes = FEMap(probeMap);
xippmex('close');
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

elecs = mappedElectrodes(stimElectrodes);
sizeE = size(elecs);
trainLength_ms=1000*(numberOfPulses/frequency_Hz);

% Generate stimulation string
elect_str = strcat('Elect=', sprintf('%d,', elecs));
tl_str = strcat('TL=', sprintf('%.3f,', repmat(trainLength_ms,sizeE)));
freq_str = strcat('Freq=', sprintf('%.0f,', repmat(frequency_Hz,sizeE)));
cDur_str = strcat('CathDur=', sprintf('%.3f,', repmat(cDur,sizeE)));
aDur_str = strcat('AnodDur=', sprintf('%.3f,', repmat(aDur,sizeE)));
cAmp_str = strcat('CathAmp=', sprintf('%d,', repmat(cAmp,sizeE)));
aAmp_str = strcat('AnodAmp=', sprintf('%d,', repmat(aAmp,sizeE)));
delay_str = strcat('TD=', sprintf('%.3f,', repmat(electrodeDelay_ms,sizeE)));
pol_str = strcat('PL=', sprintf('%d,', repmat(polarity,sizeE)));
fs_str = strcat('FS=', sprintf('%.4f,', repmat(fastSettle,sizeE)));

stimString = strcat(elect_str, ';', tl_str, ';', freq_str, ';',...
    cDur_str, ';', aDur_str, ';', cAmp_str, ';', aAmp_str, ';', delay_str, ';',...
    pol_str, ';', fs_str, ';');
% Give the NIP some time to process any commands we have sent
pause(0.5)
storeEnd = NaN;
pulseLength = 0;
trialStart = 0;
xippmex('digin');
while(true)
    [counts,ts,eventLines] = xippmex('digin');
    if(any([eventLines.sma2]>0))
        % Execute stimulation
        xippmex('stim',stimString);
        disp('Stim Delivered');
        pause(trainLength_ms/1000);
    end
    pause(0.0025);
end
