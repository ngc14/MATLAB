function badTrial = getBadTrials(eventData, segmentTimes,Fs)
    badTrial = zeros(size(eventData));
    for t=1:size(eventData,2)
        trialSig = eventData{t};
        rectified = movmean(abs(trialSig),20);
        baseline = mean(rectified(1:segmentTimes{t}(2)),'omitnan');
        Y = fft(trialSig);
        P2 = abs(Y/length(trialSig));
        P1 = P2(1:floor(length(trialSig)/2) + 1);
        P1(2:end-1) = 2*P1(2:end-1);
        freq = Fs*(0:(length(trialSig)/2))/length(trialSig);
        if(any(find(freq<15 & P1>7)))
            badTrial(t) = 1;
        elseif(range(abs(trialSig)) < 5)
            badTrial(t) = 2;
        elseif(baseline>80)
            badTrial(t) = 3;
        elseif(P1>10 & mod(freq,notch)~=0 & mod(freq,notch)<notch+.5 & mod(freq,notch)>notch-.5)
            badTrial(t) = 4;
        end
    end
end