function aligned = realign(signal, shiftVal, alignSeg)
    aligned = nan(1,alignSeg+shiftVal);
    aligned(alignSeg:-1:alignSeg-shiftVal+1) = signal(1:shiftVal);
    aligned(alignSeg+1:alignSeg+length(signal)-shiftVal) = signal(shiftVal+1:end);
end
