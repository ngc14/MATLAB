function bins = findBins(timesT,allBins)
if(~any(~isnan(timesT)))
    bins =NaN(size(timesT));
else
    if(any(timesT<min(allBins)))
        timesT(timesT<min(allBins)) = min(allBins);
    end
    bins = discretize(timesT,allBins);
end
end
