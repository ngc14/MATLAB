function bins = findBins(timesT,allBins)
if(~any(~isnan(timesT)))
    bins =NaN(size(timesT));
else
    bins = int64(discretize(timesT,allBins));
end
end