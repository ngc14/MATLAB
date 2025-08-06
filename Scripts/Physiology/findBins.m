function bins = findBins(timesT,allBins)
if(~any(~isnan(timesT)))
    bins =NaN(size(timesT));
else
    bins = discretize(timesT,allBins);
end
end
