function plotHeatMatrix(hm,hmNames,entryNames,cLimits)
for c = 1:length(hm)
    nexttile; hold on; axis tight; axis ij;
    title(hmNames(c));
    imagesc(hm{c});
    xticks(1:2:length(entryNames));xticklabels(entryNames(1:2:end));yticklabels([]);clim(cLimits);
end
end