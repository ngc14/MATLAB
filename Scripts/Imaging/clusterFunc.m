function clusterResults = clusterFunc(data,k)
optionSet = statset('UseParallel',true,'UseSubstreams',true,'Streams',RandStream('mlfg6331_64'),'Display','final');
clusterOpts= struct('Distance',@(xi,xj)dFunc(xi,xj),...
    'Algorithm','large','NumSamples',[],'PercentNeighbors',0.001,...
    'Replicates',3,'OnlinePhase','on','Options',optionSet);
startSample = randsample(length(data),k*clusterOpts.Replicates);
startMasks = permute(reshape(data(startSample,:)',[],k,clusterOpts.Replicates),[2 1 3]);
clusterOpts.Start = startMasks(1:k,:,:);
clo = namedargs2cell(clusterOpts);
clusterResults = kmedoids(data,k,clo{:});
end