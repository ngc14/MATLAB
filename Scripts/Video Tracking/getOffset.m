function offset = getOffset(prev, next)
% get information of markers in previous frame
pProps = regionprops(logical(prev(:,:,1)), 'Centroid', 'Area', 'EquivDiameter');
ind = [];
% only look at markers with area greater than 100 and equivdiameter greater
% than 8 (reference marker should always be large)
for i = 1:size(pProps,1)
    if(pProps(i).Area > 100 && pProps(i).EquivDiameter > 8)
        ind(end+1) = i;
    end
end
pProps = pProps(ind);
pVals = cat(1, pProps.Centroid);
% get the absolute minimum marker in the forward/backward direction
[~, pInd] = min(pVals(:,1));
% rounded position for indexing
loc = round(pVals(pInd,:));
% color of minimum marker
c = prev(loc(2), loc(1),:);
% if centroid happens to be black -> find nearest color value
while(~any(c))
    c = prev(loc(2), loc(1)+1,:);
    loc(1) = loc(1) + 1;
end
% match the next frame reference marker excluding markers that are not of
% the same color based on the closest distance
nProps = regionprops(logical(next(:,:,1) == c(:,:,1) & ...
    next(:,:,2) == c(:,:,2) & next(:,:,3) == c(:,:,3)), 'Centroid');
dists = cellfun(@(x) norm(x-pVals(pInd,:),2),{nProps.Centroid});
[~, nInd] = min(dists);

offset = pProps(pInd).Centroid - nProps(nInd).Centroid;
end