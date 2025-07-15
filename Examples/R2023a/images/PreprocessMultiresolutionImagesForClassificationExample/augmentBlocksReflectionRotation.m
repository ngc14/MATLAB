function out = augmentBlocksReflectionRotation(blockAndLabel)

% Copyright 2021 The MathWorks, Inc.

out = cell(4,2);

% Extract block, and spatially augment it
im = blockAndLabel{1};
out{1,1} = im;
out{2,1} = fliplr(im);
out{3,1} = flipud(out{2,1});
out{4,1} = rot90(im);

% All are the same class
out(:,2)  = deal(blockAndLabel(2));