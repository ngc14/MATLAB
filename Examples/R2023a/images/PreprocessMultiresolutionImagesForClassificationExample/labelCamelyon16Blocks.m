function [blockAndLabel, info] = labelCamelyon16Blocks(block, blockInfo, imageLabels)

% Copyright 2021 The MathWorks, Inc.

blockAndLabel = {block{1}, imageLabels(blockInfo.ImageNumber)};
info = [];

end
