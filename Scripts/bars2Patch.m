function patches = bars2Patch(yVals,groupingVar)
%       yVals: grouped barplot convention
%           N x M where N = # of observations and M = # of groups
%       patches: cell array of 1 x M patches

numGroups = size(yVals,2);
if(nargin==1)
    groupingVar = 1:numGroups;
end
xVals = 1:size(yVals,1);
xSpacing = linspace(0,.4, numGroups*2);
xSpacing = xSpacing - range(xSpacing)/2;
barWidth = mode(diff(xSpacing));
xSpacing = xSpacing(1:2:end);
pX = arrayfun(@(o) cell2mat(cellfun(@(a) a+o,num2cell(cell2mat( arrayfun(@(s) ...
     [s,s,s+barWidth,s+barWidth],xVals,'UniformOutput',false)'),2),'UniformOutput',false)),xSpacing, 'UniformOutput',false);
pY =  cellfun(@(a) cell2mat(arrayfun(@(o) [0 o o 0], a, 'UniformOutput', false)),...
    num2cell(yVals,1),'UniformOutput',false);
patches = cellfun(@(x,y,n) patch(x',y',n),pX,pY,num2cell(groupingVar), 'UniformOutput',false);
xticks(xVals);
end