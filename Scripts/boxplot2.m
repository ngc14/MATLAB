function h = boxplot2(varargin)
%BOXPLOT2 Enhanced boxplot plots
% 
% h = boxplot2(y)
% h = boxplot2(y,x)
% h = boxplot2(..., p1, v1, ...)
%
% I don't like the original boxplot function... it messes with axis
% dimensions, replaces tick labels with text (even if no rotation is
% needed), and in general hijacks the axis too much.  Plus it doesn't
% return the handles to the plotted objects in an easy-to-use way (I can
% never remember which row corresponds to which part of the boxplot).  This
% version creates more light-handed boxplots, assuming that any cosmetic
% changes (color, tick labels, line specs, etc) can be added by the user
% afterwards if necessary.  It also allows one to create clustered
% boxplots, similar to an unstacked bar graph.
%
% Input variables:
%
%   y:              either a ndata x nx array (as in boxplot) or nx x ny x
%                   ndata array where nx indicates the number of
%                   x-clusters, ny the number of boxes per cluster, and
%                   ndata the number of points per boxplot.
%
%   x:              vector of x locations for each box cluster
%
% Optional input variables (passed as parameter/value pairs)
%
%   notch:          'on' or 'off' ['off']
%
%   orientation:    'vertical' or 'horizontal' ['vertical']
%
%   barwidth:       Barwidth value used to position boxes (see bar) [0.8]
%
%   whisker:        whisker length factor (see boxplot) [1.5]
%
%   axes:           axis to plot to [gca]
%
% Output variables:
%
%   h:              structure of handles to boxplot, with the following
%                   fields: 
%                   'box':      box
%                   'ladj':     lower adjacent value
%                   'lwhis':    lower whisker
%                   'med':      median
%                   'out':      outliers
%                   'uadj':     upper adjacent value
%                   'uwhis':    upper whisker    
%

% Copyright 2012 Kelly Kearney       

% Parse input

p = inputParser;
p.addRequired('y', @isnumeric);
p.addOptional('x', [], @isnumeric);
p.addParamValue('notch', 'off', @ischar);
p.addParamValue('orientation', 'vertical', @ischar);
p.addParamValue('axes', gca, @(x) isscalar(x) && ishandle(x) && strcmp(get(x,'type'),'axes'));
p.addParamValue('barwidth', 0.8);
% p.addParamValue('boxwidth', [], @(x) isscalar(x));
p.addParamValue('whisker', 1.5, @(x) isscalar(x));


p.parse(varargin{:});

In = p.Results;
In.notch = validatestring(In.notch, {'on', 'off'});
In.orientation = validatestring(In.orientation, {'vertical', 'horizontal'});

if ndims(In.y) == 2
    In.y = permute(In.y, [2 3 1]);
end
[nx, ny, ndata] = size(In.y);

if isempty(In.x)
    In.x = 1:nx;
end

ybox = reshape(In.y, [], ndata)';

% Use bar graph to get x positions

figtmp = figure('visible', 'off');
try
    hax = axes;
    hb = bar(In.x, In.y(:,:,1), In.barwidth);
    % xb = linspace(-.5,.5,nx+2);
    % xb = arrayfun(@(a) In.x + a, xb(2:end-1),'UniformOutput',false);
    xb = reshape(categorical(["Baseline","Go","Reach","Grasp","Withdraw"]).*categorical(["Arm";"Hand"]),1,[]);
    if verLessThan('matlab', '8.4.0')
        boxwidth = diff(minmax(xbar(:,1)));
    else
        if ny > 1
            boxwidth = In.barwidth;
        else
            mindx = min(diff(In.x));
            boxwidth = mindx .* In.barwidth;
        end
    end
    delete(hb);
    yv = reshape(ybox,1,[]);
    nanY = isnan(yv);
    xv = reshape(repmat(xb,size(ybox,1),1),1,[]);
    cg = reshape(repmat(repmat(1:nx,ny,1),size(ybox,1),1),1,[]);
    bx= boxchart(hax,ybox, ...
                  'notch', In.notch, ...
                  'orientation', In.orientation, ...
                  'markerstyle', '+', ...
                  'boxwidth',  boxwidth, ...
                  'linewidth', In.whisker);
     bx   = copyobj(findall(hax,'Type','Boxchart'), In.axes);
%     h.ladj  = copyobj(findall(bx, 'tag', 'Lower Adjacent Value'), In.axes);
%     h.lwhis = copyobj(findall(hax, 'tag', 'Lower Whisker'), In.axes);
%     h.med   = copyobj(findall(hax, 'tag', 'Median'), In.axes);
%     h.out   = copyobj(findall(hax, 'tag', 'Outliers'), In.axes);
%     h.uadj  = copyobj(findall(hax, 'tag', 'Upper Adjacent Value'), In.axes);
%     h.uwhis = copyobj(findall(hax, 'tag', 'Upper Whisker'), In.axes);
%    h = structfun(@(x) reshape(flipud(x), ny, nx), h, 'uni', 0);
    
    close(figtmp);%set(figtmp,'Visible','on')
catch ME
    close(figtmp);
    rethrow(ME);
end

    h.box=bx;



