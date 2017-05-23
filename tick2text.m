function tick2text(varargin)
%TICK2TEXT Changes tick labels to text objects
%
% tick2text
% tick2text(ax)
% tick2text(ax, param1, val1, ...)
% tick2text(ax, 'undo')
%
% This function allows for more versatile labeling of axis ticks.  It
% offers two possible modifications: 1) automatic reformatting of tick
% labels based on a user-supplied format or function, and 2) conversion of
% tick labels to text objects.
%
% The first option simply offers a shortcut to manual reformatting via
% the [X/Y/Z]TickLabel properties of axes.  It also overrides the zoom/pan
% function of the axis so that ticks will update with your specified
% formatting as you zoom and pan.
%
% With the second option, the original tick labels are removed, and new
% text objects are created that mimic x, y, and/or z tick labels. The text
% object format allows more modification (such as color, rotation, etc)
% than the original tick labels, and also allows TeX strings to be utilized
% in tick labels (see xformat and yformat inputs for more details).  
%
% The handles of these new tick-like text objects are saved in the
% application data of the axis as '[X/Y/Z]TickText'.  Using these handles,
% the text objects can be modified as a group or individually (this method
% is used, rather than simply returning the handles of the new text
% objects, since zooming and panning will cause new objects to be created
% and old ones to be deleted).
%
% Note that with the conversion-to-text method, there will be a delay in
% the relabeling when panning, i.e. the old labels will drag with the mouse
% until you release the button, at which point the new labels will appear.
% Also note that the new tick label format cannot be updated if you reset
% tick locations or axis limits manually, so apply tick2text to an axis
% after these adjustments have already been made. 
%
% This function should be robust to all axis properties for 2D plots (i.e.
% log and linear, normal and reversed axis order, and left/right/top/bottom
% axis location).  This is also true in 3D, but it is not robust to
% rotation; if your axis is rotated such that the y-axis is to the left,
% ticks will not align properly.
%
% Input variables:
%
%   ax:         handle of axis to be modified (default to current axis)
%
% Optional input variables (passed as parameter/value pairs, defauts in [])
%
%   axis:       string containing characters 'x', 'y', and/or 'z',
%               indicating whether to format xticks, yticks, and/or
%               zticks, respectively. ['xy']  
%
%   convert:    logical scalar, or 1 x n array, where n is the length of
%               the axis string, indicating whether to replace ticks with
%               text objects (true) or simply reformat the text of existing
%               ticks (false).
%
%   xformat:    formatting for x tick label.  This can be either a
%               formatting string (see sprintf) or the handle to a function
%               that takes as input a scalar number and returns a character
%               array. [@num2str]
%
%   yformat:    formatting for y tick label.  This can be either a
%               formatting string (see sprintf) or the handle to a function
%               that takes as input a scalar number and returns a character
%               array. [@num2str]
%
%   zformat:    formatting for z tick label.  This can be either a
%               formatting string (see sprintf) or the handle to a function
%               that takes as input a scalar number and returns a character
%               array. [@num2str]
%
%   panon:      logical scalar.  If true, set the customized pan
%               function to update tick labels when panning.  Because
%               custom panning can only be set for a figure, not an
%               individual axis, you may want to disable this if you are
%               applying customized ticks to multiple axes in a figure.
%               [true]
%
%   'undo':     delete any added text objects associated with tick2text,
%               and reset tick labels to their defaults (auto mode).  Note
%               that this doesn't undo non-convert formatting; that can be
%               reset via set(ax, '[XYZ]TickLabelMode', 'auto');
%
% Example:
%   
%   This example creates two axes.  The top axis shows Matlab's default
%   tick labels.  The tick2text function is applied to the bottom axis,
%   changing the x axis to show ticks as multiples of pi, and the y axis to
%   label with full fixed-point values.  To demonstrate how to retrieve the
%   new tick handles, the x labels are then rotated 340 degrees (note that
%   by default x-text-ticks are centered; they may need to be left- or
%   right-aligned if rotated).
%
%     x = linspace(0,2*pi);
%     y = sin(x) + 100000;
% 
%     ax(1) = subplot(2,1,1);
%     plot(x,y);
%     title('Original tick marks');
% 
%     ax(2) = subplot(2,1,2);
%     plot(x,y);
%     title('Modified tick marks');
% 
%     set(ax, 'xlim', [0 2*pi]);
% 
%     tick2text(ax(2), 'yformat', '%.2f', ...
%                      'xformat', @(x) sprintf('%.2g\\pi', x/pi), ...
%                      'convert', [true false])
% 
%     hx = getappdata(ax(2), 'XTickText');
%     set(hx, 'Rotation', 340, 'fontsize', 12, 'horiz', 'left');

% Copyright 2009 Kelly Kearney

%----------------------------
% Parse and check input
%----------------------------

% Get axis handle

if ishandle(varargin{1})
    ax = varargin{1};
    pv = varargin(2:end);
else
    ax = gca;
    pv = varargin;
end

if ~isscalar(ax) || ~ishandle(ax) || ~strcmp(get(ax, 'type'), 'axes')
    error('Handle must be to a single axis');
end

% Check for undo option

if length(pv) == 1 && strcmp(pv{1}, 'undo')
    Ad = getappdata(ax);
    if isfield(Ad, 'YTickText')
        delete(Ad.YTickText);
        set(ax, 'yticklabelmode', 'auto');
        rmappdata(ax, 'YTickText');
    end
    if isfield(Ad, 'XTickText')
        delete(Ad.XTickText);
        set(ax, 'xticklabelmode', 'auto');
        rmappdata(ax, 'XTickText');
    end
    if isfield(Ad, 'ZTickText')
        delete(Ad.ZTickText);
        set(ax, 'zticklabelmode', 'auto');
        rmappdata(ax, 'ZTickText');
    end
    return
end

% Get parent figure (for panning)

fig = ancestor(ax, 'figure'); 

% Parse and check optional variables

Options = struct('xformat', [], 'yformat', [], 'zformat', [], ...
                 'ytickoffset', .04, 'xtickoffset', .04, ...
                 'ztickoffset', .08, 'axis', 'xy', 'panon', true, ...
                 'convert', true);
             
Options = parse_pv_pairs(Options, pv);

if ~isscalar(Options.ytickoffset) || ~isscalar(Options.xtickoffset) || ~isscalar(Options.ztickoffset)
    error('Offset values must be numerical scalars');
end

if ~ischar(Options.axis) || ~isempty(regexp(Options.axis, '[^xyz]'))
    error('Axis must be string including x, y, and/or z');
end
loc = arrayfun(@(a) strfind(Options.axis, a), 'xyz', 'uni', 0);
labelflag = ~cellfun('isempty', loc);

convertflag = false(1,3);
if isscalar(Options.convert)
    convertflag(:) = Options.convert;
elseif isequal(size(Options.axis), size(Options.convert))
    for ii = 1:length(Options.axis)
        convertflag(loc{ii}) = Options.convert(ii);
    end
else
    error('Convert array must be same length as axis string');
end
    

% Set formatting functions

if isempty(Options.xformat)
    xformatfun = @num2str;
elseif ischar(Options.xformat)
    xformatfun = @(x) num2str(x, Options.xformat);
elseif strcmp(class(Options.xformat), 'function_handle')
    xformatfun = Options.xformat;
else 
    error('xformat must be a formatting string or function handle');
end

if isempty(Options.yformat)
    yformatfun = @num2str;
elseif ischar(Options.yformat)
    yformatfun = @(x) num2str(x, Options.yformat);
elseif strcmp(class(Options.yformat), 'function_handle')
    yformatfun = Options.yformat;
else 
    error('yformat must be a formatting string or function handle');
end

if isempty(Options.zformat)
    zformatfun = @num2str;
elseif ischar(Options.zformat)
    zformatfun = @(x) num2str(x, Options.zformat);
elseif strcmp(class(Options.zformat), 'function_handle')
    zformatfun = Options.zformat;
else 
    error('zformat must be a formatting string or function handle');
end

P.ax = ax;
P.xfun = xformatfun;
P.yfun = yformatfun;
P.zfun = zformatfun;
P.lflag = labelflag;
P.cflag = convertflag;

%----------------------------
% Change exisiting tick 
% labels
%----------------------------

texttick(P);

% Add customized zoom
    
try
    hzoom = zoom(ax);
    set(hzoom, 'ActionPostCallback', {@newzoompan, P});
catch
    warning('Unable to apply customized zoom');
end
    
% Add customized pan

if Options.panon
    hpan = pan(fig);
    set(hpan, 'ActionPostCallback', {@newzoompan, P});
end

%----------------------------
% Update tick labels when
% zooming and panning
%----------------------------

function newzoompan(h, ed, P)
if P.lflag(1) && ~P.cflag(1)
    set(P.ax, 'XTickMode', 'auto', 'XTickLabelMode', 'auto');
end
if P.lflag(2) && ~P.cflag(2)
    set(P.ax, 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
end
if P.lflag(3) && ~P.cflag(3)
    set(P.ax, 'ZTickMode', 'auto', 'ZTickLabelMode', 'auto');
end

texttick(P)

%----------------------------
% Change tick labels to text 
% objects
%----------------------------

function texttick(P)

% Get current values of ticks and axis limits

A = get(P.ax);

% Create new labels

xticklab = arrayfun(P.xfun, A.XTick, 'uni', 0);
yticklab = arrayfun(P.yfun, A.YTick, 'uni', 0);
zticklab = arrayfun(P.zfun, A.ZTick, 'uni', 0);

% If no-convert option, simply replace tick labels

if P.lflag(1) && ~P.cflag(1)
    set(P.ax, 'XTick', A.XTick, 'XTickLabel', xticklab);
end
if P.lflag(2) && ~P.cflag(2)
    set(P.ax, 'YTick', A.YTick, 'YTickLabel', yticklab);
end
if P.lflag(3) && ~P.cflag(3)
    set(P.ax, 'ZTick', A.ZTick, 'ZTickLabel', zticklab);
end

% Exit if no conversion calcs necessary

if ~any(P.cflag)
    return
end

% Determine location of ticks in coordinate space

if isequal(A.View, [0 90])
    is2d = true;
elseif isequal(A.View, [-37.5 30])
    is2d = false;
else
    is2d = false;
    warning('Not robust to roation in 3D; ticks may be oddly positioned');
end

if is2d
    
    % Where in x-coords is the Y axis?
    
    islin = strcmp(A.YScale, 'linear');
    if islin
        pos = A.XLim(2) + A.TickLength(1) .* diff(A.XLim);
        neg = A.XLim(1) - A.TickLength(1) .* diff(A.XLim);
    else
        pos = 10.^(A.XLim(2) + A.TickLength(1) .* diff(log10(A.XLim)));
        neg = 10.^(A.XLim(1) - A.TickLength(1) .* diff(log10(A.XLim)));
    end
    
    if strcmp(A.XDir, 'normal') && strcmp(A.YAxisLocation, 'left')
        ytkx = neg;
        ytkh = 'right';
        ytkv = 'middle';
    elseif strcmp(A.XDir, 'normal') && strcmp(A.YAxisLocation, 'right')
        ytkx = pos;
        ytkh = 'left';
        ytkv = 'middle';
    elseif strcmp(A.XDir, 'reverse') && strcmp(A.YAxisLocation, 'left')
        ytkx = pos;
        ytkh = 'right';
        ytkv = 'middle';
    elseif strcmp(A.XDir, 'reverse') && strcmp(A.YAxisLocation, 'right')
        ytkx = neg;
        ytkh = 'left';
        ytkv = 'middle';
    end
        
    % Where in y-coords is the X axis?
    
    islin = strcmp(A.YScale, 'linear');
    if islin
        pos = A.YLim(2) + A.TickLength(1) .* diff(A.YLim);
        neg = A.YLim(1) - A.TickLength(1) .* diff(A.YLim);
    else
        pos = 10.^(log10(A.YLim(2)) + A.TickLength(1) .* diff(log10(A.YLim)));
        neg = 10.^(log10(A.YLim(1)) - A.TickLength(1) .* diff(log10(A.YLim)));
    end
    
    if strcmp(A.YDir, 'normal') && strcmp(A.XAxisLocation, 'bottom')
       	xtky = neg;
        xtkh = 'center';
        xtkv = 'top';
    elseif strcmp(A.YDir, 'normal') && strcmp(A.XAxisLocation, 'top')
        xtky = pos;
        xtkh = 'center';
        xtkv = 'bottom';
    elseif strcmp(A.YDir, 'reverse') && strcmp(A.XAxisLocation, 'bottom')
        xtky = pos;
        xtkh = 'center';
        xtkv = 'top';
    elseif strcmp(A.YDir, 'reverse') && strcmp(A.XAxisLocation, 'top')
        xtky = neg;
        xtkh = 'center';
        xtkv = 'bottom';
    end
    
    xtkz = 0;
    ytkz = 0;
    
else
    
    % Where in y-coords are the X and Z axes?
    
    islin = strcmp(A.YScale, 'linear');
    if islin
        pos = A.YLim(2) + A.TickLength(2) .* diff(A.YLim);
        neg = A.YLim(1) - A.TickLength(2) .* diff(A.YLim);
    else
        pos = 10.^(A.YLim(2) + A.TickLength(2) .* diff(log10(A.YLim)));
        neg = 10.^(A.YLim(1) - A.TickLength(2) .* diff(log10(A.YLim)));
    end
    
    if strcmp(A.YDir, 'normal')
        xtky = neg;
        ztky = pos;
    else
        xtky = pos;
        ztky = neg;
    end
    
    % Where in x-coords are the Y and Z axes?
    
    islin = strcmp(A.YScale, 'linear');
    if islin
        pos = A.XLim(2) + A.TickLength(1) .* diff(A.XLim);
        neg = A.XLim(1) - A.TickLength(1) .* diff(A.XLim);
    else
        pos = 10.^(A.XLim(2) + A.TickLength(1) .* diff(log10(A.XLim)));
        neg = 10.^(A.XLim(1) - A.TickLength(1) .* diff(log10(A.XLim)));
    end
    
    if strcmp(A.XDir, 'normal')
        ztkx = A.XLim(1);
        ytkx = neg;
    else
        ztkx = A.XLim(2);
        ytkx = pos;
    end
    
    % Where in z-coords are the X and Y axes?
    
    if strcmp(A.ZDir, 'normal')
        xtkz = A.ZLim(1);
        ytkz = A.ZLim(1);
    else
        xtkz = A.ZLim(2);
        ytkz = A.ZLim(2);
    end
    
    % Alignment (this assumes the z-axis is to the left)
    
    xtkh = 'left';
    xtkv = 'top';
    ytkh = 'right';
    ytkv = 'top';
    ztkh = 'right';
    ztkv = 'middle';
    
end
    
% Replace ticks with text

if P.lflag(1) && P.cflag(1)
    xcoord = A.XTick;
    ycoord = ones(size(A.XTick)) * xtky;
    zcoord = ones(size(A.XTick)) * xtkz;
    if isappdata(P.ax, 'XTickText')
        replaceoldticks(P.ax, 'XTickText', xcoord, ycoord, zcoord, xticklab);
    else
        set(P.ax, 'xticklabel', '');
        createnewticks(P.ax, 'XTickText', xcoord, ycoord, zcoord, xticklab, xtkh, xtkv);
    end
end

if P.lflag(2) && P.cflag(2)
    xcoord = ones(size(A.YTick)) * ytkx;
    ycoord = A.YTick;
    zcoord = ones(size(A.YTick)) * ytkz;
    set(P.ax, 'yticklabel', '');
    if isappdata(P.ax, 'YTickText')
        replaceoldticks(P.ax, 'YTickText', xcoord, ycoord, zcoord, yticklab);
    else
        createnewticks(P.ax, 'YTickText', xcoord, ycoord, zcoord, yticklab, ytkh, ytkv);
    end
end

if P.lflag(3) && P.cflag(3)
    xcoord = ones(size(A.ZTick)) * ztkx;
    ycoord = ones(size(A.ZTick)) * ztky;
    zcoord = A.ZTick;
    set(P.ax, 'zticklabel', '');
    if isappdata(P.ax, 'ZTickText')
        replaceoldticks(P.ax, 'ZTickText', xcoord, ycoord, zcoord, zticklab);
    else
        createnewticks(P.ax, 'ZTickText', xcoord, ycoord, zcoord, zticklab, ztkh, ztkv);
    end
end


%----------------------------
% Create new tick-like text
% objects
%----------------------------
    
function createnewticks(ax, appfld, xcoord, ycoord, zcoord, txt, horiz, vert);    

h = text(xcoord, ycoord, zcoord, txt, 'parent', ax, 'horizontalalignment', horiz, 'verticalalignment', vert);
setappdata(ax, appfld, h);

%----------------------------
% Replace old tick-like text
% objects
%----------------------------

function replaceoldticks(ax, appfld, xcoord, ycoord, zcoord, txt)

hxold = getappdata(ax, appfld);

% Get properties of current ticks

props = get(hxold);
params = fieldnames(props(1));
vals = struct2cell(props(1));

% Ignore properties that always change between ticks or are read only

[tf,loc] = ismember({'Extent', 'String', 'BeingDeleted', 'Children', 'Position', 'Type', 'Annotation'}, params);
params(loc) = [];
vals(loc) = [];
paramval = [params'; vals'];

% Create new ticks with same properties

delete(hxold);
hx = text(xcoord, ycoord, zcoord, txt, 'parent', ax, paramval{:});      
setappdata(ax, appfld, hx);
