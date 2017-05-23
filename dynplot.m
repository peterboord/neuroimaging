function handles = dynplot(y,e,varargin)
% FORMAT dynplot(y,e,varargin)
% 
% Draws what is commonly known as 'dynamite plot'.
% 
% 
%       |     _
%    15-|    _|_ 
%       |   |:::|  
%    10-|   |:::|           _
%       |   |:::| _        _|_
%     5-|   |:::|_|_      |:::| _
%       |   |:::|   |     |:::|_|_
%     0_|___|:::|   |_____|:::|   |___                
%       |    Group A       Group B         
% 
% 
% 
% 
% Required Input:
%   - y: n-by-k matrix of mean-values (n groups and k conditions)
%   - e: n-by-k matrix of corresponding standard errors
% 
% 
% Additional input may be passed as key-value pairs:
%   - 'bwr',[0.8]               bar-width to distance ratio
%   - 'elims',['smart']         error-limits; per default, upper limits are
%                               shown on positive bars, negative limits on 
%                               negative bars; one of 'upper','lower',
%                               'both','smart' have to be chosen
%   - 'hideaxis,[1]             hide x and y-axis
%   - 'savefig',[0]             save figure to disk
%   - 'fname',''                file name for figure
%   - 'legstr',[{{}}]           cell array of strings containing labels;
%                               a legend is drawn if this cell array is 
%                               not empty
%   - 'legloc',['EastOutside']  location for legend
%   - 'visible',['on']          per default, plots _are_ actually plotted, 
%                               but it may be convenient in certain cases 
%                               to just get the handle to a hidden figure
% 
% A structure of handles to the drawn elements is returned.
% 
% Examples:
% 
%     (1) a very basic call. Two groups, 4 conditions:
%     y = rand(2,4)-0.5
%     e = rand(2,4)/10
%     h = dynplot(y,e)
% 
%     (2) same as (1), but with a little more labeling:
%     y = rand(2,4)-0.5
%     e = rand(2,4)/10
%     h = dynplot(y,e,'xlab',{{'Group A'},{'Group B'}},'ylab','Mean Value')
% 
%     (3) same as (2), but with  a legend:
%     y = rand(2,4)-0.5
%     e = rand(2,4)/10
%     xlab = {{'Group A'},{'Group B'}};
%     ylab = 'Mean Value';
%     legstr = {'Cond I','Cond II','Cond III','Cond IV'};
%     h = dynplot(y,e,'xlab',xlab,'ylab',ylab,'legstr',legstr)
% 
%     (4) six groups (A-E), two conditions (1,2), upper and lower
%     errors, increased distance between bars, colors light gray and white: 
%     y = rand(6,2)-0.5; 
%     e = rand(6,2)/10; 
%     xlab = {{'A'},{'B'},{'C'},{'D'},{'E'},{'F'}}; 
%     ylab = 'Mean Value'; 
%     legstr = {'Condition 1','Condition 2'}; 
%     h = dynplot(y,e,'ycol',repmat({[0.8 0.8 0.8],'w'},1,prod(size(y))),...
%             'xlab',xlab,'legstr',legstr,'ylab',ylab,'elims','both',...
%             'bwr',0.5)
% -------------------------------------------------------------------------
% Apparently, SAS is unable to produce these (exact) kinds of plots. So I
% had to come up with a solution since I was annoyed by having to use MS
% Excel. Matlab needs a little hack to group the plots the way I like it,
% the bar-function by itself places the bars too far away from each other
% for my taste. Apart from that, there might be good reasons not to be too
% enthusiastic about these figures. Thus, be aware that these plots,
% although common in science, are not appreciated by all statisticians, and
% you might want to watch out of the dynamite.
% 
% See for instance:
%   - http://biostat.mc.vanderbilt.edu/wiki/Main/DynamitePlots
%   - http://biostat.mc.vanderbilt.edu/wiki/pub/Main/TatsukiKoyama/Poster3.pdf
%   - http://flowingdata.com/2008/02/15/how-to-read-and-use-a-box-and-whisker-plot/
%   - http://pablomarin-garcia.blogspot.co.nz/2010/02/why-dynamite-plots-are-bad.html
% 
% Disadvantages include (see Vanderbilt-site for more):
%   - low 'data-to-ink' ratio 
%   - raw data is hidden (only means + se of groups are shown)
%   - symmetric confidence intervals are assumed
%   - parametric assumptions about the data
% 
% However, there are more positive reviews as well:
%   - http://emdbolker.wikidot.com/blog:dynamite
% 
% -------------------------------------------------------------------------
% Credits:
%   - inspired by the solution as proposed by 'Jos' on mathworks
%     http://www.mathworks.com/matlabcentral/newsreader/view_thread/127890
%   - input parsing stuff is by 'Jonas' from stackoverflow.com
%     http://stackoverflow.com/questions/2775263/
% 
% I haven't tested most of these, but there are obviously alternatives on
% FileExchange that you might find useful:
%   - http://www.mathworks.com/matlabcentral/fileexchange/9377
%   - http://www.mathworks.com/matlabcentral/fileexchange/24718
%   - http://www.mathworks.com/matlabcentral/fileexchange/27387
%   - http://www.mathworks.com/matlabcentral/fileexchange/25613 
%   - http://www.mathworks.com/matlabcentral/fileexchange/9377
%   - http://www.mathworks.com/matlabcentral/fileexchange/9541
% 
% See also: bar, errorbar
% -------------------------------------------------------------------------
% Vers 0.11, 2013-02-28, Philipp Homan

% arguments and defaults
% test for first argument (obligatory)
if nargin < 2
  error('Need at least 2 input arguments');
end
  
% define defaults 
opts = struct(...
  'bwr',0.8,...
  'elims','smart',...
  'hideaxis',1,...
  'x',-1,...
  'ycol',{{}},...
  'xlab',{{}},...
  'ylab','Y',...
  'figname',' ',...
  'figh',-1,...
  'legstr',{{}},...
  'xcatstr',{{}},...
  'savefig',0,...
  'visible','on',...
  'fname',['dynplot_' num2str(now) '.tiff'],...
  'legloc','EastOutside' ...
  );

% read acceptable names
optkeys = fieldnames(opts);

% count arguments
nargs = length(varargin(1:end));
if round(nargs/2) ~= nargs/2
   error('need key/value pairs!')
end

for pair = reshape(varargin,2,[]) % pair is {key/value}
   inpkey = lower(pair{1}); 
   
   % overwrite options. 
   if any(strmatch(inpkey,optkeys))      
      opts.(inpkey) = pair{2};
   else
      error('%s is not a recognized parameter name',inpkey)
   end
end

% set each key to val
for i=1:numel(optkeys)
  eval([sprintf('%s = ', optkeys{i}) 'getfield(opts,optkeys{i});']);
end

% set defaults to values if necessary
if figh <= 0
  figh = figure('Visible',visible);
end

if x < 0
  x = [1:size(y,1)]';
end

if isempty(ycol)
  ycol = cell(1,size(y,2)); 
  
  % start with white per default
  cm = flipud(gray(size(y,2)));
  
  for i=1:size(cm,1)
    ycol{i} = cm(i,:);
  end
end

if isempty(xlab)
  xlab{1} = cell(1,size(y,2));
  for i=1:size(y,2)
    xlab{1}{i} = num2str(i);
  end
  for i=1:size(y,1)
    xlab{i}=xlab{1}(1:end);
  end
end

% Main programm
% return value
handles = struct;

% set XLimMode
set(gca,'XLimMode','manual');
set(gca,'XLim',[0 3]);

% cosmetic
set(figh,'Color','w');

% plot bars
ngroups = size(y,2);
intercept = x;
for i=1:ngroups
  if i == 1    
    offset = intercept - bwr/2*ngroups;    
  else
    offset = offset + bwr*1/ngroups;
  end
  bh(i) = bar(offset, y(:,i),bwr/ngroups);
  hold on
end
% bar-colors
for i=1:numel(bh)
  set(bh(i),'FaceColor',ycol{i});
end

% error limits
if ~isempty(e)
  drawerror(figh,x,y,e,bwr,elims);
end

% set up x-labels
xtick = getxtick(figh,bh,x,y,xlab);
xlab = parsexlab(xlab);
ax = drawxlab(figh,gca,xtick,xlab);

% % annotations
if ~isempty(xcatstr)
  ah = drawannots(figh,gca,xtick,xcatstr);
end

% draw legend upon request
if ~isempty(legstr)
  lh = drawleg(figh,gca,bh,legstr,legloc);
end

% y-label
yl = drawylab(figh,ax,ylab);

% hide axis lines upon request
if hideaxis
  lims = axis;
  
  % hide y-axis
  ln(1) = hideax(figh,gca,[lims(1) lims(1) lims(3) lims(4)]);
  
  % hide x-axis only if there are negative values
  if ~isempty(find(y<0))    
    ln(2) = hideax(figh,gca,[lims(1) lims(2) lims(3) lims(3)]);
  end
end

% save figure if requested
if savefig
  if exist('lh')
    saveok = savefigure(figh,fname,gca,bh,yl,lh);
  else
    saveok = savefigure(figh,fname,gca,bh,yl);
  end
end

% format return values
handles.figh = figh;
handles.bh = bh;
if exist('eh')
  handles.eh = eh;
end
if exist('lh'),
  handles.lh = lh;
end
return

function eh = drawerror(fh,x,y,e,bwr,elims)
if nargin < 6
  elims = 'smart';
end
ngroups = size(y,2);
intercept = x;
ll = e;
ul = e;
switch elims
  case 'upper'    
    ll(:) = 0;
  case 'lower';
    ul(:) = 0;
  case 'smart'    
    ll(y>0) = 0;
    ul(y<0) = 0;
  otherwise
end;
    
for i=1:ngroups
  if i == 1
    offset = intercept - bwr*0.5*ngroups;
  else
    offset = offset + bwr*1/ngroups;
  end
  eh(i) = errorbar(offset,y(:,i),ll(:,i),ul(:,i),'LineStyle','None','Color','k');
  hold on
end
return

function xtick = getxtick(figh,bh,x,y,xlab)
for i=1:numel(bh)
  xdat(i,:) = get(bh(i),'XData');
end
% xtick = [xdat(:,1)' xdat(:,2)'];
xtick = reshape(xdat,1,prod(size(xdat)));
count = 1;
ncat = numel(xlab);
ng = numel(xlab{1});
while numel(xtick) > ncat*ng
%   tmp = xtick;
  for i = 1:numel(xtick)
    if i == numel(xtick)/2 || (mod(i,2)==0 && mod(numel(xtick)/2,2)==0)
      continue
    end
    if i < numel(xtick)
      tmp(count) = (xtick(i)+xtick(i+1))/2;
      count = count + 1;
    end
  end
  xtick = tmp;
  clear tmp
  count = 1;
end
return

function xlab = parsexlab(xlab)
ncat = numel(xlab);
ng = numel(xlab{1});
count = 1;
for i=1:ncat
 for j=1:ng
   tmp{count} = xlab{i}{j};
   count = count + 1;
 end
end
xlab = tmp;
return

function ax = drawxlab(figh,ax,xtick,xlab,ticklength)
if nargin < 5
  ticklength = [0 0];
end
set(ax, ...
  'FontSize',10,...
  'FontWeight','bold',...
  'box','off',...
  'XTick',xtick,...
  'XLimMode','manual',...
  'TickLength',ticklength,...
  'XTickLabel',xlab ...
  );
return

function yl = drawylab(figh,ax,ylab)
yl = get(ax,'ylab');
set(yl,...
  'String',ylab,...
  'FontSize',10,...
  'FontWeight','bold'...
  );
return


% annotations
% -currently not implemented
function ah = drawannots(figh,ax,xtick,xcatstr)
ah = struct;
while numel(xtick) > numel(xcatstr)
  dim = numel(xtick);
  xtick = reshape(xtick,2,dim/2);
  for i=1:dim/2
    tmp(i) = sum(xtick(:,i))/2;
  end
  xtick = tmp;
  clear tmp
end
for i=1:numel(xtick)
  t1=annotation('textbox');
  pos = get(t1,'Position');
  set(t1,'Units','pixels');
  set(t1,...
    'String',xcatstr{i},...
    'Position',[xtick(i) 0.1 pos(3)*1.8 pos(4)*0.4],...
    'FontWeight','bold',...
    'FontSize',10,...
    'HorizontalAlignment','center',...
    'EdgeColor','w' ...
    );
end
return

% hide axes
% apparently, the hidden axes are not reliably reproduced if the figure is
% saved to disk in formats other than .fig
function ln = hideax(figh,ax,lims,col)
if nargin < 4
  col = 'w';
end
ln=line(lims(1:2),lims(3:4),'Color',col);
return

% resize canvas
% -currently not implemented
function figh = resize(figh,figsize,ax)
% increase canvas
set(ax,'Units','pixels');
apos = get(ax,'Position');
pos = get(figh,'Position');
set(figh,'Position',[pos(1) pos(2)*0.5 pos(3) pos(4)*1.05]);
set(ax,'Position',[apos(1) apos(2)*1.2 apos(3:4)]);
set(ax,'Units','normalized');
pos = get(figh,'Position');
set(figh,'Position',[pos(1) pos(2) pos(3)*figsize pos(4)]);
pos = get(figh,'PaperPosition');
set(figh,'PaperPosition',[pos(1) pos(2) pos(3)*0.5 pos(4)*0.8]);
return

% draw legend
function lh = drawleg(figh,ax,bh,legstr,legloc)
cc = get(bh,'Children');
if ~iscell(cc)
  tmp{1} = cc;
  cc = tmp;
  clear tmp;
end
ccv = [];
for i=1:numel(cc)
  ccv = [ccv cc{i}];
end

%-currently not implemented
% insert line breaks if string is too long
% for i=1:numel(legstr)
%   if size(legstr{i},2) > 8
% %     legstrtmp = legstr{i};
%     legstr{i} = sprintf('%s\n%s',legstr{i}(1:7),legstr{i}(8:end));
%   end
% end

% suppress warning 
warning('off');
lh = legend(ccv,legstr,'Location',legloc);
set(lh,...
  'YColor','w',...
  'XColor','w',...
  'FontSize',8 ...
  );
warning('on');

% legend cosmetic
lh = beautify(lh);
return

% legend cosmetic
function lh = beautify(lh)  
cc = get(lh,'Children');
for i=1:2:numel(cc)
  xd = get(cc(i),'XData');
  set(cc(i),'XData',[xd(1:2)' xd(3)*0.5 xd(4)*0.5]');
end
for i=2:2:numel(cc)
  pos = get(cc(i),'Position');
  set(cc(i),'Position',[pos(1)*0.6 pos(2:3)]);
end
% p = get(lh,'Position');
% set(lh,'Position',[p(1)*1 p(2)*0.98 p(3) p(4)*1.1]);
% set(lh,'Position',[p(1)*1 p(2)*1.05 p(3) p(4)*1.1]);
return

function saveok = savefigure(figh,fname,ax,bh,yl,lh)
saveok = -1;
% reformat font size
fsa = get(ax,'FontSize');
fsy = get(yl,'FontSize');
if exist('lh')
  fsl = get(lh,'FontSize');
  set(lh,'FontSize',6);
end
set(ax,'FontSize',9);
set(yl,'FontSize',9);
lh = beautify(lh);
set(figh, 'PaperPositionMode', 'auto');
try
  print(figh,'-r300','-dtiff', fname);
  saveok = 1;
catch
  saveok = 0;
end
set(ax,'FontSize',fsa);
set(yl,'FontSize',fsy);
if exist('lh')
  set(lh,'FontSize',fsl);
  lh = beautify(lh);
end
return
