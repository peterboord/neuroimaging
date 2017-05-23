function biplot(x1,y1,x2,y2,xname1,yname1,xname2,yname2,titleName)
ax1 = axes('XLimMode','manual',...
    'XLim',[min(x1),max(x1)]);
line(x1,y1,'Parent',ax1,'Color','b')
% current axes
ax1.XColor = 'r';
ax1.YColor = 'r';
if nargin > 4
    xlabel(xname1);
    ylabel(yname1);
end
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XLimMode','manual',...
    'XLim',[min(x2),max(x2)]);
line(x2,y2,'Parent',ax2,'Color','g')
if nargin > 6
    xlabel(xname2);
    ylabel(yname2);
end
if nargin > 8
    title(titleName);
end
end