function plot2dcolor(x,y)

nrX=numel(x);
cmap=colormap(hsv(nrX));
hold on
for idx=2:nrX
    plot(x(idx-1:idx),y(idx-1:idx),'Color',cmap(idx,:));
end
hold off

% 
% %x = 0:.05:2*pi;
% %y = sin(x);
% x=x(:)';
% y=y(:)';
% z = zeros(size(x));
% col = x;  % This is the color, vary with x in this case.
% surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
