function demo
% Allow a line to have its own 'ButtonDownFcn' callback.
h = zoom;
% hAx=axes;
hLine = plot(rand(1,10));
hLine.Tag = 'DoNotIgnore';
hLine.ButtonDownFcn = {@myLineButtonDownFcn,h};
hAx.ButtonDownFcn = {@myLineButtonDownFcn,h};
h.ButtonDownFilter = {@myButtonDownFilter,h};
% h.ActionPreCallback= {@myActionPreCallback,h};
% h.ActionPostCallback={@myActionPostCallback,h};
h.Enable = 'on';
% mouse click on the line
%

function myLineButtonDownFcn(obj,event_obj,h)
disp('myLineButtonDownFcn');
cP = get(gca,'Currentpoint');
x = cP(1,1);
y = cP(1,2);
disp([x,y]);
h.Enable = 'off';

% function myAxButtonDownFcn(obj,event_obj,h)
% disp('myAxButtonDownFcn');
% cP = get(gca,'Currentpoint');
% x = cP(1,1);
% y = cP(1,2);
% disp([x,y]);
% h.Enable = 'off';


function [flag] = myButtonDownFilter(obj,event_obj,h)
flag=false;
% % If the tag of the object is 'DoNotIgnore', then return true.
% disp('myButtonDownFilter');
if strcmpi(get(gcf,'SelectionType'),'alt')
    flag=true;
end
disp(flag)



% function myActionPreCallback(obj,event_obj,h)
% % disp('myActionPreCallback');
% % cP = get(gca,'Currentpoint');
% % x = cP(1,1);
% % y = cP(1,2);
% %disp([x,y]);
% %h.Enable = 'off'
% 
% function myActionPostCallback(obj,event_obj,h)
% %disp('myActionPostCallback');
% % cP = get(gca,'Currentpoint');
% % x = cP(1,1);
% % y = cP(1,2);
% %disp([x,y]);
% %h.Enable = 'off'
