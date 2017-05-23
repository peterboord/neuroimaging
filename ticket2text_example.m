x = linspace(0,2*pi);
y = sin(x) + 100000;

ax(1) = subplot(2,1,1);
plot(x,y);
title('Original tick marks');

ax(2) = subplot(2,1,2);
plot(x,y);
title('Modified tick marks');

set(ax, 'xlim', [0 2*pi]);

tick2text(ax(2), 'yformat', '%.2f', ...
    'xformat', @(x) sprintf('%.2g\\pi', x/pi), ...
    'convert', [true false])

hx = getappdata(ax(2), 'XTickText');
set(hx, 'Rotation', 340, 'fontsize', 36, 'horiz', 'left');
