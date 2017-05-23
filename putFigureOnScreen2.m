function putFigureOnScreen2
monitorPositions = get(0,'MonitorPositions');
set(gcf, 'Position', monitorPositions(2,:));
end