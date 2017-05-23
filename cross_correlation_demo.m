% function cross_correlation_demo()
% Demo to do a cross correlation to find the location of a sub-image in a larger image.
% and show how it doesn't work, while normxcorr2 does work.
clc; 
close all; 
clear all; 
workspace;
fontSize = 20;

% Read in standard MATLAB demo image.
grayImage = imread('cameraman.tif');
[rows, columns, numberOfColorChannels] = size(grayImage);
subplot(2, 3, 1);
imshow(grayImage, []);
axis on;
caption = sprintf('Original Grayscale Image\n%d rows by %d columns', rows, columns);
title(caption, 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

% Get a small sub-image around the cameraman's head
row1 = 37;
row2 = 80;
col1 = 96;
col2 = 134;
subImage = double(grayImage(row1:row2, col1:col2));
[rowsT, columnsT, numberOfColorChannels] = size(subImage);
subplot(2, 3, 4);
imshow(subImage, []);
axis on;
caption = sprintf('Template Sub-Image\n%d rows by %d columns', rowsT, columnsT);
title(caption, 'FontSize', fontSize);

% Now compute and display the regular cross correlation image.
crossCorrelationImage = xcorr2(double(grayImage), subImage);
[rows, columns, numberOfColorChannels] = size(crossCorrelationImage);
subplot(2, 3, 2);
imshow(crossCorrelationImage, []);
axis on;
caption = sprintf('"Regular" Cross Correlation Image\n%d rows by %d columns', rows, columns);
title(caption, 'FontSize', fontSize);

% Now compute and display the normalized cross correlation image.
normalizedCrossCorrelationImage = normxcorr2(subImage, grayImage);
[rows, columns, numberOfColorChannels] = size(normalizedCrossCorrelationImage);
subplot(2, 3, 5);
imshow(normalizedCrossCorrelationImage, []);
axis on;
caption = sprintf('"Normalized" Cross Correlation Image\n%d rows by %d columns', rows, columns);
title(caption, 'FontSize', fontSize);

% Find the max location of the regular correlation image and put a cross over it.
maxValue = max(crossCorrelationImage(:));
[maxRowR, maxColR] = find(crossCorrelationImage == maxValue)
% Plot the peak over the original image.
subplot(2, 3, 2);
hold on;
plot(maxColR, maxRowR, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
textSize = 14;
caption = sprintf('Max Value Location\nRow=%d, Column=%d', maxRowR, maxColR);
text(maxColR-70, maxRowR+40, caption, 'Color', 'r', 'FontSize', textSize);

% Find the max location of the normalized correlation image and put a cross over it.
maxValue = max(normalizedCrossCorrelationImage(:));
[maxRowN, maxColN] = find(normalizedCrossCorrelationImage == maxValue)
% Plot the peak over the original image.
subplot(2, 3, 5);
hold on;
plot(maxColN, maxRowN, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
textSize = 14;
caption = sprintf('Max Value Location\nRow=%d, Column=%d', maxRowN, maxColN);
text(maxColN-70, maxRowN+40, caption, 'Color', 'r', 'FontSize', textSize);

% Plot the found locations as a box over the original images.
subplot(2, 3, 3);
imshow(grayImage);
boxX2 = maxColR;
boxY2 = maxRowR;
boxX1 = boxX2 - columnsT + 1;
boxY1 = boxY2 - rowsT + 1;
boxX = [boxX1, boxX2, boxX2, boxX1, boxX1];
boxY = [boxY1, boxY1, boxY2, boxY2, boxY1];
hold on;
plot(boxX, boxY, 'r-');
title('Template "found" here.', 'FontSize', fontSize);
% Now over the normalized result.
subplot(2, 3, 6);
imshow(grayImage);
hold on;
boxX2 = maxColN;
boxY2 = maxRowN;
boxX1 = boxX2 - columnsT + 1;
boxY1 = boxY2 - rowsT + 1;
boxX = [boxX1, boxX2, boxX2, boxX1, boxX1];
boxY = [boxY1, boxY1, boxY2, boxY2, boxY1];
plot(boxX, boxY, 'r-');
title('Template found here.', 'FontSize', fontSize);


