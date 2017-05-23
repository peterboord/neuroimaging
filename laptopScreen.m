desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktopMainFrame = desktop.getMainFrame;
% % Get desktop dimensions
desktopDims = desktopMainFrame.getSize;
desktopW = desktopDims.getWidth;
desktopH = desktopDims.getHeight;
% Resize desktop to half of original size
desktopMainFrame.setSize(1310,720);
desktopMainFrame.setLocation(0,24);