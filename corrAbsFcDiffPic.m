function corrAbsFcDiffPic(fcDiffPath,picPath)
addpath('/project_space/pboord/usr/bin')
fcDiffIm=MRIread(fcDiffPath);
picIm=MRIread(picPath);
% corr abs fc with neg pic
MRIsave(fcDiffIm,reshape(corrColumns(reshape(abs(fcDiffIm.vol),[],fcDiffIm.nframes)',reshape(-picIm.vol,[],picIm.nframes)'),fcDiffIm.volsize),'corrAbsFcDiffPic.nii.gz',1);
