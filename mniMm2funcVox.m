function voxCoord = mniMm2funcVox(dataDir,funcFileName,mniMm) %#ok<INUSD>

% N.B. uses 1mm MNI

MNIcoordInMmPath = fullfile(dataDir,'tmp1_mniMm2funcVox');
save(MNIcoordInMmPath,'mniMm','-ascii');
fsl(dataDir,'std2imgcoord','-img',funcFileName,'-std','/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz','-xfm','func2mni.mat','-vox',MNIcoordInMmPath,'> ','tmp2_mniMm2funcVox');
% voxel coords, origin 1
voxCoord = load(fullfile(dataDir,'tmp2_mniMm2funcVox'),'-ascii');
% take first of weird duplicate rows produced by FSL
voxCoord = round(voxCoord(1,:)) + 1;


