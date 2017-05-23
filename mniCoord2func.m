function mniCoord2func(dataDir,funcFileName,MNIcoordInMmPath,coordName)

fsl(dataDir,'std2imgcoord','-img',funcFileName,'-std','/usr/share/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz','-xfm','func2mni.mat','-vox',MNIcoordInMmPath,'> ',['funcVoxCoord_',coordName]);
