function makeFuncVoxNrs3d(funcDir,func_mean)
vS=MRIread(fullfile(funcDir,func_mean));
MRIsave(vS,reshape(find(ones(vS.volsize)),vS.volsize),fullfile(funcDir,'funcVoxNrs3d.nii.gz'),1);
