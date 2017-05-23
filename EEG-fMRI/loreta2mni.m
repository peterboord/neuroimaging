function loreta2mni

% load loreta MNI voxel coords
% from BrainProducts Analyzer/Transformations/Special Signal
% Processing/LORETA/Export Blank ROIs file
load('loretaMniMm');
% load 1mm MNI and zero
mniS = MRIread('/usr/share/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz');
mniS.vol = zeros(mniS.volsize);
% convert MNI mm to vox
for idx = 1:size(loretaMniMm,1)
    ras = mniS.vox2ras\[loretaMniMm(idx,:),1]';
    mniS.vol(ras(2),ras(1),ras(3)) = idx; % freesurfer format is ars (not ras)
end
mniS.fspec = '/NAS_II/Home/pboord/Documents/Scripts/loretaMni.nii.gz';
MRIwrite(mniS,mniS.fspec);
% insert index of voxels into image

% save image

end
