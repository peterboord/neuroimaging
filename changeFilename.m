function changeFilename

dataPath = '/var/local/scratch/pboord/INDI/ABIDE/Caltech/Control';
dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

winLength = '11';
% RTI
funcFileName = 'rest_ssmooth.nii';
coordSuffix = 'Broca';
seedVoxCoordFileName = ['funcVoxCoord_',coordSuffix];
coordName = [coordSuffix,winLength];
for subjNr = 1:nrSubjs
    subjCode = subjDir{subjNr};
    subjPath = fullfile(dataPath,subjCode);
    movefile(fullfile(subjPath,'rest.feat','filtered_func_data.nii.gz'),fullfile(subjPath,'filtered_func_data.nii.gz'));
end