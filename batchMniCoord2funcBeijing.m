function batchMniCoord2funcBeijing(dataPath,MNIcoordInMmPath,coordName,funcFileName)

dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);
for subjNr = 1:nrSubjs
    disp(subjNr);
    subjPath = fullfile(dataPath,subjDir{subjNr});
    %funcFileName = 'rest_b0_mcf.nii';
    mniCoord2func(subjPath,funcFileName,MNIcoordInMmPath,coordName)
end