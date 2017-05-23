function pbDcm2nii(dataDir)

dbstop if error
subjDir = dir(dataDir);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

for subjNr = 1:nrSubjs
    system(['dcm2nii ',fullfile(dataDir,subjDir{subjNr})]);
end

