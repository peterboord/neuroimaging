function pbBet(dataDir)

dbstop if error
subjDir = dir(dataDir);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

for subjNr = 3:nrSubjs
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'co*.nii.gz'));
    fsl(fullfile(dataDir,subjDir{subjNr}),'bet',fullfile(dataDir,subjDir{subjNr},structFiles(1).name),'T1_brain','-m');
end

