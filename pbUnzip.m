function pbUnzip(dataDir)

dbstop if error
subjDir = dir(dataDir);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

for subjNr = 1:nrSubjs
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'*ep2*.nii.gz'));
    if numel(structFiles) > 0
        system(['gunzip ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
    
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'co*.nii.gz'));
    if numel(structFiles) > 0
        system(['gunzip ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
    
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'o*.nii.gz'));
    if numel(structFiles) > 0
        system(['rm ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
    
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'*t1mprnssags*.nii.gz'));
    if numel(structFiles) > 0
        system(['rm ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
    
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'T1_brain.nii.gz'));
    if numel(structFiles) > 0
        system(['gunzip ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
    
    structFiles = dir(fullfile(dataDir,subjDir{subjNr},'T1_brain_mask.nii.gz'));
    if numel(structFiles) > 0
        system(['gunzip ',fullfile(dataDir,subjDir{subjNr},structFiles(1).name)]);
    end
end

