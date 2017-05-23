function batchMniCoord2funcSessions(dataPath,MNIcoordInMmPath,coordName,funcFileName)

dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);
for subjNr = 1:nrSubjs
    subjSourceDir = fullfile(dataPath,subjDir{subjNr});
    sessDir = dir(subjSourceDir);
    sessDir(1:2) = [];
    sessDir = {sessDir([sessDir.isdir]').name}';
    nrSess = numel(sessDir);
    for sessNr = 1:nrSess
        sessSourceDir = fullfile(subjSourceDir,sessDir{sessNr});
        if exist(fullfile(sessSourceDir,[funcFileName,'.nii']),'file')
            mniCoord2func(sessSourceDir,funcFileName,MNIcoordInMmPath,coordName)
        else
            disp([subjDir{subjNr},' ',sessDir{sessNr},' empty']);
        end
    end
end