function copySessionDirFiles(sourceDir,targetDir,files)

%sourceDir = '/NAS_II/archive_project_space/pboord/PIC';
dbstop if error
subjDir = dir(sourceDir);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

for subjNr = 1:nrSubjs
    subjSourceDir = fullfile(sourceDir,subjDir{subjNr});
    subjTargetDir = fullfile(targetDir,subjDir{subjNr});
    if ~exist(subjTargetDir,'dir')
        mkdir(subjTargetDir);
    end
    sessDir = dir(subjSourceDir);
    sessDir(1:2) = [];
    sessDir = {sessDir([sessDir.isdir]').name}';
    nrSess = numel(sessDir);
    for sessNr = 1:nrSess
        sessSourceDir = fullfile(sourceDir,subjDir{subjNr},sessDir{sessNr});
        sessTargetDir = fullfile(targetDir,subjDir{subjNr},sessDir{sessNr});
        if ~exist(sessTargetDir,'dir')
            mkdir(sessTargetDir);
        end
        for fileNr = 1:numel(files)
            filename = files{fileNr};
            sessTargetFilePath = fullfile(sessTargetDir,filename);
            if exist(fullfile(sessSourceDir,filename),'file') && ~exist(sessTargetFilePath,'file')
                copyfile(fullfile(sessSourceDir,filename),sessTargetFilePath);
            end
        end
    end
end
