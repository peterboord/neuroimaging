function copyDirFiles(sourceDir,targetDir,files)

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
    for fileNr = 1:numel(files)
        %filename = [subjDir{subjNr},files{fileNr}];
        filename = files{fileNr};
        copyfile(fullfile(subjSourceDir,filename),fullfile(subjTargetDir,filename));
    end
end
        