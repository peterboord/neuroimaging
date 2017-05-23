function batchMni2func(dataPath,T1FileName,funcFileName,targetDir)

% T1FileName and funcFileName are names *without* extensions

dbstop if error

subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

switch nargin
    case 3
        for subjNr = 1:nrSubjs
            subjPath = fullfile(dataPath,subjDir{subjNr});
            %     subjCode = regexp(subjDir{subjNr},'_','split');
            %     T1FileName = [subjCode{1},'_avgT1'];
            %     funcFileName = [subjCode{1},'_reg_filtered_func_data'];
            %     T1FileName = 'T1';
            %     funcFileName = 'rest_b0_mcf';
            mni2func(subjPath,T1FileName,funcFileName);
        end
    case 4
        
        for subjNr = 1:nrSubjs
            disp(subjNr);
            subjSourceDir = fullfile(dataPath,subjDir{subjNr});
            subjTargetDir = fullfile(targetDir,subjDir{subjNr});
            if ~exist(subjTargetDir,'dir')
                mkdir(subjTargetDir);
            end
            sessDir = dir(subjSourceDir);
            sessDir(1:2) = [];
            sessDir = {sessDir([sessDir.isdir]').name}';
            nrSess = numel(sessDir);
            for sessNr = 1:nrSess
                tic
                disp([subjNr,sessNr]);
                sessSourceDir = fullfile(dataPath,subjDir{subjNr},sessDir{sessNr});
                sessTargetDir = fullfile(targetDir,subjDir{subjNr},sessDir{sessNr});
                mni2func(sessSourceDir,T1FileName,funcFileName,sessTargetDir);
                toc
            end
        end
    otherwise
        error('Incorrect number of arguments');
end