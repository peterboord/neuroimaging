dataPath = '/NAS_II/archive_project_space/PIC';
dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

changeList = {'picMag2_PCC11.nii','picCorrRDM2_PCC11.nii','picMeanOnly2_PCC11.nii','picMeanAndRes2_PCC11.nii','pic2_PCC11.nii','absFcSw2_PCC11.nii'};
coordSuffix = 'PCC';
seedVoxCoordFileName = ['funcVoxCoord_',coordSuffix];
winLength = '11';
coordName = [coordSuffix,winLength];
for subjNr = 1:nrSubjs
    subjCode = subjDir{subjNr};
    fileDir = dir(fullfile(dataPath,subjCode));
    fileDir(1:2) = [];
    fileDir = {fileDir(~[fileDir.isdir]').name}';
    nrFiles = numel(fileDir);
    for fileNr = 1:nrFiles
        if ismember(fileDir{fileNr},changeList)
            movefile(fullfile(dataPath,subjCode,fileDir{fileNr}),fullfile(dataPath,subjCode,regexprep(fileDir{fileNr},'PCC11','PCC15')));
        end
    end
end