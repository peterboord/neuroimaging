function batchMniCoord2func(dataPath,MNIcoordInMmPath,coordName)

dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);
for subjNr = 1:nrSubjs
    disp(subjNr);
    subjPath = fullfile(dataPath,subjDir{subjNr});
    subjCode = regexp(subjDir{subjNr},'_','split');
    funcFileName = [subjCode{1},'_reg_filtered_func_data.nii'];
    mniCoord2func(subjPath,funcFileName,MNIcoordInMmPath,coordName)
end