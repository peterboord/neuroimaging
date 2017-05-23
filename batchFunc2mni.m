function batchFunc2mni(dataPath,fileNameIn,reNameIn,res,fileNameMode)

dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);
refPath = '/usr/share/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz';
%refPath = '/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz';
for subjNr = 1:nrSubjs
    disp(subjNr);
    subjPath = fullfile(dataPath,subjDir{subjNr});
    switch fileNameMode
        case 'noPrefix'
            fileName = fileNameIn;
            reName = reNameIn;
        case 'directoryPrefix'
            fileName = [subjDir{subjNr},'_',fileNameIn];
            %reName = [subjDir{subjNr},'_',reNameIn];
            reName = reNameIn;
        otherwise
            error('unknown fileNameMode');
    end
    func2mni(subjPath,fileName,reName,refPath,res);
end