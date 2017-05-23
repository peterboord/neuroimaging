function fixInf(dataPath,zFile)

dbstop if error
subjDir = dir(dataPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);

for subjNr = 1:nrSubjs
    subjCode = subjDir{subjNr};
    zImage = load_untouch_nii(fullfile(dataPath,subjCode,zFile));
    clamp = 45 + randn;
    zImage.img(zImage.img > clamp) = clamp;
    zImage.img(zImage.img < -clamp) = -clamp;
    save_untouch_nii(zImage,fullfile(dataPath,subjCode,zFile));
end