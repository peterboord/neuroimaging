function ranksumGroups(group1Path,group2Path,fileName)

resultsDir = '/NAS_II/Projects/pboord/PIC/Results/Udall';
dbstop if error
[group1images,floatFile] = getGroupImages(group1Path,fileName);
group2images = getGroupImages(group2Path,fileName);
ranksumFile = load_untouch_nii(fullfile(resultsDir,'standard.nii'));
diffFile = floatFile;
probFile = floatFile;
imSize = size(ranksumFile.img);
for x=1:imSize(1)
    disp(x);
    for y=1:imSize(2)
        for z=1:imSize(3)
            if ranksumFile.img(x,y,z)
                group1xyz = squeeze(group1images(x,y,z,:));
                group2xyz = squeeze(group2images(x,y,z,:));
                [probFile.img(x,y,z),~,stats] = ranksum(group1xyz,group2xyz,'method','approximate');
                ranksumFile.img(x,y,z) = stats.ranksum;
                diffFile.img(x,y,z) = mean(group1xyz) - mean(group2xyz);
            end
        end
    end
end
save_untouch_nii(probFile,fullfile(resultsDir,['prob_',fileName]));
save_untouch_nii(ranksumFile,fullfile(resultsDir,['ranksum_',fileName]));
save_untouch_nii(diffFile,fullfile(resultsDir,['diff_',fileName]));
end

function [groupImages,imFile] = getGroupImages(groupPath,fileName)
subjDir = dir(groupPath);
subjDir(1:2) = [];
subjDir = {subjDir([subjDir.isdir]').name}';
nrSubjs = numel(subjDir);
groupImages = [];
for subjNr = 1:nrSubjs
    subjPath = fullfile(groupPath,subjDir{subjNr});
    imFile = load_untouch_nii(fullfile(subjPath,fileName));
    if isempty(groupImages)
        groupImages = zeros([size(imFile.img),nrSubjs]);
    end
    groupImages(:,:,:,subjNr) = imFile.img; %#ok<AGROW>
end
end
