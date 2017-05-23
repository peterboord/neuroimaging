function func2mni(dataDir,fileName,reName,refPath,res)

% [~,name] = fileparts(fileName);
% newFile = [name,'_MNIreg.nii'];
switch res
    case 4
        fsl(dataDir,'flirt','-in',fileName,'-ref',refPath,'-applyisoxfm 4','-init','func2mni.mat','-out',reName);
    case 3
        fsl(dataDir,'flirt','-in',fileName,'-ref',refPath,'-applyisoxfm 3','-init','func2mni.mat','-out',reName);
    case 1
        fsl(dataDir,'flirt','-in',fileName,'-ref',refPath,'-applyxfm','-init','func2mni.mat','-out',reName);
    otherwise
        error('resolution not supported');
end
system(['gunzip -f ',fullfile(dataDir,[reName,'.nii.gz'])]);
end
