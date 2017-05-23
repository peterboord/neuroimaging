function mni2func(dataDir,T1FileName,funcFileName,targetDir)

switch nargin
    case 3
        %dataDir = ''; % name of directory where files are
        % T1path must be prepended with'/' for fsl()
        % mean func
     %   fsl(dataDir,'fslmaths',funcFileName,'-Tmean',['mean_',funcFileName]);
        % extract func brain with BET
        fsl(dataDir,'bet',['mean_',funcFileName],['mean_',funcFileName],'-m');
        % mask func
        fsl(dataDir,'fslmaths',funcFileName,'-mas',['mean_',funcFileName,'_mask'],funcFileName);
        if ~strcmp(T1FileName,'T1_brain')
            % extract T1 brain with BET
            fsl(dataDir,'bet',T1FileName,'T1_brain','-m');
        end
        % reg T1 to MNI
        fsl(dataDir,'flirt','-in','T1_brain','-ref','/usr/share/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz','-out','T1_brain_mni','-omat','T12mni.mat','-searchrx -180 180 -searchry -180 180 -searchrz -180 180');
        %fsl(dataDir,'flirt','-in','T1_brain','-ref','/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz','-out','T1_brain_mni','-omat','T12mni.mat','-searchrx -180 180 -searchry -180 180 -searchrz -180 180');
        %reg func to T1
        fsl(dataDir,'flirt','-in',['mean_',funcFileName],'-ref','T1_brain','-out','func_T1reg','-omat','func2T1.mat','-searchrx -180 180 -searchry -180 180 -searchrz -180 180');
        % concat for FUNC to MNI
        % convert_xfm -omat <outmat_AtoC> -concat <mat_BtoC> <mat_AtoB>
        fsl(dataDir,'convert_xfm','-omat','func2mni.mat','-concat','T12mni.mat','func2T1.mat');
        % invert xfm func2mni to mni2func
        fsl(dataDir,'convert_xfm','-omat','mni2func.mat','-inverse','func2mni.mat');
        % invert xfm func2T1 to T12func
        fsl(dataDir,'convert_xfm','-omat','T12func.mat','-inverse','func2T1.mat');
        % invert xfm T12mni to T12mni
        fsl(dataDir,'convert_xfm','-omat','T12mni.mat','-inverse','T12mni.mat');
        % invert xfm T12mni to mni2T1
        fsl(dataDir,'convert_xfm','-omat','mni2T1.mat','-inverse','T12mni.mat');
    case 4
        system(['gunzip -f -c ',fullfile(dataDir,[T1FileName,'.nii.gz']),' > ',fullfile(targetDir,[T1FileName,'.nii'])]);
        system(['gunzip -f -c ',fullfile(dataDir,[funcFileName,'.nii.gz']),' > ',fullfile(targetDir,[funcFileName,'.nii'])]);
%         % reorient to MNI
%         fsl(targetDir,'fslreorient2std',fullfile(targetDir,[funcFileName,'.nii']),fullfile(targetDir,[funcFileName,'.nii']));
        % mean func
        fsl(targetDir,'fslmaths',funcFileName,'-Tmean',['mean_',funcFileName]);
        % extract func brain with BET
        fsl(targetDir,'bet',['mean_',funcFileName],['mean_',funcFileName],'-m');
        % mask func
        fsl(targetDir,'fslmaths',funcFileName,'-mas',['mean_',funcFileName,'_mask'],funcFileName);
        system(['gunzip -f ',fullfile(targetDir,['mean_',funcFileName,'_mask.nii.gz'])]);
        if ~strcmp(T1FileName,'T1_brain')
            % extract T1 brain with BET
            fsl(targetDir,'bet',T1FileName,'T1_brain','-m');
        end
        % reg T1 to MNI
        fsl(targetDir,'flirt','-in','T1_brain','-ref','/usr/share/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz','-out','T1_brain_mni','-omat','T12mni.mat','-searchrx -90 90 -searchry -90 90 -searchrz -90 90');
        %fsl(targetDir,'flirt','-in','T1_brain','-ref','/usr/share/data/fsl-mni152-templates/MNI152_T1_1mm_brain.nii.gz','-out','T1_brain_mni','-omat','T12mni.mat','-searchrx -90 90 -searchry -90 90 -searchrz -90 90');
        %reg func to T1
        fsl(targetDir,'flirt','-in',['mean_',funcFileName],'-ref','T1_brain','-out','func_T1reg','-omat','func2T1.mat','-searchrx -90 90 -searchry -90 90 -searchrz -90 90');
        % concat for FUNC to MNI
        % convert_xfm -omat <outmat_AtoC> -concat <mat_BtoC> <mat_AtoB>
        fsl(targetDir,'convert_xfm','-omat','func2mni.mat','-concat','T12mni.mat','func2T1.mat');
        % invert xfm func2mni to mni2func
        fsl(targetDir,'convert_xfm','-omat','mni2func.mat','-inverse','func2mni.mat');
    otherwise
        error('Incorrect number of arguments');
end