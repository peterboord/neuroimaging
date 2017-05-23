function segFunc(funcPath,structPath,structName)
dbstop if error

% extract T1 brain with BET
if ~exist(fullfile(structPath,'T1_brain.nii.gz'),'file')
    fsl(structPath,'bet',structName,fullfile(structPath,'T1_brain'),'-m');
end
% segment T1 with FAST
if ~exist(fullfile(structPath,'T1_brain_seg.nii.gz'),'file')
    fsl(structPath,'fast','-g','T1_brain');
end
% register func to T1 with FLIRT
if ~exist(fullfile(funcPath,'func2T1.mat'),'file')
    fsl(funcPath,'flirt','-in','mean_func','-ref',[structPath,'/T1_brain'],'-out','mean_func_reg','-omat','func2T1.mat','-searchrx -180 180 -searchry -180 180 -searchrz -180 180');
end
% xfm GM, WM, CSF to FUNC space
if ~exist(fullfile(funcPath,'mean_funcGM.nii.gz'),'file')
    % invert xfm func2T1 to T12func
    fsl(funcPath,'convert_xfm','-omat','T12func.mat','-inverse','func2T1.mat');
    % xfm GM to mean_func
    fsl(funcPath,'flirt','-in',[structPath,'/T1_brain_seg_0'],'-ref','mean_func','-out','mean_funcGM','-applyxfm','-init','T12func.mat');
    % xfm WM to mean_func
    fsl(funcPath,'flirt','-in',[structPath,'/T1_brain_seg_1'],'-ref','mean_func','-out','mean_funcWM','-applyxfm','-init','T12func.mat');
    % xfm CSF to mean_func
    fsl(funcPath,'flirt','-in',[structPath,'/T1_brain_seg_2'],'-ref','mean_func','-out','mean_funcCSF','-applyxfm','-init','T12func.mat');
end
% unzip func GM, WM, CSF
if ~exist(fullfile(funcPath,'mean_funcGm_untouch'),'file')
    system(['gunzip -c ',fullfile(funcPath,'mean_funcGM.nii.gz'),' > ',fullfile(funcPath,'mean_funcGM_untouch.nii')]);
    system(['gunzip -c ',fullfile(funcPath,'mean_funcWM.nii.gz'),' > ',fullfile(funcPath,'mean_funcWM_untouch.nii')]);
    system(['gunzip -c ',fullfile(funcPath,'mean_funcCSF.nii.gz'),' > ',fullfile(funcPath,'mean_funcCSF_untouch.nii')]);
end
% load images
mean_funcGM = load_untouch_nii(fullfile(funcPath,'mean_funcGM_untouch.nii'));
mean_funcWM = load_untouch_nii(fullfile(funcPath,'mean_funcWM_untouch.nii'));
mean_funcCSF = load_untouch_nii(fullfile(funcPath,'mean_funcCSF_untouch.nii'));
% threshold GM
mean_funcGM.img(mean_funcGM.img ~= 0 & mean_funcGM.img >= mean_funcWM.img & mean_funcGM.img >= mean_funcCSF.img) = 1;
mean_funcGM.img(mean_funcGM.img ~= 1) = 0;
% threshold WM
mean_funcWM.img(mean_funcWM.img ~= 0 & mean_funcWM.img > mean_funcGM.img & mean_funcWM.img >= mean_funcCSF.img) = 1;
mean_funcWM.img(mean_funcWM.img ~= 1 | mean_funcGM.img == 1) = 0;
% threshold CSF
mean_funcCSF.img(mean_funcCSF.img ~= 0 & mean_funcCSF.img > mean_funcGM.img & mean_funcCSF.img > mean_funcWM.img) = 1;
mean_funcCSF.img(mean_funcCSF.img ~= 1 | mean_funcGM.img == 1 | mean_funcWM.img == 1) = 0;
% save images
save_untouch_nii(mean_funcGM,fullfile(funcPath,'mean_funcGM_untouch.nii'));
save_untouch_nii(mean_funcWM,fullfile(funcPath,'mean_funcWM_untouch.nii'));
save_untouch_nii(mean_funcCSF,fullfile(funcPath,'mean_funcCSF_untouch.nii'));
