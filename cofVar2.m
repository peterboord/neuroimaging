function cofVar2(dataDir,filename,funcName)

switch nargin
    case 1
        fsl(dataDir,'fslmaths','filtered_func_data.nii.gz','-Tmean','mean_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths','filtered_func_data.nii.gz','-Tstd','std_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths','filtered_func_data.nii.gz','-bptf 4 -1','hf_filtered_func_data.nii.gz');
    case 2
        fsl(dataDir,'fslmaths',filename,'-Tmean','mean_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths',filename,'-Tstd','std_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths',filename,'-bptf 4 -1','hf_filtered_func_data.nii.gz');
    case 3
        fsl(dataDir,'fslmaths',filename,'-Tmean',['mean_',funcName]);
        fsl(dataDir,'fslmaths',filename,'-Tstd',['std_',funcName]);
        fsl(dataDir,'fslmaths',filename,'-bptf 4 -1',['hf_',funcName]);
end
switch nargin
    case {1,2}
        fsl(dataDir,'fslmaths','std_filtered_func_data.nii.gz','-div','mean_filtered_func_data.nii.gz','cv_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths','hf_filtered_func_data.nii.gz','-Tmean','mean_hf_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths','hf_filtered_func_data.nii.gz','-Tstd','std_hf_filtered_func_data.nii.gz');
        fsl(dataDir,'fslmaths','std_hf_filtered_func_data.nii.gz','-div','mean_hf_filtered_func_data.nii.gz','cv_hf_filtered_func_data.nii.gz');
    case 3
        fsl(dataDir,'fslmaths',['std_',funcName],'-div',['mean_',funcName],['cv_',funcName]);
        fsl(dataDir,'fslmaths',['hf_',funcName],'-Tmean',['mean_hf_',funcName]);
        fsl(dataDir,'fslmaths',['hf_',funcName],'-Tstd',['std_hf_',funcName]);
        fsl(dataDir,'fslmaths',['std_hf_',funcName],'-div',['mean_hf_',funcName],['cv_hf_',funcName]);
end
