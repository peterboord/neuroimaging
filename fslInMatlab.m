% Must launch Matlab with:
% LD_PRELOAD=/usr/lib64/libstdc++.so.6 /usr/local/Matlab/bin/matlab
% refs: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0902&L=FSL&D=0&P=348826
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0902&L=FSL&D=0&P=357648
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/162466
% these docs saved in pboord/Documents/FSL
fsl_path = '/usr/share/fsl/4.1/';
setenv('FSLDIR',fsl_path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));
tmp=sprintf('sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/fast -b --nopve /var/local/scratch/pboord/BreathHold_TEcheck/TE35_2_EPI_3mmISO_SENSE/TE35_2_EPI_3mmISO_SENSE/TE35_2_EPI_3mmISO_SENSE-d0130.feat/mean_func"');
system(tmp);
