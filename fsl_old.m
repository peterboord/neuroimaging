function fsl_old(sessionDir,fnName,varargin)

fsl_path = '/usr/share/fsl/4.1/';
setenv('FSLDIR',fsl_path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));

cmdString = ['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/',fnName,' '];
for i=1:numel(varargin)
    strTrimVarargin = strtrim(varargin{i});
    if ismember(strTrimVarargin(1),{'/','-','>','"'})
        cmdString = [cmdString,varargin{i},' '];
    else
        cmdString = [cmdString,fullfile(sessionDir,varargin{i}),' ']; %#ok<*AGROW>
    end
end
cmdString = [cmdString,'"'];
system(cmdString);

