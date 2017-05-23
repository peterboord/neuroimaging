function fsl(sessionDir,fnName,varargin)

% setenv('LD_PRELOAD','/usr/local/MATLAB/R2014a/sys/os/glnxa64/libstdc++.so.6');
% fsl_path = '/usr/share/fsl/5.0';
% setenv('FSLDIR',fsl_path)
% setenv('FSLOUTPUTTYPE','NIFTI_GZ')
% curpath = getenv('PATH');
% if ~strncmp(fullfile(fsl_path,'bin:'),curpath,numel(fsl_path)+1)
%     setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));
% end

cmdString = ['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/',fnName,' '];
for i=1:numel(varargin)
    strTrimVarargin = strtrim(varargin{i});
    if ismember(strTrimVarargin(1),{'/','-','>','"','~','0','1','2','3','4','5','6','7','8','9'})
        cmdString = [cmdString,varargin{i},' '];
    else
        cmdString = [cmdString,fullfile(sessionDir,varargin{i}),' ']; %#ok<*AGROW>
    end
end
cmdString = [cmdString,'"'];
[status,cmdout] = system(cmdString);
if status
    error(cmdout);
end
end
