function parrec2nii(dataDir)
dbstop if error
system(sprintf('grep -i "Protocol" %s/*.PAR > %s/ReadMe',dataDir,dataDir));
system(sprintf('grep -i "dynamics" %s/*.PAR > %s/dynamics',dataDir,dataDir));
system(sprintf('grep -i "Repetition" %s/*.PAR > %s/Repetition',dataDir,dataDir));

% options for convert_r2a.m
options.subaan= 1;
options.usealtfolder= 0;
options.altfolder= '~';
options.pathpar= [dataDir,'/'];
options.angulation= 1;
options.rescale= 1;
options.usefullprefix= 1;
options.outputformat= 1;
options.dim= 4;
options.dti_revertb0= 0;
options.prefix='';

fileProtocol = regexp(importdata(fullfile(dataDir,'ReadMe')),':','split');
parRec = cell(numel(fileProtocol),1);
niiGz = cell(numel(fileProtocol),1);
for fileIdx = 1:numel(fileProtocol)
    [~, parRec{fileIdx},ext] = fileparts(fileProtocol{fileIdx}{1});
    parRec{fileIdx} = [parRec{fileIdx},ext];
    niiGz{fileIdx} = regexprep(regexprep(fileProtocol{fileIdx}{3},'   ',''),' ','_');
    options.prefix = niiGz{fileIdx};
    convert_r2a(parRec(fileIdx),options);
end