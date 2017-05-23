function extractRefs

dbstop if error

inputPDF = '/NAS_II/Home/pboord/Documents/Dropbox/Granger/2013 Ashrafulla_Canonical granger causality between regions of interest.pdf'; 
outputfile = '/NAS_II/Home/pboord/Documents/Dropbox/Granger/2013 Ashrafulla_Canonical granger causality between regions of interest.txt'; 
cmd = ['pdftotext -raw ''',inputPDF,''' ''',outputfile,''''];
system(cmd); 
fid = fopen(outputfile); 
alltext = textscan(fid,'%s','Delimiter','\n'); 
fclose(fid);

% CM = ones(7);
% bg = biograph(CM);
% set(bg,'NodeCallback','pbNodeCallback');
% view(bg)
%% Find start of References
refIdx = [find(strncmp('REFERENCES',alltext{1},numel('REFERENCES'))),find(strncmp('References',alltext{1},numel('References')))];
if numel(refIdx) ~=1
    error('refIdx');
end
refRows = alltext{1}((refIdx+1):end);
year=regexp(refRows,'(19|20)\d\d','match');
for row = 1:numel(refRows)
    regexp
end