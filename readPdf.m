
inputPDF = '/NAS_II/Home/pboord/Documents/Dropbox/Granger/2013 Ashrafulla_Canonical granger causality between regions of interest.pdf'; 
outputfile = '/NAS_II/Home/pboord/Documents/Dropbox/Granger/2013 Ashrafulla_Canonical granger causality between regions of interest.txt'; 
cmd = ['pdftotext -raw ''',inputPDF,''' ''',outputfile,''''];
system(cmd); 
fid = fopen(outputfile); 
alltext = textscan(fid,'%s','Delimiter','\n'); 
fclose(fid);