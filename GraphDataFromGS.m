function []=GraphDataFromGS(webaddr)

dbstop if error

%level=3;
level=1;
%% out put cvs file for gephi
% open the file with write permission
fid1 = fopen('table nodes.csv', 'w');
fprintf(fid1,'%s;%s;%s;%s;%s;%s\n','Id','Label','Modularity Class','Name','Photo','Citation');
%fclose(fid);
fid2=fopen('table edges.csv','w');
fprintf(fid2,'%s;%s;%s;%s;%s;%s\n','Source','Target','Type','ID','Label','Weight');
%fclose(fid);
%%
%[a1 c1]=GSWebScarb('http://scholar.google.com/citations?user=mR1GK28AAAAJ&hl=en');
[a1 c1]=GSWebScarb(webaddr);
fprintf(fid1,'%s;%s;%s;%s;%s;%s\n','1',char(cellstr(a1{1,1})),'1',char(cellstr(a1{1,1})),char(cellstr(a1{1,2})),char(cellstr(a1{1,3})));
%fclose(fid1);
k=2;
for i=1:length(c1)-1
url=char(cellstr(c1{i,1}));
[a2 c2]=GSWebScarb(url);
fprintf(fid1,'%s;%s;%s;%s;%s;%s\n',int2str(k),char(cellstr(a2{1,1})),'1',char(cellstr(a2{1,1})),char(cellstr(a2{1,2})),char(cellstr(a2{1,3})));
fprintf(fid2,'%s;%s;%s;%s;%s;%s\n',int2str(k-1),int2str(i),'Directed','1',[char(cellstr(a1{1,1})) ' To ' char(cellstr(a2{1,1}))],'1');
k=k+1;
for j=1:length(c2)-1
url=char(cellstr(c2{j,1}));
[a3 c3]=GSWebScarb(url);
fprintf(fid1,'%s;%s;%s;%s;%s;%s\n',int2str(k),char(cellstr(a3{1,1})),'1',char(cellstr(a3{1,1})),char(cellstr(a3{1,2})),char(cellstr(a3{1,3})));
fprintf(fid2,'%s;%s;%s;%s;%s;%s\n',int2str(k-1),int2str(i),'Directed','1',[char(cellstr(a2{1,1})) ' To ' char(cellstr(a3{1,1}))],'1');
k=k+1;
%url
for p=1:length(c3)-1
url=char(cellstr(c3{p,1}));
[a4 c4]=GSWebScarb(url);
fprintf(fid1,'%s;%s;%s;%s;%s;%s\n',int2str(k),char(cellstr(a4{1,1})),'1',char(cellstr(a4{1,1})),char(cellstr(a4{1,2})),char(cellstr(a4{1,3})));
fprintf(fid2,'%s;%s;%s;%s;%s;%s\n',int2str(k-1),int2str(i),'Directed','1',[char(cellstr(a3{1,1})) ' To ' char(cellstr(a4{1,1}))],'1');
k=k+1;
end
end
end
fclose(fid1);
fclose(fid2);
end