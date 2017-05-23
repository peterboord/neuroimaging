function [author coauthor]= GSWebScarb(address)
dbstop if error
%address='http://scholar.google.com/citations?user=8jAO-t4AAAAJ&am';
html = urlread(address,'UserAgent','Mozilla 5.0');
%divopen='<div style="font-weight: bold">&nbsp;Co-authors</div>';
divopen='</a><br><span class="cit-gray">';
pstart=strfind(html,divopen);
%divend='">View all co-authors</a>';
divend='</span><br><span class="cit-gray">';
pend=strfind(html,divend);
str=html(pstart:pend);
coa='<a';
lstart=strfind(str,'href="');
tl=';hl=en" title="';
tstart=strfind(str,tl);
author=cell(1,3);
coauthor=cell(length(lstart),2);
%% Find Photo, h-index,i-index ,All citiations
tn1='<title>';
tn2= '- Google Scholar Citations</title>';
a_name1=strfind(html(1:500),tn1);
a_name2=strfind(html(1:500),tn2);
author{1,1}={html((a_name1+length(tn1)):(a_name2-2))};
%photo link lets robot-bot google too
usrInd=strfind(address,'user=');
usr=address(usrInd:end);
link=['http://scholar.google.com/citations?view_op=view_photo&' usr '&citpid=1'];
author{1,2}={link};
%h-index ,i-index
cit='<td>';
cit_2='</td';
citIdx=strfind(html,cit);
cit_str=html((citIdx(1)+length(cit)):citIdx(2));
cit=strfind(cit_str,cit_2);
author{1,3}={cit_str(1:cit-1)};
%%
if(isempty(lstart))
disp('No co-authors or Web-address not correct, please check');
else
for i=1:length(lstart)-1
%link=['http://scholar.google.com' str((lstart(i)+6) <img src="http://s0.wp.com/wp-includes/images/smilies/icon_sad.gif?m=1129645325g" alt=":(" class="wp-smiley"> tstart(i)-2))];
%coauthor(i,1)=[link];
%find fullname
nm1='">';
nm2='</a><br>';
substr=str(tstart(i):end);
namestr2=strfind(substr,nm2);
namestr1=strfind(substr,nm1);
substr2=substr(length(tl)+1:namestr1(1)-1);
coauthor{i,1}={link};coauthor{i,2}={substr2};
end
end
end