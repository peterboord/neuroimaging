function scatter(x,y,group,xName,yName)

if isdeployed, echo off
else dbstop if error; end

if nargin <2 || nargin > 5
    disp('Usage:');
    disp('scatter x y grouping [xName yName]');
    return
end

x = csvread(x);
y = csvread(y);
groupFid = fopen(group);
group = textscan(groupFid,'%s');
fclose(groupFid);
figure('Visible','on');
%gscatter(x,y,group,clr,sym,siz,doleg,xnam,ynam) 
if nargin == 5
    gscatter(x,y,group,['b','r'],['o','s'],[10,10],'on',xName,yName);
else
    gscatter(x,y,group,['b','r'],['o','s'],[10,10],'on');
end