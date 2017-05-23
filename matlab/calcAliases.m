function [fa,N,oddity]=calcAliases(f,tr)
% % fa=abs(f-N*fs)
fs=1/tr;
N=ceil(f/fs-0.5):floor(f/fs+0.5);
fa=zeros(1,numel(N));
for n=numel(N)
switch mod(N(n),2)
    case 0
        oddity=-1;
    case 1
        oddity=1;
end
fa(n)=abs(f-N(n)*fs);
end
fa=unique(fa);
% disp(fa);
% df=1/(nvol*tr);
% fScan=0:df:(nvol-1)*df;
% faIdx=round(fa/df)+1;
% faIdx=nvol/2+1;
end