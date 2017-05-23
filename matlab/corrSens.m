function corrSens(a,b)
x=a;y=b;
n=20;
s=zeros(n,5);
for i=1:n
    fc=corr(mean(x,2),mean(y,2));
    [xVpHi,xVpLo]=getVp(x);
    [yVpHi,yVpLo]=getVp(y);
    corrVpHiTs=corr(corr(xVpHi',x')',corr(yVpHi',y')');
    corrVpHiLoTsX=corr(corr(xVpHi',x')',corr(xVpLo',x')');
    corrVpHiLoTsY=corr(corr(yVpHi',y')',corr(yVpLo',y')');
    corrxy=corr(x,y);
    if rem(i,2)==1
        [r,ri]=sort(corrxy(:),'descend');
    else
        [r,ri]=sort(corrxy(:),'ascend');
    end    
    [xi,yi]=ind2sub([size(x,2),size(x,2)],ri(1));
    x(:,xi)=[];y(:,yi)=[];
    s(i,1)=fc;
    s(i,2)=r(1);
    s(i,3)=corrVpHiTs;
    s(i,4)=corrVpHiLoTsX;
    s(i,5)=corrVpHiLoTsY;
end
figure,plotyy(1:n,s(:,1:3),1:n,s(:,4:5)),legend('fc','rMax','corrVpTs','corrVpHiLoTsX','corrVpHiLoTsY');
end

function [hiVp,loVp]=getVp(x)
nf=size(x,1);
d=2;
fcTs=mean(x,2);
[~,idx]=sort(fcTs);
idxsLo=idx(1:floor(nf/d));
idxsHi=idx(end-floor(nf/d)+1:end);
anorm=x./repmat(sqrt(sum(x.^2,2)),1,size(x,2));
hiVp=sum(anorm(idxsHi,:).*repmat(fcTs(idxsHi),1,size(x,2)),1)/sum(fcTs(idxsHi));
loVp=sum(anorm(idxsLo,:).*repmat(fcTs(idxsLo),1,size(x,2)),1)/sum(fcTs(idxsLo));
hiVp=hiVp/sqrt(sum(hiVp.^2));
loVp=loVp/sqrt(sum(loVp.^2));
end