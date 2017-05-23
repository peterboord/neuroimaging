nf=50;
tr=1;
n=0.1;
x=cos(2*pi*0.05*(1:nf)'*tr);
x=x/sqrt(sum(x.^2));
y1=cos(2*pi*0.05*(1:nf)'*tr+pi/6)+n*randn(nf,1);
y2=cos(2*pi*0.05*(1:nf)'*tr+pi*8/18)+n*randn(nf,1);
y=[y1/sqrt(sum(y1.^2)),y2/sqrt(sum(y2.^2))];
w=corr(x,y);
w=w/sqrt(sum(w.^2));
wy=y*w';
wy=wy/sqrt(sum(wy.^2));
figure
%subplot(2,1,1);
plot(y1,y2)
%plot([x,wy,y])
% subplot(2,1,2);
% plot([x,t])

% w[y=mean(y.*(x>=0)-y.*(x<0));
% wz=mean(z.*(x>=0)-z.*(x<0));
