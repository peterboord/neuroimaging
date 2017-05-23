function vpfun
n=10;
vp=rand(1,27);
vp=vp/sqrt(sum(vp.^2));
ph=2*pi*(rand(1,27)-0.5);
ph=repmat((0:pi/n:2*pi)',1,27)+repmat(ph,2*n+1,1);
vpx=repmat(vp,2*n+1,1).*cos(ph);
vp2=rand(1,27);
vp2=vp2/sqrt(sum(vp2.^2));
ph2=2*pi*(rand(1,27)-0.5);
ph2=repmat((0:pi/n:2*pi)',1,27)+repmat(ph2,2*n+1,1);
vpx2=repmat(vp2,2*n+1,1).*cos(ph2);
figure,imagesc(corr(vpx,vpx2))
for y=1:21
    figure('Windowstyle','docked')
    for x=1:20
        subplot(4,5,x);
        plot(corr(vpx(y,:)',vpx')',corr(vpx2(x,:)',vpx2')');
    end
end
figure
subplot(2,1,1);plot(vpx);
subplot(2,1,2);plot(corr(vpx(1,:)',vpx')',corr(vpx(n+1,:)',vpx')')
end