function estPhase

subplot(2,1,1);
nrVols = 200;
f = 6;
t = 0:0.01:nrVols;
phi = 1*pi;
%sig = [cos(2*pi*f*t);cos(2*pi*f*t + phi)];
nf = 0.6;
sig = [(1-nf)*cos(2*pi*f*t)+nf*(rand(1,numel(t))-0.5);(1-nf)*cos(2*pi*f*t + phi)+nf*(rand(1,numel(t))-0.5)];
randIdx = randi(size(sig,2),1,nrVols);
x = sig(1,randIdx);
y = sig(2,randIdx);
plot(x,y,'.')
xlim([-1,1]);
ylim([-1,1]);
subplot(2,1,2);
rose(acos(x)-acos(y),36);
end