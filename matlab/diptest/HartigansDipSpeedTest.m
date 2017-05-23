function HartigansDipSpeedTest
%create some obviously unimodal and bimodal Gaussian distributions just to
%see what dip statistic does
% Nic Price 2006

nrRepeats = 60000;
Cent1 = ones(1,9);
Cent2 = 1:1:9;
sig = 0.5;
nboot = 500;
nrPoints = 300;
a = 1;
xpdf = repmat(sort([Cent1(a)+randn(1,floor(nrPoints/2)) Cent2(a)+randn(1,floor(nrPoints/2))]),nrRepeats,1);
tic
for a = 1:nrRepeats
    d = HartigansDipTest(xpdf(a,:));
end
% FixAxis([-2 12]);
disp(toc/nrRepeats);