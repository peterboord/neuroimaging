function picStatDsn
nrSim=10000;
nrTr=600;
nrNeigh=27;
picAngle=zeros(nrSim,2);
tic
for simNr=1:nrSim
    x=randn(nrNeigh,nrTr);
    x=x./repmat(sqrt(sum(x.^2)),nrNeigh,1);
    [eigX,~,varExp]=fastEig(x);
    rawAngle=acosd(sum(eigX)/sqrt(nrNeigh));
    picAngle(simNr,1)=abs(180*(rawAngle>90)-rawAngle);
    picAngle(simNr,2)=varExp;
end
toc
figure
subplot(1,3,1);plot(sort(picAngle(:,1)));
subplot(1,3,2);hist(picAngle(:,1));
subplot(1,3,3);hist(picAngle(:,2));
end