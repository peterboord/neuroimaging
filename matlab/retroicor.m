%%%%%%%%%%%%%%%%%%%%%%%
% retroicor
%%%%%%%%%%%%%%%%%%%%%%%
function [clean4d,noise4d,abMag,abPhase,slicePhases]=retroicor(restS,allCardioPeakTimes)
M=2;
tr=restS.tr;
restSize=size(restS.vol);
volsize=restS.volsize;
nvol=restS.nframes;
nslice=volsize(3);
rest=reshape(restS.vol,prod(volsize),nvol)';
rest=detrend(rest);
rest=reshape(rest',prod(volsize(1:2)),nslice,nvol);
slicePhases=peaksToSlicePhase(allCardioPeakTimes,tr,nslice,nvol);
a=zeros(prod(volsize(1:2)),nslice,M);
b=zeros(prod(volsize(1:2)),nslice,M);
noise4d=zeros(prod(volsize(1:2)),nslice,nvol);
abMag=zeros(prod(volsize(1:2)),nslice,M);
abPhase=zeros(prod(volsize(1:2)),nslice,M);
for sliceNr=1:nslice
    % calc coeffs
    for m=1:M
        a(:,sliceNr,m)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:))',cos(m*slicePhases(:,sliceNr))),1)/sum(cos(m*slicePhases(:,sliceNr)).^2);
        b(:,sliceNr,m)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:))',sin(m*slicePhases(:,sliceNr))),1)/sum(sin(m*slicePhases(:,sliceNr)).^2);
        abMag(:,sliceNr,m)=sqrt(b(:,sliceNr,m).^2+a(:,sliceNr,m).^2);
        abPhase(:,sliceNr,m)=atan2(b(:,sliceNr,m),a(:,sliceNr,m));
    end
    % calc noise estimate
    noise4d(:,sliceNr,:)=squeeze(a(:,sliceNr,:))*cos(slicePhases(:,sliceNr)*(1:M))'+squeeze(b(:,sliceNr,:))*sin(slicePhases(:,sliceNr)*(1:M))';
end
noise4d=reshape(noise4d,restSize);
% detrend
clean4d=reshape(detrend(reshape(reshape(rest,restSize)-noise4d,[],nvol)')',restSize);
end
