function [maxCardioPkFreq,Pxx,f]=findPeakFreq(rmsSlice,nrTr,padFactor,minCardioHz,maxCardioHz,tr,restS)
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
% find peak in cardio range
[Pxx,f]=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
dfSlice=1/(nrTr*tr);
minCardioIdx=floor(minCardioHz*padFactor/dfSlice);
maxCardioIdx=ceil(maxCardioHz*padFactor/dfSlice)+1;
cardioPeakIdxs=[0;Pxx(2:end-1)>Pxx(1:end-2) & Pxx(2:end-1)>Pxx(3:end);0];
cardioPeakIdxs(1:minCardioIdx-1)=0;
cardioPeakIdxs(maxCardioIdx+1:end)=0;
[~,maxCardioPkIdx]=max(Pxx.*cardioPeakIdxs);
% if no peak in cardio range use max boundary
if maxCardioPkIdx==1
    if Pxx(minCardioIdx)>=Pxx(maxCardioIdx)
        maxCardioPkIdx=minCardioIdx;
    else
        maxCardioPkIdx=maxCardioIdx;
    end
end
maxCardioPkFreq=(maxCardioPkIdx-1)*dfSlice/padFactor;
[aliasFreq,N,oddity]=calcAliases(maxCardioPkFreq,tr);
disp({'N','oddity'});
disp([N,oddity]);
rest=reshape(restS.vol,[],nslice,nvol);
t=reshape(0:tr/nslice:(nvol*nslice-1)*tr/nslice,[],nvol);
cosPhase=cos(2*pi*aliasFreq*t);
sinPhase=sin(2*pi*aliasFreq*t);
restInPhase=zeros(prod(volsize(1:2)),nslice);
restInQuad=zeros(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
    restInPhase(:,sliceNr)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:)),cosPhase(sliceNr,:)),2)/sum(cosPhase(sliceNr,:).^2);
    restInQuad(:,sliceNr)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:)),sinPhase(sliceNr,:)),2)/sum(sinPhase(sliceNr,:).^2);
end
restMag=sqrt(restInPhase.^2+restInQuad.^2);
restPhase=oddity*atan2(restInQuad,restInPhase);
nrPhaseBins=180;
phaseBins=(0:2*pi/nrPhaseBins:(nrPhaseBins-1)*2*pi/nrPhaseBins)-pi;
phaseHist=zeros(nslice,nrPhaseBins);
maxMag=max(restMag(:));
nrMagBins=100;
magBins=0:maxMag/nrMagBins:(nrMagBins-1)*maxMag/nrMagBins;
magHist=zeros(nslice,nrMagBins);
for sliceNr=1:nslice
    phaseHist(sliceNr,:)=hist(restPhase(restMag(:,sliceNr)>0,sliceNr),phaseBins);
    magHist(sliceNr,:)=hist(restMag(restMag(:,sliceNr)>0,sliceNr),magBins);
end
figure('WindowStyle','docked');
surf(phaseHist)
% figure
% surf(magHist)
% MRIsave(restS,reshape(restMag,restS.volsize),fullfile(fileparts(restS.fspec),'restMag.nii'),1);
% MRIsave(restS,reshape(restPhase,restS.volsize),fullfile(fileparts(restS.fspec),'restPhase.nii'),1);
% roi=[22,41,30];
% roiInd=sub2ind(restS.volsize,roi(1),roi(2),roi(3));
% fftInPhase=fft(restInPhase(roiInd,:)');
% figure
% f=0:1/tr:(nvol-1)/tr;
% fRange=1:nvol/2+1;
% plot(f(fRange),abs(fftInPhase(fRange)))

end