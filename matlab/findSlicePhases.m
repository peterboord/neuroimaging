function [maxCardioPkFreq,Pxx,f]=findSlicePhases(rmsSlice,nrTr,padFactor,minCardioHz,maxCardioHz,tr,restS,slicePhasesPulseox,cardio)
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
cosCardioEst=zeros(nslice,nvol);
sinCardioEst=zeros(nslice,nvol);
multCardioEst=zeros(nslice,nvol);
figure
axes
pkOffset=zeros(nslice,1);
for volNr=1:nvol
    for sliceNr=1:nslice
        x=(volNr-1)*nslice+sliceNr;
        t=(x-1)*tr/nslice;
        slicePhaseOffset=2*pi*maxCardioPkFreq*t;
        sortRestMagSlice=sort(restMag(:,sliceNr),'descend');
        sliceMask=restMag(:,sliceNr)>=sortRestMagSlice(50);
        nrVox=sum(sliceMask);
        sliceMag=restMag(sliceMask,sliceNr);
        restSlice=rest(sliceMask,sliceNr,volNr);
        slicePhase=restPhase(sliceMask,sliceNr);%+slicePhaseOffset;
        if volNr==1
            [~,pkOffset(sliceNr)]=max(mean(repmat(sliceMag,1,nslice).*cos(repmat(slicePhase,1,nslice)+repmat(2*pi*maxCardioPkFreq*(0:nslice-1)*tr/nslice,nrVox,1)),1));
        end
        
        % check phase offset is correct!!!
        slicePhase2d=repmat(slicePhase,1,nrPhaseBins)-repmat(phaseBins,nrVox,1);
        cosMag2d=bsxfun(@times,cos(slicePhase2d),sliceMag);
        sinMag2d=bsxfun(@times,sin(slicePhase2d),sliceMag);
        [~,maxAvIdx]=max(mean(cosMag2d,1));
        multCardioEst(sliceNr,volNr)=sum(restSlice.*cos(slicePhase+2*pi*maxCardioPkFreq*(pkOffset(sliceNr)-1)*tr/nslice))/sum(sliceMag);
        %mean((sliceMag.*cos(slicePhase)));%-2*pi*maxCardioPkFreq*(maxAvIdx-1)*tr/nslice)));
        % scale with mag!!!!
        %sliceMatch=mean(bsxfun(@times,bsxfun(@times,cos(slicePhase2d),sliceMag),restSlice),1);
        sliceMatch=restSlice'*cosMag2d;
        [~,matchIdx]=max(sliceMatch);
        cosEstPhase=phaseBins(matchIdx);
        cosCardioEst(sliceNr,volNr)=cosEstPhase;
        %sliceMatch=mean(bsxfun(@times,bsxfun(@times,cos(slicePhase2d),sliceMag),restSlice),1);
        sliceMatch=restSlice'*sinMag2d;
        [~,matchIdx]=max(sliceMatch);
        sinEstPhase=phaseBins(matchIdx);
        sinCardioEst(sliceNr,volNr)=sinEstPhase;
        %plot(cardio(1:nslice)/max(cardio));hold on;
%         plot(mean(repmat(sliceMag,1,nslice).*cos(repmat(slicePhase,1,nslice)+repmat(2*pi*maxCardioPkFreq*(0:nslice-1)*tr/nslice,nrVox,1)),1))
%         hold on
    end
    unwrapCosCardioEst=unwrap(cosCardioEst(:,volNr));
    plot(1:nslice,cardio(x-nslice+1:x)/max(cardio),1:nslice,multCardioEst(:,volNr)/max(multCardioEst(:)));
    %,1:nslice,unwrapCosCardioEst,1:nslice,cos(cosCardioEst(:,volNr)));
%     p=polyfit((1:nslice)',unwrapCosCardioEst,1);
%     pkIdxSliceFloat=((ceil((p(1)*nslice+p(2))/(2*pi)):floor((p(1)+p(2))/(2*pi)))*2*pi-p(2))/p(1);
%     pkIdxSliceFloat=((ceil((p(1)+p(2))/(2*pi)):floor((p(1)*nslice+p(2))/(2*pi)))*2*pi-p(2))/p(1);
%     pkIdxFloat=(volNr-1)*nslice+pkIdxSliceFloat;
%     pkTime=(pkIdxFloat-1)*tr/nslice;
%     hold on
%     plot(1:nslice,p(1)*(1:nslice)+p(2));
%     plot(round(pkIdxSliceFloat),cardio(round(pkIdxFloat))/max(cardio),'r*');
%     hold off
    %plot(1:nslice,unwrap(cosCardioEst(:,volNr)),1:nslice,multCardioEst(:,volNr),1:nslice,cosCardioEst(:,volNr),1:nslice,sinCardioEst(:,volNr),1:nslice,sin(sinCardioEst(:,volNr)),1:nslice,cos(cosCardioEst(:,volNr)),1:nslice,cardio(x-nslice+1:x)/max(cardio));
end
[r,lags]=xcorr(cos(cosCardioEst(:)),cardio,'coeff');
figure,plot(lags,r);
disp(corr(cosCardioEst(:),reshape(cos(slicePhasesPulseox)',[],1)));

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
