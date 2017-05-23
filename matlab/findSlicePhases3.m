function [maxCardioPkFreq,Pxx,f]=findSlicePhases3(rmsSlice,nrTr,padFactor,minCardioHz,maxCardioHz,tr,restS,slicePhasesPulseox,cardio)
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
%t=reshape(0:tr/nslice:(nvol*nslice-1)*tr/nslice,[],nvol);
t=repmat(0:tr:(nvol-1)*tr,nslice,1);

cosPhase=cos(2*pi*aliasFreq*t);
sinPhase=sin(2*pi*aliasFreq*t);
restInPhase=zeros(prod(volsize(1:2)),nslice);
restInQuad=zeros(prod(volsize(1:2)),nslice);
cosAlias=cos(2*pi*aliasFreq*(0:nvol-1)'*tr);
sinAlias=cos(2*pi*aliasFreq*(0:nvol-1)'*tr);
sumCosAliasSquared=cosAlias'*cosAlias;
sumSinAliasSquared=sinAlias'*sinAlias;
volMax=zeros(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
    restSlice=squeeze(rest(:,sliceNr,:));
    restCos=restSlice*cosAlias/sumCosAliasSquared;
    restSin=restSlice*sinAlias/sumSinAliasSquared;
    restMag=sqrt(restCos.^2+restSin.^2);
    restPhase=oddity*atan2(restSin,restCos);
    sumAngle=angle(sum(restMag.*exp(1i*restPhase)));
    volMax(:,sliceNr)=restMag.*cos(restPhase-sumAngle);
    restInPhase(:,sliceNr)=sum(bsxfun(@times,restSlice,cosPhase(sliceNr,:)),2)/sum(cosPhase(sliceNr,:).^2);
    restInQuad(:,sliceNr)=sum(bsxfun(@times,restSlice,sinPhase(sliceNr,:)),2)/sum(sinPhase(sliceNr,:).^2);
end
restMag=sqrt(restInPhase.^2+restInQuad.^2);
restPhase=oddity*atan2(restInQuad,restInPhase);
vol1phase=unwrap(angle(sum(restMag.*exp(1i*restPhase),1)))';
% volPhase=unwrap(angle(sum(restMag.*exp(1i*restPhase),1)));
%volPhase=reshape(unwrap(reshape(repmat(volPhase',1,nvol)-repmat(2*pi*maxCardioPkFreq*tr*(0:(nvol-1)),nslice,1),[],1)),nslice,nvol)';
nrPhaseBins=180;
phaseBins=(0:2*pi/nrPhaseBins:(nrPhaseBins-1)*2*pi/nrPhaseBins)-pi;
% phaseBins=phaseBins(1:end/2);
% nrPhaseBins=numel(phaseBins);
figure;
rmsCos=zeros(nrPhaseBins,nvol);
volPhase=zeros(nrPhaseBins,1);
cardio=reshape(cardio,[],nvol);
for volNr=1:nvol
%     rmsVolCos=zeros(nrPhaseBins,nslice);
%     rmsVolSin=zeros(nrPhaseBins,nslice);
    for binNr=1:nrPhaseBins
        restVol=squeeze(rest(:,:,volNr));
        %restMagCos=restMag.*cos(restPhase+phaseBins(binNr)-2*pi*maxCardioPkFreq*(volNr-1)*tr);
        %restMagCos=restMag.*cos(restPhase+phaseBins(binNr));
        phaseArg=restPhase-2*pi*maxCardioPkFreq*(volNr-1)*tr+phaseBins(binNr);
        restMagCos=restMag.*cos(phaseArg);
        volBinPhase=unwrap(angle(sum(restMag.*exp(1i*phaseArg),1)))';
        volPhase(binNr)=volBinPhase(1);
%         restMagSin=restMag.*sin(restPhase+phaseBins(binNr)-2*pi*maxCardioPkFreq*(volNr-1)*tr);
%         rmsVolCos(binNr,:)=sum((restVol.*restMagCos).^2,1)./(sum(restMagCos.^2,1).*sum(restMagCos.^2,1));
%         rmsVolSin(binNr,:)=sum((restVol.*restMagSin).^2,1)./(sum(restMagSin.^2,1).*sum(restMagSin.^2,1));
        rmsCos(binNr,volNr)=sum((restVol(:).*restMagCos(:)).^2)./(sum(restMagCos(:).^2).*sum(restMagCos(:).^2));
    end
    [pks,locs]=findpeaks(rmsCos(:,volNr));
    subplot(3,1,1);
    plot(phaseBins,rmsCos(:,volNr));
    %[~,maxBinIdx]=max(rmsCos(:,volNr));
    disp('phaseBins');
    disp(phaseBins(locs));
    disp('volPhase');
    disp(volPhase(locs));
%     sum((restVol.*restMagCos).^2,1)./(sum(restMagCos.^2,1).*sum(restMagCos.^2,1));
%     sum((restMag.*cos(restPhase+phaseBins(maxBinIdx))))./sum(restMag,1);
%     plotyy(1:nslice,cos(unwrap(angle(sum(restMag.*exp(1i*restPhase+1i*phaseBins(maxBinIdx)+1i*2*pi*maxCardioPkFreq*(volNr-1)*tr),1)))'),1:nslice,cardio(:,volNr));
%      plotyy(1:nslice,cos(unwrap(angle(sum(restMag.*exp(1i*restPhase-1i*phaseBins(maxBinIdx)-1i*2*pi*maxCardioPkFreq*(volNr-1)*tr),1)))'),1:nslice,cardio(:,volNr));
%      plotyy(1:nslice,cos(2*pi*maxCardioPkFreq*(0:nslice-1)'*tr/nslice-vol1phase(1)-phaseBins(maxBinIdx)-2*pi*maxCardioPkFreq*(volNr-1)*tr),1:nslice,cardio(:,volNr));
    subplot(3,1,2);
    plotyy(1:nslice,cos(2*pi*maxCardioPkFreq*(0:nslice-1)'*tr/nslice+phaseBins(locs(1))),1:nslice,cardio(:,volNr));
    subplot(3,1,3);
    plotyy(1:nslice,cos(2*pi*maxCardioPkFreq*(0:nslice-1)'*tr/nslice+phaseBins(locs(2))),1:nslice,cardio(:,volNr));

     %    plotyy(1:nslice,cos(unwrap(angle(sum(restMag.*exp(1i*(restPhase+repmat(sum(restPhase,1)./sum(restPhase>0,1),prod(volsize(1:2)),1))-1i*phaseBins(maxBinIdx)-1i*2*pi*maxCardioPkFreq*(volNr-1)*tr),1)))'),1:nslice,cardio(:,volNr));
%     plotyy(1:nslice,cos(unwrap(angle(sum(restMag.*exp(1i*restPhase+1i*phaseBins(maxBinIdx)-1i*2*pi*maxCardioPkFreq*(volNr-1)*tr),1)))'),1:nslice,cardio(:,volNr));
%     plotyy(1:nslice,cos(unwrap(angle(sum(restMag.*exp(1i*restPhase-1i*phaseBins(maxBinIdx)+1i*2*pi*maxCardioPkFreq*(volNr-1)*tr),1)))'),1:nslice,cardio(:,volNr));
    %subplot(5,2,[1,2,3,4]);
    %     [~,maxIdx]=max(rmsCos);
    %     subplot(5,2,[5,6]);plot(rmsCos');
    %     subplot(5,2,[7,8]);plot(cardio((volNr-1)*nslice+1:volNr*nslice))
    %     subplot(5,2,[9,10]);plot(sum((restMag.*cos(restPhase+phaseBins(maxIdx)+2*pi*maxCardioPkFreq*(volNr-1)*tr)).^2,1)./sum(restMag.^2,1))
end
surf(rmsCos)
cosCardioEst=zeros(nslice,nvol);
sinCardioEst=zeros(nslice,nvol);
multCardioEst=zeros(nslice,nvol);
volMaxMatch=zeros(nslice,nvol);
figure
axes
pkOffset=zeros(nslice,1);
for volNr=1:nvol
    for sliceNr=1:nslice
        x=(volNr-1)*nslice+sliceNr;
        t=(x-1)*tr/nslice;
        slicePhaseOffset=2*pi*maxCardioPkFreq*t;
        sortRestMagSlice=sort(restMag(:,sliceNr),'descend');
        sliceMask=restMag(:,sliceNr)>sortRestMagSlice(50);
        nrVox=sum(sliceMask);
        sliceMag=restMag(sliceMask,sliceNr);
        restSlice=rest(sliceMask,sliceNr,volNr);
        
        volMaxMatch(sliceNr,volNr)=sum((restSlice.*volMax(sliceMask,sliceNr)).^2)/(sum(volMax(sliceMask,sliceNr).^2)^2);
        
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
    plot(1:nslice,cardio(x-nslice+1:x)/max(cardio),1:nslice,volMaxMatch(:,volNr)/max(volMaxMatch(:)));
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
