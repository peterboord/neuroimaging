function [maxCardioPkFreq,Pxx,f]=findSlicePhases4(restS,nrTr,padFactor,cardio)
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
tr=restS.tr;
[maxCardioPkFreq,Pxx,f]=findPhysioPeakFreq('cardio',restS,nrTr,padFactor,nrMaskVoxPerSlice);
rest=reshape(restS.vol,[],nslice,nvol);
t=repmat(0:tr:(nvol-1)*tr,nslice,1);
bpmSideRange=10;
cardFreqRange=maxCardioPkFreq-bpmSideRange/60:1/60:maxCardioPkFreq+bpmSideRange/60;
nrCardFreq=numel(cardFreqRange);
restInPhase=zeros(prod(volsize(1:2)),nslice);
restInQuad=zeros(prod(volsize(1:2)),nslice);
restMag=zeros(prod(volsize(1:2)),nslice,nrCardFreq);
restPhase=zeros(prod(volsize(1:2)),nslice,nrCardFreq);
disp({'aliasFreq','N','oddity'});
for cardFreqIdx=1:nrCardFreq
    [aliasFreq,N,oddity]=calcAliases(cardFreqRange(cardFreqIdx),tr);
    disp([aliasFreq,N,oddity]);
    cosPhase=cos(2*pi*aliasFreq*t);
    sinPhase=sin(2*pi*aliasFreq*t);
    for sliceNr=1:nslice
        restSlice=squeeze(rest(:,sliceNr,:));
        restInPhase(:,sliceNr)=sum(bsxfun(@times,restSlice,cosPhase(sliceNr,:)),2)/sum(cosPhase(sliceNr,:).^2);
        restInQuad(:,sliceNr)=sum(bsxfun(@times,restSlice,sinPhase(sliceNr,:)),2)/sum(sinPhase(sliceNr,:).^2);
    end
    restMag(:,:,cardFreqIdx)=sqrt(restInPhase.^2+restInQuad.^2);
    restPhase(:,:,cardFreqIdx)=oddity*atan2(restInQuad,restInPhase);
end
nrPhaseBins=180;
phaseBins=(0:2*pi/nrPhaseBins:(nrPhaseBins-1)*2*pi/nrPhaseBins)-pi;
rmsCos=zeros(nrPhaseBins,nrCardFreq,nvol);
cardio=reshape(cardio,[],nvol);
for volNr=1:nvol
    restVol=squeeze(rest(:,:,volNr));
%         pkAmp=zeros(nrCardFreq,1);
%     pkPhase=zeros(nrCardFreq,1);
    pkAmp=zeros(nrCardFreq,2);
    pkPhase=zeros(nrCardFreq,2);
    for cardFreqIdx=1:nrCardFreq
        for binNr=1:nrPhaseBins
            phaseArg=restPhase(:,:,cardFreqIdx)-2*pi*cardFreqRange(cardFreqIdx)*(volNr-1)*tr+phaseBins(binNr);
            restMagCos=restMag(:,:,cardFreqIdx).*cos(phaseArg);
            restVolByRestMagCos=restVol(:).*restMagCos(:);
            rmsCos(binNr,cardFreqIdx,volNr)=sum(restVolByRestMagCos.^2)./(sum(restMagCos(:).^2).*sum(restMagCos(:).^2));
        end
%         [maxAmp,maxIdx]=max(rmsCos(:,cardFreqIdx,volNr));
%         pkAmp(cardFreqIdx)=maxAmp;
%         pkPhase(cardFreqIdx)=phaseBins(maxIdx);

        [maxAmp,maxIdx]=findpeaks(rmsCos(:,cardFreqIdx,volNr));
        [maxAmp,sortIdx]=sort(maxAmp,'descend');
        if numel(maxAmp)>1
            pkAmp(cardFreqIdx,:)=maxAmp(1:2);
            pkPhase(cardFreqIdx,:)=phaseBins(maxIdx(sortIdx(1:2)));
        else
            pkAmp(cardFreqIdx,:)=maxAmp;
            pkPhase(cardFreqIdx,:)=phaseBins(maxIdx(sortIdx));
        end
    end
    [~,maxFreqIdx]=max(pkAmp);
    cos1=cos(2*pi*cardFreqRange(maxFreqIdx(1))*(0:nslice-1)'*tr/nslice+pkPhase(maxFreqIdx(1),1));
    [~,pkIdx1]=findpeaks(cos1);
    cos2=cos(2*pi*cardFreqRange(maxFreqIdx(1))*(0:nslice-1)'*tr/nslice+pkPhase(maxFreqIdx(1),2));
    [~,pkIdx2]=findpeaks(cos2);
    disp([(sum(restVol(:,pkIdx1)).^2)./sum(restVol(:,pkIdx1)~=0),0,(sum(restVol(:,pkIdx2)).^2)./sum(restVol(:,pkIdx2)~=0)]);
    figure('WindowStyle','docked');
    plotyy(1:nslice,cos1,1:nslice,cardio(:,volNr));
    hold on
    plot(1:nslice,cos2,'g');
    slicePower=(sum(restVol.^2,1)./mean(squeeze(sum(rest.^2,1)),2)')';
    disp(corr(slicePower,cos1));
    disp(corr(slicePower,cos2));
    plot(1:nslice,slicePower,'m');
    hold off
end

end
