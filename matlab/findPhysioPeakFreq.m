function [maxPhysioPkFreq,Pxx,f]=findPhysioPeakFreq(physioType,restS,nrTr,padFactor,percentSliceVox)
% constants
minCardioHz=0.6;
maxCardioHz=1.5;
minRespHz=5/60;
maxRespHz=30/60;
% variables
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
rest=reshape(restS.vol,[],nslice,nvol);
% physioType
switch physioType
    case 'cardio'
        minPhysioHz=minCardioHz;
        maxPhysioHz=maxCardioHz;
    case 'resp'
        minPhysioHz=minRespHz;
        maxPhysioHz=maxRespHz;
    otherwise
        error('unknown physio type');
end
% rmsSlice
rmsSliceMode='getRmsSlice';
switch rmsSliceMode
    case 'getRmsSlice'
        rmsSlice=getRmsSlice2(restS.vol,percentSliceVox);
    otherwise
        sumRestSquared=squeeze(sum(rest.^2,1));
        rmsSlice=bsxfun(@rdivide,sumRestSquared,median(sumRestSquared,2));
end
% find peak in Physio range
fftRmsSlice=fft(detrend(rmsSlice,'constant'),nrTr*nslice*padFactor);
f=(0:(1/(nrTr*restS.tr))/padFactor:floor(((size(fftRmsSlice,1)-1)*(1/(nrTr*restS.tr))/padFactor)/2)+1)';
Pxx=mean(abs(fftRmsSlice(1:floor(end/2)+1,:)),2);

%[Pxx,f]=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/restS.tr);
dfSlice=1/(nrTr*restS.tr);
minPhysioIdx=floor(minPhysioHz*padFactor/dfSlice);
maxPhysioIdx=ceil(maxPhysioHz*padFactor/dfSlice)+1;
if numel(unique(Pxx))~=numel(Pxx)
    error('part of waveform could be flat, so peaks might be missed');
end
physioPeakIdxs=[0;Pxx(2:end-1)>Pxx(1:end-2) & Pxx(2:end-1)>Pxx(3:end);0];
physioPeakIdxs(1:minPhysioIdx-1)=0;
physioPeakIdxs(maxPhysioIdx+1:end)=0;
[maxPhysioPkMag,maxPhysioPkIdx]=max(Pxx.*physioPeakIdxs);
% if no peak in physio range use max boundary
if isempty(maxPhysioPkMag)
    if Pxx(minPhysioIdx)>=Pxx(maxPhysioIdx)
        maxPhysioPkIdx=minPhysioIdx;
    else
        maxPhysioPkIdx=maxPhysioIdx;
    end
end
maxPhysioPkFreq=(maxPhysioPkIdx-1)*dfSlice/padFactor;
disp([physioType,' rate is ',num2str(maxPhysioPkFreq),' Hz']);
% figure,plot(f,Pxx);
end