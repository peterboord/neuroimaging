function [allCardioPeakTimes,Pxx,f]=calcPeakTimes(restS,rmsSlice,nrTr,padFactor)
nslice=size(rmsSlice,1);
nvol=size(rmsSlice,2);
tr=restS.tr;
% find peak in cardio range
[maxCardioPkFreq,Pxx,f]=findPhysioPeakFreq('cardio',restS,nrTr,padFactor,nrMaskVoxPerSlice);
dfSlice=1/(nrTr*tr);
allCardioPeakTimes=[];
for volNr=1:nvol-nrTr+1
    x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
    rmsWave=reshape(rmsSlice(x),[],1);
    rmsSliceFft=fft(rmsWave,padFactor*numel(rmsWave));
    %rmsSliceFft=fft(rmsWave.*hann(numel(rmsWave),'periodic'),padFactor*numel(rmsWave));
    fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
    %     [Pmt,Fmt]=pmtm(rmsWave,4,padFactor*numel(x)/2,nslice/tr);
    %     subplot(2,1,2); plot(Fmt,Pmt);
    %     minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
    %     maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
    minCardio=round((maxCardioPkFreq-0.1)/(dfSlice/padFactor))+1;
    maxCardio=round((maxCardioPkFreq+0.1)/(dfSlice/padFactor))+1;
    [~,maxIdx]=max(abs(rmsSliceFft(minCardio:maxCardio)));
    maxSpecIdx=minCardio+maxIdx-1;
    maxSpecAngle=angle(rmsSliceFft(maxSpecIdx));
    % xfm (-pi,pi) to (0,2*pi)
    if maxSpecAngle < 0
        maxSpecAngle=2*pi+maxSpecAngle;
    end
    if volNr==nvol-nrTr+1
        nrTrToSearch=nrTr;
    else
        nrTrToSearch=1;
    end
    maxNrPeakInTr=floor(nrTrToSearch*tr*fSlice(maxSpecIdx))+1;
    peakTimes=(2*pi*(1:maxNrPeakInTr)'-maxSpecAngle)/(2*pi*fSlice(maxSpecIdx));
    peakTimesExtra=peakTimes(peakTimes>tr);
    peakTimes(peakTimes>=tr)=[];
    if volNr==nvol-nrTr+1 && ~isempty(peakTimesExtra)
        peakTimes=[peakTimes;peakTimesExtra]; %#ok<AGROW>
    end
    allCardioPeakTimes=cat(1,allCardioPeakTimes,(volNr-1)*tr+peakTimes);
end
% 
%figure('WindowStyle','docked');
allCardioPeakTimes=fixPhysioTimes(allCardioPeakTimes);

%diffAll=diff(allCardioPeakTimes);
%plot(diffAll)

% figure('WindowStyle','docked');
% x=1:1000;
% subplot(2,1,1);
% rr=interp1(allCardioPeakTimes(1:end-1),diff(allCardioPeakTimes),allCardioPeakTimes(1):0.1:allCardioPeakTimes(end-1));
% plot(rr(x));
% subplot(2,1,2);
% [Pxx,f]=pwelch(detrend(rr(x)),100,50,100,1/0.1);plot(f,detrend(Pxx,'constant'))
end

