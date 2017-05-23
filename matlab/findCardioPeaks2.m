function allCardioPeakTimes=findCardioPeaks2(session)

tic
dbstop if error
if nargin==0
    %session='RC4107-2';
    session='RC4103-1';
end
disp(session);
%fmri
maskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
tr=2.4;
restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],'filtered_func_data.nii'));
maskPath=fullfile(funcDir,session,maskFile);
if ~exist(maskPath,'file')
    system(['cd ',fullfile(funcDir,session),...
        ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',maskFile]);
end
maskS=MRIread(maskPath);
nvol=restS.nframes;
volsize=restS.volsize;
nslice=volsize(3);
rest=reshape(restS.vol,prod(volsize),nvol);
rest=detrend(rest')';
rest=reshape(rest,volsize(1)*volsize(2),nslice,nvol);
mask=logical(reshape(maskS.vol,volsize(1)*volsize(2),nslice));
voxInSlice=cell(nslice,1);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(mask(:,sliceNr),sliceNr,:));
    cohRestSlice=detrend(cohRestSlice')';
    voxInSlice{sliceNr}=cohRestSlice;
end
sqssSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    sqssSlice(sliceNr,:)=sqrt(mean(voxInSlice{sliceNr}.^2,1));
end
sqssSlice=detrend(sqssSlice')';
nrTr=2;
dfSlice=1/(nrTr*tr);
allCardioPeakTimes=[];
for volNr=1:nvol-nrTr+1
    x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
    sqssWave=reshape(sqssSlice(x),[],1);
    padFactor=6;
    sqssSliceFft=fft(sqssWave,padFactor*numel(sqssWave));
    sqssSliceFft=fft(sqssWave.*hann(numel(sqssWave),'periodic'),padFactor*numel(sqssWave));
    fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
    if volNr < 10
        figure('WindowStyle','docked');
        plot(fSlice(1:end/2),abs(sqssSliceFft(1:end/2)))
    else
        keyboard
    end
    minCardioHz=0.6;
    maxCardioHz=1.5;
    minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
    maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
    [~,maxIdx]=max(abs(sqssSliceFft(minCardio:maxCardio)));
    maxSpecIdx=minCardio+maxIdx-1;
    maxSpecAngle=angle(sqssSliceFft(maxSpecIdx));
    % xfm (-pi,pi) to (0,2*pi)
    if maxSpecAngle < 0
        maxSpecAngle=2*pi+maxSpecAngle;
    end
    maxNrPeaksInWindow=floor(2*tr*fSlice(maxSpecIdx))+1;
    peakTimesInWindow=(2*pi*(1:maxNrPeaksInWindow)'-maxSpecAngle)/(2*pi*fSlice(maxSpecIdx));
    nrPeaksOutsideTr=sum(peakTimesInWindow>=tr);
    if nrPeaksOutsideTr >=2
        peakTimes=peakTimesInWindow(1:end-nrPeaksOutsideTr+1);
    end
    allCardioPeakTimes=cat(1,allCardioPeakTimes,(volNr-1)*tr+peakTimes);
end
allCardioPeakTimes(allCardioPeakTimes>nvol*tr)=[];
allCardioPeakTimes=sort(allCardioPeakTimes);
% average peak times for close peaks within a tr
closePeaks=diff(allCardioPeakTimes) < median(diff(allCardioPeakTimes))/2;
avClosePeaks=mean([allCardioPeakTimes(closePeaks),allCardioPeakTimes([false;closePeaks(1:end-1)])],2);
allCardioPeakTimes(closePeaks | [false;closePeaks(1:end-1)])=[];
allCardioPeakTimes=sort([allCardioPeakTimes;avClosePeaks]);
% repeat for peaks close across tr boundaries
closePeaks=diff(allCardioPeakTimes) < median(diff(allCardioPeakTimes))/2;
avClosePeaks=mean([allCardioPeakTimes(closePeaks),allCardioPeakTimes([false;closePeaks(1:end-1)])],2);
allCardioPeakTimes(closePeaks | [false;closePeaks(1:end-1)])=[];
allCardioPeakTimes=sort([allCardioPeakTimes;avClosePeaks]);
toc
