function findCardioPeaks(session)

tic
dbstop if error
if nargin==0
    %session='RC4107-2';
    session='RC4103-1';
    %session='RC4109-1';
end
disp(session);
%fmri
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
% 
tr=2.4;
nrTr=2; % optimal for cardio
%nrTr=4; % optimal for resp
nrMaskVoxPerSlice=50;
%restS=MRIread(fullfile(funcDir,session,'prefiltered_func_data_mcf.nii.gz'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.nii.gz']));
restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],[session,'.results'],['pb01.',session,'.r01.despike.nii']));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));
%arCoeffS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_Tar1.nii'));
nvol=restS.nframes;
volsize=restS.volsize;
nslice=volsize(3);
rest=reshape(restS.vol,prod(volsize),nvol)';
%rest=[rest(1,:);rest(2:end,:)-bsxfun(@times,rest(1:end-1,:),reshape(arCoeffS.vol,1,[]))];
rest=detrend(rest);
minCardioHz=0.6;
maxCardioHz=1.5;
rest=reshape(rest',prod(volsize(1:2)),nslice,nvol);    
% make mask from 50 voxels/slice with highest Tstd
varRest=var(rest,1,3);
[varRestSize,varRestSortIdx]=sort(varRest,1,'descend');
if any(reshape(varRestSize(1:nrMaskVoxPerSlice,:),[],1)==0)
    disp('WARNING: mask has voxels with Tvar = 0');
end
topvarIdxs=varRestSortIdx(1:nrMaskVoxPerSlice,1:nslice);
varMask=false(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
   varMask(:,sliceNr)=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
end
rmsSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(varMask(:,sliceNr),sliceNr,:));
    rmsSlice(sliceNr,:)=sqrt(mean(cohRestSlice.^2,1));
    %rmsSlice(sliceNr,:)=sqrt(median(cohRestSlice.^2,1));
    % normalize 
    rmsSlice(sliceNr,:)=rmsSlice(sliceNr,:)/std(rmsSlice(sliceNr,:));
end
rmsSlice=detrend(rmsSlice')';
avSliceFft=mean(abs(fft(rmsSlice')),2);
dfSlice=1/(nrTr*tr);
padFactor=4;
disp(padFactor);
[Pxx,f]=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);

% figure
% plot(f,Pxx)
% find peak in cardio range
minCardioIdx=floor(minCardioHz*padFactor/dfSlice);
maxCardioIdx=ceil(maxCardioHz*padFactor/dfSlice)+1;
% figure
% plot(f(minCardioIdx:maxCardioIdx),Pxx(minCardioIdx:maxCardioIdx))
cardioPeakIdxs=[0;Pxx(2:end-1)>Pxx(1:end-2) & Pxx(2:end-1)>Pxx(3:end);0];
cardioPeakIdxs(1:minCardioIdx-1)=0;
cardioPeakIdxs(maxCardioIdx+1:end)=0;
[~,maxCardioPkIdx]=max(Pxx.*cardioPeakIdxs);
maxCardioPkFreq=(maxCardioPkIdx-1)*dfSlice/padFactor;
fs=1/tr;
% fa=abs(f-N*fs)
N=ceil(maxCardioPkFreq/fs-0.5):floor(maxCardioPkFreq/fs+0.5);
fa=abs(maxCardioPkFreq-N*fs);
df=1/(nvol*tr);
fScan=0:df:(nvol-1)*df;
faIdx=round(fa/df)+1;
faIdx=nvol/2+1;

if true
% make mask from nrMaskVoxPerSlice voxels with top power in cardio band
restFft=fft(reshape(rest,[],nvol)');
cardioAliasPowerUnsorted=reshape(sum(bsxfun(@times,abs(restFft),avSliceFft),1),prod(volsize(1:2)),nslice);
% kurtAbsFft=kurtosis(abs(restFft(1:end/2+1,:)));
% kurtAbsFft(isnan(kurtAbsFft))=Inf;
% cardioAliasPowerUnsorted=-1*reshape(kurtAbsFft,prod(volsize(1:2)),nslice);
% cardioAliasPowerUnsorted=reshape(std(reshape(rest,[],nvol)'),prod(volsize(1:2)),nslice);
% cardioAliasPowerUnsorted(cardioAliasPowerUnsorted==0)=-Inf;
absRestFft=reshape(abs(restFft(1:end/2+1,:))',prod(volsize(1:2)),nslice,[]);
% cardioAliasPowerUnsorted=reshape(abs(restFft(faIdx,:))',[],nslice);
[cardioAliasPower,cardioAliasSortIdx]=sort(cardioAliasPowerUnsorted,1,'descend');
topPowerIdxs=cardioAliasSortIdx(1:nrMaskVoxPerSlice,1:nslice);
cardioMask=false(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
   cardioMask(:,sliceNr)=ismember((1:prod(volsize(1:2)))',topPowerIdxs(:,sliceNr));
end
maskS=resizeFuncS(restS,3);
maskS.vol=reshape(cardioMask,restS.volsize);
maskS.fspec=fullfile(funcDir,session,'cardioMask.nii');
MRIwrite(maskS,maskS.fspec);
rmsSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(cardioMask(:,sliceNr),sliceNr,:));
    rmsSlice(sliceNr,:)=sqrt(mean(cohRestSlice.^2,1));
    % normalize 
    rmsSlice(sliceNr,:)=rmsSlice(sliceNr,:)/std(rmsSlice(sliceNr,:));
end
rmsSlice=detrend(rmsSlice')';
[Pxx_cm,f_cm]=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
% figure
% plot(f_cm,Pxx_cm)

% [Pmt,Fmt]=pmtm(rmsSlice(:),4,numel(rmsSlice(:)),nslice/tr);
% figure,
% plot(Fmt,Pmt);
end
% physio
physioDir='/projects2/udall/physio';
physioFs=100; % actually 5 replications of 100 Hz = 500 Hz, but replications removed
physioSamples=physioFs*nvol*tr;
% cardio
cardioFilePath=fullfile(physioDir,session,'rest_cardio.txt');
cardio=[];
if exist(cardioFilePath,'file')
    cardioRaw=reshape(textread(cardioFilePath),5,[]);
    cardioRaw=cardioRaw(1,:);
    cardio=resample(cardioRaw,nvol*nslice,physioSamples); %#ok<*DTXTRD>
    % % resp
    % respRaw=reshape(textread(fullfile(physioDir,session,'rest_resp.txt')),5,[]);
    % respRaw=respRaw(1,:);
    % resp=resample(respRaw,nTotalSlices,physioSamples);
end
[pks,locs]=findpeaks(cardio,'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
% figure
% plot(cardio);
% hold on
% plot(locs,pks,'r*');
% hold off
figure('WindowStyle','docked');
nrPlots=4;
subplot(nrPlots,1,1);
plot(f,Pxx);
allCardioPeakTimes=[];
%allCardioPeakTimes=(locs'-1)*tr/nslice;
for volNr=1:nvol-nrTr+1
    x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
    t=(0:nrTr*nslice-1)*tr/nslice;
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
        peakTimes=[peakTimes;peakTimesExtra];
    end
    if false%true
    % plots
    subplot(nrPlots,1,2); plot(t,rmsWave)
    cardioOffset=-4;
    cardioEst=abs(rmsSliceFft(maxSpecIdx))*cos(2*pi*fSlice(maxSpecIdx)*(cardioOffset+(0:nrTr*nslice-1))*tr/nslice+angle(rmsSliceFft(maxSpecIdx)));
    subplot(nrPlots,1,3);
    plot(t,100*cardioEst);
    if ~isempty(cardio)
        hold on; plot(t,cardio(x)); hold off
    end
    subplot(nrPlots,1,4); plot(fSlice(1:end/2),abs(rmsSliceFft(1:end/2)));
    keyboard
    end
    allCardioPeakTimes=cat(1,allCardioPeakTimes,(volNr-1)*tr+peakTimes);
end
% augment with pre- and post-beat estimate
nrToAdd=ceil(allCardioPeakTimes(1)/(allCardioPeakTimes(2)-allCardioPeakTimes(1)));
for addNr=1:nrToAdd
    allCardioPeakTimes=[2*allCardioPeakTimes(1)-allCardioPeakTimes(2);allCardioPeakTimes]; %#ok<AGROW>
end
nrToAdd=ceil((nvol*tr-allCardioPeakTimes(end))/(allCardioPeakTimes(end)-allCardioPeakTimes(end-1)));
for addNr=1:nrToAdd
    allCardioPeakTimes=[allCardioPeakTimes;2*allCardioPeakTimes(end)-allCardioPeakTimes(end-1)]; %#ok<AGROW>
end
slicePhases=zeros(nvol,nslice);
for volNr=1:nvol
    for sliceNr=1:nslice
        t=((volNr-1)*nslice+sliceNr-1)*tr/nslice;
        peakTimeBefore=allCardioPeakTimes(find(allCardioPeakTimes<t,1,'last'));
        peakTimeAfter=allCardioPeakTimes(find(allCardioPeakTimes>t,1,'first'));
        slicePhases(volNr,sliceNr)=2*pi*(t-peakTimeBefore)/(peakTimeAfter-peakTimeBefore);
    end
end
% RETROICOR
M=2;
a=zeros(prod(volsize(1:2)),nslice,M);
b=zeros(prod(volsize(1:2)),nslice,M);
noiseEst=zeros(prod(volsize(1:2)),nslice,nvol);
abPhase=zeros(prod(volsize(1:2)),nslice,M);
for sliceNr=1:nslice
    % calc coeffs
    for m=1:M
        a(:,sliceNr,m)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:))',cos(m*slicePhases(:,sliceNr))),1)/sum((cos(m*slicePhases(:,sliceNr)).^2));
        b(:,sliceNr,m)=sum(bsxfun(@times,squeeze(rest(:,sliceNr,:))',sin(m*slicePhases(:,sliceNr))),1)/sum((sin(m*slicePhases(:,sliceNr)).^2));
        abPhase(:,sliceNr,m)=atan2(b(:,sliceNr,m),a(:,sliceNr,m));
    end
    % calc noise estimate
    noiseEst(:,sliceNr,:)=squeeze(a(:,sliceNr,:))*cos(slicePhases(:,sliceNr)*(1:M))'+squeeze(b(:,sliceNr,:))*sin(slicePhases(:,sliceNr)*(1:M))';
end
% save 4D images
restS.vol=reshape(rest-noiseEst,[volsize,nvol]);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii');
MRIwrite(restS,restS.fspec);
restS.vol=reshape(noiseEst,[volsize,nvol]);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst.nii');
MRIwrite(restS,restS.fspec);
% save 3D images
restS.dim(5) = 1;
restS.nframes = 1;
restS.vol=reshape(sqrt(mean(noiseEst.^2,4)),[volsize,nvol]);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst_rms.nii');
MRIwrite(restS,restS.fspec);
varOld=var(rest,0,3);
varNew=var(rest-noiseEst,0,3);
propVar=100*(varOld-varNew)./varOld;
propVar(isnan(propVar))=0;
propVar(isinf(propVar))=0;
restS.vol=reshape(propVar,volsize);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor_propVar.nii');
MRIwrite(restS,restS.fspec);
restS.vol=reshape(abPhase(:,:,1),volsize);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_m1phase.nii');
MRIwrite(restS,restS.fspec);
restS.vol=reshape(abPhase(:,:,2),volsize);
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_m2phase.nii');
MRIwrite(restS,restS.fspec);

funcMask=~all(reshape(rest,[],nvol)==0,2);
subplot(nrPlots,1,2)
hist(propVar(funcMask))
subplot(nrPlots,1,3)
mc=load(fullfile(funcDir,session,'mocoBetHpf.feat','mc','prefiltered_func_data_mcf.par'),'-ascii');
plot(diff(mc(:,4:6),1,1))
ylabel('diff trans (mm)');
subplot(nrPlots,1,4)
plot(diff(mc(:,1:3),1,1)*180/pi)
ylabel('diff rot (deg)');
toc
