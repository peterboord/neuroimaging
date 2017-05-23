function findCardioPeaksV1_2(session)

tic
dbstop if error
if nargin==0
    session='RC4103-1';
    %session='RC4107-2';
    %session='RC4109-1';
end
disp(session);
figure('WindowStyle','docked');
nrPlots=5;
%fmri
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
% roi
switch session
    case 'RC4103-1'
        roi=[22,41,30];
    case 'RC4107-2'
        roi=[28,40,28];
    case 'RC4109-1'
        roi=[25,40,35];
    otherwise
        disp('roi not specified');
        roi=[];
end
% constants
tr=2.4;
nrTr=2; % optimal for cardio
%nrTr=4; % optimal for resp

nrMaskVoxPerSlice=50;

minCardioHz=0.6;
maxCardioHz=1.5;
padFactor=4;
physioDir='/projects2/udall/physio';
%restS=MRIread(fullfile(funcDir,session,'prefiltered_func_data_mcf.nii.gz'));
restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],[session,'.results'],['pb01.',session,'.r01.despike.nii']));
% derived constants
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
disp(padFactor);
% dirty data
[rmsSlice,varMask]=getRmsSlice(restS.vol,nrMaskVoxPerSlice);
subplot(nrPlots,1,5);
plot(rmsSlice(1:nslice*10));
[allCardioPeakTimes,Pxx,f]=calcPeakTimes(rmsSlice,nrTr,padFactor,minCardioHz,maxCardioHz,tr);
% clean data
[~,~,abMag,abPhase]=retroicor(restS.vol,allCardioPeakTimes,tr);
equalizedMask=makeEqualizedMask(restS,abMag,abPhase,nrMaskVoxPerSlice);
rmsSliceEqualized=getRmsSlice(restS.vol,nrMaskVoxPerSlice,equalizedMask);
% re-calc peaks & clean data
[allCardioPeakTimes,PxxEqualized]=calcPeakTimes(rmsSliceEqualized,nrTr,padFactor,minCardioHz,maxCardioHz,tr);
[clean4d,noise4d]=retroicor(restS.vol,allCardioPeakTimes,tr);
rmsSliceClean=getRmsSlice(clean4d,nrMaskVoxPerSlice,equalizedMask);
PxxSliceClean=pwelch(rmsSliceClean(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
% pulse ox
[~,allCardioPeakTimes,resp]=getPeaksFromPhysioData(fullfile(physioDir,session,'rest_cardio.txt'),fullfile(physioDir,session,'rest_resp.txt'),tr,nvol,nslice);
clean4dPulseox=retroicor(restS.vol,allCardioPeakTimes,tr);
rmsSlicePulseox=getRmsSlice(clean4dPulseox,nrMaskVoxPerSlice,varMask);
PxxPulseox=pwelch(rmsSlicePulseox(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
% aliased slice spectrum
PaliasDirty=zeros(nrTr*nslice/2+1,nslice);
PaliasEqualized=zeros(nrTr*nslice/2+1,nslice);
PaliasClean=zeros(nrTr*nslice/2+1,nslice);
for sliceNr=1:nslice
    [PxxAlias,falias]=pwelch(rmsSlice(sliceNr,:),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    PaliasDirty(:,sliceNr)=PxxAlias;
    PxxAliasEqualized=pwelch(rmsSliceEqualized(sliceNr,:),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    PaliasEqualized(:,sliceNr)=PxxAliasEqualized;
    PxxAliasCleaned=pwelch(rmsSliceClean(sliceNr,:),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    PaliasClean(:,sliceNr)=PxxAliasCleaned;
end
% aliased slice data
ax1=subplot(nrPlots,1,1);
plot(falias,mean(PaliasDirty,2),'b');
hold(ax1,'on');
plot(ax1,falias,mean(PaliasEqualized,2),'g');
plot(ax1,falias,mean(PaliasClean,2),'r');
hold(ax1,'off');
% unaliased slice datanoise4d
ax2=subplot(nrPlots,1,2);
plot(f,Pxx,'b');
hold(ax2,'on');
plot(ax2,f,PxxEqualized,'g');
plot(ax2,f,PxxSliceClean,'r');
plot(ax2,f,PxxPulseox,'m');
hold(ax2,'off');

% aliased roi data
if ~isempty(roi)
    subplot(nrPlots,1,3);
    rest=reshape(detrend(reshape(restS.vol,[],nvol)')',size(restS.vol));
    [Proi,froi]=pwelch(squeeze(rest(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(froi,Proi,'b');
    hold on
    Proi=pwelch(squeeze(clean4d(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(froi,Proi,'r');
    Proi=pwelch(squeeze(clean4dPulseox(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(froi,Proi,'m');
    hold off
end

% plot spec for resp peak
nrTr=4;
%rmsSliceClean=getRmsSlice(clean4d,nrMaskVoxPerSlice);
[PxxNrTr4,fNrTr4]=pwelch(rmsSliceEqualized(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
subplot(nrPlots,1,4);
plot(fNrTr4,PxxNrTr4,'b');
PxxResp=pwelch(resp(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
hold on
plot(fNrTr4,PxxResp*max(PxxNrTr4)/max(PxxResp),'m');
hold off

% calc proportion of variance
varOld=var(restS.vol,0,4);
varNew=var(clean4d,0,4);
propVar=100*(varOld-varNew)./varOld;
propVar(isnan(propVar))=0;
propVar(isinf(propVar))=0;
% save 4D images
restS.vol=clean4d;
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii');
MRIwrite(restS,restS.fspec);
restS.vol=noise4d;
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst.nii');
MRIwrite(restS,restS.fspec);
% save 3D images
restS.dim(5) = 1;
restS.nframes = 1;
restS.vol=sqrt(mean(noise4d.^2,4));
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst_rms.nii');
MRIwrite(restS,restS.fspec);
restS.vol=propVar;
restS.fspec=fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor_propVar.nii');
MRIwrite(restS,restS.fspec);
toc
end

%%%%%%%%%%%%%%%%%%%%%%%
% getRmsSlice
%%%%%%%%%%%%%%%%%%%%%%%
function [rmsSlice,voxMask]=getRmsSlice(rest,nrMaskVoxPerSlice,voxMask)
restSize=size(rest);
volsize=restSize(1:3);
nvol=restSize(4);
nslice=restSize(3);
rest=reshape(rest,prod(volsize),nvol)';
rest=detrend(rest);
rest=reshape(rest',prod(volsize(1:2)),nslice,nvol);
if nargin < 3
    % make mask from 50 voxels/slice with highest Tstd
    varRest=var(rest,1,3);
    [varRestSize,varRestSortIdx]=sort(varRest,1,'descend');
    if any(reshape(varRestSize(1:nrMaskVoxPerSlice,:),[],1)==0)
        disp('WARNING: mask has voxels with Tvar = 0');
    end
    topvarIdxs=varRestSortIdx(1:nrMaskVoxPerSlice,1:nslice);
    voxMask=false(prod(volsize(1:2)),nslice);
    for sliceNr=1:nslice
        voxMask(:,sliceNr)=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
    end
end
varRest=var(rest,1,3);
rest=rest./repmat(reshape(sqrt(varRest),prod(volsize(1:2)),nslice,1),1,1,nvol);
rest(isinf(rest))=0;
rest(isnan(rest))=0;
rmsSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(voxMask(:,sliceNr),sliceNr,:));
    rmsSlice(sliceNr,:)=sqrt(mean(cohRestSlice.^2,1));
    rmsSlice(sliceNr,:)=rmsSlice(sliceNr,:)/std(rmsSlice(sliceNr,:));
end
rmsSlice=detrend(rmsSlice')';
end
%%%%%%%%%%%%%%%%%%%%%%%
% calcPeakTimes
%%%%%%%%%%%%%%%%%%%%%%%
function [allCardioPeakTimes,Pxx,f]=calcPeakTimes(rmsSlice,nrTr,padFactor,minCardioHz,maxCardioHz,tr)
nslice=size(rmsSlice,1);
nvol=size(rmsSlice,2);
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
    if false%true
        % plots
        t=(0:nrTr*nslice-1)*tr/nslice; %#ok<UNRCH>
        subplot(nrPlots,1,2); plot(t,rmsWave)
        cardioEst=abs(rmsSliceFft(maxSpecIdx))*cos(2*pi*fSlice(maxSpecIdx)*((0:nrTr*nslice-1))*tr/nslice+angle(rmsSliceFft(maxSpecIdx)));
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
end
%%%%%%%%%%%%%%%%%%%%%%%
% retroicor
%%%%%%%%%%%%%%%%%%%%%%%
function [clean4d,noise4d,abMag,abPhase]=retroicor(rest,allCardioPeakTimes,tr)
M=2;
restSize=size(rest);
volsize=restSize(1:3);
nvol=restSize(4);
nslice=restSize(3);
rest=reshape(rest,prod(volsize),nvol)';
rest=detrend(rest);
rest=reshape(rest',prod(volsize(1:2)),nslice,nvol);
slicePhases=zeros(nvol,nslice);
% augment with pre- and post-beat estimate
nrToAdd=ceil(allCardioPeakTimes(1)/(allCardioPeakTimes(2)-allCardioPeakTimes(1)));
for addNr=1:nrToAdd
    allCardioPeakTimes=[2*allCardioPeakTimes(1)-allCardioPeakTimes(2);allCardioPeakTimes]; %#ok<AGROW>
end
nrToAdd=ceil((nvol*tr-allCardioPeakTimes(end))/(allCardioPeakTimes(end)-allCardioPeakTimes(end-1)));
for addNr=1:nrToAdd
    allCardioPeakTimes=[allCardioPeakTimes;2*allCardioPeakTimes(end)-allCardioPeakTimes(end-1)]; %#ok<AGROW>
end
for volNr=1:nvol
    for sliceNr=1:nslice
        t=((volNr-1)*nslice+sliceNr-1)*tr/nslice;
        peakTimeBefore=allCardioPeakTimes(find(allCardioPeakTimes<=t,1,'last'));
        peakTimeAfter=allCardioPeakTimes(find(allCardioPeakTimes>t,1,'first'));
        slicePhases(volNr,sliceNr)=2*pi*(t-peakTimeBefore)/(peakTimeAfter-peakTimeBefore);
    end
end
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
%%%%%%%%%%%%%%%%%%%%%%%
% getPeaksFromPhysioData
%%%%%%%%%%%%%%%%%%%%%%%
function [cardio,allCardioPeakTimes,resp,allRespPeakTimes]=getPeaksFromPhysioData(cardioFilePath,respFilePath,tr,nvol,nslice)
% physio
physioFs=100; % actually 5 replications of 100 Hz = 500 Hz, but replications removed
physioSamples=physioFs*nvol*tr;
% cardio
cardio=[];
if exist(cardioFilePath,'file')
    cardioRaw=reshape(textread(cardioFilePath),5,[]);
    cardioRaw=cardioRaw(1,:);
    cardio=resample(cardioRaw,nvol*nslice,physioSamples); %#ok<*DTXTRD>
    % resp
    respRaw=reshape(textread(respFilePath),5,[]);
    respRaw=respRaw(1,:);
    resp=resample(respRaw,nvol*nslice,physioSamples);
end
[~,locs]=findpeaks(cardio,'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
allCardioPeakTimes=(locs'-1)*tr/nslice;
[~,locs]=findpeaks(resp,'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
allRespPeakTimes=(locs'-1)*tr/nslice;
% figure
% plot(cardio);
% hold on
% plot(locs,pks,'r*');
% hold off
end

function equalizedMask=makeEqualizedMask(restS,abMag,abPhase,nrMaskVoxPerSlice)
% create new varMask with even phase dsn (or pi/2 spread) & even amplitude
% dsns across phase
% partition phase dsn into 8 bins
% fill bins in order of mag until 3 neighboring (or anti-neighboring) bins have min N voxels
% equalize mag in each phase bin
% phase dsn
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
nrBins=12;
equalizedMask=false(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
    repAbPhase=repmat(abPhase(:,sliceNr,1),1,nrBins);
    phaseRange=-pi:2*pi/nrBins:pi*(nrBins-1)/nrBins;
    phaseBins=repmat(phaseRange,prod(volsize(1:2)),1);
    % make prelim mask to use vox with max abMag (rather than max var)
    [abMagSize,abMagSortIdx]=sort(abMag(:,sliceNr,1),1,'descend');
    if any(abMagSize(1:nrMaskVoxPerSlice,:)==0)
        disp('WARNING: mask has voxels with abMag = 0');
    end
    sliceMask=ismember((1:prod(volsize(1:2)))',abMagSortIdx(1:nrMaskVoxPerSlice));
    
    slicePhases=repAbPhase>=phaseBins & repAbPhase<(phaseBins+2*pi/nrBins) & repmat(~all(reshape(restS.vol(:,:,sliceNr,:)==0,[],nvol),2),1,nrBins)...
        & repmat(sliceMask,1,nrBins);
    sliceMag=(slicePhases(:,1:nrBins/2)|slicePhases(:,nrBins/2+1:end)).*repmat(squeeze(abMag(:,sliceNr,1)),1,nrBins/2);
    [sliceMagSize,sliceMagSortIdx]=sort(sliceMag,1,'descend');
    % choose pi/2 range with max rms
    maxRms=0;
    nrVoxPerPhaseBin=5;
    while maxRms==0
        shiftIdxs=logical([1,1,1,0,0,0]);
        for shiftNr=0:nrBins/2-1
            binShift=logical(circshift(shiftIdxs,shiftNr,2));
            magBins=sliceMagSize(:,binShift);
            magBins=magBins(1:nrVoxPerPhaseBin,:);
            magBinRms=sqrt(mean(magBins(:).^2));
            if all(magBins(:)~=0) && magBinRms>maxRms
                maxRms=magBinRms;
                maxRmsBinShift=binShift;
                maxMagIdxs=sliceMagSortIdx(:,maxRmsBinShift);
                equalizedMask(:,sliceNr)=ismember((1:prod(volsize(1:2)))',reshape(maxMagIdxs(1:nrVoxPerPhaseBin,:),[],1));
            end
        end
        if maxRms==0
            nrVoxPerPhaseBin=nrVoxPerPhaseBin-1;
            if nrVoxPerPhaseBin==0
                disp('could not make equalizedMask');
                keyboard
            end
        end
    end
end
end

function calcAliases(f,tr) %#ok<DEFNU>
% % fa=abs(f-N*fs)
fs=1/tr;
N=ceil(f/fs-0.5):floor(f/fs+0.5);
fa=abs(f-N*fs);
disp(fa);
% df=1/(nvol*tr);
% fScan=0:df:(nvol-1)*df;
% faIdx=round(fa/df)+1;
% faIdx=nvol/2+1;
end