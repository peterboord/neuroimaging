function findCardioPeaksV1_4(session)

tic
dbstop if error
if nargin==0
    session='RC4103-1';
    %session='RC4107-2';
    %session='RC4109-1';
end
disp(session);
%figure('WindowStyle','docked');
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
cardioPercentSliceVox=10;
respPercentSliceVox=100;

minCardioHz=0.6;
maxCardioHz=1.5;
padFactor=4;
physioDir='/projects2/udall/physio';
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.nii.gz']));
restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,session,'prefiltered_func_data_mcf.nii.gz'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));
%restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],[session,'.results'],['pb01.',session,'.r01.despike.nii']));
restS.tr=tr;
maskPrefix=fullfile(funcDir,session,[session,'_rest']);
if ~exist([maskPrefix,'_mask.nii'],'file')
    system(['bet ',restS.fspec,' ',maskPrefix,' -f 0.2 -m -n']);
end
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
maskS=MRIread([maskPrefix,'_mask.nii']);
restS.vol=reshape(detrend(bsxfun(@times,reshape(restS.vol,[],nvol),reshape(maskS.vol,[],1))')',size(restS.vol));
% derived constants
disp(padFactor);
% pulse ox

nrTr=2;
% rest=reshape(restS.vol,[],nslice,nvol);
% % respMask=std(reshape(restS.vol,[],nslice,nvol),1,3);
% % rest=rest./repmat(reshape(respMask,prod(volsize(1:2)),nslice,1),1,1,nvol);
% % rest(isinf(rest))=0;
% % rest(isnan(rest))=0;
% sumRestSquared=squeeze(sum(rest.^2,1));
% rmsSlice=bsxfun(@rdivide,sumRestSquared,median(sumRestSquared,2));
% % [Pxx,f]=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
% % [Pxx,f]=pwelch(rmsSlice(:),nslice,0,nrTr*nslice,nslice/tr);
% fftRmsSlice=fft(detrend(rmsSlice,'constant'),nslice*padFactor);
% f=0:(1/tr)/padFactor:(size(fftRmsSlice,1)-1)*(1/tr)/padFactor;
% figure,surf(1:300,f(1:floor(end/2)+1),abs(fftRmsSlice(1:floor(end/2)+1,:)));
% figure,plot(f(1:floor(end/2)+1),mean(abs(fftRmsSlice(1:floor(end/2)+1,:)),2));

[maxCardioPkFreq,Pxx,f]=findPhysioPeakFreq('cardio',restS,nrTr,padFactor,cardioPercentSliceVox);
[maxRespPkFreq,Pxx,f]=findPhysioPeakFreq('resp',restS,nrTr,padFactor,respPercentSliceVox);

if 1
[cardio,allCardioPeakTimes,resp]=getPeaksFromPhysioData(fullfile(physioDir,session,'rest_cardio.txt'),fullfile(physioDir,session,'rest_resp.txt'),tr,nvol,nslice);
[clean4dPulseox,noise4dPulseox,~,~,slicePhasesPulseox]=retroicor(restS,allCardioPeakTimes);
end

% dirty data
rmsSlice=getRmsSlice2(restS.vol,cardioPercentSliceVox);

[maxCardioPkFreq,Pxx,f]=findSlicePhases4(restS,nrTr,padFactor,cardio);

[allCardioPeakTimes,Pxx,f]=calcPeakTimes(restS,rmsSlice,nrTr,padFactor);
% clean data
[~,noise4d,abMag,abPhase,slicePhases]=retroicor(restS,allCardioPeakTimes);
disp(corr(squeeze(noise4d(roi(1),roi(2),roi(3),:)),squeeze(noise4dPulseox(roi(1),roi(2),roi(3),:))));
clear noise4d
equalizedMask=makeEqualizedMask(restS,abMag,abPhase,nrMaskVoxPerSlice);
rmsSliceEqualized=getRmsSlice2(restS.vol,cardioPercentSliceVox,equalizedMask);
% re-calc peaks & clean data
[allCardioPeakTimes,PxxEqualized]=calcPeakTimes(restS,rmsSlice,nrTr,padFactor)
[clean4d,noise4dEqualized,~,~,slicePhasesEqualized]=retroicor(restS,allCardioPeakTimes);
disp(corr(squeeze(noise4dEqualized(roi(1),roi(2),roi(3),:)),squeeze(noise4dPulseox(roi(1),roi(2),roi(3),:))));
clear noise4dPulseox
MRIsave(restS,noise4dEqualized,fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst.nii'));
MRIsave(restS,sqrt(mean(noise4dEqualized.^2,4)),fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_noiseEst_rms.nii'),1);
clear noise4dEqualized
rmsSliceClean=getRmsSlice2(clean4d,cardioPercentSliceVox,equalizedMask);
PxxSliceClean=pwelch(rmsSliceClean(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
% figure('WindowStyle','docked');
% plot(1:nslice,diag(corr(cos(slicePhases),cos(slicePhasesPulseox))),...
%     1:nslice,diag(corr(cos(slicePhasesEqualized),cos(slicePhasesPulseox))),...
%     1:nslice,diag(corr(cos(slicePhases),cos(slicePhasesEqualized))));
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
plot(ax1,falias,mean(PaliasDirty,2),'b');
hold(ax1,'on');
plot(ax1,falias,mean(PaliasEqualized,2),'g');
plot(ax1,falias,mean(PaliasClean,2),'r');
hold(ax1,'off');
% unaliased slice data
ax2=subplot(nrPlots,1,2);
plot(ax2,f,Pxx,'b');
hold(ax2,'on');
plot(ax2,f,PxxEqualized,'g');
plot(ax2,f,PxxSliceClean,'r');
hold(ax2,'off');

% aliased roi data
if ~isempty(roi)
    ax3=subplot(nrPlots,1,3);
    rest=reshape(detrend(reshape(restS.vol,[],nvol)')',size(restS.vol));
    [Proi,froi]=pwelch(squeeze(rest(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(ax3,froi,Proi,'b');
    hold(ax3,'on');
    Proi=pwelch(squeeze(clean4d(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(ax3,froi,Proi,'r');
    Proi=pwelch(squeeze(clean4dPulseox(roi(1),roi(2),roi(3),:)),nrTr*nslice,nslice,nrTr*nslice,1/tr);
    plot(ax3,froi,Proi,'m');
    hold(ax3,'off');
end

%% CARDIAC FIGURES
figure('WindowStyle','docked');
nrPlots=2;
axWave=subplot(1,nrPlots,1);
axSpec=subplot(1,nrPlots,2);
padFactor=4;
xWave=(1:nslice*2);
tWave=(0:numel(xWave)-1)*tr/nslice;
plot(axWave,tWave,detrend(cardio(xWave),'constant')/max(detrend(cardio(xWave),'constant')),'r');
hold(axWave,'on');
plot(axWave,tWave,detrend(rmsSliceEqualized(xWave),'constant')/max(detrend(rmsSliceEqualized(xWave),'constant')),'b');
plot(axWave,tWave,detrend(rmsSlice(xWave),'constant')/max(detrend(rmsSliceEqualized(xWave),'constant')),'g');
hold(axWave,'off');
PxxCardio=pwelch(cardio(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
PxxRmsSlice=pwelch(rmsSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
plot(axSpec,f,PxxEqualized/max(PxxEqualized),'b');
hold(axSpec,'on');
plot(axSpec,f,PxxCardio/max(PxxCardio),'r');
plot(axSpec,f,PxxRmsSlice/max(PxxRmsSlice),'g');
hold(axSpec,'off');
% plots annotations
fontSize=20;
xlabel(axWave,'Time (s)','FontSize',fontSize);
ylabel(axWave,'Normalized waveform deflection','FontSize',fontSize);
xlabel(axSpec,'Frequency (Hz)','FontSize',fontSize);
ylabel(axSpec,'Normalized periodogram','FontSize',fontSize);
legend(axSpec,{'Estimated cardiac','Pulse oximeter'},'FontSize',fontSize);
xlim(axWave,[tWave(1),tWave(end)]);
xlim(axSpec,[f(1),f(end)]);
set(axWave,'FontSize',fontSize);
set(axSpec,'FontSize',fontSize);

%% RESPIRATION FIGURES
% plot spec for resp peak
restS=MRIread(fullfile(funcDir,session,[session,'_rest.nii.gz']));
restS.vol=reshape(detrend(bsxfun(@times,reshape(restS.vol,[],nvol),reshape(maskS.vol,[],1))')',size(restS.vol));
figure('WindowStyle','docked');
nrPlots=2;
axWave=subplot(1,nrPlots,1);
axSpec=subplot(1,nrPlots,2);
padFactor=8;

% n.b. resp spec better after cardio cleaned
respMask=std(reshape(restS.vol,[],nslice,nvol),1,3);
rmsSliceRespRaw=getRmsSlice2(restS.vol,respPercentSliceVox,respMask);
fftResp=fft(detrend(rmsSliceRespRaw,'constant'),nslice*padFactor);
% respMask=std(reshape(clean4d,[],nslice,nvol),1,3);
% rmsSliceRespClean=getRmsSlice2(clean4d,respPercentSliceVox,respMask);
% fftResp=fft(detrend(rmsSliceRespClean,'constant'),nslice*padFactor);

PxxResp=mean(abs(fftResp(1:round(end/2)+1,:)),2);
fResp=0:(1/tr)/padFactor:(numel(PxxResp)-1)*(1/tr)/padFactor;
plot(axSpec,fResp,PxxResp/max(PxxResp),'b');

% plot waveforms
respMask=std(reshape(restS.vol,[],nslice,nvol),1,3);
rmsSliceRespRaw=getRmsSlice2(restS.vol,respPercentSliceVox,respMask);
xWave=(1:nslice*15);
tWave=(0:numel(xWave)-1)*tr/nslice;
plot(axWave,tWave,detrend(resp(xWave),'constant')/max(detrend(resp(xWave),'constant')),'r');
hold(axWave,'on');
% plot(axWave,xWave,detrend(rmsSliceRespClean(xWave),'constant')/max(detrend(rmsSliceRespClean(xWave),'constant')),'b');
% plot(axWave,xWave,detrend(rmsSlice(xWave),'constant')/max(detrend(rmsSlice(xWave),'constant')),'b');
% plot(axWave,xWave,detrend(rmsSliceEqualized(xWave),'constant')/max(detrend(rmsSliceEqualized(xWave),'constant')),'g');
plot(axWave,tWave,detrend(rmsSliceRespRaw(xWave),'constant')/max(detrend(rmsSliceRespRaw(xWave),'constant')),'b');
hold(axWave,'off');
fftResp=fft(detrend(rmsSliceRespRaw,'constant'),nslice*padFactor);
PxxResp=mean(abs(fftResp(1:round(end/2)+1,:)),2);
fResp=0:(1/tr)/padFactor:(numel(PxxResp)-1)*(1/tr)/padFactor;
% hold(axSpec,'on');
% plot(axSpec,fResp,PxxResp/max(PxxResp),'g');
% hold(axSpec,'off');
fftResp=fft(detrend(reshape(resp,[],nvol),'constant'),nslice*padFactor);
PxxResp=mean(abs(fftResp(1:round(end/2)+1,:)),2);
hold(axSpec,'on');
plot(axSpec,fResp,PxxResp/max(PxxResp),'r');
hold(axSpec,'off');
% plots annotations
fontSize=20;
xlabel(axWave,'Time (s)','FontSize',fontSize);
ylabel(axWave,'Normalized waveform deflection','FontSize',fontSize);
xlabel(axSpec,'Frequency (Hz)','FontSize',fontSize);
ylabel(axSpec,'Normalized periodogram','FontSize',fontSize);
legend(axSpec,{'Estimated respiration','Pneumatic belt'},'FontSize',fontSize);
xlim(axWave,[tWave(1),tWave(end)]);
xlim(axSpec,[fResp(1),fResp(end)]);
set(axWave,'FontSize',fontSize);
set(axSpec,'FontSize',fontSize);

%%
% calc proportion of variance
varOld=var(restS.vol,0,4);
varNew=var(clean4d,0,4);
propVar=100*(varOld-varNew)./varOld;
propVar(isnan(propVar))=0;
propVar(isinf(propVar))=0;
% save 4D images
MRIsave(restS,clean4d,fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor.nii'));
% save 3D images
MRIsave(restS,propVar,fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_retroicor_propVar.nii'),1);
toc
end
