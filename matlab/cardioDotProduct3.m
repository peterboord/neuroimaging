function cardioDotProduct3(session)

dbstop if error
tic
if nargin==0
    session='109003';
%     session='109006';
%     session='109010';
%     session='109012';
%     session='109013';
%     session='109014';
%     session='109015';
%     session='RC4103-1';
%     session='RC4107-2';
%     session='RC4109-1';
end
switch session
    case {'RC4103-1','RC4107-2','RC4109-1'}
        tr=2.4;
        funcDir='/projects2/udall/pboord/pic/preproc/pestica';
    otherwise
        tr=2.5;
        funcDir='/project_space/pboord/act/rest';
        filePeaks=load(fullfile(funcDir,session,[session,'_rest_cardio_peaks.txt']));
end
disp(session);
%fmri
%restS=MRIread(fullfile(funcDir,session,'mocoHpf.feat','filtered_func_data.nii'));
restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,session,'rest_e002.nii.gz'));
%restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_cardioAliasBpf.nii'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));
% maskS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','mask.nii'));
% mask=maskS.vol(:);

pvesegS=MRIread(fullfile(funcDir,session,'T1_brain_pveseg_2_to_rest.nii'));
%csf=ismember(pvesegS.vol(:),1);
gm=ismember(pvesegS.vol(:),2);
wm=ismember(pvesegS.vol(:),3);
restS.vol=reshape(bsxfun(@times,reshape(restS.vol,[],restS.nframes),~(gm)),size(restS.vol));
% mask=mask & ~gm;
% MRIsave(maskS,reshape(mask,maskS.volsize),fullfile(fileparts(maskS.fspec),'mask_nogm.nii'),1);
% mask=mask & ~wm;
% MRIsave(maskS,reshape(mask,maskS.volsize),fullfile(fileparts(maskS.fspec),'mask_nogmwm.nii'),1);

%restS.vol=reshape(detrend(bsxfun(@times,reshape(restS.vol,[],restS.nframes),mask)')',size(restS.vol));

restS.tr=tr;
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
[cardio,allCardioPeakTimes]=getPeaksFromPhysioData(fullfile(funcDir,session,'rest_cardio_ds.txt'),fullfile(funcDir,session,'rest_resp_ds.txt'),tr,nvol,nslice);
allCardioPeakTimes=filePeaks;
slicePhases=peaksToSlicePhase(allCardioPeakTimes,tr,nslice,nvol);
pkSlices=round(allCardioPeakTimes*nslice/tr);
pkSlicesRound=round(pkSlices);
pkSlices=mod(pkSlices,nslice)+1;
% figure,hist(round(pkSlices),0.5:1:nslice);
% figure,hist(pkSlices,0.5:1:nslice);
percentSliceVox=1;
rest=reshape(restS.vol,[],nslice,nvol);
stdRest=std(rest,1,3);
voxMask=zeros(prod(volsize(1:2)),nslice);
% [stdRestSort,stdRestSortIdx]=sort(stdRest(:),1,'descend');
% nrVox=sum(sum(~all(rest==0,3),1));
% topvarIdxs=stdRestSortIdx(1:round(percentSliceVox*nrVox/100));
% voxMask(topvarIdxs)=1;
[~,stdRestSortIdx]=sort(stdRest,1,'descend');
nrSliceVox=sum(~all(rest==0,3),1);
for sliceNr=1:nslice
    topvarIdxs=stdRestSortIdx(1:max([1,round(percentSliceVox*nrSliceVox(sliceNr)/100)]),:);
    %topvarIdxs=stdRestSortIdx(1,:);
    voxIdxs=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
    voxMask(voxIdxs,sliceNr)=stdRest(voxIdxs,sliceNr);
end
% voxMask=stdRest;
phaseZeroAv=zeros([prod(volsize(1:2)),nslice,nvol]);
%phasePiAv=zeros([prod(volsize(1:2)),nslice,nvol]);
%restToCorr=bsxfun(@times,reshape(restS.vol,[],nvol),reshape(voxMask(:)>0,[],1));
%restToCorr=reshape(restToCorr,[],nslice,nvol);
restToCorr=reshape(restS.vol,[],nslice,nvol);
restRestPwr=zeros(nvol,nvol,nslice);
rest1Pwr=zeros(nvol,nslice);
restRestCorr=zeros(nvol,nvol,nslice);

firstSliceVoxels=restToCorr(:,1,1)~=0;
lastSliceVoxels=restToCorr(:,end,1)~=0;
commonSliceVoxels=firstSliceVoxels & lastSliceVoxels;
commonRest=restToCorr(commonSliceVoxels,:,:);
commonRest=reshape(commonRest,[],nvol);
commonRest=detrend(commonRest')';
commonRest=bsxfun(@rdivide,commonRest,std(commonRest,1,2)); commonRest(isnan(commonRest))=0;
commonRest=reshape(commonRest,[],nslice,nvol);
sliceVol=zeros(nslice,nvol);
for volNr=1:nvol
    for sliceNr=1:nslice
        sliceVol(sliceNr,volNr)=corr(commonRest(:,1,1),commonRest(:,sliceNr,volNr));
        %sliceVol(sliceNr,volNr)=commonRest(:,1,1)'*commonRest(:,sliceNr,volNr);
    end
end
sliceVol=bsxfun(@minus,sliceVol,mean(sliceVol,2));
x=(1:nslice*3)';
figure,plotyy(x,sliceVol(x),x,cardio(x));

stdSortedRest=zeros(size(restToCorr));
for sliceNr=1:nslice
    sliceData=squeeze(restToCorr(:,sliceNr,:));
    [stdRestToCorr,sortIdx]=sort(std(sliceData,0,2),1,'descend');
    stdSortedRest(:,sliceNr,:)=sliceData(sortIdx,:);
end
commonSliceVoxels=all(reshape(stdSortedRest,[],nslice*nvol)~=0,2);
commonStdSortedRest=stdSortedRest(commonSliceVoxels,:,:);
normRest=reshape(commonStdSortedRest,[],nvol);
normRest=bsxfun(@rdivide,normRest,std(normRest,1,2));
normRest=reshape(normRest,size(commonStdSortedRest));
sliceVol=zeros(nslice,nvol);
volVolSlice=zeros(nvol,nvol,nslice);
for volNr=1:nvol
    for sliceNr=1:nslice
        volVolSlice(:,:,sliceNr)=corr(squeeze(normRest(:,sliceNr,:)));
        %sliceVol(sliceNr,volNr)=corr(normRest(:,1,1),normRest(:,sliceNr,volNr));
        %sliceVol(sliceNr,volNr)=corr(normRest(:,1,volNr),normRest(:,sliceNr,volNr));
    end
end
meanVolVolSlice=squeeze(mean(volVolSlice,2))';
meanVolVolSlice=bsxfun(@minus,meanVolVolSlice,mean(meanVolVolSlice,2));
%sliceVol=bsxfun(@rdivide,sliceVol,mean(sliceVol,2));
%sliceVol=detrend(sliceVol','constant')';
x=(1:nslice*10)';
figure,plotyy(x,meanVolVolSlice(x),x,cardio(x));
[r,lags]=xcorr(meanVolVolSlice(x),cardio(x),'coeff');
%[r,lags]=xcorr(sliceVol(x),sliceVol(x),'coeff');
%[r,lags]=xcorr(sliceVol(:),cardio(:),'coeff');
figure,
plot(lags,r);

for sliceNr=1:nslice
    %restRestCorr(:,:,sliceNr)=corr(squeeze(restToCorr(:,sliceNr,:)),squeeze(restToCorr(:,sliceNr,:)));
%     b=(squeeze(restToCorr(:,sliceNr,:))'*squeeze(restToCorr(:,sliceNr,:)));
%     b=bsxfun(@rdivide,b,sum(b.^2,1));
%     
%     sliceVol=zeros(nslice,nvol);
%     phaseSteps=-pi:pi/100:pi;
%     nrSteps=numel(phaseSteps);
%     ms=zeros(nslice,nrSteps);
%     for stepNr=1:nrSteps
%         for L=1:nslice
%             for volNr=1:nvol
%                 sliceVol(:,volNr)=b(volNr,1)*cos(phaseSteps(stepNr)+2*pi*(0:nslice-1)/L);
%             end
%             ms(L,stepNr)=max(abs(diff(diff(reshape(sliceVol,[],1)))));
%         end
%     end
    a=(squeeze(restToCorr(:,sliceNr,:))'*squeeze(restToCorr(:,sliceNr,:)));
    restRestPwr(:,:,sliceNr)=a./sum(squeeze(restToCorr(:,sliceNr,1)).^2);
    rest1Pwr(:,sliceNr)=bsxfun(@rdivide,(squeeze(restToCorr(:,1,1))'*squeeze(restToCorr(:,sliceNr,:)))./sqrt(sum(squeeze(restToCorr(:,1,1)).^2)),sqrt(sum(squeeze(restToCorr(:,sliceNr,:)).^2,1)));
    restRestCorr(:,:,sliceNr)=corr(squeeze(restToCorr(:,sliceNr,:)));
%     rmsRestRestCorr=sqrt(mean(restRestCorr(:,:,sliceNr).^2,2));
%     normRmsRestRestCorr(:,sliceNr)=rmsRestRestCorr/sqrt(mean(rmsRestRestCorr.^2));
end
physioEst=restRestPwr(:,1,1);
MRIsave(restS,reshape(corr(reshape(restS.vol,[],nvol)',physioEst),volsize),fullfile(fileparts(restS.fspec),'restCorrPhysioEst.nii'),1);
save(fullfile(fileparts(restS.fspec),'physioEst.txt'),'physioEst','-ascii');

for i=1:nslice
    figure('WindowStyle','docked'),plot(restRestCorr(:,1,i));
end
% figure('WindowStyle','docked'),plot(restRestCorr(:,1,1))
figure('WindowStyle','docked'),plot(restRestCorr(:,1,end-1))
% figure('WindowStyle','docked'),plot(cos(slicePhases(1,:)))
% figure('WindowStyle','docked'),plot(cos(slicePhases(end,:)))
sliceLags=zeros(nslice,1);
sliceCorr=zeros(nslice,1);
for sliceNr=1:nslice
    [r,lags]=xcorr(slicePhases(1,:)',slicePhases(sliceNr,:)','coeff');
    [maxCorr,maxIdx]=max(r);
    sliceLags(sliceNr)=lags(maxIdx);
    sliceCorr(sliceNr)=maxCorr;
end
figure,plot(sliceLags)
figure,plot(sliceCorr)
for sliceNr=1:nslice
    figure('WindowStyle','docked'),plot(corr(restRestCorr(:,1,sliceNr),squeeze(restRestCorr(:,1,:))));
end
sliceVol=zeros(nslice,nvol);
for sliceNr=1:nslice
    sliceUpsample=reshape(interp(restRestCorr(:,1,sliceNr),nslice),nslice,[]);
    sliceVol(sliceNr,:)=sliceUpsample(sliceNr,:)./sqrt(median(sliceUpsample(sliceNr,:).^2));
end
restRestCorrShift=restRestCorr;
restRestPwrShift=restRestPwr;
for volNr=1:nvol
    restRestCorrShift(:,volNr,:)=circshift(restRestCorr(:,volNr,:),-(volNr-1),1);
    restRestPwrShift(:,volNr,:)=circshift(restRestPwr(:,volNr,:),-(volNr-1),1);
end
sliceNr=1;
volNr=3;
a=restRestPwr(:,:,sliceNr);
volSqEst=zeros(nvol*(nvol+1)/2,1);
idx=1;
for vol1=1:nvol
    for vol2=vol1:nvol
        volSqEst(idx)=a(volNr,vol1)*a(volNr,vol2)/a(vol1,vol2);
        idx=idx+1;
    end
end
disp(median(volSqEst))

figure
plot(restRestPwr(:,1,2));
hold on
volLag=zeros(nvol,1);
for volNr=2:nvol
   [r,lags]=xcorr(restRestPwr(:,volNr-1,1),restRestPwr(:,volNr,1),'coeff');
   [~,maxIdx]=max(r);
   plot(circshift(restRestPwr(:,volNr,1),-lags(maxIdx)));
   volLag(volNr)=lags(maxIdx);
end
hold off
volNr1=1;
volNr2=157;
sliceNr=1;
% volSlice=squeeze(restRestCorr(:,volNr,:));
% plot(bsxfun(@rdivide,volSlice,sqrt(mean(volSlice.^2,1))))
figure
subplot(3,1,1)
plot(restRestPwr(:,volNr1,sliceNr));
subplot(3,1,2)
plot(restRestPwr(:,volNr2,sliceNr));
subplot(3,1,3)
[r,lags]=xcorr(restRestPwr(:,volNr1,sliceNr),restRestPwr(:,volNr2,sliceNr),'coeff');
plot(lags,r);

[sortRestRest,sortRestRestIdx]=sort(restRestPwr,'descend'); %#ok<ASGLU>
sortBy='threshold';
switch sortBy
    case 'threshold'
        corMin=0.7;
        for volNr=1:nvol
            phaseZeroAv(:,:,volNr)=reshape(median(restToCorr(:,restRestPwr(:,volNr)>=corMin),2),[],nslice);
        end
    case 'nrVolToAverage'
        nrVolToAverage=round(nvol/nslice);
        for volNr=1:nvol
            phaseZeroAv(:,:,volNr)=reshape(median(restToCorr(:,sortRestRestIdx(1:nrVolToAverage,volNr)),2),[],nslice);
        end
end
sqrtSumSliceVolPwr=sqrt(squeeze(sum(reshape(phaseZeroAv.^2,[],nslice,nvol),1)));
normSliceVolPwr=bsxfun(@rdivide,sqrtSumSliceVolPwr,mean(sqrtSumSliceVolPwr,2));
x=1:nslice*20;
t=(x-1)*tr/nslice;
figure
subplot(2,1,1);
plot(t,detrend(cardio(x),'constant')/max(detrend(cardio(x),'constant')),t,detrend(normSliceVolPwr(x),'constant')/max(detrend(normSliceVolPwr(x))));
subplot(2,1,2);
[r,lags]=xcorr(cardio(:),normSliceVolPwr(:),'coeff');
[rMax,maxIdx]=max(r);
disp({'max corr ',rMax,' at ',maxIdx,lags(maxIdx),lags(maxIdx)*tr/nslice});
plot(lags,r);
% [maxPhaseZero,maxPhaseZeroIdx]=max(phaseZeroAv,[],2);
% [minPhaseZero,minPhaseZeroIdx]=min(phaseZeroAv,[],2);
toc
end