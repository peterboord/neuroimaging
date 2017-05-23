function cardioDotProduct(session)

dbstop if error
if nargin==0
    session='RC4103-1';
    %session='RC4107-2';
    %session='RC4109-1';
end
disp(session);
%fmri
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
% constants
tr=2.4;

physioDir='/projects2/udall/physio';

restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_cardioAliasBpf.nii'));
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));

restS.tr=tr;
maskPrefix=fullfile(funcDir,session,[session,'_rest']);
if ~exist([maskPrefix,'_mask.nii'],'file')
    system(['bet ',restS.fspec,' ',maskPrefix,' -f 0.2 -m -n']);
end
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
cardio=getPeaksFromPhysioData(fullfile(physioDir,session,'rest_cardio.txt'),fullfile(physioDir,session,'rest_resp.txt'),tr,nvol,nslice);
maskS=MRIread([maskPrefix,'_mask.nii']);
restS.vol=reshape(detrend(bsxfun(@times,reshape(restS.vol,[],nvol),reshape(maskS.vol,[],1))')',size(restS.vol));

% nrTr=1;
% padFactor=4;
% nrMaskVoxPerSlice=50;
% maxCardioPkFreq=findPhysioPeakFreq('cardio',restS,nrTr,padFactor,nrMaskVoxPerSlice);
% bpmSideRange=10;
% cardFreqRange=maxCardioPkFreq-bpmSideRange/60:1/60:maxCardioPkFreq+bpmSideRange/60;
% nrCardFreq=numel(cardFreqRange);
% cardAliases=zeros(nrCardFreq,3);
% for cardFreqIdx=1:nrCardFreq
%     [aliasFreq,N,oddity]=calcAliases(cardFreqRange(cardFreqIdx),tr);
%     cardAliases(cardFreqIdx,:)=[aliasFreq,N,oddity];
%     disp([aliasFreq,N,oddity]);
% end
% disp({'min alias freq is ',min(cardAliases,1)});
% disp({'fslmaths bptf hpf sigma is ',(1/min(cardAliases,1))/tr});

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

% rest=rest./repmat(reshape(voxMask,prod(volsize(1:2)),nslice,1),1,1,nvol);
% rest(isinf(rest))=0;
% rest(isnan(rest))=0;

%sliceRangeToAv={1:round(nslice/4),round(nslice/4):round(2*nslice/4),round(2*nslice/4):round(3*nslice/4),round(3*nslice/4):round(4*nslice/4)};
sliceRangeToAv={1:nslice};
phaseZeroAv=zeros([prod(volsize(1:2)),nslice,nvol]);
%phasePiAv=zeros([prod(volsize(1:2)),nslice,nvol]);
for sliceRangeNr=1:numel(sliceRangeToAv)
    %restToCorr=reshape(restS.vol(:,:,sliceRangeToAv{sliceRangeNr},:),[],nvol);

    restToCorr=bsxfun(@times,reshape(restS.vol(:,:,sliceRangeToAv{sliceRangeNr},:),[],nvol),reshape(voxMask>0,[],1));
    %restToCorr=bsxfun(@times,reshape(rest(:,sliceRangeToAv{sliceRangeNr},:),[],nvol),reshape(voxMask>0,[],1));
    
    restRestCorr=corr(restToCorr,restToCorr);
    [sortRestRest,sortRestRestIdx]=sort(restRestCorr,'descend');
    nrVolToAverage=round(nvol/nslice);
    for volNr=1:nvol
        %phaseZeroAv(:,sliceRangeToAv{sliceRangeNr},volNr)=reshape(mean(restToCorr(:,sortRestRestIdx(1:nrVolToAverage,volNr)),2),[],numel(sliceRangeToAv{sliceRangeNr}));
        phaseZeroAv(:,sliceRangeToAv{sliceRangeNr},volNr)=reshape(median(restToCorr(:,sortRestRestIdx(1:nrVolToAverage,volNr)),2),[],numel(sliceRangeToAv{sliceRangeNr}));
        %phasePiAv(:,sliceRangeToAv{sliceRangeNr},volNr)=reshape(mean(restToCorr(:,sortRestRestIdx(end-nrVolToAverage+1:end,volNr)),2),[],numel(sliceRangeToAv{sliceRangeNr}));
    end
end
sqrtSumSliceVolPwr=sqrt(squeeze(sum(reshape(phaseZeroAv.^2,[],nslice,nvol),1)));
normSliceVolPwr=bsxfun(@rdivide,sqrtSumSliceVolPwr,mean(sqrtSumSliceVolPwr,2));
x=1:nslice*10;
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
end