function allCardioPeakTimes=findCardioPeaks(session)

tic
dbstop if error
if nargin==0
    session='RC4107-2';
    %session='RC4103-1';
    %session='RC4109-1';
end
disp(session);
%fmri
%maskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
%funcDir='/projects2/udall/rest/raw';
% 
tr=2.4;
%restS=MRIread(fullfile(funcDir,session,'rest_brain_reg.nii'));
restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],'filtered_func_data.nii'));

%restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],[session,'.results'],['pb01.',session,'.r01.despike.nii']));

%restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));
% differentiate image in phase-encoding (pa) direction
% [dlrMeanFunc,dpaMeanFunc]=gradient(mean(restS.vol,4));
% dpaMeanFunc=abs(dpaMeanFunc);
% dlrMeanFunc=abs(dlrMeanFunc);
% maskS=resizeFuncS(restS,3);
% maskS.vol=reshape(dpaMeanFunc,restS.volsize);
% maskS.fspec=fullfile(funcDir,session,'dpaMeanFunc.nii');
% MRIwrite(maskS,maskS.fspec);
%restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],'filtered_func_data.nii'));
%restS=MRIread(fullfile(funcDir,[session,'_rest.nii.gz']));
%maskPath=fullfile(funcDir,session,maskFile);
% if ~exist(maskPath,'file')
%     system(['cd ',fullfile(funcDir,session),...
%         ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',maskFile]);
% end
%maskS=MRIread(maskPath);
nvol=restS.nframes;
volsize=restS.volsize;
nslice=volsize(3);
rest=reshape(restS.vol,prod(volsize),nvol)';
rest=detrend(rest);
minCardioHz=0.6;
maxCardioHz=1.5;
% restFft=abs(fft(rest));
% restFft=restFft(1:end/2+1,:);
% restFft=bsxfun(@rdivide,restFft,sum(restFft,1));
% restFft(isnan(restFft))=0;
% restFft(isinf(restFft))=0;
% df=1/(nvol*tr);
% f=0:df:(nvol*nslice-1)*df;
% minCardio=round(minCardioHz/df)+1;
% maxCardio=round(maxCardioHz/df)+1;
   

rest=reshape(rest',prod(volsize(1:2)),nslice,nvol);
    
% make mask from 50 voxels/slice with highest Tstd
nrMaskVoxPerSlice=50;
%nrMaskVoxPerSlice=40;
varRest=var(rest,1,3);
%[varRestSize,varRestSortIdx]=sort(varRest,1,'descend');
% maskCostFn=varRest;%./reshape(dpaMeanFunc,[],nslice);%;%.*reshape(dlrMeanFunc,[],nslice);%
% maskCostFn(isnan(maskCostFn))=0;
% maskCostFn(isinf(maskCostFn))=0;
[varRestSize,varRestSortIdx]=sort(varRest,1,'descend');
if any(reshape(varRestSize(1:nrMaskVoxPerSlice,:),[],1)==0)
    disp('WARNING: mask has voxels with Tvar = 0');
end
topvarIdxs=varRestSortIdx(1:nrMaskVoxPerSlice,1:nslice);
varMask=false(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
   varMask(:,sliceNr)=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
end
% maskS=resizeFuncS(restS,3);
% maskS.vol=reshape(varMask,restS.volsize);
% maskS.fspec=fullfile(funcDir,session,'varMask.nii');
% MRIwrite(maskS,maskS.fspec);
%mask=logical(reshape(maskS.vol,prod(volsize(1:2)),nslice));
sqssSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(varMask(:,sliceNr),sliceNr,:));
    sqssSlice(sliceNr,:)=sqrt(sum(cohRestSlice.^2,1));
    % normalize 
    sqssSlice(sliceNr,:)=sqssSlice(sliceNr,:)/std(sqssSlice(sliceNr,:));
end
sqssSlice=detrend(sqssSlice')';
nrTr=2;
dfSlice=1/(nrTr*tr);
padFactor=4;
disp(padFactor);
% make mask from nrMaskVoxPerSlice voxels with top power in cardio band
%[Pxx,f]=pwelch(sqssSlice(:),nrTr*nslice,nslice,minCardioHz:1/(nrTr*tr):maxCardioHz,nslice/tr);
[Pxx,f]=pwelch(sqssSlice(:),nrTr*nslice,nslice,padFactor*nrTr*nslice,nslice/tr);
%[Pxx,f]=pwelch(sqssSlice(:),nrTr*nslice,nslice,nrTr*nslice,nslice/tr);
% find peak in cardio range
minCardioIdx=floor(minCardioHz*padFactor/dfSlice);
maxCardioIdx=ceil(maxCardioHz*padFactor/dfSlice)+1;
cardioPeakIdxs=[0;Pxx(2:end-1)>Pxx(1:end-2) & Pxx(2:end-1)>Pxx(3:end);0];
cardioPeakIdxs(1:minCardioIdx-1)=0;
cardioPeakIdxs(maxCardioIdx+1:end)=0;
[~,maxCardioPkIdx]=max(Pxx.*cardioPeakIdxs);
maxCardioPkFreq=(maxCardioPkIdx-1)*dfSlice/padFactor;
fs=1/tr;
% fa=abs(f-N*fs)
% fa > 0 & fa < 0.5/tr
% abs(f-N*fs)>0
% N > f/fs
N=ceil(maxCardioPkFreq/fs-0.5):floor(maxCardioPkFreq/fs+0.5);
fa=abs(maxCardioPkFreq-N*fs);
df=1/(nvol*tr);
faIdx=round(fa/df);
restFft=fft(reshape(rest,[],nvol)');
cardioAliasPowerUnsorted=reshape(abs(restFft(faIdx,:))',[],nslice);
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

sqssSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(cardioMask(:,sliceNr),sliceNr,:));
    sqssSlice(sliceNr,:)=sqrt(sum(cohRestSlice.^2,1));
    % normalize 
    sqssSlice(sliceNr,:)=sqssSlice(sliceNr,:)/std(sqssSlice(sliceNr,:));
end
sqssSlice=detrend(sqssSlice')';

%(0.5/tr-maxCardioPkFreq)/fs
%abs(maxCardioPkFreq-
% figure
% plot(f,Pxx)

% [Pmt,Fmt]=pmtm(sqssSlice(:),4,numel(sqssSlice(:)),nslice/tr);
% figure,
% plot(Fmt,Pmt);

%figure
allCardioPeakTimes=[];
allSqssSliceFft=[];
for volNr=1:nvol-nrTr+1
    x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
    sqssWave=reshape(sqssSlice(x),[],1);
    %     figure
    %     plot(x,sqssWave);
    sqssSliceFft=fft(sqssWave,padFactor*numel(sqssWave));
    if isempty(allSqssSliceFft)
        allSqssSliceFft=sqssSliceFft;
    else
        allSqssSliceFft=[allSqssSliceFft,sqssSliceFft]; %#ok<AGROW>
    end
    %sqssSliceFft=fft(sqssWave.*hann(numel(sqssWave),'periodic'),padFactor*numel(sqssWave));
    fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
    subplot(2,1,1); plot(sqssWave);
    subplot(2,1,2); plot(fSlice(1:end/2),abs(sqssSliceFft(1:end/2)));
    keyboard
    minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
    maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
    [~,maxIdx]=max(abs(sqssSliceFft(minCardio:maxCardio)));
    maxSpecIdx=minCardio+maxIdx-1;
    maxSpecAngle=angle(sqssSliceFft(maxSpecIdx));
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
    peakTimes(peakTimes>=tr)=[];
%     if numel(allCardioPeakTimes)>1
%         curBeats=diff([allCardioPeakTimes(end);(volNr-1)*tr+peakTimes]);
%         if (any(curBeats<1) || any(curBeats>2))
%             keyboard
%         end
%     end
    allCardioPeakTimes=cat(1,allCardioPeakTimes,(volNr-1)*tr+peakTimes);
    %peakIdxs=round(peakTimes*nslice/tr)+1;
end
plot(fSlice(1:end/2),mean(abs(allSqssSliceFft(1:end/2,:)),2))
toc
