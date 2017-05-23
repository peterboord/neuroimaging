function cardioErp(session)

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

cardioPercentSliceVox=10;

padFactor=4;
physioDir='/projects2/udall/physio';
restS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data.nii'));
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
[~,allCardioPeakTimes,~]=getPeaksFromPhysioData(fullfile(physioDir,session,'rest_cardio.txt'),fullfile(physioDir,session,'rest_resp.txt'),tr,nvol,nslice);
allCardioPeakSlices=round(allCardioPeakTimes*nslice/tr+1);
allCardioPeakSlices(allCardioPeakSlices<nslice+1)=[];
allCardioPeakSlices(allCardioPeakSlices>nvol*nslice-nslice)=[];
nrBeats=numel(allCardioPeakSlices);
nrErpSamp=2*nslice+1;
pkErpVol=cell(nslice,nrErpSamp);
pkErpSlice=cell(nslice,nrErpSamp);
pkErpSamp=cell(nslice,nrErpSamp);
for beatNr=1:nrBeats
    for sampNr=1:nrErpSamp
        sliceNr=mod(allCardioPeakSlices(beatNr)-nslice+sampNr-1,nslice);
        if sliceNr==0,sliceNr=nslice;end
        volNr=floor((allCardioPeakSlices(beatNr)-nslice+sampNr-1)/nslice)+1;
        if isempty(pkErpVol{sliceNr,sampNr})
            pkErpVol{sliceNr,sampNr}=volNr;
            pkErpSlice{sliceNr,sampNr}=allCardioPeakSlices(beatNr);
            pkErpSamp{sliceNr,sampNr}=allCardioPeakSlices(beatNr)-nslice+sampNr-1;
        else
            pkErpVol{sliceNr,sampNr}=[pkErpVol{sliceNr,sampNr},volNr];
            pkErpSlice{sliceNr,sampNr}=[pkErpSlice{sliceNr,sampNr},allCardioPeakSlices(beatNr)];
            pkErpSamp{sliceNr,sampNr}=[pkErpSamp{sliceNr,sampNr},allCardioPeakSlices(beatNr)-nslice+sampNr-1];
        end
    end
end
pkErp=zeros([volsize,nrErpSamp]);
%restS.vol=reshape(bsxfun(@rdivide,reshape(restS.vol,[],nvol),std(reshape(restS.vol,[],nvol),1,2)),[volsize,nvol]);
for sliceNr=1:nslice
    for sampNr=1:nrErpSamp
        pkErp(:,:,sliceNr,sampNr)=squeeze(mean(restS.vol(:,:,sliceNr,pkErpVol{sliceNr,sampNr}),4));
    end
end
MRIsave(restS,pkErp,fullfile(fileparts(restS.fspec),'pkErp.nii'),nrErpSamp);
pkErp=reshape(pkErp,[],nrErpSamp);
maxCardioPkFreq=findPhysioPeakFreq('cardio',restS,nrTr,padFactor,cardioPercentSliceVox);
t=(-nslice:nslice)*tr/nslice;
cosPhase=cos(2*pi*maxCardioPkFreq*t);
sinPhase=sin(2*pi*maxCardioPkFreq*t);
pkErpInPhase=sum(bsxfun(@times,pkErp,cosPhase),2)/sum(cosPhase.^2);
pkErpInQuad=sum(bsxfun(@times,pkErp,sinPhase),2)/sum(sinPhase.^2);
pkErpMag=sqrt(pkErpInPhase.^2+pkErpInQuad.^2);
pkErpPhase=atan2(pkErpInQuad,pkErpInPhase);
MRIsave(restS,reshape(pkErpPhase,volsize),fullfile(fileparts(restS.fspec),'pkErpPhase.nii'),nrErpSamp);
MRIsave(restS,reshape(pkErpMag,volsize),fullfile(fileparts(restS.fspec),'pkErpMag.nii'),nrErpSamp);
end
