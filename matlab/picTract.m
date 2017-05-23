function picTract

if ~isdeployed
    dbstop if error
end
%% files & directories
%filePrefix='rest_medn_40_10';
filePrefix='rest_medn';
maskName='mask.nii.gz';
subjCodeNrStr='109003';
dirName='/project_space/pboord/act/rest';
%% constants
tr=2.5;
sz = 1;
minNeighNr=3;
%% seeds
seed(1).name = 'WM';
seed(2).name = 'WM2';
seed(3).name = 'WM3';
seed(4).name = 'WM4';
%% load data
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
funcS.vol=reshape(funcS.vol,[],funcS.nframes);
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));
maskInd=find(maskS.vol(:));
wmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_wm.nii.gz'));

for seedNr=1:numel(seed)
    seedFileName=['func_',seed(seedNr).name];
    seed(seedNr).vox=round(load(fullfile(dirName,subjCodeNrStr,seedFileName),'-ascii'));
    y=seed(seedNr).vox(1); % n.b. swapped
    x=seed(seedNr).vox(2); % x/y for FSL->FS
    z=seed(seedNr).vox(3);
    seed(seedNr).ind=sub2ind(funcS.volsize,x,y,z);
    seed(seedNr).mind=find(maskInd==seed(seedNr).ind);
end
%
funcS.vol=reshape(funcS.vol,[],funcS.nframes);
[pathstr,name,ext]=fileparts(funcS.fspec);
name=strtok([name,ext],'.');
filePrefixDir=fullfile(pathstr,name);
if ~exist(filePrefixDir,'dir')
    mkdir(filePrefixDir);
end
picDir=fullfile(filePrefixDir,'pic');
if ~exist(picDir,'dir')
    mkdir(picDir);
end
nf=funcS.nframes;
nrSeed=numel(seed);
maskInd=find(maskS.vol(:));
nrMaskVox=numel(maskInd);
funcInMask=reshape(funcS.vol,[],funcS.nframes);
funcInMask=funcInMask(maskInd,:)';
funcInMask=bsxfun(@minus,funcInMask,mean(funcInMask,1));
% for seedNr=1:numel(seed)
%     seedTs=funcInMask(:,seed(seedNr).mind);
%     lsPicPcc=seedTs\funcInMask;
%     tmpFuncS=funcS;
%     tmpFuncS.vol=zeros(prod(funcS.volsize),1);
%     tmpFuncS.vol(maskInd,:)=lsPicPcc;
%     tmpFuncS.vol=reshape(tmpFuncS.vol,[tmpFuncS.volsize,1]);
%     MRIsave(tmpFuncS,tmpFuncS.vol,fullfile(picDir,['tract_',seed(seedNr).name,'.nii.gz']),1);
%     clear tmpFuncS
% end

% mask subscripts
[x,y,z]=ind2sub(maskS.volsize,maskInd);
maskSub=[x,y,z];
[roiShape,nrNeighMax] = getRoiShape(sz);
% neighbor subscripts: [nrNeigh*nrVox,3] (includes out-of-range neighbors)
neighSub=reshape(repmat(reshape(maskSub,1,nrMaskVox,3),nrNeighMax,1,1),[],3) + repmat(roiShape,nrMaskVox,1);
neighOutOfRange=any(neighSub<1,2) | neighSub(:,1)>funcS.volsize(1) | neighSub(:,2)>funcS.volsize(2) | neighSub(:,3)>funcS.volsize(3);
neighSub(neighOutOfRange,:)=1;
% neighbor indices
neighInd=sub2ind(funcS.volsize,neighSub(:,1),neighSub(:,2),neighSub(:,3));
[inMask,neighMind]=ismember(neighInd,maskInd);
neighInd(~inMask)=1;
% neighbor func data
neighFunc=funcS.vol(neighInd,:);
neighFunc(neighInd==1,:)=NaN;
neighInd(neighInd==1)=NaN;
neighInd=reshape(neighInd,nrNeighMax,[]);
neighMind(neighMind==0)=NaN;
neighMind=reshape(neighMind,nrNeighMax,[]);
nrIc=3;
tic
for seedNr=1:nrSeed
    seed(seedNr).neighInd=neighInd(:,seed(seedNr).mind);
    seed(seedNr).neighMind=neighMind(:,seed(seedNr).mind);
    [ic, A, W] = fastica(funcInMask(:,seed(seedNr).neighMind)','numOfIC',nrIc);
    seed(seedNr).ic=ic;
    seed(seedNr).A=A;
end
for f=0:2,figure('WindowStyle','docked'),imagesc(corr(seed(1+f).ic',seed(2+f).ic'));end
figure
for seedNr=1:nrSeed, for icNr=1:nrIc, for z=1:3,
    subplot(nrSeed,nrIc*3,(seedNr-1)*nrIc*3+mod(icNr,3)*nrIc+z);
    imagesc(reshape(seed(seedNr).A(roiShape(:,3)==z-2,icNr),3,3));
end, end, end
wmS=MRIread('/project_space/pboord/act/rest/109003/func_wm.nii.gz');
sum(wmS.vol(:)>0)*toc/nrSeed/60
for icNr=1:nrIc
    tmpFuncS=funcS;
    tmpFuncS.vol=zeros(prod(funcS.volsize),1);
    tmpFuncS.vol(maskInd)=corr(funcInMask,icasig(icNr,:)');
    tmpFuncS.vol=reshape(tmpFuncS.vol,[tmpFuncS.volsize,1]);
    MRIsave(tmpFuncS,tmpFuncS.vol,fullfile(picDir,['icFc_',num2str(icNr),'_',seed(seedNr).name,'.nii.gz']),1);
end

nrNeigh=sum(~isnan(neighMind));
neighFunc=reshape(neighFunc,nrNeighMax,[],nf);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
varExp=zeros(nrMaskVox,1);
timePerVp=0;
for neighNr=uniqueNrNeigh
    tic
    neighNrVox=nrNeigh==neighNr;
    % nr voxel patterns
    nrVp=sum(neighNrVox);
    if neighNr==27
        % conserve memory
        neighFunc(:,~neighNrVox,:)=[];
        neighNrFunc=neighFunc;
        clear neighFunc
    else
        neighNrFunc=neighFunc(:,neighNrVox,:);
    end
    neighNrFunc(isnan(neighNrFunc))=[];
    neighNrFunc=reshape(neighNrFunc,neighNr,[],nf);
    neighNrFunc=neighNrFunc-repmat(mean(neighNrFunc,3),1,1,nf);
    neighNrInd=find(neighNrVox);
    
    neighNrMind=neighMind(:,neighNrInd);
    lsVec=reshape(lsPicPcc(neighNrMind(~isnan(neighNrMind))),neighNr,nrVp);
    lsVec=lsVec./repmat(sqrt(sum(lsVec.^2)),neighNr,1);
    varExp(neighNrInd,:)=100*sum(squeeze(sum(neighNrFunc.*repmat(lsVec,1,1,nf))).^2,2)./sum(sum(neighNrFunc.^2,3))';
    disp([neighNr,nrVp,nrVp*timePerVp]);
    clear neighNrFunc2d
    timePerVp=toc/nrVp;
    
    clear neighNrFunc
    
end
tmpFuncS=funcS;
tmpFuncS.vol=zeros(prod(funcS.volsize),1);
tmpFuncS.vol(maskInd,:)=varExp;
tmpFuncS.vol=reshape(tmpFuncS.vol,[tmpFuncS.volsize,1]);
MRIsave(tmpFuncS,tmpFuncS.vol,fullfile(picDir,['varExp_',seed(seedNr).name,'.nii.gz']),1);
end

function [roiShape,nrVox] = getRoiShape(sz)
roiShape = zeros((2*sz+1)^3,3);
nrVox = 0;
for x = -1:1
    for y = -1:1
        for z = -1:1
            nrVox = nrVox + 1;
            roiShape(nrVox,:) = [x,y,z];
        end
    end
end
end

