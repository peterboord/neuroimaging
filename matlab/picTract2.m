function picTract2

if ~isdeployed
    dbstop if error
end
%% files & directories
%filePrefix='rest_medn_40_10';
filePrefix='rest_medn';
%filePrefix='rest_tsoc';
maskName='mask.nii.gz';
subjCodeNrStr='109003';
dirName='/project_space/pboord/act/rest';
%% constants
tr=2.5;
sz = 1;
minNeighNr=3;
%% load data
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
nf=funcS.nframes;
funcS.vol=reshape(funcS.vol,[],nf);
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));
maskInd=find(maskS.vol(:));
wmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_wm.nii.gz'));
f2S=MRIread(fullfile(dirName,subjCodeNrStr,'inv_median_func_4d.nii.gz'));
f2d=reshape(f2S.vol,[],nf);
f2d=f2d(maskInd,:);

funcS.vol=reshape(funcS.vol,[],nf);
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
maskInd=find(maskS.vol(:));
nrMaskVox=numel(maskInd);
funcInMask=reshape(funcS.vol,[],nf);
funcInMask=funcInMask(maskInd,:)';
funcInMask=bsxfun(@minus,funcInMask,mean(funcInMask,1));
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
neighIc=nan(nrNeighMax,nf,nrMaskVox);
neighA=nan(nrNeighMax,nrNeighMax,nrMaskVox);
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
    disp([neighNr,nrVp,nrVp*timePerVp]);
    for vpNr=1:nrVp
        %rng('default');
        D=squeeze(neighNrFunc(:,vpNr,:));
        D=bsxfun(@minus,D,mean(D,2));
        [ic, A, W] = fastica(D,'approach','symm','verbose','off','stabilization','on');% defl or symm
        [propVar,sortIdx]=sort(sum((A*ic).^2,2)/sum(D(:).^2),'descend');
        nrIc=size(ic,1);
        D2=f2d(neighNrMind(~isnan(neighNrMind(:,vpNr)),vpNr),:);
        D2=bsxfun(@minus,D2,mean(D2,2));
        [ic2, ~,~] = fastica(D2,'approach','symm','verbose','off','stabilization','on');% defl or symm
        corrIc=corr(ic',ic2');
        if any(abs(corrIc(:))>0.8)
            keyboard
        end
%         neighIc(1:nrIc,:,neighNrInd(vpNr))=ic(sortIdx,:);
%         neighA(~isnan(neighNrMind(:,vpNr)),1:nrIc,neighNrInd(vpNr))=A(:,sortIdx);
    end
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

