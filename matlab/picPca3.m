function picPca3(dirName,subjCodeNrStr,filePrefix,maskName)
% from SGE_roiCorr52
if ~isdeployed
    dbstop if error
end
if nargin < 3
    % 0.01-0.1 Hz
    tr=2.5;
    for fNr=0:5
        d=(0.1-0.01)/3;
        lo=0.01+fNr*d;
        hi=lo+d;
        disp([1/(tr*lo) 1/(tr*hi)]);
    end
    %filePrefix='rest_medn_40_10';
    %filePrefix='rest_medn_10_5.7';
    %filePrefix='rest_medn_5.7_4';
    %filePrefix='rest_medn_4_3.08';
    %filePrefix='rest_medn_3.08_2.5';
    %filePrefix='rest_medn_2.5_2.1';
    %filePrefix='regress';
    %filePrefix='aicor_brain';
    filePrefix='raw_resting';
end
if nargin < 4
    maskName='mask.nii.gz';
    %maskName='func_gm.nii.gz';
end
if nargin < 2
    subjCodeNrStr='109003';
    %subjCodeNrStr='RC4101-1';
end
if nargin < 1
    dirName='/project_space/pboord/act/rest';
    %dirName='/mnt/pboord_5TB/pic/act';
    %nasDir='/project_space/pboord/act/rest';
    %dirName='/mnt/pboord_5TB/udall';
end
%% constants
sz = 1;
minNeighNr=3;
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
funcS.vol=reshape(funcS.vol,[],funcS.nframes);
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));

[pathstr,name,ext]=fileparts(funcS.fspec);
if true
    %volToHist='sumFreq.nii.gz';
    volToHist='eigVarExp.nii.gz';
    filePrefix='raw_resting';
    wmVolToHistS=MRIread(fullfile(dirName,subjCodeNrStr,'mask',filePrefix,'pic',volToHist));
    gmVolToHistS=MRIread(fullfile(dirName,subjCodeNrStr,'func_gm',filePrefix,'pic',volToHist));
    x=1:10:100;
    %x=0:0.1:1;
    wmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_wm.nii.gz'));
    gmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_gm.nii.gz'));
    wm=hist(wmVolToHistS.vol(wmS.vol>0),x);gm=hist(gmVolToHistS.vol(gmS.vol>0),x);figure,bar(x,[wm/sum(wm);gm/sum(gm)]')
    % threshold
    wmSort=sort(wmVolToHistS.vol(wmS.vol>0));
    p=0.05;
    t=wmSort(round((1-p)*numel(wmSort)));
    disp(t);
end
if false
    vecS=MRIread('/project_space/pboord/act/rest/109003/eigVec6.nii.gz');
    %vecS=MRIread('/project_space/pboord/act/rest/109003/eigVec6mask.nii.gz');
    vecVol=reshape(vecS.vol,vecS.volsize(1),vecS.volsize(2),vecS.volsize(3),27,[]);
    nrFreq=size(vecVol,5);
    vecVol=vecVol./repmat(sqrt(sum(vecVol.^2,4)),1,1,1,27,1);
    sumFreq=sqrt(sum(sum(vecVol,5).^2,4))/nrFreq;
    MRIsave(vecS,sumFreq,'/project_space/pboord/act/rest/109003/sumFreq.nii.gz',1);
    gmVolToHistS=MRIread(fullfile(dirName,subjCodeNrStr,'sumFreq.nii.gz'));
    wmVolToHistS=MRIread(fullfile(dirName,subjCodeNrStr,'mask','sumFreq.nii.gz'));
    x=0:0.1:1;
    wmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_wm_ero.nii.gz'));
    gmS=MRIread(fullfile(dirName,subjCodeNrStr,'func_gm.nii.gz'));
    wm=hist(wmVolToHistS.vol(wmS.vol>0),x);gm=hist(gmVolToHistS.vol(gmS.vol>0),x);figure,bar(x,[wm/sum(wm);gm/sum(gm)]')
    % threshold
    wmSort=sort(wmVolToHistS.vol(wmS.vol>0));
    p=0.05;
    t=wmSort(round((1-p)*numel(wmSort)));
    disp(t);
end
%picDsn(funcS,maskS,wmS,gmS);

name=strtok([name,ext],'.');
fastPic(funcS,maskS,sz,minNeighNr,pathstr,name);

end
function picDsn(funcS,maskS,wmS,gmS)
nrSim=10000;
picAngle=zeros(nrSim,2);
nf=funcS.nframes;
func2d=reshape(funcS.vol,[],nf);
func2d=func2d(maskS.vol(:)>0,:)';
nvox=size(func2d,2);
func2deq=ifft(exp(1i*angle(fft(func2d))),'symmetric');
nrNeigh=27;
for simNr=1:nrSim
    x=func2deq(:,randperm(nvox,nrNeigh));
    x=x./repmat(sqrt(sum(x.^2,2)),1,nrNeigh);
    [eigX,~,varExp]=fastEig(x');
    rawAngle=acosd(sum(eigX)/sqrt(nrNeigh));
    picAngle(simNr,1)=abs(180*(rawAngle>90)-rawAngle);
    picAngle(simNr,2)=varExp;
end
figure,hist(picAngle(:,2))
end

function fastPic(funcS,maskS,sz,minNeighNr,pathstr,name)
nf=funcS.nframes;
maskInd=find(maskS.vol(:));
nrMaskVox=numel(maskInd);
filePrefixDir=fullfile(pathstr,name);
if ~exist(filePrefixDir,'dir')
    mkdir(filePrefixDir);
end
picDir=fullfile(filePrefixDir,'pic');
if ~exist(picDir,'dir')
    mkdir(picDir);
end
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
% neighbor func dataeigVarExp
neighFunc=funcS.vol(neighInd,:);
neighFunc(neighInd==1,:)=NaN;
neighMind(neighMind==0)=NaN;
neighMind=reshape(neighMind,nrNeighMax,[]);
nrNeigh=sum(~isnan(neighMind));
neighFunc=reshape(neighFunc,nrNeighMax,[],nf);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
eigVec=zeros(nrMaskVox,27);
eigVals=zeros(nrMaskVox,1);
eigVarExp=zeros(nrMaskVox,1);
eigActAngle=zeros(nrMaskVox,1);
nrNeighVol=zeros(nrMaskVox,1);
timePerVp=0;
for neighNr=uniqueNrNeigh
    neighNrVox=nrNeigh==neighNr;
    nrNeighVol(neighNrVox)=neighNr;
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
    disp([neighNr,nrVp,nrVp*timePerVp]);
    clear neighNrFunc2d
    tic
    for vpNr=1:nrVp
        X=squeeze(neighNrFunc(:,vpNr,:));
        % ref: pcasecon.m
        mu = mean(X,2);
        X = bsxfun(@minus,X,mu);        
        % equalize
        X=ifft(exp(1i*angle(fft(X'))),'symmetric')';
        C = X*X';
        [U,D] = eigs(C,1);
        diagD=diag(D);
        clear C;
        eigVals(neighNrInd(vpNr))=diagD;
        eigVarExp(neighNrInd(vpNr))=100*diagD/sum(X(:).^2);
        actAngle=acosd(sum(U)/sqrt(neighNr));
        % IS THIS RIGHT???... needed for sign ambiguity of PCA
        flipAngle=actAngle>90;
        eigActAngle(neighNrInd(vpNr),:)=abs(180*flipAngle-actAngle);
        U=U.*repmat(~flipAngle-flipAngle,neighNr,1);
        eigVec(neighNrInd(vpNr),1:neighNr)=U;
    end
    timePerVp=toc/nrVp;
    
    clear neighNrFunc

end
% nrNeighVol=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=nrNeighVol;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'nrNeighVol.nii.gz'),1);
% eigActAngle=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=eigActAngle;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigActAngle_equalized.nii.gz'),1);
% eigVec=zeros(nrMaskVox,27);
funcS.vol=zeros(prod(funcS.volsize),27);
funcS.vol(maskInd,:)=eigVec;
funcS.vol=reshape(funcS.vol,[funcS.volsize,27]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVec.nii.gz'),27);
% eigVals=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=eigVals;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVals.nii.gz'),1);
% eigVarExp=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=eigVarExp;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVarExp.nii.gz'),1);
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
