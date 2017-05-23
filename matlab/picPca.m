function picPca(dirName,subjCodeNrStr,filePrefix,maskName)
% from SGE_roiCorr52
if ~isdeployed
    dbstop if error
end
if nargin < 3
    filePrefix='regress';
end
if nargin < 4
    maskName='mask.nii.gz';
end
if nargin < 2
    subjCodeNrStr='109003';
end
if nargin < 1
    dirName='/mnt/pboord_5TB/pic/act';
end
%% constants
sz = 1;
minNeighNr=3;
%% seeds
seed(1).name = 'PCC';
seed(1).mni = [-6,-58,28];
seed(2).name = 'WM';
seed(2).mni = [-26,-6,28];%[30,-46,28];
seed(3).name = 'LpIPS';
seed(3).mni = [-26,-65,52];
seed(4).name = 'LaIPS';
seed(4).mni = [-45,-37,48];
seed(5).name = 'RpIPS';
seed(5).mni = [28,-65,51];
seed(6).name = 'RaIPS';
seed(6).mni = [43,-36,46];
seed(7).name = 'LFEF';
seed(7).mni = [-29,-5,55];
seed(8).name = 'RFEF';
seed(8).mni = [31,-5,54];
seed(9).name = 'RDLPFC';
seed(9).mni = [0,0,0];
seed(10).name = 'LDLPFC';
seed(10).mni = [0,0,0];
seed(11).name = 'Rfrontal';
seed(11).mni = [0,0,0];
seed(12).name = 'Lfrontal';
seed(12).mni = [0,0,0];
seed(13).name = 'RIPL';
seed(13).mni = [0,0,0];
seed(14).name = 'LIPL';
seed(14).mni = [0,0,0];
seed(15).name = 'RIPS';
seed(15).mni = [0,0,0];
seed(16).name = 'LIPS';
seed(16).mni = [0,0,0];
seed(17).name = 'mPFC';
seed(17).mni = [-2,50,17];
seed(18).name = 'RAG';
seed(18).mni = [0,0,0];
seed(19).name = 'LBG';
seed(20).name = 'RBG';
seed(21).name = 'LHIP';
seed(22).name = 'RHIP';
seed(23).name = 'RTP';
% seed(24).name = 'Rfinger';
% seed(25).name = 'Lfinger';
% seed(26).name = 'LpIPS1';
% seed(27).name = 'LpIPS2';
% seed(28).name = 'LpIPS3';
% seed(29).name = 'LpIPS4';
% seed(30).name = 'LpIPS5';
% seed(31).name = 'LpIPS6';
% seed(32).name = 'lowCorrPicAct';
%
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));
maskInd=find(maskS.vol(:));
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
fastPic(funcS,maskS,sz,minNeighNr,pathstr,name,seed);

end

function fastPic(funcS,maskS,sz,minNeighNr,pathstr,name,seed)
nf=funcS.nframes;
nrSeed=numel(seed);
maskInd=find(maskS.vol(:));
nrMaskVox=numel(maskInd);
funcInMask=reshape(funcS.vol,[],funcS.nframes);
funcInMask=funcInMask(maskInd,:);
funcInMask=bsxfun(@minus,funcInMask,mean(funcInMask,2));

% funcInMask=funcInMask';
% FF=funcInMask'*funcInMask;
% % FF=zeros(nrMaskVox,nrMaskVox);
% % for vNr=1:nrMaskVox
% %     tic
% %     %FF(:,vNr)=mean(funcInMask.*repmat(funcInMask(vNr,:)>=0,nrMaskVox,1)-funcInMask.*repmat(funcInMask(vNr,:)<0,nrMaskVox,1),2);
% %     FF(:,vNr)=corr(funcInMask(:,vNr),funcInMask);
% %     (nrMaskVox-vNr)*toc/60%disp(toc*nrMaskVox/n/60)
% % end


filePrefixDir=fullfile(pathstr,name);
if ~exist(filePrefixDir,'dir')
    mkdir(filePrefixDir);
end
picFcCorrDir=fullfile(filePrefixDir,'picFcCorr');
if ~exist(picFcCorrDir,'dir')
    mkdir(picFcCorrDir);
end
picDir=fullfile(filePrefixDir,'pic');
if ~exist(picDir,'dir')
    mkdir(picDir);
end
fcDir=fullfile(filePrefixDir,'fc');
if ~exist(fcDir,'dir')
    mkdir(fcDir);
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
% neighbor func data
neighFunc=funcS.vol(neighInd,:);
neighFunc(neighInd==1,:)=NaN;
neighInd(neighInd==1)=NaN;
neighInd=reshape(neighInd,nrNeighMax,[]);
neighMind(neighMind==0)=NaN;
neighMind=reshape(neighMind,nrNeighMax,[]);
for seedNr=1:nrSeed
    seed(seedNr).neighInd=neighInd(:,seed(seedNr).mind);
    seed(seedNr).neighMind=neighMind(:,seed(seedNr).mind);
end
nrNeigh=sum(~isnan(neighMind));
neighFunc=reshape(neighFunc,nrNeighMax,[],nf);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
mnVec=zeros(nrMaskVox,27);
mnVarExp=zeros(nrMaskVox,1);
mnTs=zeros(nrMaskVox,nf);
mnActCorr=zeros(nrMaskVox,1);
mnVoxCorr=zeros(nrMaskVox,1);
mnActAngle=zeros(nrMaskVox,1);
ctrActCorr=zeros(nrMaskVox,1);
corrVpTs=zeros(nrMaskVox,nf);
eig1vec=zeros(nrMaskVox,27);
eig2vec=zeros(nrMaskVox,27);
nrPc=2;
eigVals=zeros(nrMaskVox,nrPc);
eigVarExp=zeros(nrMaskVox,nrPc);
eig1ts=zeros(nrMaskVox,nf);
eig2ts=zeros(nrMaskVox,nf);
eigActCorr=zeros(nrMaskVox,nrPc);
eigVoxCorr=zeros(nrMaskVox,nrPc);
eigActAngle=zeros(nrMaskVox,nrPc);
act=zeros(nrMaskVox,nf);
timePerVp=0;
for neighNr=uniqueNrNeigh
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
    ctrVoxIdx=sum(~isnan(neighNrFunc(1:14,:,1)));
    neighNrFunc(isnan(neighNrFunc))=[];
    neighNrFunc=reshape(neighNrFunc,neighNr,[],nf);
    neighNrFunc=neighNrFunc-repmat(mean(neighNrFunc,3),1,1,nf);
    %
    ROpermute=reshape(permute(neighNrFunc,[1,3,2]),neighNr,[]);
    corrWithExtraSample=reshape([corrColumns(ROpermute(:,1:end-1),ROpermute(:,2:end)),1],[],nrVp)';
    clear ROpermute
    corrVpTs(neighNrVox,:)=[corrWithExtraSample(:,1:end-1),mean(abs(corrWithExtraSample(:,1:end-1)),2)];

    neighNrInd=find(neighNrVox);
    act(neighNrInd,:)=squeeze(mean(neighNrFunc,1));
    disp([neighNr,nrVp,nrVp*timePerVp]);
    
    posInd=repmat(mean(neighNrFunc,1)>=0,neighNr,1,1);
    posnegVp=mean((neighNrFunc.*posInd)-(neighNrFunc.*~posInd),3);
    posnegVp=posnegVp./repmat(sqrt(sum(posnegVp.^2,1)),neighNr,1);
    actAngle=acosd(sum(posnegVp)/sqrt(neighNr));
    % don't flip bec no ambiguity for mean, cf PCA!
    %     flipAngle=actAngle>90;
    %     mnActAngle(neighNrInd,:)=abs(180*flipAngle-actAngle);
    %     posnegVp=posnegVp.*repmat(~flipAngle-flipAngle,neighNr,1);
    mnActAngle(neighNrInd,:)=actAngle;
    mnVec(neighNrInd,1:neighNr)=posnegVp';
    ts=squeeze(sum(repmat(posnegVp,1,1,nf).*neighNrFunc));
    mnTs(neighNrInd,:)=ts;
    mnVarExp(neighNrInd,:)=100*sum(ts.^2,2)./squeeze(sum(sum(neighNrFunc.^2,1),3))';
    actCorr=corrColumns(squeeze(mean(neighNrFunc))',ts');
    mnActCorr(neighNrInd,:)=actCorr;
    neighNrFunc2d=reshape(neighNrFunc,[],nf);
    voxCorr=corrColumns(neighNrFunc2d(ctrVoxIdx+(0:neighNr:(nrVp-1)*neighNr),:)',ts');
    mnVoxCorr(neighNrInd,:)=voxCorr;
    ctrActCorr(neighNrInd,:)=corrColumns(neighNrFunc2d(ctrVoxIdx+(0:neighNr:(nrVp-1)*neighNr),:)',squeeze(mean(neighNrFunc))');
    clear neighNrFunc2d
%     tic
%     for vpNr=1:nrVp
%         X=squeeze(neighNrFunc(:,vpNr,:));
%         % ref: pcasecon.m
%         mu = mean(X,2);
%         X = bsxfun(@minus,X,mu);
%         C = X'*X;
%         [V,D] = eigs(C,nrPc);
%         diagD=diag(D);
%         clear C;
%         U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
%         s = sqrt(abs(diagD));
%         U = bsxfun(@(x,c)x./c, U, s');
%         eigVals(neighNrInd(vpNr),:)=diagD;
%         eigVarExp(neighNrInd(vpNr),:)=100*diagD/sum(X(:).^2);
%         actAngle=[acosd(sum(U(:,1))/sqrt(neighNr)),acosd(sum(U(:,2))/sqrt(neighNr))];
%         % IS THIS RIGHT???... needed for sign ambiguity of PCA
%         flipAngle=actAngle>90;
%         eigActAngle(neighNrInd(vpNr),:)=abs(180*flipAngle-actAngle);
%         U=U.*repmat(~flipAngle-flipAngle,neighNr,1);
%         eig1vec(neighNrInd(vpNr),1:neighNr)=U(:,1);
%         eig2vec(neighNrInd(vpNr),1:neighNr)=U(:,2);
%         ts=U'*X;
%         eig1ts(neighNrInd(vpNr),:)=ts(1,:);
%         eig2ts(neighNrInd(vpNr),:)=ts(2,:);
%         actCorr=[corr(mean(X,1)',ts(1,:)'),corr(mean(X,1)',ts(2,:)')];
%         eigActCorr(neighNrInd(vpNr),:)=actCorr;
%         voxCorr=[corr(X(ctrVoxIdx(1),:)',ts(1,:)'),corr(X(ctrVoxIdx(2),:)',ts(2,:)')];
%         eigVoxCorr(neighNrInd(vpNr),:)=voxCorr;
%     end
%     timePerVp=toc/nrVp;
    
    clear neighNrFunc

end

for seedNr=1:nrSeed
    seed(seedNr).seedTs=funcInMask(seed(seedNr).neighMind,:);
    seed(seedNr).corrVpTs=corrVpTs(seed(seedNr).mind,:);
    seed(seedNr).mnVec=mnVec(seed(seedNr).mind,:);
    seed(seedNr).mnVarExp=mnVarExp(seed(seedNr).mind,:);
    seed(seedNr).mnTs=mnTs(seed(seedNr).mind,:);
    seed(seedNr).mnVoxCorr=mnVoxCorr(seed(seedNr).mind,:);
    seed(seedNr).mnActAngle=mnActAngle(seed(seedNr).mind,:);
    seed(seedNr).act=act(seed(seedNr).mind,:);
    seed(seedNr).posnegVp=zeros(27,nrSeed);
    seed(seedNr).posnegVpNonNorm=zeros(27,nrSeed);
    seed(seedNr).vpSe=zeros(27,nrSeed);
    seed(seedNr).projTs=zeros(nf,nrSeed);
    seed(seedNr).relProjTs=zeros(nf,nrSeed);
end
for seedNrLoc=1:nrSeed
    locSeedTs=seed(seedNrLoc).seedTs;
    for seedNrRem=1:nrSeed
        posInd=repmat(seed(seedNrRem).act>=0,27,1);
        %posInd=repmat(seed(seedNrRem).mnTs>=0,27,1);
        %absRemWeight=abs(repmat(seed(seedNrRem).mnTs,27,1));
        unsignedLocSeedTs=(locSeedTs.*posInd)-(locSeedTs.*~posInd);
        %unsignedLocSeedTs=(locSeedTs.*posInd.*absRemWeight)-(locSeedTs.*~posInd.*absRemWeight);
        %normUnsignedLocSeedTs=unsignedLocSeedTs./repmat(sqrt(sum(unsignedLocSeedTs.^2)),27,1);
        posnegVpNonNorm=mean(unsignedLocSeedTs,2);
        %posnegVp=mean(normUnsignedLocSeedTs,2);
        seed(seedNrLoc).vpSe(:,seedNrRem)=std(unsignedLocSeedTs,0,2)/sqrt(nf);
        seed(seedNrLoc).posnegVpNonNorm(:,seedNrRem)=posnegVpNonNorm;
        seed(seedNrLoc).posnegVp(:,seedNrRem)=posnegVpNonNorm/sqrt(sum(posnegVpNonNorm.^2));
        seed(seedNrLoc).projTs(:,seedNrRem)...
            =(seed(seedNrLoc).posnegVp(:,seedNrRem)'*(locSeedTs./repmat(sqrt(sum(locSeedTs.^2)),27,1)));
    end
end
relSeed=1;
for seedNrLoc=1:nrSeed
    locSeedTs=seed(seedNrLoc).seedTs./repmat(sqrt(sum(locSeedTs.^2)),27,1);
    for seedNrRem=1:nrSeed
        seed(seedNrLoc).relProjTs(:,seedNrRem)=seed(seedNrLoc).posnegVp(:,seedNrRem)'*locSeedTs;
    end
end
% % Figures
% set(0,'DefaultAxesFontSize',16);
% figure;s=6;t=1;
% subplot(3,2,1);plot(seed(s).seedTs');
% xlim([1,nf]);title([seed(s).name,' neighborhood voxels']);
% subplot(3,2,[2,4,6]);
% plot(seed(s).posnegVpNonNorm(:,s),'r-','Linewidth',3);hold on;
% plot(seed(s).posnegVpNonNorm(:,t),'b-','Linewidth',3);
% plot(seed(s).posnegVpNonNorm(:,s)+1.96*seed(s).vpSe(:,s),'k:');plot(seed(s).posnegVpNonNorm(:,s)-1.96*seed(s).vpSe(:,s),'k:');
% plot(seed(s).posnegVpNonNorm(:,t)+1.96*seed(s).vpSe(:,t),'k:');plot(seed(s).posnegVpNonNorm(:,t)-1.96*seed(s).vpSe(:,t),'k:');
% legend(seed(s).name,seed(t).name,'Location','Northwest');
% %plot(eig1vec(seed(s).mind,:));
% hold off;
% xlim([1,27]);xlabel('Neighborhood voxel');ylabel('Pattern amplitude (a.u.)');
% title(['Sign weighted patterns in ',seed(s).name]);
% disp(corr(seed(s).posnegVpNonNorm(:,s),seed(s).posnegVpNonNorm(:,t)));
% subplot(3,2,3);
% plot(seed(t).seedTs');hold on;
% title([seed(t).name,' neighborhood voxels']);xlim([1,nf]);
% hold off
% subplot(3,2,5);
% hold on;plot(mean(seed(s).seedTs),'r-','Linewidth',3);plot(mean(seed(t).seedTs),'b-','Linewidth',3);hold off
% xlim([1,nf]);xlabel('fMRI volume');
% title(['Mean ',seed(s).name,' and ',seed(t).name,' timeseries'])
% disp(corr(mean(seed(s).seedTs)',mean(seed(t).seedTs)'));
% %disp(corr(seed(s).act',seed(t).act'));
% 
% subplot(2,2,4);plot(seed(s).posnegVpNonNorm(:,t));hold on;
% %plot(seed(t).posnegVpNonNorm(:,t)+1.96*seed(t).vpSe(:,t),'k:');plot(seed(t).posnegVpNonNorm(:,t)-1.96*seed(t).vpSe(:,t),'k:');
% hold off;
% xlim([1,27]);xlabel('Neighborhood voxel');ylabel('Pattern vector (a.u.)');

%
% corrVpTs=zeros(nrMaskVox,nf);
funcS.vol=zeros(prod(funcS.volsize),nf);
funcS.vol(maskInd,:)=corrVpTs;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'corrVpTs.nii.gz'),nf);
%
s=1;t=4;
posInd=repmat(reshape(act>=0,1,nrMaskVox,nf),27,1,1);
seedTsS=repmat(reshape(seed(s).seedTs,27,1,nf),1,nrMaskVox,1);
seedTsT=repmat(reshape(seed(t).seedTs,27,1,nf),1,nrMaskVox,1);
posnegVpS=mean((seedTsS.*posInd)-(seedTsS.*~posInd),3); clear seedTsS
posnegVpT=mean((seedTsT.*posInd)-(seedTsT.*~posInd),3); clear seedTsT
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=corrColumns(posnegVpS,posnegVpT);
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,['corr_',seed(s).name,'_',seed(t).name,'.nii.gz']),1);
nrRand=100;
for seedNr=1:nrSeed
    seedTs=repmat(reshape(seed(seedNr).seedTs,27,1,nf),1,nrMaskVox,1);
    tic
    posnegVpSrssRand=zeros(nrRand,nrMaskVox);
    for randNr=1:nrRand
        posIndRand=repmat(reshape(circshift(flip(act>=0,2),randi(nf),2),1,nrMaskVox,nf),27,1,1);
        %posIndRand=repmat(reshape(circshift(flip(act>=0,2),0,2),1,nrMaskVox,nf),27,1,1);
        posnegVpSrssRand(randNr,:)=sqrt(sum((mean((seedTs.*posIndRand)-(seedTs.*~posIndRand),3)).^2,1));
    end
    posnegVpSrssRand=sort(posnegVpSrssRand);
    toc
    posInd=repmat(reshape(act>=0,1,nrMaskVox,nf),27,1,1);
    %posInd=repmat(reshape(mnTs>=0,1,nrMaskVox,nf),27,1,1);
    posnegVp=mean((seedTs.*posInd)-(seedTs.*~posInd),3);
%     posnegVp=mean((seedTs(:,:,1:nf/2).*posInd(:,:,1:nf/2))-(seedTs(:,:,1:nf/2).*~posInd(:,:,1:nf/2)),3);
%     posnegVp=mean((seedTs(:,:,nf/2+1:end).*posInd(:,:,nf/2+1:end))-(seedTs(:,:,nf/2+1:end).*~posInd(:,:,nf/2+1:end)),3);
    posnegVpSrss=sqrt(sum(posnegVp.^2,1));
    sigMask=posnegVpSrss>=posnegVpSrssRand(99,:);
    posnegVp=posnegVp./repmat(posnegVpSrss,neighNr,1);
%     funcS.vol=zeros(prod(funcS.volsize),1);
%     funcS.vol(maskInd,:)=sigMask.*acosd(sum(repmat(seed(seedNr).mnVec',1,nrMaskVox).*posnegVp));
%     funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
%     MRIsave(funcS,funcS.vol,fullfile(picDir,['angleMnTs_',seed(seedNr).name,'.nii.gz']),1);
    %figure('WindowStyle','docked');plotyy(1:nf,seed(seedNr).mnTs',1:nf,sum(posnegVp(:,sigMask)*act(sigMask,:))')
    %disp(corr(seed(seedNr).mnTs',sum(posnegVp(:,sigMask)*act(sigMask,:))'));
    a=(corr(posnegVp(:,sigMask)));
    figure('WindowStyle','docked');hist(a(triu(true(sum(sigMask)),1)))
end
for seedNr=1:nrSeed
    funcS.vol=zeros(prod(funcS.volsize),1);
    funcS.vol(maskInd,:)=corr(act',act(seed(seedNr).mind,:)');
    funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
    MRIsave(funcS,funcS.vol,fullfile(picDir,['corrAct_',seed(seedNr).name,'.nii.gz']),1);   
end

% save
% mnVec=zeros(nrMaskVox,27);
funcS.vol=zeros(prod(funcS.volsize),27);Ts
funcS.vol(maskInd,:)=mnVec;
funcS.vol=reshape(funcS.vol,[funcS.volsize,27]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnVec.nii.gz'),27);
% mnVarExp=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd,:)=mnVarExp;
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnVarExp.nii.gz'),1);
% mnTs=zeros(nrMaskVox,nf);
funcS.vol=zeros(prod(funcS.volsize),nf);
funcS.vol(maskInd,:)=mnTs;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnTs.nii.gz'),nf);
% mnActCorr=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd,:)=mnActCorr;
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnActCorr.nii.gz'),1);
% mnVoxCorr=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd,:)=mnVoxCorr;
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnVoxCorr.nii.gz'),1);
% mnActAngle=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd,:)=mnActAngle;
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'mnActAngle.nii.gz'),1);
% ctrActCorr=zeros(nrMaskVox,1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd,:)=ctrActCorr;
funcS.vol=reshape(funcS.vol,[funcS.volsize,1]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'ctrActCorr.nii.gz'),1);
% act=zeros(prod(funcS.volsize),nf);
funcS.vol=zeros(prod(funcS.volsize),nf);
funcS.vol(maskInd,:)=act;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'act.nii.gz'),nf);
% eig1vec=zeros(nrMaskVox,27);
funcS.vol=zeros(prod(funcS.volsize),27);
funcS.vol(maskInd,:)=eig1vec;
funcS.vol=reshape(funcS.vol,[funcS.volsize,27]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eig1vec.nii.gz'),27);
% eig2vec=zeros(nrMaskVox,27);
funcS.vol=zeros(prod(funcS.volsize),27);
funcS.vol(maskInd,:)=eig2vec;
funcS.vol=reshape(funcS.vol,[funcS.volsize,27]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eig2vec.nii.gz'),27);
% eigVals=zeros(nrMaskVox,nrPc);
funcS.vol=zeros(prod(funcS.volsize),nrPc);
funcS.vol(maskInd,:)=eigVals;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nrPc]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVals.nii.gz'),nrPc);
% eigVarExp=zeros(nrMaskVox,nrPc);
funcS.vol=zeros(prod(funcS.volsize),nrPc);
funcS.vol(maskInd,:)=eigVarExp;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nrPc]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVarExp.nii.gz'),nrPc);
% eig1ts=zeros(nrMaskVox,nf);
funcS.vol=zeros(prod(funcS.volsize),nf);
funcS.vol(maskInd,:)=eig1ts;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eig1ts.nii.gz'),nf);
% eig2ts=zeros(nrMaskVox,nf);
funcS.vol=zeros(prod(funcS.volsize),nf);
funcS.vol(maskInd,:)=eig2ts;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eig2ts.nii.gz'),nf);
% eigActCorr=zeros(nrMaskVox,nrPc);
funcS.vol=zeros(prod(funcS.volsize),nrPc);
funcS.vol(maskInd,:)=eigActCorr;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nrPc]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigActCorr.nii.gz'),nrPc);
% eigVoxCorr=zeros(nrMaskVox,nrPc);
funcS.vol=zeros(prod(funcS.volsize),nrPc);
funcS.vol(maskInd,:)=eigVoxCorr;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nrPc]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigVoxCorr.nii.gz'),nrPc);
% eigActAngle=zeros(nrMaskVox,nrPc);
funcS.vol=zeros(prod(funcS.volsize),nrPc);
funcS.vol(maskInd,:)=eigActAngle;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nrPc]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'eigActAngle.nii.gz'),nrPc);
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
