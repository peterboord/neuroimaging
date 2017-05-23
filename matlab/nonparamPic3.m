function nonparamPic(dirName,subjCodeNrStr,filePrefix,maskName,nrRand)
% from SGE_roiCorr52
if ~isdeployed
    dbstop if error
end
if nargin < 5
    nrRand='0';
end
if nargin < 3
    filePrefix='regress';
    %filePrefix='regress.hpf';
    %filePrefix='volreg';
    %filePrefix='tshift';
    %filePrefix='ricor';
    %filePrefix='despike';
    %filePrefix='raw';
end
if nargin < 4
    maskName='mask.nii.gz';
%     switch filePrefix
%         case {'raw','despike','ricor','tshift'}
%             %maskName='raw_mask.nii.gz';
%             maskName='mask.nii.gz';
%         case {'volreg','regress'}
%             maskName = 'mask.nii.gz';
%         otherwise
%             error('unassigned mask');
%     end
end
if nargin < 2
    %subjCodeNrStr='3211';
    subjCodeNrStr='109003';
end
if nargin < 1
    %dirName='/mnt/pboord_5TB/pic/iowa/raw';
    %dirName='/project_space/pboord/PIC/iowa/raw';
    %dirName='/project_space/pboord/act/rest';
    dirName='/mnt/pboord_5TB/pic/act';
end
%% constants
nrRand=str2double(nrRand);
sz = 1;
%mniRef = '/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz';
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
seed(21).name = 'Rfinger';
seed(22).name = 'Lfinger';
seed(23).name = 'LHIP';
seed(24).name = 'RHIP';
seed(25).name = 'RTP';
seed(26).name = 'LpIPS1';
seed(27).name = 'LpIPS2';
seed(28).name = 'LpIPS3';
seed(29).name = 'LpIPS4';
seed(30).name = 'LpIPS5';
seed(31).name = 'LpIPS6';
% seed(9).name = 'salACC';
% seed(9).mni = [7,33,19];
% seed(10).name = 'salRFIC';
% seed(10).mni = [34,26,-6];
% seed(11).name = 'salLFIC';
% seed(11).mni = [-32,25,-10];
%%
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
for seedNr=1:numel(seed)
    seedFileName=['func_',seed(seedNr).name];
%     switch filePrefix
%         case {'raw','despike','ricor','tshift'}
%             %seedFileName=['raw_func_',seed(seedNr).name];
%             seedFileName=['func_',seed(seedNr).name];
%         case {'volreg','regress'}
%             seedFileName = ['func_',seed(seedNr).name];
%         otherwise
%             error('unassigned mask');
%     end
    seed(seedNr).vox=round(load(fullfile(dirName,subjCodeNrStr,seedFileName),'-ascii'));
    y=seed(seedNr).vox(1); % n.b. swapped
    x=seed(seedNr).vox(2); % x/y for FSL->FS
    z=seed(seedNr).vox(3);
    seed(seedNr).ind=sub2ind(funcS.volsize,x,y,z);
end
nrVox=prod(funcS.volsize);
halfUniqueFftLength = (funcS.nframes/2)-1;
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));
funcS.vol=reshape(funcS.vol,[],funcS.nframes);

[pathstr,name,ext]=fileparts(funcS.fspec);
name=strtok([name,ext],'.');
fastPic(0,funcS,maskS,[],sz,minNeighNr,seed,pathstr,name);

absFftFunc = abs(fft(funcS.vol'));
%angleFftFunc = angle(fft(funcS.vol'));
funcS.vol=[];
randNr=0;
while randNr<nrRand
    randNr=randNr+1;
    rng(str2double(subjCodeNrStr)+randNr,'twister');
    randPhase2d = reshape(pi*(2*rand(halfUniqueFftLength*nrVox,1) - 1),halfUniqueFftLength,nrVox);
    randPhase2d=[pi*(rand(1,nrVox)>0.5);randPhase2d;pi*(rand(1,nrVox)>0.5);-1*flipud(randPhase2d)]; %#ok<AGROW>
    tic
    
    %randNr=fastPic(randNr,funcS,maskS,...
    
    fastPic(randNr,funcS,maskS,...
        ifft(absFftFunc.*exp(1i*randPhase2d),'symmetric')',...
        sz,minNeighNr,seed,pathstr,name);
%     % randomize, maintaining covariance structure
%     % ref: 2012 Allen. Tracking Whole-Brain Connectivity Dynamics in the Resting State
%     randPhase2d=angleFftFunc+repmat(randPhase2d(:,1),1,nrVox);
%     fastPic(randNr,funcS,maskS,...
%         ifft(absFftFunc.*exp(1i*randPhase2d),'symmetric')',...
%         sz,minNeighNr,seed,pathstr,['cov_',name]);
    toc
end
end

function randNr=fastPic(randNr,funcS,maskS,func,sz,minNeighNr,seed,pathstr,name)
s=1;t=18;u=9;v=8;w=11;
normMode={'none'};%{'none','space','spaceTime','time','timeSpace'}
nf=funcS.nframes;
maskInd=find(maskS.vol(:));
nrMaskVox=numel(maskInd);
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
neighInd(~ismember(neighInd,maskInd))=1;
% neighbor func data
if isempty(func)
    neighFunc=funcS.vol(neighInd,:);
else
    neighFunc=func(neighInd,:);
end
neighFunc(neighInd==1,:)=NaN;
nrNeigh=sum(reshape(neighInd~=1,nrNeighMax,[]));
neighFunc=reshape(neighFunc,nrNeighMax,[],nf);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
picSlideInMask=zeros(nrMaskVox,nf);
fcInMask=zeros(nrMaskVox,nf);
ctrInMask=zeros(nrMaskVox,nf);
nzpInMaskLo=zeros(nrMaskVox,1);
nzpInMaskHi=zeros(nrMaskVox,1);
papInMask=zeros(nrMaskVox,1);
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
    ctrInMask(neighNrVox,:)=reshape(squeeze(neighNrFunc(14,:,:)),nrVp,[]);
    neighNrFunc(isnan(neighNrFunc))=[];
    neighNrFunc=reshape(neighNrFunc,neighNr,[],nf);
    fcInMask(neighNrVox,:)=reshape(mean(neighNrFunc,1),nrVp,[]);
    % 
    d=2;
    curFc=fcInMask(neighNrVox,:)';
    actVoxCorr=reshape(corrColumns(reshape(repmat(curFc,neighNr,1),nf,[]),reshape(neighNrFunc,[],nf)'),neighNr,nrVp);
    papInMask(neighNrVox)=sum(actVoxCorr<0)*100/neighNr;
    [sortedFc,idx]=sort(curFc,1);
    % convert row indices to matrix indices, so that sortedFc=curFc(idx)
    idx=idx+nf*repmat((0:nrVp-1),nf,1);
    % 2 May 16 begin
    % detrend
    neighNrFunc=neighNrFunc-repmat(mean(neighNrFunc,3),1,1,nf);
    % normalize
    neighNrFuncNorm=neighNrFunc./repmat(sqrt(sum(neighNrFunc.^2,3)),1,1,nf);
    actWeights=repmat(reshape(fcInMask(neighNrVox,:),1,nrVp,nf),neighNr,1,1);
    aFcWeightHi=squeeze(sum(neighNrFunc.*actWeights,3)./sum(actWeights,3));
    aFcWeight1=squeeze(sum(neighNrFunc(:,1,:).*actWeights(:,1,:),3)./sum(actWeights(:,1,:),3));
    aFcWeight2=squeeze(sum(neighNrFunc(:,2,:).*actWeights(:,2,:),3)./sum(actWeights(:,2,:),3));
    [sortedFc,idx]=sort(curFc,1);
    idx1Hi=idx(end-floor(nf/d)+1:end,1);idx2Hi=idx(end-floor(nf/d)+1:end,2);
    aFcWeightHi1=squeeze(sum(neighNrFunc(:,1,idx1Hi).*actWeights(:,1,idx1Hi),3)./sum(actWeights(:,1,idx1Hi),3));
    aFcWeightHi2=squeeze(sum(neighNrFunc(:,2,idx2Hi).*actWeights(:,2,idx2Hi),3)./sum(actWeights(:,2,idx2Hi),3));
    aFcWeightHi=aFcWeightHi./repmat(sqrt(sum(aFcWeightHi.^2,1)),neighNr,1);
    neighNrFuncNorm=permute(neighNrFuncNorm,[1,3,2]);
    % 2 May 16 end
    %neighNrFuncNorm=permute(neighNrFunc./repmat(sqrt(sum(neighNrFunc.^2,1)),neighNr,1,1),[1,3,2]);
    neighNrFuncNorm=reshape(neighNrFuncNorm,neighNr,[]);
    neighNrFuncNormSorted=reshape(neighNrFuncNorm(:,idx(:)),neighNr,nf,nrVp);
    loFc=repmat(reshape(sortedFc(1:floor(nf/d),:),1,[],nrVp),neighNr,1,1);
    aFcWeightLo=squeeze(sum(neighNrFuncNormSorted(:,1:floor(nf/d),:).*loFc,2)./sum(loFc,2));
    hiFc=repmat(reshape(sortedFc(end-floor(nf/d)+1:end,:),1,[],nrVp),neighNr,1,1);
    curFc=repmat(reshape(curFc,1,[],nrVp),neighNr,1,1);
    aFcWeightHiOld=squeeze(sum(neighNrFuncNormSorted(:,end-floor(nf/d)+1:end,:).*hiFc,2)./sum(hiFc,2));
    hiFc=repmat(reshape(sortedFc,1,[],nrVp),neighNr,1,1);
    aFcWeightHiNew=squeeze(sum(neighNrFuncNormSorted.*hiFc,2)./sum(hiFc,2));
    aFcWeightHiOld=aFcWeightHiOld./repmat(sqrt(sum(aFcWeightHiOld.^2,1)),neighNr,1);
    aFcWeightHiNew=aFcWeightHiNew./repmat(sqrt(sum(aFcWeightHiNew.^2,1)),neighNr,1);
    % weight with ALL volumes
    aFcWeightLo=aFcWeightLo./repmat(sqrt(sum(aFcWeightLo.^2,1)),neighNr,1);
    %aFcWeightHi=squeeze(sum(reshape(neighNrFuncNorm,neighNr,nf,nrVp).*curFc,2)./sum(curFc,2));
    % non-zero phase angle
    nzpInMaskLo(neighNrVox)=acosd(sum(aFcWeightLo,1)/sqrt(neighNr));
    nzpInMaskHi(neighNrVox)=acosd(sum(aFcWeightHi,1)/sqrt(neighNr));
    
    %
    %RO=reshape(neighNrFunc,neighNr,[]);
    clear neighNrFunc
    %bo=ones(1,neighNr)/sqrt(neighNr);
    % wrong! picSlideInMask(neighNrVox,:)=reshape(sqrt(sum(RO-bo'*(bo*RO).^2,1)),nrVp,[]);
picSlideInMask(neighNrVox,:)=reshape(corrColumns(neighNrFuncNorm,reshape(repmat(aFcWeightHi,nf,1),neighNr,[])),nf,nrVp)';
clear neighNrFuncNorm
%% magnitude orthogonal to bo projection:
%  rejRO=RO - repmat(mean(RO,1),neighNr,1); % rejection is orthogonal to projection
%  clear RO
%  picSlideInMask(neighNrVox,:)=reshape(sqrt(sum((rejRO).^2,1)),nrVp,[]);
%  clear rejRO
%% Correlation of succesive differences
%ROpermute=reshape(permute(neighNrFunc,[1,3,2]),neighNr,[]);
%clear RO neighNrFunc
%corrWithExtraSamples=reshape([corrColumns(ROpermute(:,1:end-2)-ROpermute(:,2:end-1),ROpermute(:,2:end-1)-ROpermute(:,3:end)),1,1],[],nrVp)';
%picSlideInMask(neighNrVox,:)=abs([corrWithExtraSamples(:,1:end-2),repmat(mean(abs(corrWithExtraSamples(:,1:end-2)),2),1,2)]);
%% Correlation of succesive position vectors - Not used bec doesn't oscillate through origin
%corrWithExtraSample=reshape([corrColumns(ROpermute(:,1:end-1),ROpermute(:,2:end)),1],[],nrVp)';
%picSlideInMask(neighNrVox,:)=abs([corrWithExtraSample(:,1:end-1),mean(abs(corrWithExtraSample(:,1:end-1)),2)]);
%% Correlation of succesive residuals
%rejRO=ROpermute - repmat(mean(ROpermute,1),neighNr,1);
%corrWithExtraSamples=reshape([corrColumns(rejRO(:,1:end-1),rejRO(:,2:end)),1],[],nrVp)';
% equivalent as corr independent of mean, bec corr first zero-centers variables
%corrWithExtraSamples=reshape([corrColumns(ROpermute(:,1:end-1),ROpermute(:,2:end)),1],[],nrVp)';
% abs version
%picSlideInMask(neighNrVox,:)=abs([corrWithExtraSamples(:,1:end-1),mean(abs(corrWithExtraSamples(:,1:end-1)),2)]);
% non abs version
%picSlideInMask(neighNrVox,:)=[corrWithExtraSamples(:,1:end-1),mean(corrWithExtraSamples(:,1:end-1),2)];
 
  %bo=repmat(ones(neighNr,1)/sqrt(neighNr),1,size(RO,2));
    % cosAngle between RO and bo:
  %picSlideInMask(neighNrVox,:)=reshape(sum(RO.*bo,1)./(sqrt(sum(RO.^2,1))),nrVp,[]);
    
end
% seedNr=1;
% [~,seedInd]=ismember(seed(seedNr).ind,maskInd);
% localIdx=sum(neighNrVox(1:seedInd));
% rx=radius2d(localIdx,:);
% pcc=squeeze(residual3d(:,localIdx,:))';
% [~,score]=pca(pcc);
% figure,plot2dcolor(score(:,1),score(:,2))

clear neighNrFunc
% in-volume results
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
ctrDir=fullfile(filePrefixDir,'ctr');
if ~exist(ctrDir,'dir')
    mkdir(ctrDir);
end

func2d=funcS.vol;
func3d=reshape(func2d,[funcS.volsize,nf]);

% save percentage anti-phase
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=papInMask;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'pap.nii.gz'),1);
% save non-zero phase
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=nzpInMaskLo;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'nzpLo.nii.gz'),1);
funcS.vol=zeros(prod(funcS.volsize),1);
funcS.vol(maskInd)=nzpInMaskHi;
funcS.vol=reshape(funcS.vol,funcS.volsize);
MRIsave(funcS,funcS.vol,fullfile(picDir,'nzpHi.nii.gz'),1);
% save picSlideInMask
funcS.vol=zeros([prod(funcS.volsize),nf]);
funcS.vol(maskInd,:)=picSlideInMask;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(picDir,'pic.nii.gz'),nf);
% save fcInMask
funcS.vol=zeros([prod(funcS.volsize),nf]);
funcS.vol(maskInd,:)=fcInMask;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(fcDir,'fc.nii.gz'),nf);
funcS.vol=[];
% save ctrInMask
funcS.vol=zeros([prod(funcS.volsize),nf]);
funcS.vol(maskInd,:)=ctrInMask;
funcS.vol=reshape(funcS.vol,[funcS.volsize,nf]);
MRIsave(funcS,funcS.vol,fullfile(fcDir,'ctr.nii.gz'),nf);
funcS.vol=[];
% picFcCorr
picFcCorr=zeros(prod(funcS.volsize),1);
picSlideInMask=picSlideInMask';
fcInMask=fcInMask';
ctrInMask=ctrInMask';
% % abs of FC, as pic >=0
picFcCorr(maskS.vol(:)>0,:)=corrColumns(picSlideInMask,abs(fcInMask));
%picFcCorr(maskS.vol(:)>0,:)=corrColumns(picSlideInMask,fcInMask);
MRIsave(funcS,reshape(picFcCorr,funcS.volsize),fullfile(picFcCorrDir,[num2str(randNr),'.nii.gz']),1);
% seeds
picTsAll=zeros(nf,numel(seed));
fcTsAll=zeros(nf,numel(seed));
aall=zeros(nf,27,numel(seed));
ball=zeros(nf,27,numel(seed));
aFcWeightLo=zeros(27,numel(seed));
aFcWeightHi=zeros(27,numel(seed));
anormAll=zeros(nf,27,numel(seed));
bFcWeightLo=zeros(27,numel(seed));
bFcWeightHi=zeros(27,numel(seed));
bnormAll=zeros(nf,27,numel(seed));
nzpAngle=zeros(numel(seed),1);
idxsHiAll=zeros(floor(nf/2),numel(seed));
for seedNr=1:numel(seed)
    [~,seedInd]=ismember(seed(seedNr).ind,maskInd);
    if seedInd~=0
        picSeedDir=fullfile(picDir,seed(seedNr).name);
        fcSeedDir=fullfile(fcDir,seed(seedNr).name);
        ctrSeedDir=fullfile(ctrDir,seed(seedNr).name);
        if ~exist(picSeedDir,'dir')
            mkdir(picSeedDir);
        end
        if ~exist(fcSeedDir,'dir')
            mkdir(fcSeedDir);
        end
        if ~exist(ctrSeedDir,'dir')
            mkdir(ctrSeedDir);
        end
        pic=zeros(prod(funcS.volsize),1);
        fc=pic;
        ctr=pic;
        
        %seed(seedNr).name
        %figure,plot([fcInMask(:,seedInd),detrend(picSlideInMask(:,seedInd),'constant')])
        %figure,plotyy(1:nf,fcInMask(:,seedInd),1:nf,picSlideInMask(:,seedInd))
        %corr(fcInMask(:,seedInd),picSlideInMask(:,seedInd))
        picTs=picSlideInMask(:,seedInd);
        pic(maskS.vol(:)>0,:)=corr(picSlideInMask,picTs);
        MRIsave(funcS,reshape(pic,funcS.volsize),fullfile(picSeedDir,[num2str(randNr),'.nii.gz']),1);
        save(fullfile(picSeedDir,[num2str(randNr),'.txt']),'picTs','-ascii');
        fcTs=fcInMask(:,seedInd);
        fc(maskS.vol(:)>0,:)=corr(fcInMask,fcTs);
        MRIsave(funcS,reshape(fc,funcS.volsize),fullfile(fcSeedDir,[num2str(randNr),'.nii.gz']),1);
        save(fullfile(fcSeedDir,[num2str(randNr),'.txt']),'fcTs','-ascii');
        ctrTs=ctrInMask(:,seedInd);
        ctr(maskS.vol(:)>0,:)=corr(ctrInMask,ctrTs);
        MRIsave(funcS,reshape(ctr,funcS.volsize),fullfile(ctrSeedDir,[num2str(randNr),'.nii.gz']),1);
        save(fullfile(ctrSeedDir,[num2str(randNr),'.txt']),'ctrTs','-ascii');
        fcTsAll(:,seedNr)=fcTs;
        picTsAll(:,seedNr)=picTs;
        
        p=repmat(seed(seedNr).vox,27,1)+roiShape;
        a=zeros(nf,27);for i=1:27,a(:,i)=squeeze(func3d(p(i,2),p(i,1),p(i,3),:));end
        
        a=a-repmat(mean(a,1),nf,1);
        %normMode='none';
        switch normMode{1}
            case 'none'
                anorm=a;
            case 'regress'
                for i=1:27
                    [~,~,r] = regress(a(:,i),fcspaceTimeTs);
                    anorm(:,i)=r;
                end
            case 'space'
                anorm=a-repmat(mean(a,2),1,27);
            case 'spaceTime'
                anorm=a-repmat(mean(a,2),1,27);
                anorm=anorm./repmat(sqrt(sum(anorm.^2,1)),nf,1);
            case 'time'
                anorm=a./repmat(sqrt(sum(a.^2,1)),nf,1);
            case 'timeSpace'
                anorm=a./repmat(sqrt(sum(a.^2,1)),nf,1);
                anorm=anorm-repmat(mean(anorm,2),1,27);
        end
        anormAll(:,:,seedNr)=anorm;
        b=a-repmat(mean(a,2),1,27);
        aall(:,:,seedNr)=a;
        ball(:,:,seedNr)=b;
        
%         figure('WindowStyle','docked')
%         subplot(2,1,1);plot(1:288,a),subplot(2,1,2);plot(1:288,b)
        
        [~,idx]=sort(fcTs);
%                 figure('Windowstyle','docked');
%         %subplot(10,1,1:2);plot(b);xlim([1,nf]);
%         %subplot(10,1,3:5);
%         subplot(2,1,1);
%         imagesc(corr(a',a'));
%         %subplot(10,1,6:8);
%         subplot(2,1,2);
%         imagesc(corr(a(idx,:)',a(idx,:)'))
%         %subplot(10,1,9:10);
        %figure
        d=2;
        idxsLo=idx(1:floor(nf/d));
        idxsHi=idx(end-floor(nf/d)+1:end);
        idxsHiAll(:,seedNr)=idxsHi;
        % 2 May 2016
%         if seedNr==s
%             figure('Windowstyle','docked');
%             subplot(3,1,1);
%             imagesc(corr(anorm',anorm'));
%             colorbar; zlim([-1,1]);
%             title(['Voxel pattern autocorrelation at ',seed(seedNr).name]);
%             xlabel('Volumes');ylabel('Volumes');
%             subplot(3,1,2);
%             imagesc(corr(anorm(idx,:)',anorm(idx,:)'))
%             colorbar; zlim([-1,1]);
%             title('Autocorrelations sorted from low to high activation levels');
%             xlabel('Sorted volumes');ylabel('Sorted volumes');
%             subplot(3,1,3);
%             plot(anorm);xlim([1,nf]);colorbar
%         end
        
        aFcWeightLo(:,seedNr)=sum(anorm(idxsLo,:).*repmat(fcTs(idxsLo),1,27),1)/sum(fcTs(idxsLo));
        aFcWeightHi(:,seedNr)=sum(anorm(idxsHi,:).*repmat(fcTs(idxsHi),1,27),1)/sum(fcTs(idxsHi));
        aFcWeightLo(:,seedNr)=aFcWeightLo(:,seedNr)/sqrt(sum(aFcWeightLo(:,seedNr).^2));
        aFcWeightHi(:,seedNr)=aFcWeightHi(:,seedNr)/sqrt(sum(aFcWeightHi(:,seedNr).^2));
        % non-zero phase angle
        nzpAngle(seedNr)=acosd(sum(aFcWeightHi(:,seedNr))/sqrt(27));
        % 2 May 2016
        bnorm=b./repmat(sqrt(sum(b.^2,1)),nf,1);
        %bnorm=b./repmat(sqrt(sum(b.^2,2)),1,27);
        bFcWeightLo(:,seedNr)=sum(bnorm(idxsLo,:).*repmat(fcTs(idxsLo),1,27),1)/sum(fcTs(idxsLo));
        bFcWeightHi(:,seedNr)=sum(bnorm(idxsHi,:).*repmat(fcTs(idxsHi),1,27),1)/sum(fcTs(idxsHi));
        bFcWeightLo(:,seedNr)=bFcWeightLo(:,seedNr)/sqrt(sum(bFcWeightLo(:,seedNr).^2));
        bFcWeightHi(:,seedNr)=bFcWeightHi(:,seedNr)/sqrt(sum(bFcWeightHi(:,seedNr).^2));
        bnormAll(:,:,seedNr)=bnorm;
        %         plot(reshape(repmat(1:27,floor(nf/d),1)+0.125*(rand(floor(nf/d),27)-1)-0.125,[],1),reshape(bnorm(idxLo,:),[],1),'b.')
%         hold on
%         plot(reshape(repmat(1:27,floor(nf/d),1)+0.125*(rand(floor(nf/d),27))+0.125,[],1),reshape(bnorm(idxHi,:),[],1),'r.')
%         plot(1:27,bFcWeightLo(:,seedNr),'bo',1:27,bFcWeightHi(:,seedNr),'ro');
%         hold off
%         %figure
%         plot(corr(bFcWeightLo(:,seedNr),bnorm'),corr(bFcWeightHi(:,seedNr),bnorm'));
%         disp([corr(fcTs,corr(bFcWeightLo(:,seedNr),bnorm')'),corr(fcTs,corr(bFcWeightHi(:,seedNr),bnorm')')]);
%         disp([num2str(seedNr),' ',seed(seedNr).name]);
% if seedNr==s
%    figure('WindowStyle','docked')
%    plot([aFcWeightLo(:,seedNr),aFcWeightHi(:,seedNr)])
%    xlim([1,27]);xlabel('Voxels')
%    title('Normalized activation weighted patterns');
%    legend('Low activation weighted pattern','High activation weighted pattern');
% end
% figure('WindowStyle','docked')
% plot([aFcWeightLo(:,seedNr),aFcWeightHi(:,seedNr)])
% xlim([1,27]);xlabel('Voxels')
% title('Normalized activation weighted patterns');
% legend('Low activation weighted pattern','High activation weighted pattern');
    end
end
for ss=1:numel(seed)
    aa=zeros(27,numel(seed));
    figure('Windowstyle','docked','NumberTitle','off','Name',num2str(ss))
    for tt=1:numel(seed)
        d=2;
        %[~,idxSs]=sort(fcTsAll(:,ss));idxsHiSs=idxSs(end-floor(nf/d)+1:end);
        [~,idxTt]=sort(fcTsAll(:,tt));idxsHiTt=idxTt(end-floor(nf/d)+1:end);
        aa(:,tt)=sum(anormAll(idxsHiTt,:,ss).*repmat(fcTsAll(idxsHiTt,tt),1,27),1)...
            /sum(fcTsAll(idxsHiTt,tt));
        aa(:,tt)=aa(:,tt)/sqrt(sum(aa(:,tt).^2));
    end;
    imagesc(corr(aa));colorbar;
end

for seedNr=1:numel(seed)
    figure('Windowstyle','docked');
    for i=1:27
        subplot(5,6,i);
        h=plotyy(aall(:,i,seedNr),fcTsAll(:,seedNr),aall(:,i,seedNr),corr(aFcWeightHi(:,seedNr),anormAll(:,:,seedNr)')');
        axis(h(1),'off');axis(h(2),'off');
    end
end
for ss=1:numel(seed)
    figure('WindowStyle','docked');title(normMode);
    subplot(1,3,1);plotyy(1:288,corr(aFcWeightHi(:,ss),anormAll(:,:,ss)')',1:288,fcTsAll(:,ss));
    subplot(1,3,2);plot([aFcWeightLo(:,ss),aFcWeightHi(:,ss)]);xlim([1,27]);
    subplot(1,3,3);plot(corr(aFcWeightLo(:,ss),anormAll(:,:,ss)')',corr(aFcWeightHi(:,ss),anormAll(:,:,ss)')');
    xlim([-1,1]); ylim([-1,1]);
end
if false
    bFcWeightHi=ones(27,numel(seed))/sqrt(27);
    his=sum(repmat(bFcWeightHi(:,s),1,nf).*bnormAll(:,:,s)',1);
    hit=sum(repmat(bFcWeightHi(:,t),1,nf).*bnormAll(:,:,t)',1);
    hiu=sum(repmat(bFcWeightHi(:,u),1,nf).*bnormAll(:,:,u)',1);
    hiv=sum(repmat(bFcWeightHi(:,v),1,nf).*bnormAll(:,:,v)',1);
    hiw=sum(repmat(bFcWeightHi(:,w),1,nf).*bnormAll(:,:,w)',1);
else
    his=corr(aFcWeightHi(:,s),anormAll(:,:,s)');
    hit=corr(aFcWeightHi(:,t),anormAll(:,:,t)');
    hiu=corr(aFcWeightHi(:,u),anormAll(:,:,u)');
    hiv=corr(aFcWeightHi(:,v),anormAll(:,:,v)');
    hiw=corr(aFcWeightHi(:,w),anormAll(:,:,w)');
    lis=corr(aFcWeightLo(:,s),anormAll(:,:,s)');
    lit=corr(aFcWeightLo(:,t),anormAll(:,:,t)');
    liu=corr(aFcWeightLo(:,u),anormAll(:,:,u)');
    liv=corr(aFcWeightLo(:,v),anormAll(:,:,v)');
    liw=corr(aFcWeightLo(:,w),anormAll(:,:,w)');
    figure
    plot(lis,his)
    xlabel('Correlation with low activation pattern');
    ylabel('Correlation with high activation pattern');
    title('Trajectory of RpIPS voxel pattern across time')
    figure
    plot(sum(anormAll(:,:,s)',1)./(sqrt(27)*sqrt(sum(anormAll(:,:,s)'.^2,1))))
    
    xx=find(fcTsAll(1:end-1,s)>0 & fcTsAll(2:end,s)<0);
    %[~,xx]=findpeaks(fcTsAll(:,s),'minpeakprominence',100);
    figure
    plot(1:288,fcTsAll(:,s),'b-',xx,fcTsAll(xx,s),'r*')
    title('Peak activation at RpIPS')
    xlabel('Volumes')
    ylabel('Activation')
    
    
    %[~,xx]=findpeaks(his,'minpeakprominence',.2);
    xx1=[xx(1:end-1),xx(2:end)];
    ndp=9;
    vq=zeros(ndp+1,size(xx1,1),27);
    figure
    for i=1:27
        subplot(5,6,i);
        for xidx=1:size(xx1,1)
            np=xx1(xidx,2)-xx1(xidx,1);
            vv=anormAll(xx1(xidx,1):xx1(xidx,2),i,s);
            %vv=detrend(anormAll(xx1(xidx,1):xx1(xidx,2),i,s),'constant');
            %vv=vv/sqrt(sum(vv.^2));
            vq(:,xidx,i)=interp1(xx1(xidx,1):xx1(xidx,2),vv,xx1(xidx,1):np/ndp:xx1(xidx,2));
            hold on;
            %plot(vq(:,xidx,i));%
            plot(vq(:,xidx,i)/sqrt(sum(vq(:,xidx,i).^2)));
            hold off;
        end;
        hold on
        meanVp=mean(vq(:,:,i),2);
        %plot(1:ndp+1,meanVp,'k*-','linewidth',5);%
        plot(1:ndp+1,meanVp/sqrt(sum(meanVp.^2)),'k*-','linewidth',5);
        ylim([-1,1]);
        hold off
    end
    figure
    i=8;
    subplot(2,1,1);
    for xidx=1:(numel(xx)-1)
        hold on;
        plot(aall(xx(xidx):xx(xidx+1),i,s)/sqrt(sum(aall(xx(xidx):xx(xidx+1),i,s).^2)));
        hold off;
    end;
    subplot(2,1,2);
    for xidx=1:(numel(xx)-1)
        np=xx(xidx+1)-xx(xidx);
        vq(:,xidx,i)=interp1(xx(xidx):xx(xidx+1),anormAll(xx(xidx):xx(xidx+1),i,s),xx(xidx):np/ndp:xx(xidx+1));
        hold on;
        %plot(vq(:,xidx,i));%
        plot(vq(:,xidx,i)/sqrt(sum(vq(:,xidx,i).^2)));
        hold off;
    end;
    hold on
    meanVp=mean(vq(:,:,i),2);
    %plot(1:ndp+1,meanVp,'k*-','linewidth',5);%
    plot(1:ndp+1,meanVp/sqrt(sum(meanVp.^2)),'k*-','linewidth',5);
    ylim([-1,1]);
    hold off
  
    %subplot(5,6,28);
    dhis=diff(his);
    trp=find(dhis(1:end-1)>0 & dhis(2:end)<0);
    trn=find(dhis(1:end-1)<0 & dhis(2:end)>0);
    trp=trp(1:min([numel(trp),numel(trn)]));
    trn=trn(1:min([numel(trp),numel(trn)]));
    nrTr=numel(trp);
    sh=floor(sqrt(nrTr));wh=ceil(sqrt(nrTr))+1;
    figure
    for trNr=1:nrTr-1
        subplot(sh,wh,trNr);plot(lis(trp(trNr)+1:trp(trNr+1)+1)',his(trp(trNr)+1:trp(trNr+1)+1)');xlim([-1,1]);ylim([-1,1]);
    end
%     figure
%     for trNr=1:nrTr-1
%         subplot(sh,wh,trNr);plot(lit(trp(trNr)+1:trp(trNr+1)+1)',hit(trp(trNr)+1:trp(trNr+1)+1)');xlim([-1,1]);ylim([-1,1]);
%     end
%     figure
%     for trNr=1:nrTr-1
%         subplot(sh,wh,trNr);plot(liu(trp(trNr)+1:trp(trNr+1)+1)',hiu(trp(trNr)+1:trp(trNr+1)+1)');xlim([-1,1]);ylim([-1,1]);
%     end
%     figure
%     for trNr=1:nrTr-1
%         subplot(sh,wh,trNr);plot(liv(trp(trNr)+1:trp(trNr+1)+1)',hiv(trp(trNr)+1:trp(trNr+1)+1)');xlim([-1,1]);ylim([-1,1]);
%     end
%     figure
%     for trNr=1:nrTr-1
%         subplot(sh,wh,trNr);plot(liw(trp(trNr)+1:trp(trNr+1)+1)',hiw(trp(trNr)+1:trp(trNr+1)+1)');xlim([-1,1]);ylim([-1,1]);
%     end
%     figure
%     for trNr=1:nrTr-1
%         %subplot(sh,wh,trNr);plot(aall(tr(trNr):tr(trNr+1),:,s));xlim([1,tr(trNr+1)-tr(trNr)+1]);
%         subplot(sh,wh,trNr);plot(aall(trp(trNr)+1:trp(trNr+1)+1,:,s)');
%         hold on;plot(1:27,aall(trp(trNr)+1,:,s)','b','linewidth',2);
%         plot(1:27,aall(trn(trNr+1)+1,:,s)','r','linewidth',2);hold off;
%     end

    figure;plot(squeeze(mean(vq,2)))

    
%     figure
%     subplot(4,5,1),plot(lis,his,'.-');
%     subplot(4,5,2),plot(lit,hit,'.-');
%     subplot(4,5,3),plot(liu,hiu,'.-');
%     subplot(4,5,4),plot(liv,hiv,'.-');
%     subplot(4,5,5),plot(liw,hiw,'.-');
%     subplot(4,5,6),plot([lis',his']);
%     subplot(4,5,7),plot([lit',hit']);
%     subplot(4,5,8),plot([liu',hiu']);
%     subplot(4,5,9),plot([liv',hiv']);
%     subplot(4,5,10),plot([liw',hiw']);
%     subplot(4,5,11),plot(aall(:,:,s));
%     subplot(4,5,12),plot(aall(:,:,t));
%     subplot(4,5,13),plot(aall(:,:,u));
%     subplot(4,5,14),plot(aall(:,:,v));
%     subplot(4,5,15),plot(aall(:,:,w));
%     k=2;
%     subplot(4,5,16),plot(his,aall(:,k,s),'.');
%     subplot(4,5,17),plot(hit,aall(:,k,t),'.');
%     subplot(4,5,18),plot(hiu,aall(:,k,u),'.');
%     subplot(4,5,19),plot(hiv,aall(:,k,v),'.');
%     subplot(4,5,20),plot(hiw,aall(:,k,w),'.');
end
disp([seed(s).name,' ',seed(t).name,' ',seed(u).name,' ',seed(v).name,' ',seed(w).name]);
[N,edges,bin]=histcounts(lis);
vp=zeros(numel(N),27);
vpstd=zeros(numel(N),27);
for vpNr=1:numel(N)
    vp(vpNr,:)=mean(anormAll(bin==vpNr,:,s),1);
    vpstd(vpNr,:)=std(anormAll(bin==vpNr,:,s),0,1);
end
% efig=figure;
% subplot(2,1,1);plot(vp)
% subplot(2,1,2);plot(vpstd)
% for vNr=1:27
%     plot(edges(2:end),vp(:,vNr));
%     ylim([-1,1]);
%     E(vNr)=getframe(efig);
% end
%movie(efig,E,10,1)
fig=figure;
for x=1:288
    subplot(2,4,1);plot(his,hit,'b.-',his(x),hit(x),'ro');
    title(num2str(nzpAngle(s),2));
    xlabel([num2str(corr(his',hit'),2),' ',seed(t).name,' ',num2str(nzpAngle(t),2)]);
    subplot(2,4,3);plot(his,hiu,'b.-',his(x),hiu(x),'ro');
    xlabel([num2str(corr(his',hiu'),2),' ',seed(u).name,' ',num2str(nzpAngle(u),2)]);
    subplot(2,4,5);plot(his,hiv,'b.-',his(x),hiv(x),'ro');
    xlabel([num2str(corr(his',hiv'),2),' ',seed(v).name,' ',num2str(nzpAngle(v),2)]);
    subplot(2,4,7);plot(his,hiw,'b.-',his(x),hiw(x),'ro');
    xlabel([num2str(corr(his',hiw'),2),' ',seed(w).name,' ',num2str(nzpAngle(w),2)]);
    subplot(2,4,2);plot(fcTsAll(:,s),fcTsAll(:,t),'b.-',fcTsAll(x,s),fcTsAll(x,t),'ro');xlim([min(fcTsAll(:,s)),max(fcTsAll(:,s))]);ylim([min(fcTsAll(:,t)),max(fcTsAll(:,t))]);
    xlabel(num2str(corr(fcTsAll(:,s),fcTsAll(:,t)),2));
    subplot(2,4,4);plot(fcTsAll(:,s),fcTsAll(:,u),'b.-',fcTsAll(x,s),fcTsAll(x,u),'ro');xlim([min(fcTsAll(:,s)),max(fcTsAll(:,s))]);ylim([min(fcTsAll(:,u)),max(fcTsAll(:,u))]);
    xlabel(num2str(corr(fcTsAll(:,s),fcTsAll(:,u)),2));
    subplot(2,4,6);plot(fcTsAll(:,s),fcTsAll(:,v),'b.-',fcTsAll(x,s),fcTsAll(x,v),'ro');xlim([min(fcTsAll(:,s)),max(fcTsAll(:,s))]);ylim([min(fcTsAll(:,v)),max(fcTsAll(:,v))]);
    xlabel(num2str(corr(fcTsAll(:,s),fcTsAll(:,v)),2));
    subplot(2,4,8);plot(fcTsAll(:,s),fcTsAll(:,w),'b.-',fcTsAll(x,s),fcTsAll(x,w),'ro');xlim([min(fcTsAll(:,s)),max(fcTsAll(:,s))]);ylim([min(fcTsAll(:,w)),max(fcTsAll(:,w))]);
    xlabel(num2str(corr(fcTsAll(:,s),fcTsAll(:,w)),2));
    F(x)=getframe(fig);
    keyboard
end
movie(fig,F,10,6)
randNr=numel(dir(fullfile(fcSeedDir,'*.nii.gz')))-1;
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

% function [funcRoiLag,funcRoiLagDiff]=lag1func(func,minNeighNr)
% sz = 1;
% sizeFunc = size(func);
% funcRoiLag = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
% funcRoiLagDiff = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
% nrVoxels = (2*sz+1)^3;
% for x=1:size(func,1)
%     for y=1:size(func,2)
%         for z=1:size(func,3)
%             if all(func(x,y,z,:) ~= 0)
%                 roi = zeros(size(func,4),nrVoxels);
%                 idx = 0;
%                 for i=x-sz:x+sz
%                     for j=y-sz:y+sz
%                         for k=z-sz:z+sz
%                             if i>0 && i<=sizeFunc(1) && j>0 && j<=sizeFunc(2) && k>0 && k<=sizeFunc(3)
%                                 if all(func(i,j,k,:) ~= 0)
%                                     idx = idx + 1;
%                                     roi(:,idx) = squeeze(func(i,j,k,:));
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 if idx >= minNeighNr
%                     roi(:,idx+1:end) = [];
%                     funcRoiLag(x,y,z,:) = corrColumns(roi(1:end-1,:)',roi(2:end,:)');
%                     funcRoiLagDiff(x,y,z,:) = mean(diff(roi),2);
%                 end
%             end
%         end
%     end
% end
% end

% function [funcRoiLag,funcRoiLagDiff] = lag1func(func)
% sz = 1;
% sizeFunc = size(func);
% funcRoiLag = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
% funcRoiLagDiff = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
% nrVoxels = (2*sz+1)^3;
% for x=2:size(func,1)-1
%     for y=2:size(func,2)-1
%         for z=2:size(func,3)-1
%             if all(func(x,y,z,:) ~= 0)
%                 roi = zeros(size(func,4),nrVoxels);
%                 idx = 0;
%                 for i=x-sz:x+sz
%                     for j=y-sz:y+sz
%                         for k=z-sz:z+sz
%                             if all(func(i,j,k,:) ~= 0)
%                                 idx = idx + 1;
%                                 roi(:,idx) = squeeze(func(i,j,k,:));
%                             end
%                         end
%                     end
%                 end
%                 if idx >= 2
%                     roi(:,idx+1:end) = [];
%                     funcRoiLag(x,y,z,:) = corrColumns(roi(1:end-1,:)',roi(2:end,:)');
%                     roi = bsxfun(@minus,roi,mean(roi,1)); % remove temporal mean
%                     roiSpatialMean = mean(roi,2);
%                     funcRoiLagDiff(x,y,z,:) = diff(roiSpatialMean);
%                 end
%             end
%         end
%     end
% end
% end

% function funcRoiLag1 = lag1func(func)
% sz = 1;
% sizeFunc = size(func);
% func2d = reshape(func,[],sizeFunc(4))';
% %reshape(permute(func,[4,1,2,3]),[sizeFunc(4),sizeFunc(1:3)]);
% funcRoiLag1 = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
% [roiShape,nrVox] = getRoiShape(sz);
% for x=2:size(func,1)-1
%     for y=2:size(func,2)-1
%         for z=2:size(func,3)-1
%             if all(func(x,y,z,:) ~= 0)
%                 roiCoord = roiShape + repmat([x,y,z],[nrVox,1]);
%                 roiLinInd = sub2ind(sizeFunc(1:3),roiCoord(:,1),roiCoord(:,2),roiCoord(:,3));
%                 roi = func2d(:,roiLinInd);
%                 roi(:,all(roi == 0,1)) = [];
%                 roiRes = roiResidual(roi);
%                 funcRoiLag1(x,y,z,:) = corrColumns(roiRes(1:end-1,:)',roiRes(2:end,:)');
%             end
%         end
%     end
% end
% end
% function roiOut = roiResidual(roi)
% roi = roi - repmat(mean(roi,1),[size(roi,1),1]); % remove temporal mean
% roiOut = roi - repmat(mean(roi,2),[1,size(roi,2)]); % remove spatial mean
% end