function nonparamPic(dirName,subjCodeNrStr,filePrefix,maskName,nrRand)
% from SGE_roiCorr52
if ~isdeployed
    dbstop if error
end
if nargin < 5
    nrRand='0';
end
if nargin < 4
    maskName='mask.nii.gz';
end
if nargin < 3
    filePrefix='regress.hpf';
    %filePrefix='volreg';
    %filePrefix='tshift';
    %filePrefix='ricor';
    %filePrefix='despike';
    %filePrefix='tcat';
end
if nargin < 2
    %subjCodeNrStr='3211';
    subjCodeNrStr='109003';
end
if nargin < 1
    % dirName='/var/local/pboord_5TB/pic/iowa/raw';
    %dirName='/project_space/pboord/PIC/iowa/raw';
    dirName='/project_space/pboord/act/rest';
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
seed(9).name = 'salACC';
seed(9).mni = [7,33,19];
seed(10).name = 'salRFIC';
seed(10).mni = [34,26,-6];
seed(11).name = 'salLFIC';
seed(11).mni = [-32,25,-10];
seed(12).name = 'mPFC';
seed(12).mni = [-2,50,17];
%%
funcS=MRIread(fullfile(dirName,subjCodeNrStr,[filePrefix,'.nii.gz']));
for seedNr=1:numel(seed)
	seed(seedNr).vox=round(load(fullfile(dirName,subjCodeNrStr,['func_',seed(seedNr).name]),'-ascii'));
    y=seed(seedNr).vox(1); % n.b. swapped
    x=seed(seedNr).vox(2); % x/y for FSL->FS
    z=seed(seedNr).vox(3);
    seed(seedNr).ind=sub2ind(funcS.volsize,x,y,z);
end
nrVox=prod(funcS.volsize);
halfUniqueFftLength = (funcS.nframes/2)-1;
maskS=MRIread(fullfile(dirName,subjCodeNrStr,maskName));
funcS.vol=reshape(funcS.vol,[],funcS.nframes);

[pathstr,name]=fileparts(funcS.fspec);
fastPic(0,funcS,maskS,[],sz,minNeighNr,seed,pathstr,name);

absFftFunc = abs(fft(funcS.vol'));
%angleFftFunc = angle(fft(funcS.vol'));
funcS.vol=[];
randNr=0;
while randNr<nrRand
    randNr=randNr+1;
    rng(str2double(subjCodeNrStr)+randNr);
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
neighFunc=reshape(neighFunc,nrNeighMax,[],funcS.nframes);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
picSlideInMask=zeros(nrMaskVox,funcS.nframes);
fcInMask=picSlideInMask;
for neighNr=uniqueNrNeigh
    neighNrVox=nrNeigh==neighNr;
    if neighNr==27
        % conserve memory
        neighFunc(:,~neighNrVox,:)=[];
        neighNrFunc=neighFunc;
        clear neighFunc
    else
        neighNrFunc=neighFunc(:,neighNrVox,:);
    end
    neighNrFunc(isnan(neighNrFunc))=[];
    neighNrFunc=reshape(neighNrFunc,neighNr,[],funcS.nframes);
    fcInMask(neighNrVox,:)=reshape(mean(neighNrFunc,1),sum(neighNrVox),[]);
    RO=reshape(neighNrFunc,neighNr,[]);
    clear neighNrFunc
    projRO=repmat(mean(RO,1),neighNr,1);
    %bo=ones(1,neighNr)/sqrt(neighNr);
    % wrong! picSlideInMask(neighNrVox,:)=reshape(sqrt(sum(RO-bo'*(bo*RO).^2,1)),sum(neighNrVox),[]);
    picSlideInMask(neighNrVox,:)=reshape(sqrt(sum((RO-projRO).^2,1)),sum(neighNrVox),[]);
    %%picSlideInMask(neighNrVox,:)=reshape(sign(projRO(1,:)).*sum(RO.*projRO,1)./(sqrt(sum(RO.^2,1)).*sqrt(sum(projRO.^2,1))),sum(neighNrVox),[]);
    %bo=repmat(ones(neighNr,1)/sqrt(neighNr),1,size(RO,2));
    %picSlideInMask(neighNrVox,:)=reshape(sum(RO.*bo,1)./(sqrt(sum(RO.^2,1))),sum(neighNrVox),[]);
    
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
% picFcCorr
picFcCorr=zeros(prod(funcS.volsize),1);
picSlideInMask=picSlideInMask';
fcInMask=fcInMask';
picFcCorr(maskS.vol(:)>0,:)=corrColumns(picSlideInMask,fcInMask);
MRIsave(funcS,reshape(picFcCorr,funcS.volsize),fullfile(picFcCorrDir,[num2str(randNr),'.nii.gz']),1);
% seeds
for seedNr=1:numel(seed)
    [~,seedInd]=ismember(seed(seedNr).ind,maskInd);
    picSeedDir=fullfile(picDir,seed(seedNr).name);
    fcSeedDir=fullfile(fcDir,seed(seedNr).name);
    if ~exist(picSeedDir,'dir')
        mkdir(picSeedDir);
    end
    if ~exist(fcSeedDir,'dir')
        mkdir(fcSeedDir);
    end
    pic=zeros(prod(funcS.volsize),1);
    fc=pic;
    
    %seed(seedNr).name
    %figure,plot([fcInMask(:,seedInd),detrend(picSlideInMask(:,seedInd),'constant')])
    %figure,plotyy(1:funcS.nframes,fcInMask(:,seedInd),1:funcS.nframes,picSlideInMask(:,seedInd))
    %corr(fcInMask(:,seedInd),picSlideInMask(:,seedInd))
    picTs=picSlideInMask(:,seedInd);
    pic(maskS.vol(:)>0,:)=corr(picSlideInMask,picTs);
    MRIsave(funcS,reshape(pic,funcS.volsize),fullfile(picSeedDir,[num2str(randNr),'.nii.gz']),1);
    save(fullfile(picSeedDir,[num2str(randNr),'.txt']),'picTs','-ascii');
    fcTs=fcInMask(:,seedInd);
    fc(maskS.vol(:)>0,:)=corr(fcInMask,fcTs);
    MRIsave(funcS,reshape(fc,funcS.volsize),fullfile(fcSeedDir,[num2str(randNr),'.nii.gz']),1);
    save(fullfile(fcSeedDir,[num2str(randNr),'.txt']),'fcTs','-ascii');
end
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