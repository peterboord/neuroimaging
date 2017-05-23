function randPhase(dirName,subjCodeNrStr,maskName,nrRand)
% from SGE_roiCorr52
if ~isdeployed
    dbstop if error
end
if nargin < 4
    nrRand=10;
end
if nargin < 3
    maskName='mask.nii';
end
if nargin < 2
    subjCodeNrStr='3211';
end
if nargin < 1
    dirName='/project_space/pboord/PIC/iowa/raw/3211';
end
%% constants
sz = 1;
%mniRef = '/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz';
minNeighNr=3;
%% seeds
seed(1).name = 'PCC';
seed(1).mni = [-6,-58,28];
seed(2).name = 'mPFC';
seed(2).mni = [-2,50,17];
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
seed(12).name = 'WM';
seed(12).mni = [30,-46,28];
%%
fileName=[subjCodeNrStr,'_reg_filtered_func_data.nii'];
funcS=MRIread(fullfile(dirName,fileName));
nrVox=prod(funcS.volsize);
halfUniqueFftLength = (funcS.nframes/2)-1;
maskS=MRIread(fullfile(dirName,maskName));
funcS.vol=reshape(funcS.vol,[],funcS.nframes);
fastPic(0,funcS,maskS,func2d,sz,minNeighNr,seed);
absFftFunc = abs(fft(func2d'));
clear func2d
for randNr=1:nrRand
    rng(str2double(subjCodeNrStr)+randNr);
    randPhase2d = reshape(pi*(2*rand(halfUniqueFftLength*nrVox,1) - 1),halfUniqueFftLength,nrVox);
    tic
    fastPic(randNr,funcS,maskS,...
        ifft(absFftFunc.*exp(1i*[zeros(1,nrVox);randPhase2d;pi*(2*rand(1,nrVox)-1);-1*flipud(randPhase2d)]),'symmetric')',...
        sz,minNeighNr,seed);
    toc
end
end

function fastPic(randNr,funcS,maskS,func,sz,minNeighNr,seed)
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
neighFunc=func(neighInd,:);
neighFunc(neighInd==1,:)=NaN;
nrNeigh=sum(reshape(neighInd~=1,nrNeighMax,[]));
neighFunc=reshape(neighFunc,nrNeighMax,[],funcS.nframes);
uniqueNrNeigh=unique(nrNeigh);
uniqueNrNeigh=uniqueNrNeigh(uniqueNrNeigh>=minNeighNr);
% in-mask results
picInMask=zeros(nrMaskVox,funcS.nframes-1);
fcDiffInMask=zeros(nrMaskVox,funcS.nframes-1);
for neighNr=uniqueNrNeigh
    neighNrVox=nrNeigh==neighNr;
    neighNrFunc=neighFunc(:,neighNrVox,:);
    neighNrFunc=reshape(neighNrFunc(~isnan(neighNrFunc)),neighNr,[],funcS.nframes);
    picInMask(neighNrVox,:)=0.5-0.5*reshape(corrColumns(reshape(neighNrFunc(:,:,1:end-1),neighNr,[]),reshape(neighNrFunc(:,:,2:end),neighNr,[])),sum(neighNrVox),[]);
    fcDiffInMask(neighNrVox,:)=reshape(mean(neighNrFunc(:,:,2:end)-neighNrFunc(:,:,1:end-1),1),sum(neighNrVox),[]);
end
% in-volume results
% pic=zeros(prod(funcS.volsize),funcS.nframes-1);
% pic(maskS.vol(:)>0,:)=picInMask;
% pic=reshape(pic,[funcS.volsize,funcS.nframes-1]);
% fcDiff=zeros(prod(funcS.volsize),funcS.nframes-1);
% fcDiff(maskS.vol(:)>0,:)=fcDiffInMask;
% fcDiff=reshape(fcDiff,[funcS.volsize,funcS.nframes-1]);
picFcCorr=zeros(prod(funcS.volsize),1);
picInMask=picInMask';
fcDiffInMask=fcDiffInMask';
picFcCorr(maskS.vol(:)>0,:)=corrColumns(picInMask,fcDiffInMask);
picFcCorr=reshape(picFcCorr,funcS.volsize);
validatePicFc=0;
if validatePicFc
    tic %#ok<UNRCH>
    func(~maskS.vol(:),:)=0;
    [funcRoiLag,funcRoiLagDiff]=lag1func(reshape(func,[funcS.volsize,funcS.nframes]),minNeighNr);
    funcRoiLag=reshape(funcRoiLag,[],funcS.nframes-1);
    all(pic(:)==funcRoiLag(:))
    all(fcDiff(:)==funcRoiLagDiff(:))
    toc
end
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

function [funcRoiLag,funcRoiLagDiff]=lag1func(func,minNeighNr)
sz = 1;
sizeFunc = size(func);
funcRoiLag = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
funcRoiLagDiff = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
nrVoxels = (2*sz+1)^3;
for x=1:size(func,1)
    for y=1:size(func,2)
        for z=1:size(func,3)
            if all(func(x,y,z,:) ~= 0)
                roi = zeros(size(func,4),nrVoxels);
                idx = 0;
                for i=x-sz:x+sz
                    for j=y-sz:y+sz
                        for k=z-sz:z+sz
                            if i>0 && i<=sizeFunc(1) && j>0 && j<=sizeFunc(2) && k>0 && k<=sizeFunc(3)
                                if all(func(i,j,k,:) ~= 0)
                                    idx = idx + 1;
                                    roi(:,idx) = squeeze(func(i,j,k,:));
                                end
                            end
                        end
                    end
                end
                if idx >= minNeighNr
                    roi(:,idx+1:end) = [];
                    funcRoiLag(x,y,z,:) = corrColumns(roi(1:end-1,:)',roi(2:end,:)');
                    funcRoiLagDiff(x,y,z,:) = mean(diff(roi),2);
                end
            end
        end
    end
end
end

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