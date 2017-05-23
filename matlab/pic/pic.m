function pic(dirName,fileName,maskName)
tic
if nargin < 2
    disp('Usage: run_pic ${MCRROOT} dirName fileName maskName');
else
    imInS = MRIread(fullfile(dirName,fileName));
    if nargin == 2
        [funcRoiLag,funcRoiLagDiff,funcRoiLagRes,funcRoiLagResSlice] = lag1func(imInS.vol);
    else
        maskS = MRIread(fullfile(dirName,maskName));
        [funcRoiLag,funcRoiLagDiff,funcRoiLagRes,funcRoiLagResSlice]=lag1func(imInS.vol.*repmat(reshape((maskS.vol>0)*1,[maskS.volsize,1]),[1,1,1,imInS.nframes]));
    end
    imInS.nframes = imInS.nframes - 1;
    imInS.dim(5) = imInS.nframes;
    imInS.vol = funcRoiLag;
    MRIwrite(imInS,fullfile(dirName,['pic_',fileName]),'float');
    imInS.vol = funcRoiLagDiff;
    MRIwrite(imInS,fullfile(dirName,['fcDiff_',fileName]),'float');
%     imInS.vol = funcRoiLagRes;
%     MRIwrite(imInS,fullfile(dirName,['picRes_',fileName]),'float');
%     imInS.vol = funcRoiLagResSlice;
%     MRIwrite(imInS,fullfile(dirName,['picResSlice_',fileName]),'float');
    % Fisher xfm
    imInS.vol = atanh(funcRoiLag);
    MRIwrite(imInS,fullfile(dirName,['z_pic_',fileName]),'float');
    imInS.vol = atanh(funcRoiLagDiff);
    MRIwrite(imInS,fullfile(dirName,['z_fcDiff_',fileName]),'float');
%     imInS.vol = atanh(funcRoiLagRes);
%     MRIwrite(imInS,fullfile(dirName,['z_picRes_',fileName]),'float');
%     imInS.vol = atanh(funcRoiLagResSlice);
%     MRIwrite(imInS,fullfile(dirName,['z_picResSlice_',fileName]),'float');
    toc
end
end
fastPic
function [funcRoiLag,funcRoiLagDiff,funcRoiLagRes,funcRoiLagResSlice] = lag1func(func)
sz = 1;
sizeFunc = size(func);
funcRoiLag = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
funcRoiLagDiff = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
funcRoiLagRes = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
funcRoiLagResSlice = zeros([sizeFunc(1:3),sizeFunc(4)-1]);
nrVoxels = (2*sz+1)^3;
nrVoxelsSlice = (2*sz+1)^2;
for x=2:size(func,1)-1
    for y=2:size(func,2)-1
        for z=2:size(func,3)-1
            if all(func(x,y,z,:) ~= 0)
                roi = zeros(size(func,4),nrVoxels);
                roiSlice = zeros(size(func,4),nrVoxelsSlice);
                idx = 0;
                idxSlice = 0;
                for i=x-sz:x+sz
                    for j=y-sz:y+sz
                        for k=z-sz:z+sz
                            if all(func(i,j,k,:) ~= 0)
                                idx = idx + 1;
                                roi(:,idx) = squeeze(func(i,j,k,:));
                                if k == z
                                    idxSlice = idxSlice + 1;
                                    roiSlice(:,idxSlice) = squeeze(func(i,j,k,:));
                                end
                            end
                        end
                    end
                end
                if idx >= 2
                    roi(:,idx+1:end) = [];
                    funcRoiLag(x,y,z,:) = corrColumns(roi(1:end-1,:)',roi(2:end,:)');
                    roi = bsxfun(@minus,roi,mean(roi,1)); % remove temporal mean
                    roiSpatialMean = mean(roi,2);
                    roiRes = bsxfun(@minus,roi,roiSpatialMean); % remove spatial mean
                    funcRoiLagDiff(x,y,z,:) = diff(roiSpatialMean);
                    funcRoiLagRes(x,y,z,:) = corrColumns(roiRes(1:end-1,:)',roiRes(2:end,:)');
                end
                if idxSlice >= 2
                    roiSlice(:,idxSlice+1:end) = [];
                    roiSlice = bsxfun(@minus,roiSlice,mean(roiSlice,1)); % remove temporal mean
                    roiSliceSpatialMean = mean(roiSlice,2);
                    roiResSlice = bsxfun(@minus,roiSlice,roiSliceSpatialMean); % remove spatial mean
                    funcRoiLagResSlice(x,y,z,:) = corrColumns(roiResSlice(1:end-1,:)',roiResSlice(2:end,:)');
                end
            end
        end
    end
end
end

