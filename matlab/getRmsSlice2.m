function [rmsSlice,voxMask]=getRmsSlice2(rest,percentSliceVox,voxMask)
restSize=size(rest);
volsize=restSize(1:3);
nvol=restSize(4);
nslice=restSize(3);
rest=reshape(rest,prod(volsize(1:2)),nslice,nvol);
if nargin<3
    nrSliceVox=sum(~all(rest==0,3),1);
    stdRest=std(rest,1,3);
    [~,stdRestSortIdx]=sort(stdRest,1,'descend');
    voxMask=zeros(prod(volsize(1:2)),nslice);
    for sliceNr=1:nslice
        topvarIdxs=stdRestSortIdx(1:max([1,round(percentSliceVox*nrSliceVox(sliceNr)/100)]),:);
        voxIdxs=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
        voxMask(voxIdxs,sliceNr)=stdRest(voxIdxs,sliceNr);
    end
end
rest=rest./repmat(reshape(voxMask,prod(volsize(1:2)),nslice,1),1,1,nvol);
rest(isinf(rest))=0;
rest(isnan(rest))=0;
rmsSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    sumSliceSquared=sum(squeeze(rest(voxMask(:,sliceNr)>0,sliceNr,:).^2),1);
    rmsSlice(sliceNr,:)=bsxfun(@rdivide,sumSliceSquared,median(sumSliceSquared,2));

%     cohRestSlice=squeeze(rest(voxMask(:,sliceNr)>0,sliceNr,:));
%     
%     % rms calculation includes zero voxels!!!
%     %rmsSlice(sliceNr,:)=sqrt(mean(cohRestSlice.^2,1));
%     rmsSlice(sliceNr,:)=sqrt(sum(cohRestSlice(~all(cohRestSlice==0,2),:).^2,1)/sum(~all(cohRestSlice==0,2)));
% 
% % ????? necessary?
% %rmsSlice(sliceNr,:)=detrend(rmsSlice(sliceNr,:));
%     rmsSlice(sliceNr,:)=rmsSlice(sliceNr,:)/std(rmsSlice(sliceNr,:));

end

% % ???? detrend before normalization!?!
% rmsSlice=detrend(rmsSlice')';

end
