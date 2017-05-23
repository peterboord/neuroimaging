function [rmsSlice,voxMask]=getRmsSlice(rest,nrMaskVoxPerSlice,voxMask)
restSize=size(rest);
volsize=restSize(1:3);
nvol=restSize(4);
nslice=restSize(3);
rest=reshape(rest,prod(volsize(1:2)),nslice,nvol);
if nargin < 3
    nrSliceVox=sum(~all(rest==0,3),1);
    percentNrSliceVox=(100*50)./nrSliceVox;
    % make mask from 50 voxels/slice with highest Tstd
    stdRest=std(rest,1,3);
    [stdRestSize,stdRestSortIdx]=sort(stdRest,1,'descend');
    if any(reshape(stdRestSize(1:nrMaskVoxPerSlice,:),[],1)==0)
        disp('WARNING: mask has voxels with Tvar = 0');
    end
    topvarIdxs=stdRestSortIdx(1:nrMaskVoxPerSlice,:);
    voxMask=zeros(prod(volsize(1:2)),nslice);
    for sliceNr=1:nslice
        voxIdxs=ismember((1:prod(volsize(1:2)))',topvarIdxs(:,sliceNr));
        voxMask(voxIdxs,sliceNr)=stdRest(voxIdxs,sliceNr);
    end
end

%stdRest=std(rest,1,3);
%rest=rest./repmat(reshape(stdRest.*(voxMask>0),prod(volsize(1:2)),nslice,1),1,1,nvol);
rest=rest./repmat(reshape(voxMask,prod(volsize(1:2)),nslice,1),1,1,nvol);

rest(isinf(rest))=0;
rest(isnan(rest))=0;
rmsSlice=zeros(nslice,nvol);
for sliceNr=1:nslice
    cohRestSlice=squeeze(rest(voxMask(:,sliceNr)>0,sliceNr,:));
    
    % rms calculation includes zero voxels!!!
    %rmsSlice(sliceNr,:)=sqrt(mean(cohRestSlice.^2,1));
    rmsSlice(sliceNr,:)=sqrt(sum(cohRestSlice(~all(cohRestSlice==0,2),:).^2,1)/sum(~all(cohRestSlice==0,2)));

% ????? necessary?
%rmsSlice(sliceNr,:)=detrend(rmsSlice(sliceNr,:));
    rmsSlice(sliceNr,:)=rmsSlice(sliceNr,:)/std(rmsSlice(sliceNr,:));

end

% ???? detrend before normalization!?!
rmsSlice=detrend(rmsSlice')';

end
