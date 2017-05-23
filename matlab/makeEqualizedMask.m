function equalizedMask=makeEqualizedMask(restS,abMag,abPhase,nrMaskVoxPerSlice)
disp('fix makeEqualizedMask to use percentSliceVox instead of nrMaskVoxPerSlice');
% create new varMask with even phase dsn (or pi/2 spread) & even amplitude
% dsns across phase
% partition phase dsn into 8 bins
% fill bins in order of mag until 3 neighboring (or anti-neighboring) bins have min N voxels
% equalize mag in each phase bin
% phase dsn
volsize=restS.volsize;
nslice=volsize(3);
nvol=restS.nframes;
nrBins=12;
equalizedMask=zeros(prod(volsize(1:2)),nslice);
for sliceNr=1:nslice
    repAbPhase=repmat(abPhase(:,sliceNr,1),1,nrBins);
    phaseRange=-pi:2*pi/nrBins:pi*(nrBins-1)/nrBins;
    phaseBins=repmat(phaseRange,prod(volsize(1:2)),1);
    % make prelim mask to use vox with max abMag (rather than max var)
    [abMagSize,abMagSortIdx]=sort(abMag(:,sliceNr,1),1,'descend');
    if any(abMagSize(1:nrMaskVoxPerSlice,:)==0)
        disp('WARNING: mask has voxels with abMag = 0');
    end
    sliceMask=ismember((1:prod(volsize(1:2)))',abMagSortIdx(1:nrMaskVoxPerSlice));
    
    slicePhases=repAbPhase>=phaseBins & repAbPhase<(phaseBins+2*pi/nrBins) & repmat(~all(reshape(restS.vol(:,:,sliceNr,:)==0,[],nvol),2),1,nrBins)...
        & repmat(sliceMask,1,nrBins);
    sliceMag=(slicePhases(:,1:nrBins/2)|slicePhases(:,nrBins/2+1:end)).*repmat(squeeze(abMag(:,sliceNr,1)),1,nrBins/2);
    [sliceMagSortSize,sliceMagSortIdx]=sort(sliceMag,1,'descend');
    % choose pi/2 range with max rms
    maxRms=0;
    nrVoxPerPhaseBin=5;
    while maxRms==0 && nrVoxPerPhaseBin>0
        shiftIdxs=logical([1,1,1,0,0,0]);
        
        shiftIdxs=logical([1,1,1,1,1,1]);
        
        for shiftNr=0:nrBins/2-1
            binShift=logical(circshift(shiftIdxs,shiftNr,2));
            magBins=sliceMagSortSize(:,binShift);
            magBins=magBins(1:nrVoxPerPhaseBin,:);
            magBinRms=sqrt(mean(magBins(:).^2));
            
            if 1%all(magBins(:)~=0) && magBinRms>maxRms
                
                maxRms=magBinRms;
                maxRmsBinShift=binShift;
                maxMagIdxs=sliceMagSortIdx(1:nrVoxPerPhaseBin,maxRmsBinShift).*(magBins>0);
                equalizedMask(:,sliceNr)=1*ismember((1:prod(volsize(1:2)))',reshape(maxMagIdxs,[],1));
                normalizedMagBins=magBins(:)./reshape(repmat(sqrt(sum(magBins.^2,1)./sum(magBins>0,1)),nrVoxPerPhaseBin,1),[],1);
                equalizedMask(equalizedMask(:,sliceNr)>0,sliceNr)=normalizedMagBins(normalizedMagBins>0);
                %equalizedMask(equalizedMask(:,sliceNr)>0,sliceNr)=magBins(magBins>0);
            end
        end
        if maxRms==0
            nrVoxPerPhaseBin=nrVoxPerPhaseBin-1;
            if nrVoxPerPhaseBin==0
                disp('could not make equalizedMask');
                keyboard
            end
        end
    end
end
end
