function hippoTheta

subjPath = '/var/local/scratch/pboord/INDI/ABIDE/SocialBrainLab_descending/SBL/0051557/session_1/rest_1/rest.feat';
lHippoXYZ = [49,29,18;49,30,18;49,31,18;49,32,17;49,33,17;49,34,16;49,35,15;49,36,15;49,37,15;49,38,15;49,39,14;49,40,14;49,41,13;49,42,12];
funcFile = load_untouch_nii(fullfile(subjPath,'filtered_func_data.nii'));
roiFile = load_untouch_nii(fullfile(subjPath,'ROI_MNI_V4_4mm_funcReg.nii'));
funcImg = double(funcFile.img);
nrVols = size(funcImg,4);
roiImg = roiFile.img;
hippoImg = double(roiImg == 4101);
hippoZ = sum(squeeze(sum(hippoImg,1)),1);
hippoSliceFirst = find(hippoZ,1,'first');
hippoSliceLast = find(hippoZ,1,'last');
nrHippoSlices = hippoSliceLast-hippoSliceFirst+1;
hippoActivity = zeros(nrHippoSlices,nrVols);
for volNr = 1:nrVols
    for zIdx = 1:nrHippoSlices
        z = hippoSliceFirst+zIdx-1;
        funcSlice = funcImg(:,:,z,volNr);
        hippoSlice = hippoImg(:,:,z);
        [i,j] = ind2sub(size(hippoSlice),find(hippoSlice(:)));
        hippoActivity(zIdx,volNr) = sum(funcSlice(:).*hippoSlice(:))/sum(hippoSlice(:));
        %hippoActivity(zIdx,volNr) = funcSlice(i(end),j(end));
    end
end
hippoActivity = hippoActivity - repmat(mean(hippoActivity,2),1,nrVols);
hippoActivity = hippoActivity./repmat(std(hippoActivity,0,2),1,nrVols);
surf(hippoActivity);
view(0,90);
end