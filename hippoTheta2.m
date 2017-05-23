function hippoTheta2

subjPath = '/var/local/scratch/pboord/INDI/ABIDE/SocialBrainLab_descending/SBL/0051557/session_1/rest_1/rest.feat';
%lHippoXYZ = [49,29,18;49,32,17;49,34,16;49,35,15;49,39,14;49,41,13;49,42,12];
lHippoXYZ = [49,29,17;49,31,16;49,33,15;49,35,14;49,37,13;49,39,12];
funcFile = load_untouch_nii(fullfile(subjPath,'filtered_func_data.nii'));
funcImg = double(funcFile.img);
nrVols = size(funcImg,4);
nrHippoSamples = size(lHippoXYZ,1);
hippoActivity = zeros(nrHippoSamples,nrVols);
for volNr = 1:nrVols
    for idx = 1:nrHippoSamples
        hippoActivity(idx,volNr) = funcImg(lHippoXYZ(idx,1)+1,lHippoXYZ(idx,2),lHippoXYZ(idx,3),volNr);
    end
end
hippoActivity = hippoActivity - repmat(mean(hippoActivity,2),1,nrVols);
hippoActivity = hippoActivity./repmat(std(hippoActivity,0,2),1,nrVols);
surf(hippoActivity);view(0,90)
end