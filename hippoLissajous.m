function hippoLissajous

dbstop if error

subjPath = '/var/local/scratch/pboord/INDI/ABIDE/SocialBrainLab_descending/SBL/0051557/session_1/rest_1/rest.feat';
funcFile = load_untouch_nii(fullfile(subjPath,'filtered_func_data.nii'));
roiFile = load_untouch_nii(fullfile(subjPath,'ROI_MNI_V4_4mm_funcReg.nii'));
funcImg = double(funcFile.img);
roiImg = roiFile.img;
hippoImg = double(roiImg == 4101);
hippoIndices = find(hippoImg);
idx = 200;
[i1,j1,k1] = ind2sub(size(hippoImg),hippoIndices(idx));
for idx2 = 1:25%numel(hippoIndices);
    subplot(5,5,idx2);
    [i2,j2,k2] = ind2sub(size(hippoImg),hippoIndices(idx2));
    vox1 = detrend(squeeze(funcImg(i1,j1,k1,:)));
    vox2 = detrend(squeeze(funcImg(i2,j2,k2,:)));
%     Fs = 1/2.185;
%     h=fdesign.lowpass('N,F3dB',12,0.01,Fs);
%     d1 = design(h,'butter');
%     vox1 = filtfilt(d1.sosMatrix,d1.ScaleValues,vox1);
%     vox2 = filtfilt(d1.sosMatrix,d1.ScaleValues,vox2);
    sgOrder = 5;
    vox1 = vox1 - sgolayfilt(vox1,sgOrder,41);
    vox2 = vox2 - sgolayfilt(vox2,sgOrder,41);
    vox1 = vox1/std(vox1);
    vox2 = vox2/std(vox2);
    
    plot([vox1,vox2]);
    plot(vox1,vox2,'.');
%     xlim([-1000,1000]);
%     ylim([-1000,1000]);
    %xlim([-5,5]);
    ylim([-5,5]);
end
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
figure,surf(hippoActivity);view(0,90)
end