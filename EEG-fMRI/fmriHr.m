function fmriHr

dbstop if error

nrTr = 120; % tr
isi = 20; % sec
tr = 2.52; % sec
hrLen = round(isi/tr);
nrEventPerSess = floor(nrTr*tr/isi);
nrRestTr = 240;
fs = 250; % Hz

baseDir = '/NAS_II/Home/pboord/Documents/EEG_FMRI_002/EEG_FMRI_002/PART1';
sessDir = {'WIP_ME-VisReaction1_SENSE','WIP_ME-VisReaction2_SENSE','WIP_ME-VisReaction3_SENSE','WIP_ME-VisReaction4_SENSE'};
nrSess = numel(sessDir);
roiList = {'roi1.txt','roi2.txt','roi3.txt','roi4.txt','roi5.txt','roi6.txt','roi7.txt'};
nrRoi = numel(roiList);
trSessRoi = zeros(nrTr,nrSess,nrRoi);
trEventRoi = zeros(hrLen,nrSess*nrEventPerSess,nrRoi);
sessTime = 0:tr:tr*(nrTr-1);
periTime = 0:tr:tr*(hrLen-1);
eventTime = 0:isi:isi*(nrEventPerSess-1);
periEventTime = repmat(eventTime,[hrLen,1]) + repmat(periTime(:),[1,nrEventPerSess]);
for sessNr = 1:nrSess
    eventIdxs = (sessNr-1)*nrEventPerSess + (1:nrEventPerSess);
    for roiNr = 1:nrRoi
         roiTr = load(fullfile(baseDir,sessDir{sessNr},roiList{roiNr}),'-ascii');
         trEventRoi(:,eventIdxs,roiNr) = interp1(sessTime,roiTr,periEventTime,'spline');
         trSessRoi(:,sessNr,roiNr) = roiTr;
    end
end
trRoi = squeeze(mean(trEventRoi,2));
win = repmat(hann(hrLen),[1,nrSess*nrEventPerSess,nrRoi]);
winTrEventRoi = win.*trEventRoi;
fftTrEventRoi = fft(trEventRoi);
absFftTrEventRoi = abs(fftTrEventRoi);
angleFftTrEventRoi = angle(fftTrEventRoi);
% figure('WindowStyle','docked'),plot(trSessRoi(:,:,6));legend('show')
% figure('WindowStyle','docked'),plot(trEventRoi(:,:,6))
% figure('WindowStyle','docked'),plot(trRoi);legend('show')
% figure('WindowStyle','docked'),plot(mean(absFftTrEventRoi(:,:,6),2))
% figure('WindowStyle','docked'),plot(angleFftTrEventRoi(:,:,6))
% figure('WindowStyle','docked'),for i=1:5,subplot(2,3,i);rose(angleFftTrEventRoi(i,:,6));end
restTrRoi = zeros(nrRestTr,nrRoi);
for roiNr = 1:nrRoi
    restTrRoi(:,roiNr) = load(fullfile(baseDir,'WIP_ME-RS_SENSE',roiList{roiNr}),'-ascii');
end
% figure('WindowStyle','docked'),plot(restTrRoi(:,6))
posRoi = cell(nrRoi,1);
negRoi = cell(nrRoi,1);
for roiNr = 1:nrRoi
    posSeq = zeros(nrRestTr,2); % [,1] = start tr; [,2] = length
    negSeq = zeros(nrRestTr,2);
    posIdx = 1;
    negIdx = 1;
    for trNr = 1:nrRestTr
        if restTrRoi(trNr,roiNr) >= 0
            if posSeq(posIdx,2) == 0
                posSeq(posIdx,1) = trNr;
            end
            posSeq(posIdx,2) = posSeq(posIdx,2) + 1;
            if negSeq(negIdx,2) ~= 0
                negIdx = negIdx + 1;
            end
        elseif restTrRoi(trNr,roiNr) < 0
            if posSeq(posIdx,2) ~= 0
                posIdx = posIdx + 1;
            end
            if negSeq(negIdx,2) == 0
                negSeq(negIdx,1) = trNr;
            end
            negSeq(negIdx,2) = negSeq(negIdx,2) + 1;
        end
    end
    if posSeq(posIdx,1) == 0
        posIdx = posIdx - 1;
    end
    if negSeq(negIdx,1) == 0
        negIdx = negIdx - 1;
    end
    posRoi{roiNr} = posSeq(1:posIdx,:);
    negRoi{roiNr} = negSeq(1:negIdx,:);
end
% figure,subplot(2,1,1);hist(posRoi{6}(:,2));subplot(2,1,2);hist(negRoi{6}(:,2))
clusterS = MRIread('/NAS_II/Home/pboord/Documents/Scripts/cluster_mask_zstat1_1mm.nii.gz');
%cluster2d = reshape(clusterS.vol,[1,prod(clusterS.volsize)]);
loretaS = MRIread('/NAS_II/Home/pboord/Documents/Scripts/loretaMni.nii.gz');
%loreta2d = reshape(loretaS.vol,[1,prod(loretaS.volsize)]);
loretaVox = cell(nrRoi,1);
[~,~,xls] = xlsread('/NAS_II/Home/pboord/Documents/EEG_FMRI_002/EEG/export/visualROIfile.xls');
for roiNr = 1:nrRoi
    loretaVox{roiNr} = unique(loretaS.vol(clusterS.vol == roiNr));
    loretaVox{roiNr}(1) = [];
    for voxNr = 1:numel(loretaVox{roiNr})
        xls{3+loretaVox{roiNr}(voxNr),7} = roiNr;
    end
    roi2save = loretaVox{roiNr};
    save(fullfile('/NAS_II/Home/pboord/Documents/EEG_FMRI_002/EEG/export',['roi',num2str(roiNr)]),'roi2save','-ascii');
end
fid = fopen('/NAS_II/Home/pboord/Documents/EEG_FMRI_002/EEG/export/visualROIfileOut.csv','w');
fprintf(fid,'%s\n',xls{1,1});
fprintf(fid,'%s\n',xls{2,1});
fprintf(fid,'%s;%s;%s;%s;%s;%s;%s\n',xls{3,:});
for rowNr = 1:2394
    fprintf(fid,'%d;%d;%d;%s;%s;%s;%d\n',xls{3+rowNr,:});
end
fclose(fid);
eegRoi = load('/NAS_II/Home/pboord/Documents/EEG_FMRI_002/EEG/export/PRB_131107_Loreta_7roi_rest.csv');
hrDelay = 0;
[posEpochPsd,f] = getEpochPsd(eegRoi,fs,nrRoi,posRoi,hrDelay);
negEpochPsd = getEpochPsd(eegRoi,fs,nrRoi,negRoi,hrDelay);
figure('WindowStyle','docked'),plot(f(3:end),mean(posEpochPsd(3:end,:),2),f(3:end),mean(negEpochPsd(3:end,:),2))
end
function [epochPsd,f] = getEpochPsd(eegRoi,fs,nrRoi,signRoi,hrDelay)
for roiNr = 6%1:nrRoi
    epochPsd = zeros(fs/2+1,sum(signRoi{roiNr}(:,2)));
    epochIdx = 0;
    for epochNr = 1:size(signRoi{roiNr},1)
        for sampNr = 1:signRoi{roiNr}(epochNr,2)
            epochIdx = epochIdx + 1;
            epochStartSamp = (signRoi{roiNr}(epochNr,1)-1)*fs+sampNr;
            epochEndSamp = epochStartSamp + fs - 1;
            epochEeg = eegRoi(epochStartSamp:epochEndSamp,roiNr);
            [epochPsd(:,epochIdx),f] = pwelch(epochEeg,fs,[],fs,fs);
        end
    end
end
end