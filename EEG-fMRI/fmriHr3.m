function fmriHr3

dbstop if error

nrTr = 120; % tr
isi = 20; % sec
tr = 2.52; % sec
hrLen = round(isi/tr);
nrEventPerSess = floor(nrTr*tr/isi);
nrRestTr = 240;
%fs = 250; % Hz
fs = 256; % Hz

baseDir = '/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/MRI/PART1';
sessDir = {'WIP_ME-VisReaction1_SENSE','WIP_ME-VisReaction2_SENSE','WIP_ME-VisReaction3_SENSE','WIP_ME-VisReaction4_SENSE'};
nrSess = numel(sessDir);
roiList = {'roi1.txt','roi2.txt','roi3.txt','roi4.txt','roi5.txt','roi6.txt','roi7.txt'};
nrRoi = numel(roiList);
% trSessRoi = zeros(nrTr,nrSess,nrRoi);
% trEventRoi = zeros(hrLen,nrSess*nrEventPerSess,nrRoi);
% sessTime = 0:tr:tr*(nrTr-1);
% periTime = 0:tr:tr*(hrLen-1);
% eventTime = 0:isi:isi*(nrEventPerSess-1);
% periEventTime = repmat(eventTime,[hrLen,1]) + repmat(periTime(:),[1,nrEventPerSess]);
% for sessNr = 1:nrSess
%     eventIdxs = (sessNr-1)*nrEventPerSess + (1:nrEventPerSess);
%     for roiNr = 1:nrRoi
%          roiTr = load(fullfile(baseDir,sessDir{sessNr},roiList{roiNr}),'-ascii');
%          trEventRoi(:,eventIdxs,roiNr) = interp1(sessTime,roiTr,periEventTime,'spline');
%          trSessRoi(:,sessNr,roiNr) = roiTr;
%     end
% end
% trRoi = squeeze(mean(trEventRoi,2));
% win = repmat(hann(hrLen),[1,nrSess*nrEventPerSess,nrRoi]);
% winTrEventRoi = win.*trEventRoi;
% fftTrEventRoi = fft(trEventRoi);
% absFftTrEventRoi = abs(fftTrEventRoi);
% angleFftTrEventRoi = angle(fftTrEventRoi);
% % figure('WindowStyle','docked'),plot(trSessRoi(:,:,6));legend('show')
% % figure('WindowStyle','docked'),plot(trEventRoi(:,:,6))
% % figure('WindowStyle','docked'),plot(trRoi);legend('show')
% % figure('WindowStyle','docked'),plot(mean(absFftTrEventRoi(:,:,6),2))
% % figure('WindowStyle','docked'),plot(angleFftTrEventRoi(:,:,6))
% % figure('WindowStyle','docked'),for i=1:5,subplot(2,3,i);rose(angleFftTrEventRoi(i,:,6));end
restTrRoi = zeros(nrRestTr,nrRoi);
for roiNr = 1:nrRoi
    restTrRoi(:,roiNr) = load(fullfile(baseDir,'WIP_ME-RS_SENSE',roiList{roiNr}),'-ascii');
end
% % figure('WindowStyle','docked'),plot(restTrRoi(:,6))
% % figure,subplot(2,1,1);hist(posRoi{6}(:,2));subplot(2,1,2);hist(negRoi{6}(:,2))
% clusterS = MRIread('/NAS_II/Home/pboord/Documents/Scripts/cluster_mask_zstat1_1mm.nii.gz');
% %cluster2d = reshape(clusterS.vol,[1,prod(clusterS.volsize)]);
% loretaS = MRIread('/NAS_II/Home/pboord/Documents/Scripts/loretaMni.nii.gz');
% %loreta2d = reshape(loretaS.vol,[1,prod(loretaS.volsize)]);
% loretaVox = cell(nrRoi,1);
% [~,~,xls] = xlsread('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/EEG_FMRI_002/EEG/export/visualROIfile.xls');
% for roiNr = 1:nrRoi
%     loretaVox{roiNr} = unique(loretaS.vol(clusterS.vol == roiNr));
%     loretaVox{roiNr}(1) = [];
%     for voxNr = 1:numel(loretaVox{roiNr})
%         xls{3+loretaVox{roiNr}(voxNr),7} = roiNr;
%     end
%     roi2save = loretaVox{roiNr};
%     save(fullfile('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/EEG_FMRI_002/EEG/export',['roi',num2str(roiNr)]),'roi2save','-ascii');
% end
% fid = fopen('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/EEG_FMRI_002/EEG/export/visualROIfileOut.csv','w');
% fprintf(fid,'%s\n',xls{1,1});
% fprintf(fid,'%s\n',xls{2,1});
% fprintf(fid,'%s;%s;%s;%s;%s;%s;%s\n',xls{3,:});
% for rowNr = 1:2394
%     fprintf(fid,'%d;%d;%d;%s;%s;%s;%d\n',xls{3+rowNr,:});
% end
% fclose(fid);
%eegRoi = load('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/EEG_FMRI_002/EEG/export/PRB_131107_Loreta_7roi_rest.csv');
%eegRoi = load('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/export/PRB_131107_Loreta_7roi_rest.csv');
eegRoi = load('/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/EEG/part1/export/PRB_131107_Loreta_7roi_rest.csv');
delaySeconds = 0:1:8; % sec
winLen = round(fs*tr);
absSpecRatio = zeros(winLen,numel(delaySeconds));
binRange = 2:126;
freq = (binRange/tr)';
posRoi = restTrRoi>=0;
negRoi = restTrRoi<0;
win = hamming(winLen,'periodic');
for ch = 1:7
    for delayIdx = [1,5]%1:numel(delaySeconds)
        hrDelaySec = delaySeconds(delayIdx); % sec
        hrDelayTr = ceil(hrDelaySec/tr);
        posRoi = restTrRoi>=0;
        posRoi(1:hrDelayTr,:) = false;
        negRoi = restTrRoi<0;
        negRoi(1:hrDelayTr,:) = false;
        delayedEegRoi = [zeros(floor(hrDelaySec*fs),nrRoi);eegRoi(1:end-floor(hrDelaySec*fs),:)];
        posFft = getEpochFft(delayedEegRoi,winLen,nrRestTr,nrRoi,posRoi,win,ch);
        negFft = getEpochFft(delayedEegRoi,winLen,nrRestTr,nrRoi,negRoi,win,ch);
        absSpecRatio(:,delayIdx) = mean(abs(posFft),2)./mean(abs(negFft),2);
        figure('WindowStyle','docked'),plot(freq,mean(abs(posFft(binRange,:)),2),freq,mean(abs(negFft(binRange,:)),2))
        %figure('WindowStyle','docked'),plot(freq,log10(absSpecRatio(binRange,delayIdx)));
    end
    %figure('WindowStyle','docked'),surf(log10(absSpecRatio(binRange,:)),'XData',1:9,'YData',freq)
    %xdata = [1,5];figure('WindowStyle','docked'),plot(freq,log10(absSpecRatio(binRange,xdata)));
end
end

function roiFft = getEpochFft(eegRoi,winLen,nrRestTr,nrRoi,signRoi,win,ch)
signStartSamp = cell(nrRoi,1);
for roiNr = ch%1:nrRoi
    signStartSamp{roiNr} = round(1:winLen:nrRestTr*winLen)'.*signRoi(:,roiNr);
    signStartSamp{roiNr} = unique(signStartSamp{roiNr});
    signStartSamp{roiNr}(1) = [];
    nrEpoch = numel(signStartSamp{roiNr});
    roiFft = complex(zeros(winLen,nrEpoch),zeros(winLen,nrEpoch));
    for epochNr = 1:nrEpoch
        startEegSamp = signStartSamp{roiNr}(epochNr);
        endEegSamp = startEegSamp + winLen - 1;
        if endEegSamp > size(eegRoi,1)
            roiFft(:,epochNr:end) = [];
            break
        else
            epochEeg = eegRoi(startEegSamp:endEegSamp,roiNr);
            roiFft(:,epochNr) = fft(epochEeg).*win;
        end
    end
end
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