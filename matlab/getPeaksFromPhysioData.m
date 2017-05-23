%%%%%%%%%%%%%%%%%%%%%%%
% getPeaksFromPhysioData
%%%%%%%%%%%%%%%%%%%%%%%
function [cardio,allCardioPeakTimes,resp,allRespPeakTimes]=getPeaksFromPhysioData(cardioFilePath,respFilePath,tr,nvol,nslice)
% physio
physioFs=100; % actually 5 replications of 100 Hz = 500 Hz, but replications removed
physioSamples=physioFs*nvol*tr;
% cardio
cardio=[];
if exist(cardioFilePath,'file')
%     cardioRaw=reshape(textread(cardioFilePath),5,[]);
    cardioRaw=reshape(textread(cardioFilePath),1,[]); % assumes file already has replications removed
    cardioRaw=detrend(cardioRaw(1,:)');
    cardioRaw=cardioRaw/std(cardioRaw);
    cardio=resample(cardioRaw,nvol*nslice,physioSamples); %#ok<*DTXTRD>
    % resp
    %respRaw=reshape(textread(respFilePath),5,[]);
    respRaw=reshape(textread(respFilePath),1,[]); % assumes file already has replications removed
    respRaw=respRaw(1,:);
    resp=resample(respRaw,nvol*nslice,physioSamples);
end
maxCardioHz=1.5;
[~,locs]=findpeaks(cardioRaw,'MinPeakHeight',0,...
    'MinPeakDistance',round(0.5*physioFs/maxCardioHz),...
    'MinPeakProminence',1);
allCardioPeakTimes=(locs'-1)/physioFs;
% [~,locs]=findpeaks(cardio,'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
% allCardioPeakTimes_ds=(locs'-1)*tr/nslice;
% figure
% tMin=0;
% tMax=tr*10;
% t=tMin:1/physioFs:tMax;
% subplot(2,1,1)
% x=round(t*physioFs)+1;
% cardioPeakTimes=allCardioPeakTimes(allCardioPeakTimes>=tMin & allCardioPeakTimes<=tMax);
% plot(t,cardioRaw(x),'b',cardioPeakTimes,cardioRaw(round(cardioPeakTimes*physioFs)+1),'r*');
[~,locs]=findpeaks(resp,'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
allRespPeakTimes=(locs'-1)*tr/nslice;
% figure
% plot(cardio);
% hold on
% plot(locs,pks,'r*');
% hold off
end
