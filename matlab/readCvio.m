function readCvio
subjects={ 3211  3402  3407  3422  3424  3484  3485 3486  3487  3488  3492  3496  3497  3503  3504  3505};
% no file: 3425
% other dir: 3498
% short file: 3414
Hd = mainsFilter;
for subjNr=1:numel(subjects)
    subject=num2str(subjects{subjNr});
    dataDir=fullfile('/NAS_II/Projects/mprojects/wfe_subject_data',subject,subject,'resting_state');
    fidCvioTtl=fopen(fullfile(dataDir,'acquire_ttl_1.cvio.tim'));
    % EKG
    fidCvioTim=fopen(fullfile(dataDir,'acquire_EKG_1.cvio.tim'));
    fidCvioDat=fopen(fullfile(dataDir,'acquire_EKG_1.cvio.dat'));
    cvioTtl=fread(fidCvioTtl,Inf,'uint32',0,'b');
    cvioTim=fread(fidCvioTim,Inf,'uint32',0,'b');
    cvioDat=fread(fidCvioDat,Inf,'uint32',0,'b');
    cvioDat=filter(Hd,cvioDat);
    fclose(fidCvioTim);
    fclose(fidCvioDat);
    figure('WindowStyle','docked');
    plot(cvioDat(16000:16000*11))
    % respiration
    fidCvioTim=fopen(fullfile(dataDir,'acquire_respiration_1.cvio.tim'));
    fidCvioDat=fopen(fullfile(dataDir,'acquire_respiration_1.cvio.dat'));
    cvioTtl=fread(fidCvioTtl,Inf,'uint32',0,'b');
    cvioTim=fread(fidCvioTim,Inf,'uint32',0,'b');
    cvioDat=fread(fidCvioDat,Inf,'uint32',0,'b');
    cvioDat=filter(Hd,cvioDat);
    fclose(fidCvioTim);
    fclose(fidCvioDat);
    figure('WindowStyle','docked');
    plot(cvioDat(16000:16000*11))
    fclose(fidCvioTtl);
end
end