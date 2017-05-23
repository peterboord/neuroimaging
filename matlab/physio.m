clear
dbstop if error
physioDir='/projects2/udall/physio';
subjDirs=dir(fullfile(physioDir,'RC4*'));
subjDirs={subjDirs.name};
for sessionCell=subjDirs
    session=sessionCell{1};
    if exist(fullfile(physioDir,session,'rest_cardio.txt'),'file') && exist(fullfile(physioDir,session,'rest_resp.txt'),'file')
        tic
        disp(session);
        maskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
        funcDir='/projects2/udall/pboord/pic/preproc/pestica';
        tr=2.4;
        physioFs=500;
        nfft=32;
        nyqIdx=nfft/2+1;
        restS=MRIread(fullfile(funcDir,session,'rest_brain_reg.nii'));
        maskPath=fullfile(funcDir,session,maskFile);
        if ~exist(maskPath,'file')
            system(['cd ',fullfile(funcDir,session),...
                ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',maskFile]);
        end
        maskS=MRIread(maskPath);
        nvol=restS.nframes;
        volsize=restS.volsize;
        nslice=volsize(3);
        % physio variables
        physioSamples=physioFs*nvol*tr;
        nTotalSlices=nvol*nslice;
        nPhysioFft=nslice*nfft;
        physioFftIdx=1:nPhysioFft/2:nTotalSlices-nPhysioFft;
        nPhysioEpoch=numel(physioFftIdx);
        physioFftIdx=repmat(physioFftIdx,nPhysioFft,1) + repmat((0:nPhysioFft-1)',1,nPhysioEpoch);
        physioNyqIdx=nPhysioFft/2+1;
        dfPhysio=1/(nPhysioFft*tr/nslice);
        xPhysio=0:dfPhysio:dfPhysio*(physioNyqIdx-1);
        % cardio
        cardio=resample(textread(fullfile(physioDir,session,'rest_cardio.txt')),nTotalSlices,physioSamples); %#ok<*DTXTRD>
        fftCardio=fft(detrend(cardio(reshape(physioFftIdx,nPhysioFft,[]))));
        fftCardio=mean(abs(fftCardio(1:physioNyqIdx,:)),2);
        fftCardioAliased=fft(detrend(cardio(reshape(physioFftIdx(1:nslice:end,:),nfft,[]))));
        fftCardioAliased=mean(abs(fftCardioAliased(1:nyqIdx,:)),2);
        % resp
        resp=resample(textread(fullfile(physioDir,session,'rest_resp.txt')),nTotalSlices,physioSamples);
        fftResp=fft(detrend(resp(reshape(physioFftIdx,nPhysioFft,[]))));
        fftResp=mean(abs(fftResp(1:physioNyqIdx,:)),2);
        fftRespAliased=fft(detrend(resp(reshape(physioFftIdx(1:nslice:end,:),nfft,[]))));
        fftRespAliased=mean(abs(fftRespAliased(1:nyqIdx,:)),2);
        % fmri
        rest=detrend(reshape(restS.vol,prod(restS.volsize),nvol)');
        mask=logical(reshape(maskS.vol,volsize(1)*volsize(2),nslice));
        nvox=size(rest,2);
        fftIdx=1:nfft/2:nvol-nfft;
        nepoch=numel(fftIdx);
        fftIdx=repmat(fftIdx,nfft,1) + repmat((0:nfft-1)',1,nepoch);
        fftRest=fft(reshape(rest(fftIdx,:),nfft,[],nvox));
        fftRest=reshape(fftRest(1:nyqIdx,:,:),nyqIdx,[],volsize(1)*volsize(2),nslice);
        cohFft=zeros(nyqIdx,nslice);
        noncohFft=zeros(nyqIdx,nslice);
        for sliceNr=1:nslice
            cohFft(:,sliceNr)=squeeze(mean(abs(mean(squeeze(fftRest(:,:,mask(:,sliceNr),sliceNr)),3)),2));
            noncohFft(:,sliceNr)=squeeze(mean(abs(mean(squeeze(abs(fftRest(:,:,mask(:,sliceNr),sliceNr))),3)),2));
        end
        df=1/(nfft*tr);
        x=0:df:df*(nyqIdx-1);
        figH=figure;%('Visible','off');
        xmax=1;
        subplot(7,1,1);plot(x,mean(cohFft,2));xlim([0,xmax]);
        subplot(7,1,2);plot(x,mean(noncohFft,2));xlim([0,xmax]);
        subplot(7,1,3);plot(x,mean(cohFft,2)./mean(noncohFft,2));xlim([0,xmax]);
        subplot(7,1,4);plot(xPhysio,fftCardio);xlim([0,xmax]);
        subplot(7,1,5);plot(xPhysio,fftResp);xlim([0,xmax]);
        subplot(7,1,6);plot(x,fftCardioAliased);xlim([0,xmax]);
        subplot(7,1,7);plot(x,fftRespAliased);xlim([0,xmax]);
        print(figH,'-dpng','-r300',fullfile(funcDir,session,'physioFft.png'));
        copyfile(fullfile(funcDir,session,'physioFft.png'),fullfile(physioDir,['physioFft_',session,'.png']),'f');
        clearvars -except physioDir subjDirs sessionCell
        toc
    end
end