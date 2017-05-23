function physioSge(session)
dbstop if error
if nargin==0
    %session='RC4107-2';
    session='RC4103-1';
end
physioDir='/projects2/udall/physio';
if exist(fullfile(physioDir,session,'rest_cardio.txt'),'file') && exist(fullfile(physioDir,session,'rest_resp.txt'),'file')
    tic
    disp(session);
    maskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
    funcDir='/projects2/udall/pboord/pic/preproc/pestica';
    tr=2.4;
    physioFs=100; % actually 5 replications of 100 Hz = 500 Hz, but replications removed
    
    nfft=30;
    
    nyqIdx=nfft/2+1;
    %restS=MRIread(fullfile(funcDir,session,'rest_brain_reg.nii'));
    restS=MRIread(fullfile(funcDir,session,'RC4103-1_rest.feat','filtered_func_data.nii'));
    % slomoco
    %restS=MRIread(fullfile(funcDir,session,'RC4103-1_rest.slicemocoxy_afni.zalg_moco2.nii'));
    %sloTxt=load(fullfile(funcDir,session,[session,'_rest_pestica'],'slomoco.combined.reg6dof.txt'),'-ascii');
    
    maskPath=fullfile(funcDir,session,maskFile);
%     if ~exist(maskPath,'file')
%         system(['cd ',fullfile(funcDir,session),...
%             ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',maskFile]);
%     end
    maskS=MRIread(maskPath);
    nvol=restS.nframes;
    volsize=restS.volsize;
    nvox=prod(volsize(1:2));
    nslice=volsize(3);
    % physio variables
    physioSamples=physioFs*nvol*tr;
    nTotalSlices=nvol*nslice;
    nPhysioFft=nslice*nfft;
    physioFftIdx=1:nPhysioFft:nTotalSlices-nPhysioFft;
    nPhysioEpoch=numel(physioFftIdx);
    physioFftIdx=repmat(physioFftIdx,nPhysioFft,1) + repmat((0:nPhysioFft-1)',1,nPhysioEpoch);
    physioNyqIdx=nPhysioFft/2+1;
    dfPhysio=1/(nPhysioFft*tr/nslice);
    xPhysio=0:dfPhysio:dfPhysio*(physioNyqIdx-1);
    % cardio
    cardioRaw=reshape(textread(fullfile(physioDir,session,'rest_cardio.txt')),5,[]);
    cardioRaw=cardioRaw(1,:);
    cardio=resample(cardioRaw,nTotalSlices,physioSamples); %#ok<*DTXTRD>
    fftCardio=fft(detrend(cardio(reshape(physioFftIdx,nPhysioFft,[]))));
    fftCardio=mean(abs(fftCardio(1:physioNyqIdx,:)),2);
    fftCardioAliased=fft(detrend(cardio(reshape(physioFftIdx(1:nslice:end,:),nfft,[]))));
    fftCardioAliased=mean(abs(fftCardioAliased(1:nyqIdx,:)),2);
    % resp
    respRaw=reshape(textread(fullfile(physioDir,session,'rest_resp.txt')),5,[]);
    respRaw=respRaw(1,:);
    resp=resample(respRaw,nTotalSlices,physioSamples);
    fftResp=fft(detrend(resp(reshape(physioFftIdx,nPhysioFft,[]))));
    fftResp=mean(abs(fftResp(1:physioNyqIdx,:)),2);
    fftRespAliased=fft(detrend(resp(reshape(physioFftIdx(1:nslice:end,:),nfft,[]))));
    fftRespAliased=mean(abs(fftRespAliased(1:nyqIdx,:)),2);
    % fmri
    rest=reshape(restS.vol,prod(volsize),nvol);
    rest=detrend(rest')';
    rest=reshape(rest,volsize(1)*volsize(2),nslice,nvol);
    mask=logical(reshape(maskS.vol,volsize(1)*volsize(2),nslice));
    fftIdx=1:nfft:nvol;
    nepoch=numel(fftIdx);
    fftIdx=repmat(fftIdx,nfft,1) + repmat((0:nfft-1)',1,nepoch);
    rest4d=reshape(reshape(rest(:,:,fftIdx),nvox*nslice,[])',nfft,[],volsize(1)*volsize(2),nslice);
    disp(size(rest4d));
    fftRest=fft(rest4d);
    fftRest=reshape(fftRest(1:nyqIdx,:,:),nyqIdx,[],volsize(1)*volsize(2),nslice);
    fftRest=fftRest(1:nyqIdx,:,:,:);
    disp(size(fftRest));
    cohFft=zeros(nyqIdx,nslice);
    cohFft2=zeros(nyqIdx,nslice);
    noncohFft=zeros(nyqIdx,nslice);
    voxInSlice=cell(nslice,1);
    % retroicor
    oba=load('oba.slibase.1D_matlab','-ascii');
    nrReg=8;
    reg=zeros(nvol,nslice,nrReg);
    for regNr=1:nrReg
        reg(:,:,regNr)=oba(:,5+regNr:13:end);
        regX=reg(:,1,regNr);
        save(['reg',num2str(regNr)],'regX','-ascii');
    end
    maxSliceCor=zeros(nslice,1);
    maxSliceLag=zeros(nslice,1);
    for sliceNr=1:nslice
        cohFft(:,sliceNr)=squeeze(mean(abs(mean(squeeze(fftRest(:,:,mask(:,sliceNr),sliceNr)),3)),2));
        restSlice=squeeze(rest4d(:,:,mask(:,sliceNr),sliceNr));
        cohRestSlice=squeeze(rest(mask(:,sliceNr),sliceNr,:));
        cohRestSlice=detrend(cohRestSlice')';
        nrVox=size(cohRestSlice,1);
        regSlice=squeeze(reg(:,sliceNr,:));
        regCor=zeros(nrVox,1);
        regLag=zeros(nrVox,1);
        
        regNr = 1;
        
        for voxNr=1:nrVox
           [r,lag]=xcorr(regSlice(:,regNr),cohRestSlice(voxNr,:)',10,'coeff');
           [maxR,maxIdx]=max(r);
           regCor(voxNr)=maxR;
           regLag(voxNr)=lag(maxIdx);
        end
        [maxR,maxIdx]=max(regCor);
        maxSliceCor(sliceNr)=maxR;
        maxSliceLag(sliceNr)=regLag(maxIdx);
        %cohRestSlice=bsxfun(@rdivide,cohRestSlice,std(cohRestSlice,1,2));
        %voxInSlice(sliceNr,:)=mean(cohRestSlice,1);
        voxInSlice{sliceNr}=cohRestSlice;
        
%         figure
        df=1/(nvol*tr);
        fLo=1+round(0.01/df);
        %fLo=0;
        restAliasFft=fft(mean(cohRestSlice,1)');
        restAliasFftFilt=restAliasFft;
        restAliasFftFilt(1:fLo)=0;
        % reflect about Nyquist and ifft
        respUnaliasFft=[0;restAliasFft(2:fLo);zeros(nvol/2-fLo,1);conj(flipud(restAliasFftFilt(2:nvol/2+1)));0;restAliasFftFilt(2:nvol/2+1);zeros(nvol/2-fLo,1);conj(flipud(restAliasFft(2:fLo)))];
        respUnaliasComplex=ifft(respUnaliasFft);
        respUnalias=real(respUnaliasComplex);
%         x1=1:1290; x2=1:60;
%         %x1=1:12900; x2=1:600;
%         biplot(x1,resp(x1),x2,respUnalias(x2))
%         figure
%         plot(xcorr(resample(resp',600,12900),respUnalias,'coeff'))
%         
%         figure
%         x1=2:300;
%         y2=abs(fft(resp));
%         biplot(x1,abs(respUnaliasFft(x1)),x1,y2(x1));
%         
%         figure
%         df=1/(nvol*tr);
%         f=0:df:(nvol-1)*df;
%         plot(f(2:end/2),abs(restAliasFft(2:end/2)))
        
        
        restSliceMinusCoh=restSlice-repmat(mean(restSlice,3),1,1,size(restSlice,3));
        cohFft21d=mean(abs(fft(mean(restSliceMinusCoh,3))),2);
        cohFft2(:,sliceNr)=cohFft21d(1:nyqIdx);
        noncohFft(:,sliceNr)=squeeze(mean(mean(squeeze(abs(fftRest(:,:,mask(:,sliceNr),sliceNr))),3),2));
    end
    %figure,plot(maxSliceCor)
    %figure,plot(maxSliceLag)
    x1=1000:4000;%nTotalSlices-1000;%6000;
    %figure
    sqssSlice=zeros(nslice,nvol);
    for sliceNr=1:nslice
        sqssSlice(sliceNr,:)=sqrt(mean(voxInSlice{sliceNr}.^2,1));
    end
    windowSize=16;
    sqssSlice=detrend(sqssSlice')';
    sqssSliceSmooth=filtfilt((1/windowSize)*ones(1,windowSize),1,sqssSlice(:));
    sqssSliceSmooth=sqssSlice;
    sqssSliceFft=zeros(nslice,nvol);
    nrTr=2;
    dfSlice=1/(nrTr*tr);
    for volNr=1:nrTr:nrTr*10
        figure('WindowStyle','docked');
        nrPlot=5;
        x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
        subplot(nrPlot,1,1);
        sqssWave=reshape(sqssSlice(x),[],1);
        padFactor=6;
        sqssSliceFft=fft(sqssWave,padFactor*numel(sqssWave));
        %sqssSliceFft=fft(sqssWave.*hann(numel(sqssWave),'periodic'),padFactor*numel(sqssWave));
        plot(sqssWave);
        subplot(nrPlot,1,2);
        fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
        minCardioHz=0.6;
        maxCardioHz=1.5;
        minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
        maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
        [maxSpec,maxIdx]=max(abs(sqssSliceFft(minCardio:maxCardio)));
        maxSpecIdx=minCardio+maxIdx-1;
        plot(fSlice(1:floor(end/2)),abs(sqssSliceFft(1:floor(end/2))),'b',fSlice(maxSpecIdx),maxSpec,'r*');
        nrX=numel(x);
        cardioEst=abs(sqssSliceFft(maxSpecIdx))*cos(2*pi*fSlice(maxSpecIdx)*(0:nrX-1)*tr/nslice+angle(sqssSliceFft(maxSpecIdx)));
        subplot(nrPlot,1,3);
        plot(resp(x));
        subplot(nrPlot,1,4);
        [pks,locs]=findpeaks(cardio(x),'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
        plot(cardio(x));
        xCardioPlot=xlim;
        hold on
        plot(locs,pks,'r*');
        plot(1:nrX,cardioEst,'g');
        hold off
        subplot(nrPlot,1,5);
        cardioWidth=round(0.2/(dfSlice/padFactor))+1;
        cardioIdxs=(maxSpecIdx-cardioWidth:maxSpecIdx+cardioWidth)';
        sizePadFft=size(sqssSliceFft,1);
        cardioSliceFft=complex(zeros(sizePadFft,1),zeros(sizePadFft,1));
        cardioSliceFft(cardioIdxs)=sqssSliceFft(cardioIdxs);
        imCardioIdxs=sizePadFft-cardioIdxs+2;
        cardioSliceFft(imCardioIdxs)=sqssSliceFft(imCardioIdxs);
        ifftCardioSliceFft=ifft(cardioSliceFft);
        ifftCardioSliceFft=ifftCardioSliceFft(1:numel(sqssWave));
        plot(1:nrX,ifftCardioSliceFft,'m');
        %plot(locs(2:end),(diff(locs*tr)/nslice).^-1);
        xlim(xCardioPlot);
    end
   keyboard
    sqssSliceSmooth=sqssSliceSmooth(:);
    figure,biplot(x1,resp(x1),x1,sqssSliceSmooth(x1))
    sqssSliceSmoothFft=fft(hann(nTotalSlices,'periodic').*sqssSliceSmooth);
    df=1/(nvol*tr);
    f=0:df:df*(nTotalSlices-1);
    figure,plot(f(2:end/2),abs(sqssSliceSmoothFft(2:end/2)));
    % tr harmonics at nvol+1 & 2*nvol+1
    deleteRange=-1:1;
    sqssSliceSmoothFftFilt=sqssSliceSmoothFft;
%     sqssSliceSmoothFftFilt=sqssSliceSmoothFft;
%     for harmonicNr=1:floor(nslice/2)
%         sqssSliceSmoothFftFilt(harmonicNr*nvol+1+deleteRange)=0;
%         sqssSliceSmoothFftFilt(nTotalSlices-harmonicNr*nvol+1+deleteRange)=0;
%     end
    figure,plot(f(2:end/2),abs(sqssSliceSmoothFftFilt(2:end/2)));
    %sqssSliceClean=ifft(sqssSliceSmoothFftFilt)./hann(nTotalSlices,'periodic');
    sqssSliceClean=ifft(sqssSliceSmoothFftFilt);
    %sqssSliceClean=filtfilt((1/windowSize)*ones(1,windowSize),1,sqssSliceClean(:));
%     lpfCutoffIdx=round(0.4/df)-1;
%     sqssSliceSmoothFftFilt=sqssSliceSmoothFft.*[0;ones(lpfCutoffIdx,1);zeros(nTotalSlices/2-lpfCutoffIdx-1,1);0;zeros(nTotalSlices/2-lpfCutoffIdx-1,1);ones(lpfCutoffIdx,1)];
%     sqssSliceClean=ifft(sqssSliceSmoothFftFilt,'symmetric');
    figure,biplot(x1,resp(x1),x1,sqssSliceClean(x1))
    
    figure,plot(xcorr(resp',sqssSliceSmooth,'coeff'))
    df=1/(nfft*tr);
    x=0:df:df*(nyqIdx-1);
    figH=figure;%('Visible','off');
    xmax=1;
    subplot(8,1,1);plot(x,mean(cohFft,2));xlim([0,xmax]);
    subplot(8,1,2);plot(x,mean(noncohFft,2));xlim([0,xmax]);
    subplot(8,1,3);plot(x,mean(cohFft,2)./mean(noncohFft,2));xlim([0,xmax]);
    subplot(8,1,4);plot(xPhysio,fftCardio);xlim([0,xmax]);
    subplot(8,1,5);plot(xPhysio,fftResp);xlim([0,xmax]);
    subplot(8,1,6);plot(x,fftCardioAliased);xlim([0,xmax]);
    subplot(8,1,7);plot(x,fftRespAliased);xlim([0,xmax]);
    subplot(8,1,8);plot(x,mean(cohFft2,2));xlim([0,xmax]);
    %print(figH,'-dpng','-r300',fullfile(funcDir,session,'physioFft.png'));
    copyfile(fullfile(funcDir,session,'physioFft.png'),fullfile(physioDir,['physioFft_',session,'.png']),'f');
end