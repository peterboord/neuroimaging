function makeCardioRegressors(session)
dbstop if error
if nargin==0
    %session='RC4107-2';
    session='RC4103-1';
end
physioDir='/projects2/udall/physio';
if exist(fullfile(physioDir,session,'rest_cardio.txt'),'file') && exist(fullfile(physioDir,session,'rest_resp.txt'),'file')
    tic
    disp(session);
    %fmri
    maskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
    funcDir='/projects2/udall/pboord/pic/preproc/pestica';
    tr=2.4;
    %restS=MRIread(fullfile(funcDir,session,'rest_brain_reg.nii'));
    restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],'filtered_func_data.nii'));
    maskPath=fullfile(funcDir,session,maskFile);
    if ~exist(maskPath,'file')
        system(['cd ',fullfile(funcDir,session),...
            ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',maskFile]);
    end
    maskS=MRIread(maskPath);
    nvol=restS.nframes;
    volsize=restS.volsize;
    nvox=prod(volsize(1:2));
    nslice=volsize(3);
    nTotalSlices=nvol*nslice;
    nfft=30;
    nyqIdx=nfft/2+1;
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
    voxInSlice=cell(nslice,1);
    for sliceNr=1:nslice
        cohRestSlice=squeeze(rest(mask(:,sliceNr),sliceNr,:));
        cohRestSlice=detrend(cohRestSlice')';
        voxInSlice{sliceNr}=cohRestSlice;
    end
    sqssSlice=zeros(nslice,nvol);
    for sliceNr=1:nslice
        sqssSlice(sliceNr,:)=sqrt(mean(voxInSlice{sliceNr}.^2,1));
    end
    sqssSlice=detrend(sqssSlice')';
    nrTr=2;
    dfSlice=1/(nrTr*tr);
    
    % physio
    physioFs=100; % actually 5 replications of 100 Hz = 500 Hz, but replications removed
    physioSamples=physioFs*nvol*tr;
    % cardio
    cardioRaw=reshape(textread(fullfile(physioDir,session,'rest_cardio.txt')),5,[]);
    cardioRaw=cardioRaw(1,:);
    cardio=resample(cardioRaw,nTotalSlices,physioSamples); %#ok<*DTXTRD>
    % resp
    respRaw=reshape(textread(fullfile(physioDir,session,'rest_resp.txt')),5,[]);
    respRaw=respRaw(1,:);
    resp=resample(respRaw,nTotalSlices,physioSamples);
    
    for volNr=1:10
        figure('WindowStyle','docked');
        nrPlot=5;
        x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
        t=(x-1)*tr/nslice;
        nrX=numel(x);
        subplot(nrPlot,1,1);
        sqssWave=reshape(sqssSlice(x),[],1);
        padFactor=6;
        sqssSliceFft=fft(sqssWave,padFactor*numel(sqssWave));
        %sqssSliceFft=fft(sqssWave.*hann(numel(sqssWave),'periodic'),padFactor*numel(sqssWave));
        plot(t,sqssWave);
        subplot(nrPlot,1,2);
        fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
        minCardioHz=0.6;
        maxCardioHz=1.5;
        minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
        maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
        [maxSpec,maxIdx]=max(abs(sqssSliceFft(minCardio:maxCardio)));
        maxSpecIdx=minCardio+maxIdx-1;
        cardioEst=abs(sqssSliceFft(maxSpecIdx))*cos(2*pi*fSlice(maxSpecIdx)*(0:nrX-1)*tr/nslice+angle(sqssSliceFft(maxSpecIdx)));
        plot(fSlice(1:floor(end/2)),abs(sqssSliceFft(1:floor(end/2))),'b',fSlice(maxSpecIdx),maxSpec,'r*');
%         % sidebands
%         hold on
%         sbWidth=round(0.5/(dfSlice/padFactor))+1;
%         [minSb,minIdx]=findpeaks(abs(sqssSliceFft(maxSpecIdx-sbWidth:maxSpecIdx)),'NPeaks',2,'SortStr','descend');
%         [maxSb,maxIdx]=findpeaks(abs(sqssSliceFft(maxSpecIdx:maxSpecIdx+sbWidth)),'NPeaks',2,'SortStr','descend');
%         minSbIdx=maxSpecIdx-sbWidth+minIdx(1)-1;
%         maxSbIdx=maxSpecIdx+maxIdx(1)-1;
%         %plot(fSlice(1:floor(end/2)),abs(sqssSliceFft(1:floor(end/2))),'b',fSlice(maxSpecIdx),maxSpec,'r*',fSlice(minSbIdx),minSb(1),'g*',fSlice(maxSbIdx),maxSb(1),'g*');
%         plot(fSlice(minSbIdx),minSb(1),'g*',fSlice(maxSbIdx),maxSb(1),'g*');
%         hold off
%         respEst=abs(sqssSliceFft(maxSpecIdx))*cos(2*pi*fSlice(maxSbIdx-maxSpecIdx)*(0:nrX-1)*tr/nslice+angle(sqssSliceFft(maxSbIdx-maxSpecIdx)));
        
        subplot(nrPlot,1,3);
        plot(t,resp(x));
        subplot(nrPlot,1,4);
        [pks,locs]=findpeaks(cardio(x),'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
        plot(t,cardio(x));
        hold on
        plot(t(locs),pks,'r*');
        plot(t,cardioEst,'g');
        maxSpecAngle=angle(sqssSliceFft(maxSpecIdx));
        % xfm (-pi,pi) to (0,2*pi)
        if maxSpecAngle < 0
            maxSpecAngle=2*pi+maxSpecAngle;
        end
        nrPeakInTr=floor(tr*fSlice(maxSpecIdx)+1);
        peakTimes=(2*pi*(1:nrPeakInTr)-maxSpecAngle)/(2*pi*fSlice(maxSpecIdx));
        peakIdxs=round(peakTimes*nslice/tr)+1;
        plot(t(1)+(peakIdxs-1)*tr/nslice,cardioEst(peakIdxs),'m*');
        hold off
%         subplot(nrPlot,1,5);
%         sbRange=maxSpecIdx-sbWidth:maxSpecIdx+sbWidth;
%         unwrappedPhase=unwrap(angle(sqssSliceFft));
%         unwrappedPhase=unwrappedPhase-unwrappedPhase(maxSpecIdx);
%         plot(fSlice(sbRange),unwrappedPhase(sbRange),'b',fSlice(maxSpecIdx),unwrappedPhase(maxSpecIdx),'r*',fSlice(minSbIdx),unwrappedPhase(minSbIdx),'g*',fSlice(maxSbIdx),unwrappedPhase(maxSbIdx),'g*');
%         phaseDiff=mod(unwrappedPhase(minSbIdx)-unwrappedPhase(maxSbIdx),2*pi)/pi

        %         N = 80;                                                             % Order of Filter
%         Wn = 0.1;                                                    % Pass Band Edge Frequency.
%         a = fir1(N,Wn);                                             % Return Numerator of Low Pass FIR filter
%         b = 1;                                                                % Denominator of Low Pass FIR Filter
%         filtDemod = filtfilt(a,b,abs([0;diff(sqssWave)]));
%         plot(1:nrX,filtDemod,'m');
        
%         cardioWidth=round(0.2/(dfSlice/padFactor))+1;
%         cardioIdxs=(maxSpecIdx-cardioWidth:maxSpecIdx+cardioWidth)';
%         sizePadFft=size(sqssSliceFft,1);
%         cardioSliceFft=complex(zeros(sizePadFft,1),zeros(sizePadFft,1));
%         cardioSliceFft(cardioIdxs)=sqssSliceFft(cardioIdxs);
%         imCardioIdxs=sizePadFft-cardioIdxs+2;
%         cardioSliceFft(imCardioIdxs)=sqssSliceFft(imCardioIdxs);
%         ifftCardioSliceFft=ifft(cardioSliceFft);
%         ifftCardioSliceFft=ifftCardioSliceFft(1:numel(sqssWave));
%         plot(1:nrX,ifftCardioSliceFft,'m');
%         %plot(locs(2:end),(diff(locs*tr)/nslice).^-1);
    end
end