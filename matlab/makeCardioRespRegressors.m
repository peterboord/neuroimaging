function makeCardioRespRegressors(session)
dbstop if error
if nargin==0
    %session='RC4107-2';
    session='RC4103-1';
    %session='RC4206-1';
end
physioDir='/projects2/udall/physio';
if exist(fullfile(physioDir,session,'rest_cardio.txt'),'file') && exist(fullfile(physioDir,session,'rest_resp.txt'),'file')
    tic
    disp(session);
    %fmri
    %brainMaskS=MRIread(fullfile('/projects2/udall/rest/fast',
    varMaskFile='rest_brain_reg_hpf_std_thrP95_bin.nii';
    funcDir='/projects2/udall/pboord/pic/preproc/pestica';
    tr=2.4;
    %restS=MRIread(fullfile(funcDir,session,'rest_brain_reg.nii'));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest.feat'],'filtered_func_data.nii'));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_brain.nii']));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest.nii.gz']));RC4103-1_rest.slicemocoxy_afni.zalg_moco2.nii
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_despike.nii']));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_hpf.nii']));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_skinMask.nii']));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_skinMask.feat'],'filtered_func_data.nii'));
    %restS=MRIread(fullfile(funcDir,session,[session,'_rest_ribbon.nii']));
    restS=MRIread(fullfile(funcDir,session,[session,'_rest.slicemocoxy_afni.zalg_moco2.nii']));
    maskPath=fullfile(funcDir,session,varMaskFile);
    if ~exist(maskPath,'file')
        system(['cd ',fullfile(funcDir,session),...
            ';fslmaths rest_brain_reg -bptf 41.7 -1 rest_brain_reg_hpf;fslmaths rest_brain_reg_hpf -Tstd rest_brain_reg_hpf_std;fslmaths rest_brain_reg_hpf_std -thrP 95 -bin ',varMaskFile]);
    end
    maskS=MRIread(maskPath);
    nvol=restS.nframes;
    volsize=restS.volsize;
    nslice=volsize(3);
    nTotalSlices=nvol*nslice;
    rest=reshape(restS.vol,prod(volsize),nvol);
    rest=detrend(rest')';
    rest=reshape(rest,volsize(1)*volsize(2),nslice,nvol);
    mask=logical(reshape(maskS.vol,volsize(1)*volsize(2),nslice));
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
    % calc center of mass
    % ref: 2001 Raj_Respiratory effects in human functional magnetic resonance imaging due to bulk susceptibility changes
    paIdx=repmat((1:volsize(1))',1,volsize(2));
    lrIdx=repmat(1:volsize(2),volsize(1),1);
    paCom=zeros(nvol,nslice);
    lrCom=zeros(nvol,nslice);
    sumSliceImg=zeros(nvol,nslice);
    for volNr=1:nvol
        for sliceNr=1:nslice
            sliceImg=restS.vol(:,:,sliceNr,volNr);
            paCom(volNr,sliceNr)=sum(paIdx(:).*sliceImg(:))/sum(sliceImg(:));
            lrCom(volNr,sliceNr)=sum(lrIdx(:).*sliceImg(:))/sum(sliceImg(:));
            sumSliceImg(volNr,sliceNr)=sum(sliceImg(:));
        end
    end
    slo=load('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/RC4103-1_rest_pestica/slomoco.combined.reg6dof.txt');
    slo=reshape(slo,nslice,nvol,6);
    paCom=paCom';
    sumSliceImg=sumSliceImg';
    x=1:430;
    % remove linear trend
    paComDetrend=detrend(paCom')';
    paComDetrendNorm=bsxfun(@rdivide,paComDetrend,std(paComDetrend,1,2));
    sumSliceDetrend=detrend(sumSliceImg')';
    sumSliceDetrendNorm=bsxfun(@rdivide,sumSliceDetrend,std(sumSliceDetrend,1,2));
    figure
    biplot(x,sumSliceDetrend(x),x,resp(x),'slice number','A-P center of mass','slice number','Respiration');
    figure
    biplot(x,sumSliceDetrendNorm(x),x,resp(x),'slice number','A-P center of mass minus mean over scan','slice number','Respiration');
    keyboard
    lrCom=lrCom';
    [b,bint,r,rint,stats] = regress(lrCom(:),[ones(nvol*nslice,1),slo(:,1:6)]);
    sumSliceImg=bsxfun(@minus,sumSliceImg,mean(sumSliceImg,1));
    sumSliceImg=bsxfun(@rdivide,sumSliceImg,std(sumSliceImg,1,1));
    sumSliceImg=bsxfun(@rdivide,sumSliceImg,mean(sumSliceImg,1));
    dfPaComFft=1/(nvol*tr);
    fPaComFft=(0:dfPaComFft:(nvol-1)*dfPaComFft)';
    paComFft=complex(zeros(nvol,nslice),zeros(nvol,nslice));
    sumSliceImgFft=complex(zeros(nvol,nslice),zeros(nvol,nslice));
    for sliceNr=1:nslice
        paComFft(:,sliceNr)=fft(detrend(paCom(sliceNr,:)'));
        sumSliceImgFft(:,sliceNr)=fft(detrend(sumSliceImg(sliceNr,:)'));
    end
    paComFftVol=complex(zeros(nvol,nslice),zeros(nvol,nslice));
    df=1/tr;
    f=0:df:(nslice-1)*df;
    for volNr=1:nvol
        paComFftVol(volNr,:)=fft(detrend(paComDetrendNorm(:,volNr)));
        sumSliceImgFft(volNr,:)=fft(detrend(sumSliceImg(:,volNr)));
    end
    figure,plot(f(1:22),abs(paComFftVol(:,1:22))')
%     winLenInTr=20;% must be an even number, unless code modified
%     avIdx=repmat((0:winLenInTr-1)',1,nvol-winLenInTr+1)+repmat(1:nvol-winLenInTr+1,winLenInTr,1);
%     avIdx=[repmat(avIdx(:,1),1,winLenInTr/2-1),avIdx,repmat(avIdx(:,end),1,winLenInTr/2)];
%     paCom=paCom-squeeze(mean(reshape(paCom(avIdx(:),:),winLenInTr,nvol,nslice)));
%     lrCom=lrCom-squeeze(mean(reshape(lrCom(avIdx(:),:),winLenInTr,nvol,nslice)));
%     
%     paCom=bsxfun(@rdivide,paCom,std(paCom,1,1));
%     lrCom=bsxfun(@rdivide,lrCom,std(lrCom,1,1));

    nrTr=5;
    cardioDelayPhase=-0.12*2*pi;
    respDelayPhase=0;%0.4*2*pi;
    dfSlice=1/(nrTr*tr);
    %padFactor=6;
    padFactor=1;
    fSlice=0:dfSlice/padFactor:dfSlice*(nrTr*nslice);fSlice(end)=[];
    for volNr=1:10
        figure('WindowStyle','docked');
        nrPlot=6;
        x=(volNr-1)*nslice+1:(volNr+nrTr-1)*nslice;
        t=(x-1)*tr/nslice;
        nrX=numel(x);
        
        subplot(nrPlot,1,1);
        sqssWave=reshape(sqssSlice(x),[],1);
        sqssSliceFft=fft(detrend(sqssWave,'constant'),padFactor*numel(sqssWave));
        %sqssSliceFft=fft(sqssWave.*hann(numel(sqssWave),'periodic'),padFactor*numel(sqssWave));
        plot(t,sqssWave);
        
        subplot(nrPlot,1,2);
        minCardioHz=0.6;
        maxCardioHz=1.5;
        minCardio=round(minCardioHz/(dfSlice/padFactor))+1;
        maxCardio=round(maxCardioHz/(dfSlice/padFactor))+1;
        [maxCardioSpec,maxCardioIdx]=max(abs(sqssSliceFft(minCardio:maxCardio)));
        maxCardioSpecIdx=minCardio+maxCardioIdx-1;
        cardioEst=abs(sqssSliceFft(maxCardioSpecIdx))*cos(2*pi*fSlice(maxCardioSpecIdx)*(0:nrX-1)*tr/nslice+angle(sqssSliceFft(maxCardioSpecIdx))+cardioDelayPhase);
        plot(fSlice(1:floor(end/2)),abs(sqssSliceFft(1:floor(end/2))),'b',fSlice(maxCardioSpecIdx),maxCardioSpec,'r*');
        
        subplot(nrPlot,1,3);
        [pks,locs]=findpeaks(cardio(x),'MinPeakHeight',0,'MinPeakDistance',nslice/tr);
        plot(t,cardio(x));
        hold on
        plot(t(locs),pks,'r*');
        plot(t,cardioEst,'g');
        maxCardioSpecAngle=angle(sqssSliceFft(maxCardioSpecIdx))+cardioDelayPhase;
        % xfm (-pi,pi) to (0,2*pi)
        if maxCardioSpecAngle < 0
            maxCardioSpecAngle=2*pi+maxCardioSpecAngle;
        end
        nrPeakInWin=floor(nrTr*tr*fSlice(maxCardioSpecIdx)+1);
        peakTimes=(2*pi*(1:nrPeakInWin)-maxCardioSpecAngle)/(2*pi*fSlice(maxCardioSpecIdx));
        peakIdxs=round(peakTimes*nslice/tr)+1;
        peakIdxs(peakIdxs>numel(cardioEst))=[];
        plot(t(1)+(peakIdxs-1)*tr/nslice,cardioEst(peakIdxs),'m*');
        hold off
        
        subplot(nrPlot,1,4);
        respCom=paCom;%lrCom;%
        respComWave=reshape(respCom(x(:)),[],1);
        respComFft=fft(detrend(respComWave,'constant'),padFactor*numel(respComWave));
        %respComFft=fft(respComWave.*hann(numel(respComWave),'periodic'),padFactor*numel(respComWave));
        plot(t,respComWave);
        
        subplot(nrPlot,1,5);
        minRespHz=0.15;
        maxRespHz=0.4;
        minResp=round(minRespHz/(dfSlice/padFactor))+1;
        maxResp=round(maxRespHz/(dfSlice/padFactor))+1;
        [maxRespSpec,maxRespIdx]=max(abs(respComFft(minResp:maxResp)));
        maxRespSpecIdx=minResp+maxRespIdx-1;
        respEst=abs(respComFft(maxRespSpecIdx))*cos(2*pi*fSlice(maxRespSpecIdx)*(0:nrX-1)*tr/nslice+angle(respComFft(maxRespSpecIdx))+respDelayPhase);
        plot(fSlice(1:floor(end/2)),abs(respComFft(1:floor(end/2))),'b',fSlice(maxRespSpecIdx),maxRespSpec,'r*');       
        
        subplot(nrPlot,1,6);
        plot(t,resp(x)*std(respEst)/std(resp(x)));
        hold on
        plot(t,respEst,'g');
        maxRespSpecAngle=angle(respComFft(maxRespSpecIdx))+respDelayPhase;
        % xfm (-pi,pi) to (0,2*pi)
        if maxRespSpecAngle < 0
            maxRespSpecAngle=2*pi+maxRespSpecAngle;
        end
        nrPeakInWin=floor(nrTr*tr*fSlice(maxRespSpecIdx)+1);
        peakTimes=(2*pi*(1:nrPeakInWin)-maxRespSpecAngle)/(2*pi*fSlice(maxRespSpecIdx));
        peakIdxs=round(peakTimes*nslice/tr)+1;
        peakIdxs(peakIdxs>numel(respEst))=[];
        plot(t(1)+(peakIdxs-1)*tr/nslice,respEst(peakIdxs),'m*');
        hold off

    end
end