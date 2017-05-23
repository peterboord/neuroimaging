function csfThreshold

dataDir = '/projects2/udall/pboord/pic/preproc/pestica/RC4103-1';
adcS = MRIread(fullfile(dataDir,'dti_MD_in_T1.nii'));
gm_maskS = MRIread(fullfile(dataDir,'gm_mask.nii'));
csf_maskS = MRIread(fullfile(dataDir,'csf_mask.nii'));
skin_maskS=MRIread(fullfile(dataDir,'T1_brain_outskin_mask.nii'));
skull_maskS=MRIread(fullfile(dataDir,'T1_brain_outskull_mask.nii'));
adcS.vol(~skin_maskS.vol)=0;
adc=adcS.vol(:);
gm=logical(gm_maskS.vol(:));
csf=logical(csf_maskS.vol(:));

x=linspace(min(adc),max(adc),1000);
gmHist=hist(adc(gm),x);
gmHist=cumsum(gmHist/sum(gmHist));
csfHist=hist(adc(csf),x);
csfHist=cumsum(csfHist/sum(csfHist));
figure
plot(x,csfHist,x,gmHist);
legend('csf','gmHist')
% set threshold for 95% of gm
gmThreshPercent=0.95;
[~,gmThreshIdx] = min(abs(gmHist-gmThreshPercent));
gmThreshVal=x(gmThreshIdx);
% threshold ADC to make CSF image
adcS.vol(adcS.vol<gmThreshVal)=0;
% discard anterior & lateral adc values outside skull mask
nonCsfS=adcS;
% find image boundaries
apBoundaries=any(skin_maskS.vol(:,:,floor(end/2)),2);
apBoundaries=[find(apBoundaries,1,'first'),find(apBoundaries,1,'last')];
lrBoundaries=any(skin_maskS.vol(:,:,floor(end/2)),1);
lrBoundaries=[find(lrBoundaries,1,'first'),find(lrBoundaries,1,'last')];
nonCsfS.vol(skull_maskS.vol>0)=0;
nonCsfS.vol(1:floor(end/2),floor(end/4):floor(3*end/4),:)=0;
adcS.vol(nonCsfS.vol>0)=0;
nonCsfS.fspec=fullfile(dataDir,'nonCsf.nii');
MRIwrite(nonCsfS,nonCsfS.fspec);
% % get eyeballs
% eyesS = adcS;
% eyesS.vol(skull_maskS.vol>0)=0;
% eyeClusters = bwconncomp(eyesS.vol);
% numPixels = cellfun(@numel,eyeClusters.PixelIdxList);
% [clusterSize,clusterIdx]=sort(numPixels,'descend');
% eyePixels = [eyeClusters.PixelIdxList{clusterIdx(1)};eyeClusters.PixelIdxList{clusterIdx(2)}];
% eyesS.vol(:)=0;
% eyesS.vol(eyePixels)=1;
% 
% eyesS.vol(:)=0;
% eyesS.vol(floor(end/2):end,:,:)=1;
% eyesS.fspec=fullfile(dataDir,'eyes.nii');
% MRIwrite(eyesS,eyesS.fspec);
adcS.fspec=fullfile(dataDir,'adc_csf.nii');
MRIwrite(adcS,adcS.fspec);
disp([num2str(gmHist(gmThreshIdx)*100),'% of gm retained']);
disp([num2str(csfHist(gmThreshIdx)*100),'% of csf discarded']);
end
