% fluidS=MRIread('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/rest_brain_reg_hpf_masStd.nii');
fluidS=MRIread('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/rest_brain_reg_hpf_s24mm.nii');
fluid3d=reshape(fluidS.vol,fluidS.volsize(1)*fluidS.volsize(2),[],fluidS.nframes);
fluidMaskS=MRIread('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/rest_brain_reg_hpf_std_thr40_bin.nii');
mask=logical(fluidMaskS.vol);
% % minus voxel mean
% fluid=bsxfun(@minus,fluid,mean(fluid,4));
% average masked voxels in slice
%fluid=bsxfun(@rdivide,squeeze(sum(reshape(fluid,[80*80,43,300]),1)),sum(reshape(fluidMaskS.vol>0,[80*80,43]),1)');
% find voxels (x,y,z,:) with contiguous voxels (x,y,z+1,:) - wrap around
% from end to beginning slice
contigMask = reshape(cat(3,mask(:,:,1:end-1) & mask(:,:,2:end),mask(:,:,end) & mask(:,:,1)),fluidS.volsize(1)*fluidS.volsize(2),[]);
% calc slice phase sum
slicePhases = zeros(fluidS.volsize(3),fluidS.nframes-1);
for sliceNr = 1:fluidS.volsize(3)
   slice1 = squeeze(fluid3d(contigMask(:,sliceNr),sliceNr,1:end-1));
   if sliceNr < fluidS.volsize(3)
       slice2 = squeeze(fluid3d(contigMask(:,sliceNr),sliceNr+1,2:end));
   else
       slice2 = squeeze(fluid3d(contigMask(:,sliceNr),1,2:end));
   end
%    phaseQuad = (((slice1 > 0))*1 + ((slice1 > slice2))*2)*pi/2;
%    slicePhases(sliceNr,:) = squeeze(angle(sum(exp(i*phaseQuad),1)));
   phaseQuad = slice2-slice1;
   slicePhases(sliceNr,:) = mean(phaseQuad,1);
end
cardio=resample(textread('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/rest_cardio.txt'),12900,360000);
slicePhases=slicePhases(:);
x=1:numel(slicePhases);
figure,plotyy(x,slicePhases,x,cardio(x))
for i=1:10
    figure('WindowStyle','docked');
    plot(slicePhases(i));
end
corr(slicePhases,cardio(x))
disp('ignores phase shift between pulse ox and fmri');
contigT1 = reshape(fluidS.vol,fluidS.volsize(1)*fluidS.volsize(2),[],fluidS.nframes);
contigT2 = contigT1;
contigT1(~repmat(cat(2,contigMask,false(fluidS.volsize(1)*fluidS.volsize(2),1)),1,1,fluidS.nframes)) = 0;
contigT1 = contigT1(:,1:end-1,:);
contigT2(~repmat(cat(2,false(fluidS.volsize(1)*fluidS.volsize(2),1),contigMask),1,1,fluidS.nframes)) = 0;
contigT2 = contigT2(:,1:end-1,:);
% estimate phase quadrant
phaseQuad = (contigT1(repmat(contigMask,1,1,300)) > contigT2(repmat(contigMask,1,1,300)))*1 & (contigT1(repmat(contigMask,1,1,300)) > 0)*2;
phaseQuad = (contigT1 > contigT2)*1 & (contigT1 > 0)*2;
phaseQuad = reshape(phaseQuad,fluidS.volsize(1)*fluidS.volsize(2),[],fluidS.nframes);
contigT1 = reshape(fluid(repmat(cat(3,contigMask,false(fluidMaskS.volsize(1:2))),1,1,1,fluidS.nframes)),[],fluidS.volsize(3),fluidS.nframes);
contigT2 = fluid(repmat(cat(3,false(fluidMaskS.volsize(1:2)),contigMask),1,1,1,fluidS.nframes));
fluid=squeeze(sum(reshape(fluid,[80*80,43,300]),1))./repmat(sum(reshape(mask,[80*80,43]),1)',1,300);
fluid=fluid(:);
if any(isnan(fluid(:))), error('nan present'); end
% % interpolate missing samples
% badLogical=isnan(fluid) | isinf(fluid);
% goodDataIndices=find(~badLogical);
% badLogical([1:(min(goodDataIndices)-1) (max(goodDataIndices)+1):end])=0;
% fluid(badLogical)=interp1(goodDataIndices,fluid(goodDataIndices),find(badLogical));
% clear badLogical goodDataIndices
% if isnan(fluid(end)), fluid(end)=fluid(end-1); end
% normalize slice variance
fluid=reshape(fluid,43,300);
fluid=bsxfun(@rdivide,fluid,var(fluid')');
cardio=resample(textread('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/rest_cardio.txt'),12900,360000);
fluidS.vol = repmat(reshape(cardio,[1,1,43,300]),80,80,1,1);
fluidS.fspec='/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/fluid_cardio.nii';
MRIwrite(fluidS,fluidS.fspec);
% normalize
fluid=(fluid-mean(fluid))/(max(fluid)-min(fluid));
cardio=(cardio-mean(cardio))/(max(cardio)-min(cardio));
% volume marker
volMarker=repmat([1;NaN;zeros(41,1)],300,1);
fluid=reshape(fluid,43,300);
fluidFft=fft(bsxfun(@times,fluid,ones(43,1)));%hann(43,'periodic')));
figure
plot(fluid(:,1));
x=zeros(43);
realFluid=zeros(43);
for j=1:43
    for k=1:43
        x(k,j)=(1/43)*fluidFft(k,1)*power(exp(-2*pi*1i/43),-(j-1)*(k-1));
    end
end
hold on
plot(real(sum(x,2)));
hold on
for j=1:43
    plot(real(x(:,j)))
end
figure
plot(sum(real(x),2))
for i=1:43
    realFluid(:,i)=abs(fluidFft(2,1))*cos(2*pi*(1:43)*(2.4/43)*(1/2.4)-angle(fluidFft(2,1)));
end
plot(sum(realFluid,2))
figure
plot([volMarker,cardio,fluid])
keyboard