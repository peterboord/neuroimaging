function splitPhaseVol(session)

dbstop if error
if nargin==0
    session='RC4103-1';
    %session='RC4107-2';
    %session='RC4109-1';
end
disp(session);
%fmri
funcDir='/projects2/udall/pboord/pic/preproc/pestica';
% constants
phase1S=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','pkErpPhase.nii'));
maskS=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','mask.nii'));
maskVol=reshape(maskS.vol,[],1);
phaseVol=reshape(phase1S.vol,[],1);
phaseVol(phaseVol==0)=Inf;
volsize=phase1S.volsize;
nrPhaseBins=8;
phaseRes=2*pi/nrPhaseBins;
phaseVol4d=zeros([prod(volsize),nrPhaseBins]);
for binNr=1:nrPhaseBins
    hiPhase=pi-(binNr-1)*phaseRes;
    loPhase=hiPhase-phaseRes;
    phaseVol4d(:,binNr)=(phaseVol>=loPhase & phaseVol<hiPhase)*binNr;
end
phaseVol4d(phaseVol4d<0)=phaseVol4d(phaseVol4d<0)+2*pi;
phaseVol4d(isinf(phaseVol4d))=0;
%phaseVol4d=bsxfun(@times,phaseVol4d,maskVol);
MRIsave(phase1S,reshape(phaseVol4d,[volsize,nrPhaseBins]),fullfile(fileparts(phase1S.fspec),'phaseVol.nii'),nrPhaseBins);
end

function stuff
curZ=29;
phaseSlice=phase1S.vol(:,:,curZ);
dxPhase=angle(exp(1i*(phaseSlice(2:end,:)-phaseSlice(1:end-1,:))));
dxPhase=dxPhase(:,1:end-1);
dyPhase=angle(exp(1i*(phaseSlice(:,2:end)-phaseSlice(:,1:end-1))));
dyPhase=dyPhase(1:end-1,:);
[x,y]=ind2sub([volsize(1)-1,volsize(2)-1],1:((volsize(1)-1)*(volsize(2)-1)));
xy=1:((volsize(1)-1)*(volsize(2)-1));
figure
quiver(x,y,cos(dxPhase(xy)),cos(dyPhase(xy)));


end