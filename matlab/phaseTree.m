function phaseTree(session)

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
mag1S=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_m1mag.nii'));
phase1S=MRIread(fullfile(funcDir,session,'mocoBetHpf.feat','filtered_func_data_m1phase.nii'));
volsize=mag1S.volsize;
nslice=volsize(3);
nrVol=10;
ves=zeros([volsize,nrVol]);
shiftVec=zeros(27,3);
vecIdx=1;
for ap=-1:1
    for xy=-1:1
        for z=-1:1
            shiftVec(vecIdx,:)=[ap,xy,z];
            vecIdx=vecIdx+1;
        end
    end
end
nrNear=27;
shiftVec=[shiftVec(14,:);shiftVec(1:13,:);shiftVec(15:nrNear,:)];
magPad=cat(3,zeros(volsize(1:2)),mag1S.vol,zeros(volsize(1:2)));
phasePad=cat(3,zeros(volsize(1:2)),phase1S.vol,zeros(volsize(1:2)));
padSize=size(magPad);
[ap,xy,z]=ind2sub(padSize,1:prod(padSize));
nearIdx=reshape(repmat(reshape([ap',xy',z'],1,[],3),nrNear,1,1)+repmat(reshape(shiftVec,[nrNear,1,3]),1,prod(padSize),1),[],3);
nearIdx(any(nearIdx==0,2),:)=1;
nearIdx(any(nearIdx>repmat(padSize,size(nearIdx,1),1),2),:)=1;
nearIdx=reshape(sub2ind(padSize,nearIdx(:,1),nearIdx(:,2),nearIdx(:,3)),nrNear,[]);
for volNr=1:nrVol
    disp(volNr);
    tic
    phaseRes=volNr*(pi/100);
    colorVol=zeros(volsize);
    colorVol(mag1S.vol(:)>0)=-sum(mag1S.vol(:)>0):-1;
    colorPad=cat(3,zeros(volsize(1:2)),colorVol,zeros(volsize(1:2)));
    phaseNear=phasePad(nearIdx);
    magNear=magPad(nearIdx);
    voxNear=abs(angle(repmat(exp(1i*phaseNear(1,:)),nrNear-1,1).*exp(-1i*phaseNear(2:end,:))));
    voxNear(magNear(2:end,:)==0)=Inf;
    voxNear(:,magPad(:)==0)=Inf;
    closestVox=min(voxNear,[],1);
    closestVox(isinf(closestVox))=0;
    voxNear=voxNear<=phaseRes;
    voxNear=voxNear & magNear(2:end,:)>0;
    anyVoxNear=any(voxNear,1).*reshape(magPad>0,1,[]);
    colorNear=colorPad(nearIdx).*cat(1,anyVoxNear,voxNear);
    anyVoxNearIdx=find(anyVoxNear);
    colorNearShort=colorNear(:,anyVoxNearIdx);
    for voxNr=1:numel(anyVoxNearIdx)
        voxColorNear=colorNearShort(:,voxNr);
        voxColorNear(voxColorNear==0)=[];
        [voxCol,minIdx]=min(voxColorNear);
        voxColorNear(minIdx)=[];
        for colChange=voxColorNear'
            colorNearShort(colorNearShort==colChange)=voxCol;
        end
    end
    voxColor=zeros(1,prod(padSize));
    voxColor(anyVoxNearIdx)=min(colorNearShort,[],1);   
    voxColor=reshape(voxColor,padSize);
    voxColor=voxColor(:,:,2:end-1);
    [C,~,ic]=unique(voxColor(:));
    C=0:numel(C)-1;
    voxColor=reshape(C(ic),volsize);
    ves(:,:,:,volNr)=voxColor;
    toc
end
MRIsave(mag1S,ves,fullfile(fileparts(mag1S.fspec),'ves.nii'),nrVol);
end
