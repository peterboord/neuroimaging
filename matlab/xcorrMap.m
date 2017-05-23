function xcorrMap(funcFile)

mc=load('/projects2/udall/pboord/pic/preproc/pestica/RC4103-1/RC4103-1_rest+.feat/mc/prefiltered_func_data_mcf.par','-ascii');
physioDir='/projects2/udall/physio';
session='RC4103-1';
respRaw=reshape(textread(fullfile(physioDir,session,'rest_resp.txt')),5,[]);
respRaw=respRaw(1,:);
resp=resample(respRaw,300,numel(respRaw))';

dbstop if error
if nargin < 1
    funcFile='RC4103-1_rest_despike_hpf.nii';
end
dataDir='/projects2/udall/pboord/pic/preproc/pestica/RC4103-1';
eyesS=MRIread(fullfile(dataDir,'eyes.nii'));
stdMaskS=MRIread(fullfile(dataDir,'RC4103-1_rest_hpf_std_thrP99_bin.nii'));
funcS=MRIread(fullfile(dataDir,funcFile));
eyes=logical(reshape(eyesS.vol,prod(eyesS.volsize(1:2)),[]));
stdMask=logical(reshape(stdMaskS.vol,prod(eyesS.volsize(1:2)),[]));
func=reshape(funcS.vol,prod(eyesS.volsize(1:2)),[],funcS.nframes);
eyeZ=find(any(eyes,1));
nrEyeZ=numel(eyeZ);
eyeZtime=zeros(funcS.nframes,nrEyeZ);
for zNr=1:nrEyeZ
    eyeZtime(:,zNr)=squeeze(mean(func(eyes(:,eyeZ(zNr)),eyeZ(zNr),:),1));
end
end