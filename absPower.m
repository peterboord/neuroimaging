function absPower(dataDir,niiName,TR)

switch nargin
    case 1
        nii = dir(fullfile(dataDir,'*.nii'));
        niiName = nii.name;
        fileOutName = 'absPower.nii';
    case {2,3}
        fileOutName = ['absFft_',niiName];
end
func = load_untouch_nii(fullfile(dataDir,niiName));
T = TR*size(func.img,4);
df = 1/T;
fftFunc = fft(double(permute(func.img,[4,1,2,3])));
fftFunc(round(end/2)+1:end,:,:,:) = [];
fftFunc(1,:,:,:) = 0;
func.img = permute(single(abs(fftFunc)),[2,3,4,1]);
func.hdr.dime.dim(5) = size(func.img,4);
save_untouch_nii(func,fullfile(dataDir,fileOutName));





