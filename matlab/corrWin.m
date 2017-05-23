function corrWin(funcPath,fslCoord,idxs,outDir,prefix)
% fslCoord is FSL order, but 1-based
x=fslCoord(2);
y=fslCoord(1);
z=fslCoord(3);
funcS=MRIread(funcPath);
func2d=reshape(funcS.vol,[],funcS.nframes)';
for k=1:numel(idxs)-1
    win=idxs(k):idxs(k+1);
    vol=reshape(corr(func2d(win,:),func2d(win,sub2ind(funcS.volsize,x,y,z))),funcS.volsize);
    MRIsave(funcS,vol,fullfile(outDir,[prefix,'_',num2str(idxs(k),'%03u'),'_',num2str(idxs(k+1),'%03u'),'.nii.gz']),1);
end