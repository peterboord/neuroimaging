function sW = slideCorr(wav1,wav2,w)
n=numel(wav1);
nW=n-w+1;
wI=repmat(1:nW,w,1)+repmat(0:w-1,nW,1)';
sW=[zeros(1,ceil(w/2)-1),corrColumns(reshape(wav1(wI),w,nW),reshape(wav2(wI),w,nW)),zeros(1,floor(w/2))]';
end