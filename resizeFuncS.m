function funcS = resizeFuncS(funcS,dim,nframes)
funcS.dim(1) = dim;
if nargin > 2
    if nframes < funcS.nframes
        funcS.vol(:,:,:,nframes+1:end) = [];
    else
        funcS.vol = repmat(funcS.vol(:,:,:,1),[1,1,1,nframes]);
    end
    funcS.nframes = nframes;
    funcS.dim(5) = nframes;
else
    funcS.dim(5) = 1;
    funcS.nframes = 1;
end
end
