function MRIsave(funcS,vol,fspec,nvol)

funcS.vol=vol;
clear vol
funcS.fspec=fspec;
if nargin<4
    MRIwrite(funcS,fspec);
else
    %funcS.dim(5) = nvol;
    funcS.dim(4) = nvol;
    funcS.nframes = nvol;
    MRIwrite(funcS,fspec);
end