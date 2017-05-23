function mAtan2(yFile,xFile,outFile)

if nargin < 3
    disp('Usage:');
    disp('mAtan2 $MCRROOT yFile xFile outFile');
    disp('outFile = atan2(yFile,xFile)');
else
    yS=MRIread(yFile);
    xS=MRIread(xFile);
    oS=yS;
    oS.vol=reshape(atan2(reshape(yS.vol,[],yS.nframes)',reshape(xS.vol,[],xS.nframes)'),[oS.volsize,oS.nframes]);
    oS.fspec=outFile;
    MRIwrite(oS,oS.fspec);
end