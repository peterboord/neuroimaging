function fcMap(funcPath,roiPath,fcMapPath,mefc)

if nargin < 4
    mefc = 0;
end

switch class(mefc)
    case 'char'
        mefc = str2double(mefc);
end

funcStruct = MRIread(funcPath);
if ~isempty(funcStruct)
    roiStruct = MRIread(roiPath);
    nrColumns = prod(roiStruct.volsize);
    roi2d  = reshape(roiStruct.vol,[1,nrColumns]);
    func2d  = reshape(permute(funcStruct.vol,[4,1,2,3]),[funcStruct.nframes,nrColumns]);
    func2d = detrend(func2d);
    seedTs = mean(func2d(:,roi2d>0),2);
    roiStruct.vol = reshape(corr2z(corr(seedTs,func2d),funcStruct.nframes,mefc),roiStruct.volsize);
    roiStruct.fspec = fcMapPath;
    MRIwrite(roiStruct,roiStruct.fspec);
end
end

function z = corr2z(corrMap,N,mefc)
sizeCorrMap = size(corrMap);
corrMap(isnan(corrMap)) = 0;
clamp = 0.999;
corrMap(corrMap > clamp) = clamp;
corrMap(isnan(corrMap)) = 0;
if mefc
    % From README.mefc:
    % Denoised time series may be used for task fMRI analysis or seed-based functional connectivity analysis.
    % Hint: try computing seed-based Pearson correlation using the rest30_mefc.nii.gz dataset,
    % then compute statistical significance using the standard Fisher Z-transform
    % with the number of components as the degrees of freedom.
    % This is already implemented in AFNI InstaCorr.
    %
    % From http://users.fmrib.ox.ac.uk/~stuart/thesis/chapter_6/section6_4.html
    % The correlation coefficient, r, can be transformed so that it has a Z distribution
    % (that is a Gaussian distribution with zero mean and unit variance), by applying the Fisher Z transform
    %
    % For a time course of n independent time points, there are n-2 degrees of freedom.
    % So in Fisher's Z-transform n-2 is replaced by nframes for mefc
    z = atanh(corrMap)*sqrt(N - 1);
    %z = (sqrt(N-1)/2)*log((ones(sizeCorrMap)+corrMap)./(ones(sizeCorrMap)-corrMap));
else
    % From http://users.fmrib.ox.ac.uk/~stuart/thesis/chapter_6/section6_4.html
    % The correlation coefficient, r, can be transformed so that it has a Z distribution
    % (that is a Gaussian distribution with zero mean and unit variance), by applying the Fisher Z transform
    z = atanh(corrMap)*sqrt(N - 3);
    %z = (sqrt(N-3)/2)*log((ones(sizeCorrMap)+corrMap)./(ones(sizeCorrMap)-corrMap));
end
end