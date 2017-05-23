function [peakTimes,peakLocs,normPhysio,sampleTimes]=getPhysioPeaks(physioFile,maxCardioPeriod,minPeakProminence)

physioFs=100; % Hz
try
    normPhysio=textread(physioFile);
catch
    normPhysio=[];
end
if mod(numel(normPhysio),5)~=0
    normPhysio=[];
end
if isempty(normPhysio)
    peakTimes=[];
    peakLocs=[];
    sampleTimes=[];
else
    normPhysio=reshape(normPhysio,5,[]);
    normPhysio=detrend(normPhysio(1,:)');
    fn=physioFs/2;
    pbf=2;
    sbf=12;
    d = designfilt('lowpassfir', ...
        'PassbandFrequency',pbf/fn,'StopbandFrequency',sbf/fn, ...
        'PassbandRipple',1,'StopbandAttenuation',60, ...
        'DesignMethod','equiripple');
    normPhysio = filtfilt(d,normPhysio);
    normPhysio=normPhysio/std(normPhysio);
    [~,peakLocs]=findpeaks(normPhysio,...
        'MinPeakDistance',round(0.5*physioFs*maxCardioPeriod),...
        'MinPeakProminence',minPeakProminence);
    peakTimes=(peakLocs'-1)/physioFs;
    sampleTimes=((1:numel(normPhysio))-1)/physioFs;
end
%figure, plot(sampleTimes,normPhysio,'b',peakTimes,normPhysio(peakLocs),'r*');
end