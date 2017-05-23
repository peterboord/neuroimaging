function physioQA

dbstop if error
physioDir='/projects2/act-plus/physio';

physioFiles=dir(fullfile(physioDir,'*_rest_cardio.txt'));
physioFiles={physioFiles.name}';
nrFiles=numel(physioFiles);
for fileNr=1:nrFiles
    [~,fileName]=fileparts(physioFiles{fileNr});
    tuneParamPath=fullfile(physioDir,[fileName,'_tuneParam.mat']);
    [peakTimes,peakLocs,normPhysio,sampleTimes]=getPhysioPeaks(fullfile(physioDir,physioFiles{fileNr}),tuneParamPath);
    [fixedPeakTimes,fixedPeakLocs]=fixPhysioTimes(peakTimes,normPhysio,peakLocs,sampleTimes);
    %[fixedPeakTimes,fixedPeakLocs]=fixPhysioTimes(peakTimes);
end
end