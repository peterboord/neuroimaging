function savePulse(dataDir,fileName)

physioFile = load(fullfile(dataDir,fileName));
if ~strcmp(physioFile.output.physio.names{2},'Pulse'), break; end
pulse = physioFile.output.physio.data(:,2);
pulseResamp = resample(pulse,22,500);
save(fullfile(dataDir,'pulse22Hz'),pulseResamp);