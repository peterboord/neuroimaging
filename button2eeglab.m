dbstop if error


EEGsrate = 5000;
%eegDir = '/project_space/pboord/EEG-fMRI/data/Ojemann/test2/EEG/export';
eegDir = '/project_space/pboord/EEG-fMRI/data/Ojemann/EEG/export';
filename = {'tr3p1','tr3p2','tr3p3'};%'tr0p1','tr0p2','tr1p1','tr1p2',
for fileNr = 1:numel(filename)
%     % buttons
%     fid1 = fopen(fullfile(eegDir,[filename{fileNr},'press.txt']), 'r');
%     D = textscan(fid1, '%s %s %d %d %s', 'Delimiter', ',', 'HeaderLines', 2);
%     fclose(fid1);
%     points = double(D{:,3});
%     rowSelect = [true;diff(points)>=EEGsrate];
%     type = repmat({'B'},sum(rowSelect),1);
%     latency = (points(rowSelect)-1)/EEGsrate;
type=[];
latency=[];
    % R128
    if ~strcmp(filename{fileNr}(3),'0')
        fid1 = fopen(fullfile(eegDir,[filename{fileNr},'R128.txt']), 'r');
        D = textscan(fid1, '%s %s %d %d %s', 'Delimiter', ',', 'HeaderLines', 2);
        fclose(fid1);
        points = double(D{:,3});
        type = [type;repmat({'R128'},numel(points),1)];
        latency = [latency;(points-1)/EEGsrate];
    end
    % R
    fid1 = fopen(fullfile(eegDir,[filename{fileNr},'R.txt']), 'r');
    D = textscan(fid1, '%s %s %d %d %s', 'Delimiter', ',', 'HeaderLines', 2);
    fclose(fid1);
    points = double(D{:,3});
    type = [type;repmat({'R'},numel(points),1)];
    latency = [latency;(points-1)/EEGsrate];
    OT = table(type,latency);
    % figure
    % plot(diff(T{:,3}))
    writetable(OT,fullfile(eegDir,[filename{fileNr},'_EEGLAB.txt']));
end
