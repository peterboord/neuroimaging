function varargout = eegFmri3(varargin)
% EEGFMRI3 MATLAB code for eegFmri3.fig
%      EEGFMRI3, by itself, creates a new EEGFMRI3 or raises the existing
%      singleton*.
%
%      H = EEGFMRI3 returns the handle to a new EEGFMRI3 or the handle to
%      the existing singleton*.
%
%      EEGFMRI3('CALLBACK',hObj,eventData,h,...) calls the local
%      function named CALLBACK in EEGFMRI3.M with the given input arguments.
%
%      EEGFMRI3('Property','Value',...) creates a new EEGFMRI3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eegFmri3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eegFmri3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eegFmri3

% Last Modified by GUIDE v2.5 04-Feb-2014 09:19:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eegFmri3_OpeningFcn, ...
                   'gui_OutputFcn',  @eegFmri3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

% --- Outputs from this function are returned to the command line.
function varargout = eegFmri3_OutputFcn(hObj, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = h.output;

end


function hiFreq_Callback(hObj, eventdata, h)
% hObj    handle to hiFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of hiFreq as text
%        str2double(get(hObj,'String')) returns contents of hiFreq as a double

h.ud.hiFreq = str2double(get(hObj,'String'));
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function hiFreq_CreateFcn(hObj, eventdata, h)
% hObj    handle to hiFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end


function lowFreq_Callback(hObj, eventdata, h)
% hObj    handle to lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of lowFreq as text
%        str2double(get(hObj,'String')) returns contents of lowFreq as a double

h.ud.lowFreq = str2double(get(hObj,'String'));
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function lowFreq_CreateFcn(hObj, eventdata, h)
% hObj    handle to lowFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end

% --- Executes on selection change in taskSelect.
function taskSelect_Callback(hObj, eventdata, h)
% hObj    handle to taskSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObj,'String')) returns taskSelect contents as cell array
%        contents{get(hObj,'Value')} returns selected item from taskSelect

contents = cellstr(get(hObj,'String'));
h.ud.task = contents{get(hObj,'Value')};
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function taskSelect_CreateFcn(hObj, eventdata, h)
% hObj    handle to taskSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end

% --- Executes on selection change in roiSelect.
function roiSelect_Callback(hObj, eventdata, h)
% hObj    handle to roiSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObj,'String')) returns roiSelect contents as cell array
%        contents{get(hObj,'Value')} returns selected item from roiSelect

h.ud.roiNr = get(hObj,'Value');
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function roiSelect_CreateFcn(hObj, eventdata, h)
% hObj    handle to roiSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end

% --- Executes just before eegFmri3 is made visible.
function eegFmri3_OpeningFcn(hObj, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eegFmri3 (see VARARGIN)

% Choose default command line output for eegFmri3
h.output = hObj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize user data (ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defaults
def.lowFreq = 0.01;
def.hiFreq = 0.1;
% constants
h.ud.nrTr = 120; % tr
h.ud.isi = 20; % sec
h.ud.tr = 2.52; % sec
h.ud.nrRestTr = 240;
h.ud.fs = 256.4102564103; % Hz
h.ud.mriDir = '/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/MRI/PART1';
h.ud.taskEegDir = '/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/EEG/part2/export';
h.ud.restEegDir = '/NAS_II/Projects/pboord/Projects/EEG-fMRI/Data/Test2/EEG/part1/export';
h.ud.taskName = {'visual','auditory'};
h.ud.sessDirs = {'WIP_ME-VisReaction1_SENSE','WIP_ME-VisReaction2_SENSE','WIP_ME-VisReaction3_SENSE','WIP_ME-VisReaction4_SENSE';...
    'WIP_ME-AudReaction1_SENSE','WIP_ME-AudReaction2_SENSE','WIP_ME-AudReaction3_SENSE','WIP_ME-AudReaction4_SENSE'};
h.ud.roiList = {'roi1.txt','roi2.txt','roi3.txt','roi4.txt','roi5.txt','roi6.txt','roi7.txt'};
% derived from constants
hrLen = round(h.ud.isi/h.ud.tr);
nrEventPerSess = floor(h.ud.nrTr*h.ud.tr/h.ud.isi);
nrTask = size(h.ud.sessDirs,1);
nrSess = size(h.ud.sessDirs,2);
nrRoi = numel(h.ud.roiList);
%% update from gui
h.ud.hiFreq = str2double(get(h.hiFreq,'String'));
h.ud.lowFreq = str2double(get(h.lowFreq,'String'));
contents = cellstr(get(h.taskSelect,'String'));
h.ud.taskNr = get(h.taskSelect,'Value');
h.ud.roiNr = get(h.roiSelect,'Value');
%% update to gui
if isnan(h.ud.lowFreq)
    h.ud.lowFreq = def.lowFreq;
    set(h.lowFreq,'String',num2str(h.ud.lowFreq));
end
if isnan(h.ud.hiFreq)
    h.ud.hiFreq = def.hiFreq;
    set(h.hiFreq,'String',num2str(h.ud.hiFreq));
end
%% load task fMRI roi data
h.ud.roiFunc = cell(nrTask,nrSess,nrRoi);
tic
for taskNr = 1:nrTask
    for sessNr = 1:nrSess
        sessDir = fullfile(h.ud.mriDir,h.ud.sessDirs{taskNr,sessNr});
        funcS = MRIread(fullfile(sessDir,[h.ud.taskName{taskNr},'ME_medn.nii.gz']));
        for roiNr = 1:nrRoi
            roiS = MRIread(fullfile(sessDir,[h.ud.taskName{taskNr},'ME_tsoc_ro_withFILM.feat'],['roi',num2str(roiNr),'.nii.gz']));
            h.ud.roiFunc{taskNr,sessNr,roiNr} = getRoiInFunc(funcS,roiS);
        end
    end
end
% Load Laplacian EEG for roi
h.ud.eeg = load(fullfile(h.ud.taskEegDir,'visualALL_MR_sr256_300s_CB_rej3_C3_ascii'),'-ascii');
toc
% Update handles structure
guidata(hObj, h);
% update GUI
updateGui(h);

% UIWAIT makes eegFmri3 wait for user response (see UIRESUME)
% uiwait(h.figure1);

end

function roi = getRoiInFunc(funcS,roiS)
funcVoxT = reshape(funcS.vol,[],funcS.nframes);
roi = funcVoxT(reshape(roiS.vol,[],1)>0,:);
end

function updateGui(h)

% derived from constants
hrLen = round(h.ud.isi/h.ud.tr);
nrEventPerSess = floor(h.ud.nrTr*h.ud.tr/h.ud.isi);
nrSess = size(h.ud.sessDirs,2);
nrRoi = numel(h.ud.roiList);
%% BOLD
% plot taskBoldTimeseries
nrEvents = floor(h.ud.nrTr*h.ud.tr/h.ud.isi);
isiTr = ceil(h.ud.isi/h.ud.tr);
eventIdxs = reshape(repmat(round((0:nrEvents-1)*h.ud.isi/h.ud.tr),isiTr+1,1),1,[]);
eventIdxs = eventIdxs + repmat(1:isiTr+1,1,nrEvents);
boldTimeSeries = [];
eventData = [];
for sessNr = 1:nrSess
    sessData = h.ud.roiFunc{h.ud.taskNr,sessNr,h.ud.roiNr};
    avSessData = mean(sessData,1);
    sessEventData = reshape(avSessData(eventIdxs),isiTr+1,nrEvents);
    if isempty(boldTimeSeries)
        boldTimeSeries = sessData;
        eventData = sessEventData;
    else
        boldTimeSeries = cat(1,boldTimeSeries,sessData);
        eventData = cat(2,eventData,sessEventData);
    end
end
plot(h.ax_taskBoldTimeseries,boldTimeSeries');
% average BOLD time series in roi
plot(h.ax_taskBoldHr,eventData);
hold(h.ax_taskBoldHr,'on');
plot(h.ax_taskBoldHr,mean(eventData,2),'Color','k','LineWidth',3);
hold(h.ax_taskBoldHr,'off');
%% EEG
%hpf
d = fdesign.lowpass('Fp,Fst,Ap,Ast',50,60,3,20,h.ud.fs);
Hd = design(d);
eegLpf = filtfilt(Hd.Numerator,1,h.ud.eeg);
fsInt = 96;
eegLpf = resample(eegLpf,fsInt,floor(h.ud.fs));
fs = h.ud.fs*fsInt/floor(h.ud.fs);
% EEG time-frequency maps
[stft,tfFreq,tfTime] = spectrogram(eegLpf,hamming(fsInt),fsInt-1,fsInt,fs);
stft = stft(2:end,:);
% % smooth & downsample frequencies to 2 Hz resolution
% stft = squeeze(mean(reshape(stft,2,size(stft,1)/2,[]),1));
% seperate real and imag for k-means
% limit size bec of memory
% limit frequencies
freqLimit = 10;
stft = stft(1:freqLimit,:);
stft = [real(stft);imag(stft)];
nrFreqVars = size(stft,1);
stftLen = size(stft,2);
% limit time
timeLimit = 40;
memoryDivsor = floor((size(stft,2)/h.ud.fs)/timeLimit); % 60 = 5 secs
stftLen = floor(stftLen/memoryDivsor);
winLen = fsInt;
kMatrixIdx = repmat(1:stftLen-winLen,[winLen,1]) + repmat((0:winLen-1)',[1,stftLen-winLen]);
% offset to avoid startup/filter transients
offsetSeconds = 60;
kMatrixIdx = kMatrixIdx + offsetSeconds*fsInt*ones(size(kMatrixIdx));
kMatrix = reshape(stft(:,kMatrixIdx),nrFreqVars*winLen,[])';
nrCluster = 2;
tic
[IDX,C,sumd,D] = kmeans(kMatrix,nrCluster,'display','iter');
toc
% reassemble real and imag into complex number
clusters = zeros(freqLimit,winLen,nrCluster);
clusters = complex(clusters,clusters);
for clusterNr = 1:nrCluster
    clusterRealImag = reshape(C(clusterNr,:),[],2,winLen);
    clusters(:,:,clusterNr) = complex(squeeze(clusterRealImag(:,1,:)),squeeze(clusterRealImag(:,2,:)));
    figure, surf(abs(clusters(2:end,:,clusterNr))); view(-90,90);
    % figure, surf(angle(clusters(2:end,:,clusterNr))); view(-90,90);
end
figure, plot(D);
disp(['corr between distance series = ',num2str(corr(D(:,1),D(:,2)))]);
figure
tic
[silh2,h] = silhouette(kMatrix,IDX);
toc
set(get(gca,'Children'),'FaceColor',[.8 .8 1])

% smooth distance measures for downsampling
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.3/h.ud.tr,0.5/h.ud.tr,3,20,h.ud.fs);
tic
Hd = design(d);
toc
dist2clust1 = filtfilt(Hd.Numerator,1,D(:,1));
dist2clust2 = filtfilt(Hd.Numerator,1,D(:,2));
figure, plot([dist2clust1,dist2clust2]);
% downsample to BOLD sampling rate
dist2clust1_ds = resample(dist2clust1,size(boldTimeSeries,2),numel(eegLpf));
dist2clust2_ds = resample(dist2clust2,size(boldTimeSeries,2),numel(eegLpf));
meanBold = mean(boldTimeSeries,1)';
dist2clust1_ds_norm = detrend(dist2clust1_ds/std(dist2clust1_ds),'constant');
dist2clust2_ds_norm = detrend(dist2clust2_ds/std(dist2clust2_ds),'constant');
meanBold_norm = detrend(meanBold(1:numel(dist2clust1_ds))/std(meanBold(1:numel(dist2clust1_ds))),'constant');
figure, plot([dist2clust1_ds_norm,dist2clust2_ds_norm,meanBold_norm]);
disp(['corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold(1:numel(dist2clust1_ds)),[dist2clust1_ds,dist2clust2_ds]))]);
% convolve with HRF
hrf = mean(eventData,2);
dist2clust1_ds_conv = conv(dist2clust1_ds,hrf);
dist2clust2_ds_conv = conv(dist2clust2_ds,hrf);
dist2clust1_ds_conv = dist2clust1_ds_conv(1:numel(dist2clust1_ds));
dist2clust2_ds_conv = dist2clust2_ds_conv(1:numel(dist2clust2_ds));
dist2clust1_ds_conv_norm = detrend(dist2clust1_ds_conv/std(dist2clust1_ds_conv),'constant');
dist2clust2_ds_conv_norm = detrend(dist2clust2_ds_conv/std(dist2clust2_ds_conv),'constant');
figure, plot([dist2clust1_ds_conv_norm,dist2clust2_ds_conv_norm,meanBold_norm]);
disp(['corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold(1:numel(dist2clust1_ds_conv)),[dist2clust1_ds_conv,dist2clust2_ds_conv]))]);
clear kMatrix
% stft = repmat(stft,[1,1,size(sampEvents,2)]);
% for eventNr = 2:size(sampEvents,2)
%     stft(:,:,eventNr) = spectrogram(sampEvents(:,eventNr),hamming(128),120,128,256);
% end
% tfBaseline = mean(stft(:,floor((0)*numel(tfTime)):floor((0.25)*numel(tfTime)),:),2);
% tfPostStim = mean(stft(:,floor((0.25)*numel(tfTime)):floor((1)*numel(tfTime)),:),2);
% tfBaseline3d = repmat(tfBaseline,[1,numel(tfTime),1]);
% erdTfPsd = (stft - tfBaseline3d)./tfBaseline3d;
% avTfPsd = mean(erdTfPsd,3);
% surf(h.ax_taskEegTimeFreq,tfTime,tfFreq,avTfPsd,'edgecolor','none'); axis tight; 
% %surf(tfTime,tfFreq,avTfPsd,'edgecolor','none'); axis tight; 
% view(0,90);
% xlabel('Time (Seconds)'); ylabel('Hz');
% % post alpha ERS
% plot(h.ax_restBoldTimeseries,tfFreq,mean(squeeze(tfPostStim)./squeeze(tfBaseline),2));
% xlabel('Hz');
% % Scatterplot band power vs BOLD for Task and Rest

end
