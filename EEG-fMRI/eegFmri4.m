function varargout = eegFmri4(varargin)
% mod of eegFmri3 for shift invariant TF maps

% EEGFMRI4 MATLAB code for eegFmri4.fig
%      EEGFMRI4, by itself, creates a new EEGFMRI4 or raises the existing
%      singleton*.
%
%      H = EEGFMRI4 returns the handle to a new EEGFMRI4 or the handle to
%      the existing singleton*.
%
%      EEGFMRI4('CALLBACK',hObj,eventData,h,...) calls the local
%      function named CALLBACK in EEGFMRI4.M with the given input arguments.
%
%      EEGFMRI4('Property','Value',...) creates a new EEGFMRI4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eegFmri4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eegFmri4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eegFmri4

% Last Modified by GUIDE v2.5 06-Feb-2014 22:33:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eegFmri4_OpeningFcn, ...
                   'gui_OutputFcn',  @eegFmri4_OutputFcn, ...
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
function varargout = eegFmri4_OutputFcn(hObj, eventdata, h) 
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

function eegLength_Callback(hObj, eventdata, h)
% hObj    handle to eegLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of eegLength as text
%        str2double(get(hObj,'String')) returns contents of eegLength as a double

h.ud.eegLength = str2double(get(hObj,'String'));
% Update h structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function eegLength_CreateFcn(hObj, eventdata, h)
% hObj    handle to eegLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end

function offset_Callback(hObj, eventdata, h)
% hObj    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of offset as text
%        str2double(get(hObj,'String')) returns contents of offset as a double

h.ud.offset = str2double(get(hObj,'String'));
% Update h structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObj, eventdata, h)
% hObj    handle to offset (see GCBO)
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

%%%%%%%%%%%%%
% INITIALIZE GUI
%%%%%%%%%%%%%%%%%%%
% --- Executes just before eegFmri4 is made visible.
function eegFmri4_OpeningFcn(hObj, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eegFmri4 (see VARARGIN)

% Choose default command line output for eegFmri4
h.output = hObj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize user data (ud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defaults
def.lowFreq = 0.01;
def.hiFreq = 0.1;
def.freqLimit = 10;
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
h.ud.offset = str2double(get(h.offset,'String'));
h.ud.eegLength = str2double(get(h.eegLength,'String'));
contents = cellstr(get(h.taskSelect,'String'));
h.ud.taskNr = get(h.taskSelect,'Value');
contents = cellstr(get(h.distanceMetric,'String'));
h.ud.distanceMetric = contents{get(h.distanceMetric,'Value')};
h.ud.freqLimit = str2double(get(h.freqLimit,'String'));
h.ud.winLen = str2double(get(h.winLen,'String'));
contents = cellstr(get(h.patternType,'String'));
h.ud.patternType = contents{get(h.patternType,'Value')};
%h.ud.roiNr = get(h.roiSelect,'Value');
%% update to gui
% set roi to left M1
h.ud.roiNr = 5;
set(h.roiSelect,'Value',h.ud.roiNr);
if isnan(h.ud.lowFreq)
    h.ud.lowFreq = def.lowFreq;
    set(h.lowFreq,'String',num2str(h.ud.lowFreq));
end
if isnan(h.ud.hiFreq)
    h.ud.hiFreq = def.hiFreq;
    set(h.hiFreq,'String',num2str(h.ud.hiFreq));
end
if isnan(h.ud.freqLimit)
    h.ud.freqLimit = def.freqLimit;
    set(h.freqLimit,'String',num2str(h.ud.freqLimit));
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
%h.ud.eeg = load(fullfile(h.ud.taskEegDir,'visualALL_MR_sr256_300s_CB_rej3_C3_ascii'),'-ascii');
h.ud.eeg = load(fullfile(h.ud.taskEegDir,'visualALL_MR_sr256_300s_CB_rej3_lapRefC3_ascii'),'-ascii');
toc
% Update handles structure
guidata(hObj, h);
% update GUI
updateGui(h);

% UIWAIT makes eegFmri4 wait for user response (see UIRESUME)
% uiwait(h.figure1);

end

function roi = getRoiInFunc(funcS,roiS)
funcVoxT = reshape(funcS.vol,[],funcS.nframes);
roi = funcVoxT(reshape(roiS.vol,[],1)>0,:);
end

%%%%%%%%%%%%%
%% UPDATE GUI
%%%%%%%%%%%%%%%%%%%%%%%%%
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
        boldTimeSeries = avSessData;
        eventData = sessEventData;
    else
        boldTimeSeries = cat(2,boldTimeSeries,avSessData);
        eventData = cat(2,eventData,sessEventData);
    end
end
% HPF BOLD
d = fdesign.highpass('Fst,Fp,Ast,Ap' ,0.01,0.02,10,3,1/h.ud.tr);
Hd = design(d);
boldTimeSeriesHpf = filtfilt(Hd.Numerator,1,boldTimeSeries);
plot(h.ax_taskBoldTimeseries,boldTimeSeries');
% average BOLD time series in roi
plot(h.ax_taskBoldHr,eventData);
hold(h.ax_taskBoldHr,'on');
plot(h.ax_taskBoldHr,mean(eventData,2),'Color','k','LineWidth',3);
hold(h.ax_taskBoldHr,'off');
%% EEG
% LPF EEG
d = fdesign.lowpass('Fp,Fst,Ap,Ast',50,60,3,20,h.ud.fs);
Hd = design(d);
eegLpf = filtfilt(Hd.Numerator,1,h.ud.eeg);
fsInt = 96;
eegLpf_ds = resample(eegLpf,fsInt,floor(h.ud.fs));
fs = h.ud.fs*fsInt/floor(h.ud.fs);
% EEG time-frequency maps
[stft,tfFreq,tfTime] = spectrogram(eegLpf_ds,hamming(fsInt),fsInt-1,fsInt,fs);
stft = stft(2:end,:);
%% limit size bec of memory
% limit frequencies
freqLimit = h.ud.freqLimit;
stft = stft(1:freqLimit,:);
% offset to avoid startup/filter transients
offsetSeconds = h.ud.offset;
% limit time
timeLimit = h.ud.eegLength; % length in seconds
stft = stft(:,1+floor(offsetSeconds*fs):floor(offsetSeconds*fs)+floor(timeLimit*fs));
% BOLD section for this time window
meanBold = mean(boldTimeSeries((1+floor(offsetSeconds/h.ud.tr)):(floor(offsetSeconds/h.ud.tr)+floor(timeLimit/h.ud.tr))),1)';
stftLen = size(stft,2);
winLen = floor(h.ud.winLen*fsInt);
kMatrixIdx = repmat(1:stftLen-winLen,[winLen,1]) + repmat((0:winLen-1)',[1,stftLen-winLen]);
kMatrix = stft(:,kMatrixIdx);
% shape to TF vs time
%kMatrixComplex = reshape(kMatrix,freqLimit,winLen,[]);
kMatrix = reshape(kMatrix,freqLimit,winLen,[]);
switch h.ud.patternType
    case 'abs'
        % look only at magnitude (discard phase info)
        kMatrix = abs(kMatrix);
        % normalize
        kMatrix = kMatrix./repmat(mean(mean(kMatrix,1),2),[freqLimit,winLen,1]);
    case 'angle'
        kMatrix = angle(kMatrix);
    case 'complex'
end
% make TF patterns shift-invariant
kMatrix = abs(fft2(kMatrix));
nrFreqVars = size(kMatrix,1);
kMatrix = reshape(kMatrix,nrFreqVars*winLen,[])';
nrCluster = 2;
tic
[IDX,C,sumd,D] = kmeans(kMatrix,nrCluster,'distance',h.ud.distanceMetric,'display','iter');
toc
%% convert distance D to activation
switch h.ud.distanceMetric
    case 'sqEuclidean'
        D = -1*D;
    case {'cosine','correlation'}
        D = 1 - D;
end
disp(['corr between distance series = ',num2str(corr(D(:,1),D(:,2)))]);
figure('WindowStyle','docked','NumberTitle','off','Name','D'), plot(D);

% smooth distance measures for downsampling
% % FFT LPF design
% fNyquistSamples = floor(0.5*(1/h.ud.tr)/(fs/size(D,1)));
% filterMask = [ones(fNyquistSamples,2);zeros(floor((size(D,1)-1)/2)-fNyquistSamples,2)];
% if rem(size(D,1),2)
%     % [0,0] = remove dc
%     filterMask = [[0,0];filterMask;flipud(filterMask)];
% else
%     filterMask = [[0,0];filterMask;[0,0];flipud(filterMask)];
% end
% D_lpf = ifft(fft(D).*filterMask,'symmetric');
% FIR LPF design
% remove dc
D = detrend(D,'constant');
% d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.4*(1/h.ud.tr),0.5*(1/h.ud.tr),3,10,fs);
% Hd = design(d);
% save('/NAS_II/Home/pboord/Documents/MATLAB/EEG-fMRI/Hd_D_lpf.mat','Hd');
load('/NAS_II/Home/pboord/Documents/MATLAB/EEG-fMRI/Hd_D_lpf.mat');
D_lpf = filtfilt(Hd.Numerator,1,D);
figure('WindowStyle','docked','NumberTitle','off','Name','D_lpf'), plot(D_lpf);
% downsample to BOLD sampling rate
dist2clust1_ds = resample(D_lpf(:,1),numel(meanBold),size(D_lpf,1));
dist2clust2_ds = resample(D_lpf(:,2),numel(meanBold),size(D_lpf,1));
dist2clust1_ds_norm = detrend(dist2clust1_ds,'constant')/std(dist2clust1_ds);
dist2clust2_ds_norm = detrend(dist2clust2_ds,'constant')/std(dist2clust2_ds);
meanBold_norm = detrend(meanBold,'constant')/std(meanBold);
figure('WindowStyle','docked','NumberTitle','off','Name','D_ds'), plot([dist2clust1_ds_norm,dist2clust2_ds_norm,meanBold_norm]);
disp(['no conv: corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold,[dist2clust1_ds,dist2clust2_ds]))]);
% convolve with HRF
hrf = mean(eventData,2);
dist2clust1_ds_conv = conv(dist2clust1_ds,hrf);
dist2clust2_ds_conv = conv(dist2clust2_ds,hrf);
dist2clust1_ds_conv = dist2clust1_ds_conv(1:numel(dist2clust1_ds));
dist2clust2_ds_conv = dist2clust2_ds_conv(1:numel(dist2clust2_ds));
dist2clust1_ds_conv_norm = detrend(dist2clust1_ds_conv,'constant')/std(dist2clust1_ds_conv);
dist2clust2_ds_conv_norm = detrend(dist2clust2_ds_conv,'constant')/std(dist2clust2_ds_conv);
figure('WindowStyle','docked','NumberTitle','off','Name','hrfAv'), plot([dist2clust1_ds_conv_norm,dist2clust2_ds_conv_norm,meanBold_norm]);
disp(['hrfAv: corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold,[dist2clust1_ds_conv,dist2clust2_ds_conv]))]);
% convolve with local HRF
hrf_local_events = eventData(:,(floor(h.ud.offset/h.ud.isi)+1):floor((h.ud.offset+h.ud.eegLength)/h.ud.isi));
hrf_local = mean(hrf_local_events,2);
dist2clust1_ds_conv = conv(dist2clust1_ds,hrf_local);
dist2clust2_ds_conv = conv(dist2clust2_ds,hrf_local);
dist2clust1_ds_conv = dist2clust1_ds_conv(1:numel(dist2clust1_ds));
dist2clust2_ds_conv = dist2clust2_ds_conv(1:numel(dist2clust2_ds));
dist2clust1_ds_conv_norm = detrend(dist2clust1_ds_conv,'constant')/std(dist2clust1_ds_conv);
dist2clust2_ds_conv_norm = detrend(dist2clust2_ds_conv,'constant')/std(dist2clust2_ds_conv);
figure('WindowStyle','docked','NumberTitle','off','Name','cf hrf')
plot(1:numel(hrf),hrf,'Color','k','LineWidth',3)
hold('on');
plot(hrf_local_events);
plot(1:numel(hrf),hrf_local,'Color','m','LineWidth',3);
hold('off');
figure('WindowStyle','docked','NumberTitle','off','Name','hrfLocal'), plot([dist2clust1_ds_conv_norm,dist2clust2_ds_conv_norm,meanBold_norm]);
disp(['hrfLocal: corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold,[dist2clust1_ds_conv,dist2clust2_ds_conv]))]);

%% deconvolve
hrfLen = 9;
% y = dist2d*hrfEst; hrfEst = dist2d\y
% B = A*x; Matlab version. i.e. x == hrfEst; Ax = B; x = A\B
distPad1 = [zeros(hrfLen-1,1);dist2clust1_ds;zeros(hrfLen-1,1)];
distPad2 = [zeros(hrfLen-1,1);dist2clust2_ds;zeros(hrfLen-1,1)];
nrRows = hrfLen + numel(dist2clust1_ds) - 1;
distIdx = fliplr(repmat(1:hrfLen,[nrRows,1]) + repmat((0:nrRows-1)',[1,hrfLen]));
dist2d1 = distPad1(distIdx);
dist2d2 = distPad2(distIdx);
y = [meanBold;zeros(nrRows-numel(meanBold),1)];
hrfEst1 = dist2d1\y;
hrfEst2 = dist2d2\y;
figure('WindowStyle','docked','NumberTitle','off','Name','hrfEst'),plot([hrfEst1,hrfEst2])
% convolve with estimated HRF
dist2clust1_ds_conv = conv(dist2clust1_ds,hrfEst1);
dist2clust2_ds_conv = conv(dist2clust2_ds,hrfEst2);
dist2clust1_ds_conv = dist2clust1_ds_conv(1:numel(dist2clust1_ds));
dist2clust2_ds_conv = dist2clust2_ds_conv(1:numel(dist2clust2_ds));
dist2clust1_ds_conv_norm = detrend(dist2clust1_ds_conv,'constant')/std(dist2clust1_ds_conv);
dist2clust2_ds_conv_norm = detrend(dist2clust2_ds_conv,'constant')/std(dist2clust2_ds_conv);
% dist2clust1_ds_conv_norm = dist2clust1_ds_conv;
% dist2clust2_ds_conv_norm = dist2clust2_ds_conv;
figure('WindowStyle','docked','NumberTitle','off','Name','estBOLD'), plot([dist2clust1_ds_conv_norm,dist2clust2_ds_conv_norm,meanBold_norm]);
disp(['hrfEst: corr between BOLD and cluster 1, cluster2: ',num2str(corr(meanBold,[dist2clust1_ds_conv,dist2clust2_ds_conv]))]);

% 
% clusters = zeros(freqLimit,winLen,nrCluster);
% clusters = complex(clusters,clusters);
% for clusterNr = 1:nrCluster
%     clusterRealImag = reshape(C(clusterNr,:),[],2,winLen);
%     clusters(:,:,clusterNr) = complex(squeeze(clusterRealImag(:,1,:)),squeeze(clusterRealImag(:,2,:)));
%     figure('WindowStyle','docked'), surf(abs(clusters(2:end,:,clusterNr))); view(-90,90);
%     % figure('WindowStyle','docked'), surf(angle(clusters(2:end,:,clusterNr))); view(-90,90);
% end
% 
% % reassemble real and imag into complex number
% clusters = zeros(freqLimit,winLen,nrCluster);
% clusters = complex(clusters,clusters);
% for clusterNr = 1:nrCluster
%     clusterRealImag = reshape(C(clusterNr,:),[],2,winLen);
%     clusters(:,:,clusterNr) = complex(squeeze(clusterRealImag(:,1,:)),squeeze(clusterRealImag(:,2,:)));
%     figure('WindowStyle','docked'), surf(abs(clusters(2:end,:,clusterNr))); view(-90,90);
%     % figure('WindowStyle','docked'), surf(angle(clusters(2:end,:,clusterNr))); view(-90,90);
% end
% figure('WindowStyle','docked')
% tic
% [silh2,h] = silhouette(kMatrix,IDX);
% toc
% set(get(gca,'Children'),'FaceColor',[.8 .8 1])

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


% --- Executes on selection change in distanceMetric.
function distanceMetric_Callback(hObj, eventdata, h)
% hObj    handle to distanceMetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObj,'String')) returns distanceMetric contents as cell array
%        contents{get(hObj,'Value')} returns selected item from distanceMetric

contents = cellstr(get(hObj,'String'));
h.ud.distanceMetric = contents{get(hObj,'Value')};
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function distanceMetric_CreateFcn(hObj, eventdata, h)
% hObj    handle to distanceMetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end



function freqLimit_Callback(hObj, eventdata, h)
% hObj    handle to freqLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of freqLimit as text
%        str2double(get(hObj,'String')) returns contents of freqLimit as a double

h.ud.freqLimit = str2double(get(hObj,'String'));
% Update h structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function freqLimit_CreateFcn(hObj, eventdata, h)
% hObj    handle to freqLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end



function winLen_Callback(hObj, eventdata, h)
% hObj    handle to winLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: get(hObj,'String') returns contents of winLen as text
%        str2double(get(hObj,'String')) returns contents of winLen as a double

h.ud.winLen = str2double(get(hObj,'String'));
% Update h structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function winLen_CreateFcn(hObj, eventdata, h)
% hObj    handle to winLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end

end


% --- Executes on selection change in patternType.
function patternType_Callback(hObj, eventdata, h)
% hObj    handle to patternType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObj,'String')) returns patternType contents as cell array
%        contents{get(hObj,'Value')} returns selected item from patternType
contents = cellstr(get(hObj,'String'));
h.ud.patternType = contents{get(hObj,'Value')};
% Update handles structure
guidata(h.figure1, h);
updateGui(h);
end

% --- Executes during object creation, after setting all properties.
function patternType_CreateFcn(hObj, eventdata, h)
% hObj    handle to patternType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObj,'BackgroundColor','white');
end
end
