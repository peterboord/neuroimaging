function varargout = eegFmri2(varargin)
% EEGFMRI2 MATLAB code for eegFmri2.fig
%      EEGFMRI2, by itself, creates a new EEGFMRI2 or raises the existing
%      singleton*.
%
%      H = EEGFMRI2 returns the handle to a new EEGFMRI2 or the handle to
%      the existing singleton*.
%
%      EEGFMRI2('CALLBACK',hObj,eventData,h,...) calls the local
%      function named CALLBACK in EEGFMRI2.M with the given input arguments.
%
%      EEGFMRI2('Property','Value',...) creates a new EEGFMRI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eegFmri2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eegFmri2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eegFmri2

% Last Modified by GUIDE v2.5 31-Jan-2014 16:25:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eegFmri2_OpeningFcn, ...
                   'gui_OutputFcn',  @eegFmri2_OutputFcn, ...
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
function varargout = eegFmri2_OutputFcn(hObj, eventdata, h) 
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

% --- Executes just before eegFmri2 is made visible.
function eegFmri2_OpeningFcn(hObj, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eegFmri2 (see VARARGIN)

% Choose default command line output for eegFmri2
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

% UIWAIT makes eegFmri2 wait for user response (see UIRESUME)
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
d = fdesign.highpass('N,Fc',floor(100*h.ud.fs),0.1,h.ud.fs);
Hd = design(d);
eegHpf = h.ud.eeg;%filter(Hd,h.ud.eeg);
nrEvents = nrSess*floor(h.ud.nrTr*h.ud.tr/h.ud.isi);
preStim = floor(-5*h.ud.fs):0;
postStim = 1:floor(2*h.ud.fs);
periSamples = [preStim,postStim];
eventIdxs = reshape(repmat(round((1:nrEvents-1)*h.ud.isi*h.ud.fs),numel(periSamples),1),1,[]);
eventIdxs = eventIdxs + repmat(periSamples,1,nrEvents-1);
sampEvents = reshape(eegHpf(eventIdxs),numel(periSamples),nrEvents-1);
sampEvents = detrend(sampEvents);
% ERP
plot(h.ax_taskBoldTimeseries,sampEvents);
hold(h.ax_taskBoldTimeseries,'on');
plot(h.ax_taskBoldTimeseries,mean(sampEvents,2),'Color','k','LineWidth',3);
hold(h.ax_taskBoldTimeseries,'off');
% Event-average EEG time-frequency maps
[~,tfFreq,tfTime,tfPsd] = spectrogram(sampEvents(:,1),hamming(128),120,128,256);
tfPsd = repmat(tfPsd,[1,1,size(sampEvents,2)]);
for eventNr = 2:size(sampEvents,2)
    [~,~,~,tfPsd(:,:,eventNr)] = spectrogram(sampEvents(:,eventNr),hamming(128),120,128,256);
end
tfPostStim = mean(tfPsd(:,floor((5.5/7)*numel(tfTime)):floor((6/7)*numel(tfTime)),:),2);
tfBaseline = mean(tfPsd(:,floor((4/7)*numel(tfTime)):floor((5/7)*numel(tfTime)),:),2);
tfBaseline3d = repmat(tfBaseline,[1,numel(tfTime),1]);
erdTfPsd = (tfPsd - tfBaseline3d)./tfBaseline3d;
avTfPsd = mean(erdTfPsd,3);
surf(h.ax_taskEegTimeFreq,tfTime,tfFreq,avTfPsd,'edgecolor','none'); axis tight; 
%surf(tfTime,tfFreq,avTfPsd,'edgecolor','none'); axis tight; 
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
% post alpha ERS
plot(h.ax_restBoldTimeseries,tfFreq,mean(squeeze(tfPostStim)./squeeze(tfBaseline),2));
xlabel('Hz');
% Scatterplot band power vs BOLD for Task and Rest

end

