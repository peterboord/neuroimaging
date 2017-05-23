function varargout = editPeaks(varargin)
% EDITPEAKS MATLAB code for editPeaks.fig
%      EDITPEAKS, by itself, creates a new EDITPEAKS or raises the existing
%      singleton*.
%
%      H = EDITPEAKS returns the handle to a new EDITPEAKS or the handle to
%      the existing singleton*.
%
%      EDITPEAKS('CALLBACK',hObject,eventData,h,...) calls the local
%      function named CALLBACK in EDITPEAKS.M with the given input arguments.
%
%      EDITPEAKS('Property','Value',...) creates a new EDITPEAKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before editPeaks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to editPeaks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help editPeaks

% Last Modified by GUIDE v2.5 28-May-2015 10:54:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editPeaks_OpeningFcn, ...
                   'gui_OutputFcn',  @editPeaks_OutputFcn, ...
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


% --- Executes just before editPeaks is made visible.
function editPeaks_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% varargin   command line arguments to editPeaks (see VARARGIN)

% Choose default command line output for editPeaks
h.output = hObject;

% Update h structure
guidata(hObject, h);

% UIWAIT makes editPeaks wait for user response (see UIRESUME)
% uiwait(h.figure1);
set(hObject,'Toolbar','figure');
physioDir='/projects2/act-plus/physio';
physioFiles=dir(fullfile(physioDir,'*_rest_cardio.txt'));
physioFiles={physioFiles.name}';
h.physioFs=100; % Hz
fileNr=3;
[~,fileName]=fileparts(physioFiles{fileNr});
set(h.text2,'String',fileName);
tuneParamPath=fullfile(physioDir,[fileName,'_tuneParam.mat']);
[peakTimes,peakLocs,normPhysio,sampleTimes]=getPhysioPeaks(fullfile(physioDir,physioFiles{fileNr}),tuneParamPath);
pkFilePath=fullfile(physioDir,[fileName,'_peaks.txt']);
h.pkFilePath=pkFilePath;
if exist(pkFilePath,'file')
    fixedPeakTimes=load(pkFilePath,'-ascii');
    fixedPeakLocs=round(fixedPeakTimes*h.physioFs+1);
else
    [fixedPeakTimes,fixedPeakLocs]=fixPhysioTimes(peakTimes,normPhysio,peakLocs,sampleTimes);
end
h.physioDir=physioDir;
h.physioFiles=physioFiles;
h.fileNr=fileNr;
h.fileName=fileName;
h.peakTimes=peakTimes;
h.peakLocs=peakLocs;
h.normPhysio=normPhysio;
h.sampleTimes=sampleTimes;
h.fixedPeakTimes=fixedPeakTimes;
h.fixedPeakLocs=fixedPeakLocs;
plot(h.axes1,peakTimes(1:end-1),diff(peakTimes));
plot(h.axes2,h.sampleTimes,h.normPhysio,'b',h.peakTimes,h.normPhysio(h.peakLocs),'g*',h.fixedPeakTimes,h.normPhysio(h.fixedPeakLocs),'r*');
plot(h.axes3,h.fixedPeakTimes(1:end-1),diff(h.fixedPeakTimes));
updateAxes(h)
guidata(hObject, h);

% --- Outputs from this function are returned to the command line.
function varargout = editPeaks_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Get default command line output from h structure
varargout{1} = h.output;

% --- Executes on button press in axes1.
function axes1_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
disp('axes1');

function assignAxesFcn(h)
set(h.axes1,'ButtonDownFcn', @axes1_Callback)
set(get(h.axes1,'Children'),'ButtonDownFcn', @axes1_Callback)
set(h.axes2,'ButtonDownFcn', {@axes2_Callback,h})
set(get(h.axes2,'Children'),'ButtonDownFcn', {@axes2_Callback,h})
set(h.axes3,'ButtonDownFcn', @axes3_Callback)
set(get(h.axes3,'Children'),'ButtonDownFcn', @axes3_Callback)

function updateAxes(h)
fixedPeakLocs=round(h.fixedPeakTimes*h.physioFs+1);
ax2Xlim=xlim(h.axes2);
ax2Ylim=ylim(h.axes2);
plot(h.axes2,h.sampleTimes,h.normPhysio,'b',h.peakTimes,h.normPhysio(h.peakLocs),'g*',h.fixedPeakTimes,h.normPhysio(fixedPeakLocs),'r*');
xlim(h.axes2,ax2Xlim);
ylim(h.axes2,ax2Ylim);
plot(h.axes3,h.fixedPeakTimes(1:end-1),diff(h.fixedPeakTimes));
assignAxesFcn(h);

% --- Executes on button press in axes2.
function axes2_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
cP = get(gca,'Currentpoint');
x = cP(1,1);
y = cP(1,2);
disp([x,y]);
switch get(gcf,'SelectionType')
    case 'normal' % Click left mouse button.
        % add peak
        disp('normal');
        h.fixedPeakTimes=sort([h.fixedPeakTimes;x]);
        updateAxes(h);
        %s = sprintf('left: (%1.4g, %1.4g) level = %1.4g',x,y, x.*exp(-x.^2-y.^2));
    case 'alt'    % Control - click left mouse button or click right mouse button.
        % delete peak
        disp('alt');
%         s = sprintf('right: (%1.4g, %1.4g level = %1.4g)',x,y, x.*exp(-x.^2-y.^2));
    case 'extend' % Shift - click left mouse button or click both left and right mouse buttons.
        disp('extend');
%         s = sprintf('2-click: (%1.4g, %1.4g level = %1.4g)',x,y, x.*exp(-x.^2-y.^2));
    case 'open'   % Double-click any mouse button.
        disp('open');
%         s = sprintf('double click: (%1.4g, %1.4g) level = %1.4g',x,y, x.*exp(-x.^2-y.^2));
end
disp('axes2');
guidata(h.figure1, h);

% --- Executes on button press in axes3.
function axes3_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
disp('axes3');

% --- Executes on button press in pushbutton1 'Save'.
function pushbutton1_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
disp('pushbutton1');
fixedPeakTimes=h.fixedPeakTimes; %#ok<NASGU>
save(h.pkFilePath,'fixedPeakTimes','-ascii');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
disp('pushbutton2');
