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
% See also: GUIDE, GUIDATA, GUIh

% Edit the above text to modify the response to help editPeaks

% Last Modified by GUIDE v2.5 01-Jun-2015 15:04:26

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
h.zoom=zoom;
h.zoom.Enable = 'off';
h.physioFs=100; % Hz
h.physioDir='/projects2/act-plus/physio';
h.maxCardioInterval=0.67; % sec
h.minPeakProminence=1.5; % standard deviations
h.keepXlim=false;
h.fileName=[];
set(h.correctPk,'Value',1);
h.viewWindow=30; % sec
set(h.subjectIndicator,'String','');

guidata(hObject, h);

function loadData(hObject, eventdata, h)
set(h.subjectIndicator,'String',h.fileName);
set(h.messageText,'String','');
[peakTimes,peakLocs,normPhysio,sampleTimes]=getPhysioPeaks(fullfile(h.physioDir,[h.fileName,'.txt']),h.maxCardioInterval,h.minPeakProminence);
pkFilePath=fullfile(h.physioDir,[h.fileName,'_peaks.txt']);
badFilePath=fullfile(h.physioDir,[h.fileName,'_peaks_bad.txt']);
h.fixedPeakTimes=[];
if exist(pkFilePath,'file')
    h.fixedPeakTimes=load(pkFilePath,'-ascii');
end
if exist(badFilePath,'file')
    set(h.messageText,'String',['bad peaks file ',badFilePath,' exists']);
end
if isempty(normPhysio)
    disp('physio file is empty');
    badButton_Callback(hObject, eventdata, h);
    cla(h.axes1);
    cla(h.axes2);
    cla(h.axes3);
else
    if isempty(h.fixedPeakTimes)
        h.fixedPeakTimes=fixPhysioTimes(peakTimes,normPhysio,peakLocs,sampleTimes);
    elseif logical(get(h.correctPk,'Value'))
        for peakNr=1:numel(h.fixedPeakTimes)
            searchWidth=0.15; %sec
            minIdx=floor(h.physioFs*max([0,h.fixedPeakTimes(peakNr)-searchWidth])+1);
            maxIdx=ceil(h.physioFs*min([h.fixedPeakTimes(peakNr)+searchWidth,(numel(normPhysio)-1)/h.physioFs])+1);
            [~,pkIdx]=max(normPhysio(minIdx:maxIdx));
            h.fixedPeakTimes(peakNr)=(minIdx+pkIdx-2)/h.physioFs;
        end
    end
    h.fixedPeakLocs=round(h.fixedPeakTimes*h.physioFs+1);
    h.peakTimes=peakTimes;
    h.peakLocs=peakLocs;
    h.normPhysio=normPhysio;
    h.sampleTimes=sampleTimes;
%     plot(h.axes1,peakTimes(1:end-1),diff(peakTimes));
%     plot(h.axes2,h.sampleTimes,h.normPhysio,'b',h.peakTimes,h.normPhysio(h.peakLocs),'g*',h.fixedPeakTimes,h.normPhysio(h.fixedPeakLocs),'r*');
%     plot(h.axes3,h.fixedPeakTimes(1:end-1),diff(h.fixedPeakTimes));
    updateAxes(h,h.keepXlim);
end
guidata(h.figure1, h);

% --- Outputs from this function are returned to the command line.
function varargout = editPeaks_OutputFcn(hObject, eventdata, h) 
varargout{1} = h.output;

function assignCallbackFcns(h)
set(h.zoom,'ButtonDownFilter', @zoomButtonDownFilter)
set(h.axes2,'ButtonDownFcn', {@axes2_ButtonDownFcn,h})
set(get(h.axes2,'Children'),'ButtonDownFcn', {@axes2_ButtonDownFcn,h})
set(h.axes3,'ButtonDownFcn', {@axes3_ButtonDownFcn,h})
set(get(h.axes3,'Children'),'ButtonDownFcn', {@axes3_ButtonDownFcn,h})

function updateAxes(h,keepXlim)
set(h.medianPeakToPeak,'String',num2str(median(diff(h.fixedPeakTimes))));
plot(h.axes1,h.peakTimes(1:end-1),diff(h.peakTimes));
fixedPeakLocs=round(h.fixedPeakTimes*h.physioFs+1);
if islogical(keepXlim) && keepXlim
    ax2Xlim=xlim(h.axes2);
    ax2Ylim=ylim(h.axes2);
end
plot(h.axes2,h.sampleTimes,h.normPhysio,'b',h.peakTimes,h.normPhysio(h.peakLocs),'g*',h.fixedPeakTimes,h.normPhysio(fixedPeakLocs),'r*');
if isnumeric(keepXlim)
    xlim(h.axes2,[keepXlim-h.viewWindow,keepXlim+h.viewWindow]);
    ylim(h.axes2,[min(h.normPhysio),max(h.normPhysio)]);
elseif keepXlim
    xlim(h.axes2,ax2Xlim);
    ylim(h.axes2,ax2Ylim);
end
plot(h.axes3,sort(h.fixedPeakTimes(1:end-1)),diff(sort(h.fixedPeakTimes)));
if logical(get(h.respBox,'Value'))
    showRespPhase(h);
end
assignCallbackFcns(h);

function [flag] = zoomButtonDownFilter(hObject, eventdata)
% ref: help/zoom, Example 4 — Coding a ButtonDown Callback
switch get(gcf,'SelectionType')
    case {'alt','extend'}
        flag=true;
    otherwise
        flag=false;
end

function axes3_ButtonDownFcn(hObject, eventdata, h)
% ref: help/zoom, Example 4 — Coding a ButtonDown Callback
cP = get(gca,'Currentpoint');
t = cP(1,1);
switch get(gcf,'SelectionType')
    case 'alt'
        h.zoom.Enable='off';
end
updateAxes(h,t);

function axes2_ButtonDownFcn(hObject, eventdata, h)
% ref: help/zoom, Example 4 — Coding a ButtonDown Callback
cP = get(gca,'Currentpoint');
t = cP(1,1);
switch get(gcf,'SelectionType')
    case 'normal' % Click left mouse button.
        % add peak
        searchHalfWidth=0.1; % sec
        xMin=max([1,round((t-searchHalfWidth)*h.physioFs)+1]);
        xMax=round((t+searchHalfWidth)*h.physioFs)+1;
        [~,xPkMax]=max(h.normPhysio(xMin:xMax));
        h.fixedPeakTimes=[h.fixedPeakTimes;(xMin+xPkMax-2)/h.physioFs];
        updateAxes(h,true);
    case 'alt'    % Control - click left mouse button or click right mouse button.
        switch h.zoom.Enable
            case 'on'
                h.zoom.Enable='off';
            case 'off'
                h.zoom.Enable='on';
        end
     case 'extend' % Shift - click left mouse button or click both left and right mouse buttons.
        % delete nearest peak
        [~,minIdx]=min(abs(h.fixedPeakTimes-t));
        h.fixedPeakTimes(minIdx)=[];
        updateAxes(h,true);
    case 'open'   % Double-click any mouse button.
        h.fixedPeakTimes(end)=[]; % bec double-click does 'normal'+'open'
        updateAxes(h,false);
        h.zoom.Enable='on';
end
guidata(h.figure1, h);

% Save button
function saveButton_Callback(hObject, eventdata, h)
fixedPeakTimes=sort(h.fixedPeakTimes); %#ok<NASGU>
pkFilePath=fullfile(h.physioDir,[h.fileName,'_peaks.txt']);
save(pkFilePath,'fixedPeakTimes','-ascii');
set(h.messageText,'String',['peaks saved to ',pkFilePath]);

% Load button
function loadButton_Callback(hObject, eventdata, h)
[FileName,PathName]=uigetfile(fullfile(h.physioDir,'*.txt'));
if ischar(FileName)
    h.fileName=FileName(1:end-4);
    h.physioDir=PathName;
    loadData(hObject, eventdata, h);
    h.zoom.Enable = 'off';
end

% Bad button
function badButton_Callback(hObject, eventdata, h)
badPeaksPath=fullfile(h.physioDir,[h.fileName,'_bad.txt']);
pkFilePath=fullfile(h.physioDir,[h.fileName,'_peaks.txt']);
if exist(pkFilePath,'file')
    movefile(pkFilePath,badPeaksPath);
else
    fixedPeakTimes=sort(h.fixedPeakTimes); %#ok<NASGU>
    save(badPeaksPath,'fixedPeakTimes','-ascii');
    set(h.messageText,'String',['bad peaks saved to ',badPeaksPath]);
end

% Undo
function undoButton_Callback(hObject, eventdata, h)
h.fixedPeakTimes(end)=[];
updateAxes(h,true);

function maxCardioInterval_Callback(hObject, eventdata, h)
% hObject    handle to maxCardioInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxCardioInterval as text
%        str2double(get(hObject,'String')) returns contents of maxCardioInterval as a double
h.maxCardioInterval=str2double(get(hObject,'String'));
h.keepXlim=true;
guidata(h.figure1, h);
loadData(hObject, eventdata, h)

% --- Executes during object creation, after setting all properties.
function maxCardioInterval_CreateFcn(hObject, eventdata, h)
% hObject    handle to maxCardioInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',num2str(get(hObject,'Value')));

function minPeakProminence_Callback(hObject, eventdata, h)
% hObject    handle to minPeakProminence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minPeakProminence as text
%        str2double(get(hObject,'String')) returns contents of minPeakProminence as a double
h.minPeakProminence=str2double(get(hObject,'String'));
h.keepXlim=true;
guidata(h.figure1, h);
loadData(hObject, eventdata, h)

% --- Executes during object creation, after setting all properties.
function minPeakProminence_CreateFcn(hObject, eventdata, h)
% hObject    handle to minPeakProminence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',num2str(get(hObject,'Value')));


% --- Executes on button press in prevFile.
function prevFile_Callback(hObject, eventdata, h)
% hObject    handle to prevFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(h.fileName)
    splitExpr=regexp(h.fileName,'\d*','split');
    physioFiles=dir(fullfile(h.physioDir,[splitExpr{1},'*',splitExpr{2},'.txt']));
    physioFiles={physioFiles.name}';
    [~,fileNr]=ismember([h.fileName,'.txt'],physioFiles);
    if fileNr==1
        msgbox(['this is the first file of type ',splitExpr{1},'*',splitExpr{2},'.txt'],'modal');
    else
        [~,fileName]=fileparts(physioFiles{fileNr-1});
        h.fileName=fileName;
        guidata(h.figure1, h);
        loadData(hObject, eventdata, h);
        h.zoom.Enable = 'off';
    end
end

function nextFile_Callback(hObject, eventdata, h)
% hObject    handle to nextFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
if ~isempty(h.fileName)
    splitExpr=regexp(h.fileName,'\d*','split');
    physioFiles=dir(fullfile(h.physioDir,[splitExpr{1},'*',splitExpr{2},'.txt']));
    physioFiles={physioFiles.name}';
    [~,fileNr]=ismember([h.fileName,'.txt'],physioFiles);
    if fileNr==numel(physioFiles)
        msgbox(['no more files of type ',splitExpr{1},'*',splitExpr{2},'.txt'],'modal');
    else
        [~,fileName]=fileparts(physioFiles{fileNr+1});
        h.fileName=fileName;
        guidata(h.figure1, h);
        loadData(hObject, eventdata, h);
        h.zoom.Enable = 'off';
    end
end

% --- Executes on button press in correctPk.
function correctPk_Callback(hObject, eventdata, h)
% hObject    handle to correctPk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of correctPk

function showRespPhase(h)
r=h.normPhysio-min(h.normPhysio);
drDt=filtfilt(ones(h.physioFs,1)/h.physioFs,1,r);
rMax=max(r);
rHist=hist(r,rMax/100:rMax/100:rMax);
rPhase=zeros(numel(r),1);
for sampNr=1:numel(r)
    rPhase(sampNr)=sign(drDt(sampNr))*pi*sum(rHist(1:round(100*r(sampNr)./rMax)))/sum(rHist);
end
hold(h.axes2,'on');
plot(h.axes2,h.sampleTimes,-cos(rPhase));
hold(h.axes2,'off');

% --- Executes on button press in respBox.
function respBox_Callback(hObject, eventdata, h)
% hObject    handle to respBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of respBox
if ~isempty(h.subjectIndicator)
    updateAxes(h,true);
end

