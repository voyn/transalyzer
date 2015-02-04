function varargout = GUI_transalyzer(varargin)
% GUI_TRANSALYZER M-file for GUI_transalyzer.fig
%      GUI_TRANSALYZER, by itself, creates a new GUI_TRANSALYZER or raises the existing
%      singleton*.
%
%      H = GUI_TRANSALYZER returns the handle to a new GUI_TRANSALYZER or the handle to
%      the existing singleton*.
%
%      GUI_TRANSALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRANSALYZER.M with the given input arguments.
%
%      GUI_TRANSALYZER('Property','Value',...) creates a new GUI_TRANSALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_transalyzer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_transalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_transalyzer

% Last Modified by GUIDE v2.5 12-Jan-2011 22:44:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_transalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_transalyzer_OutputFcn, ...
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


% --- Executes just before GUI_transalyzer is made visible.
function GUI_transalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_transalyzer (see VARARGIN)

% Choose default command line output for GUI_transalyzer
handles.output = hObject;

%addpath('K:\ns\mb\mb-shared\Calin\Scripts\Matlab_Transalyzer')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_transalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_transalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttonexit.
function buttonexit_Callback(hObject, eventdata, handles)
% hObject    handle to buttonexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)




% --- Executes on button press in buttonsort.
function buttonsort_Callback(hObject, eventdata, handles)
% hObject    handle to buttonsort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_events



% --- Executes on button press in buttondetect.
function buttondetect_Callback(hObject, eventdata, handles)
% hObject    handle to buttondetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_detect



% --- Executes on button press in buttonnoise.
function buttonnoise_Callback(hObject, eventdata, handles)
% hObject    handle to buttonnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_noise



% --- Executes on button press in buttonlocal.
function buttonlocal_Callback(hObject, eventdata, handles)
% hObject    handle to buttonlocal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_localstructures

