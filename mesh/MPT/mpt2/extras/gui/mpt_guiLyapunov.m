function varargout = mpt_guiLyapunov(varargin)
% MPT_GUILYAPUNOV M-file for mpt_guiLyapunov.fig
%      MPT_GUILYAPUNOV, by itself, creates a new MPT_GUILYAPUNOV or raises the existing
%      singleton*.
%
%      H = MPT_GUILYAPUNOV returns the handle to a new MPT_GUILYAPUNOV or the handle to
%      the existing singleton*.
%
%      MPT_GUILYAPUNOV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUILYAPUNOV.M with the given input arguments.
%
%      MPT_GUILYAPUNOV('Property','Value',...) creates a new MPT_GUILYAPUNOV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiLyapunov_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiLyapunov_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiLyapunov

% Last Modified by GUIDE v2.5 20-Feb-2005 20:52:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiLyapunov_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiLyapunov_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mpt_guiLyapunov is made visible.
function mpt_guiLyapunov_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiLyapunov (see VARARGIN)

global mpt__ctrl

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiLyapunov
handles.output = 'pwq';

% Update handles structure
guidata(hObject, handles);

if isa(mpt__ctrl, 'mptctrl'),
    norm = mpt__ctrl.probStruct.norm;
    if norm==2,
        set(handles.radio_pwa, 'Enable', 'off');
    end
end

% UIWAIT makes mpt_guiLyapunov wait for user response (see UIRESUME)
uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiLyapunov_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles),
    varargout{1} = 'pwq';
else
    varargout{1} = handles.output;
    % The figure can be deleted now
    delete(handles.figure1);
end


% --- Executes on button press in radio_pwq.
function radio_pwq_Callback(hObject, eventdata, handles)
% hObject    handle to radio_pwq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_pwq

set(handles.radio_pwa, 'Value', 0);
set(handles.radio_common, 'Value', 0);

% --- Executes on button press in radio_pwa.
function radio_pwa_Callback(hObject, eventdata, handles)
% hObject    handle to radio_pwa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_pwa

set(handles.radio_pwq, 'Value', 0);
set(handles.radio_common, 'Value', 0);


% --- Executes on button press in radio_common.
function radio_common_Callback(hObject, eventdata, handles)
% hObject    handle to radio_common (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_common

set(handles.radio_pwq, 'Value', 0);
set(handles.radio_pwa, 'Value', 0);



% --- Executes on button press in button_return.
function button_return_Callback(hObject, eventdata, handles)
% hObject    handle to button_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radio_pwq, 'Value'),
    handles.output = 'pwq';
elseif get(handles.radio_pwa, 'Value'),
    handles.output = 'pwa';
else
    handles.output = 'quad';
end

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
