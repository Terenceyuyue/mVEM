function varargout = mpt_guiSimulate(varargin)
% MPT_GUISIMULATE M-file for mpt_guiSimulate.fig
%      MPT_GUISIMULATE, by itself, creates a new MPT_GUISIMULATE or raises the existing
%      singleton*.
%
%      H = MPT_GUISIMULATE returns the handle to a new MPT_GUISIMULATE or the handle to
%      the existing singleton*.
%
%      MPT_GUISIMULATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUISIMULATE.M with the given input arguments.
%
%      MPT_GUISIMULATE('Property','Value',...) creates a new MPT_GUISIMULATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiSimulate_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiSimulate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiSimulate

% Last Modified by GUIDE v2.5 19-Feb-2005 23:17:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiSimulate_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiSimulate_OutputFcn, ...
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


% --- Executes just before mpt_guiSimulate is made visible.
function mpt_guiSimulate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiSimulate (see VARARGIN)

clear global mpt__sim_ctrl

global mpt__sim_ctrl

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiSimulate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiSimulate wait for user response (see UIRESUME)
% uiwait(handles.figure1);

mpt__sim_ctrl = [];
if length(varargin)>0,
    ctrl = varargin{1};
    
    if isa(ctrl, 'mptctrl'),
        if ctrl.probStruct.tracking==0,
            set(handles.textfield_reference, 'Enable', 'off');
            set(handles.text_reference, 'Enable', 'off');
            if isfield(ctrl.probStruct, 'yref'),
                set(handles.textfield_reference, 'String', mat2str(ctrl.probStruct.yref));
            end
        end
        mpt__sim_ctrl = ctrl;
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiSimulate_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function textfield_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

global mpt__sim_options

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

nx = mpt_setSysStruct('nx');
if nx>1,
    set(hObject, 'String', mat2str(zeros(nx,1)));
    mpt__sim_options.x0 = zeros(nx,1);
end


function textfield_x0_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_x0 as text
%        str2double(get(hObject,'String')) returns contents of textfield_x0 as a double

global mpt__sim_options

try
    value = evalin('base', get(hObject, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

mpt__sim_options.x0 = value;


% --- Executes on button press in button_simulate.
function button_simulate_Callback(hObject, eventdata, handles)
% hObject    handle to button_simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__sim_options mpt__sim_ctrl

opt = mpt__sim_options;
ctrl = mpt__sim_ctrl;

if isempty(ctrl)
    errordlg('No controller loaded!', 'Error', 'modal');
    return
elseif ~isa(ctrl, 'mptctrl'),
    errordlg('Unknown controller type!', 'Error', 'modal');
end

Options = [];
N = [];
x0 = [];

if ~isempty(opt),
    if isfield(opt, 'N'),
        N = opt.N;
    end
    if isfield(opt, 'x0')
        x0 = opt.x0;
    end
    if isfield(opt, 'reference'),
        Options.reference = opt.reference;
    end
end

if isempty(x0),
    errordlg('Initial state must be given!', 'Error', 'modal');
    return
end
    
Options.newfigure = 1;
Options.guierrors = 1;
Options.openloop = get(handles.radio_openloop, 'Value');
x0 = x0(:);
try
    T = evalc('mpt_plotTimeTrajectory(ctrl, x0, N, Options);');
catch
    sub_errordlg(lasterr);
end

% --- Executes on button press in radio_closedloop.
function radio_closedloop_Callback(hObject, eventdata, handles)
% hObject    handle to radio_closedloop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_closedloop

global mpt__sim_options
mpt__sim_options.openloop = 0;

set(handles.radio_openloop, 'Value', 0);


% --- Executes on button press in radio_openloop.
function radio_openloop_Callback(hObject, eventdata, handles)
% hObject    handle to radio_openloop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_openloop

global mpt__sim_options
mpt__sim_options.openloop = 1;

set(handles.radio_closedloop, 'Value', 0);

% --- Executes during object creation, after setting all properties.
function textfield_reference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_reference_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_reference as text
%        str2double(get(hObject,'String')) returns contents of textfield_reference as a double

global mpt__sim_options

try
    value = evalin('base', get(hObject, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

mpt__sim_options.reference = value;

% --- Executes during object creation, after setting all properties.
function textfield_horizon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_horizon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_horizon_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_horizon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_horizon as text
%        str2double(get(hObject,'String')) returns contents of textfield_horizon as a double


global mpt__sim_options

try
    value = evalin('base', get(hObject, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

if ~isa(value, 'double')
    errordlg('Number of simulation steps must be an integer!', 'Error', 'modal');
end
if ~isempty(value),
    if value ~= round(value),
        errordlg('Number of simulation steps must be an integer!', 'Error', 'modal');
    end
    if value<1,
        errordlg('Number of simulation steps must be greater zero!', 'Error', 'modal');
    end
end
mpt__sim_options.N = value;



%=============================================================
function sub_errordlg(msg)

origmsg = msg;

try
    newmsg = msg;
    errorpos = findstr(msg, 'Error using ==>');
    if ~isempty(errorpos),
        newmsg = msg(length('Error using ==>')+2:end);
        newlinepos = findstr(newmsg, 10);
        spacepos = findstr(newmsg, ' ');
        if ~isempty(newlinepos) & ~isempty(spacepos),
            breakpos = min(newlinepos(1), spacepos(1));
            newmsg = newmsg(breakpos+1:end);
        elseif isempty(newlinepos),
            newmsg = newmsg(spacepos(1)+1:end);
        elseif isempty(spacepos),
            newmsg = newmsg(newlinepos(1)+1:end);
        end
    end
    errordlg(newmsg, 'Error', 'modal');
catch
    errordlg(origmsg, 'Error', 'modal');
end