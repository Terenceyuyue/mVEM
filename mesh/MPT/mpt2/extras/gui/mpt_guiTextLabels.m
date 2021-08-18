function varargout = mpt_guiTextLabels(varargin)
% MPT_GUITEXTLABELS M-file for mpt_guiTextLabels.fig
%      MPT_GUITEXTLABELS, by itself, creates a new MPT_GUITEXTLABELS or raises the existing
%      singleton*.
%
%      H = MPT_GUITEXTLABELS returns the handle to a new MPT_GUITEXTLABELS or the handle to
%      the existing singleton*.
%
%      MPT_GUITEXTLABELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUITEXTLABELS.M with the given input arguments.
%
%      MPT_GUITEXTLABELS('Property','Value',...) creates a new MPT_GUITEXTLABELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiTextLabels_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiTextLabels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to button_help mpt_guiTextLabels

% Last Modified by GUIDE v2.5 22-Feb-2005 13:37:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiTextLabels_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiTextLabels_OutputFcn, ...
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


% --- Executes just before mpt_guiTextLabels is made visible.
function mpt_guiTextLabels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiTextLabels (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiTextLabels
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiTextLabels wait for user response (see UIRESUME)


sysStruct = mpt_setSysStruct;

if ~isstruct(sysStruct),
    return
end

nx = mpt_setSysStruct('nx');
nu = mpt_setSysStruct('nu');
ny = mpt_setSysStruct('ny');

xnames = {};
if ~isfield(sysStruct, 'StateName'),
    for ii = 1:nx,
        xnames{end+1} = sprintf('x_%d', ii);
    end
    sysStruct.StateName = xnames;
end
unames = {};
if ~isfield(sysStruct, 'InputName'),
    for ii = 1:nu,
        unames{end+1} = sprintf('u_%d', ii);
    end
    sysStruct.InputName = unames;
end
ynames = {};
if ~isfield(sysStruct, 'OutputName'),
    for ii = 1:ny,
        ynames{end+1} = sprintf('y_%d', ii);
    end
    sysStruct.OutputName = ynames;
end
mpt_setSysStruct(sysStruct);

%uiwait(handles.TextLabels);

% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiTextLabels_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = [];


% --- Executes during object creation, after setting all properties.
function SelectState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

if ~isstruct(sysStruct),
    return
end

nx = mpt_setSysStruct('nx');

xstrings = '';
for ii=1:nx
    xstrings = strvcat(xstrings, sprintf('State %d',ii));
end
set(hObject, 'String', xstrings);
set(hObject, 'Value', 1);

% --- Executes on selection change in SelectState.
function SelectState_Callback(hObject, eventdata, handles)
% hObject    handle to SelectState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SelectState contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectState

sysStruct = mpt_setSysStruct;
index = get(hObject, 'Value');

string = sprintf('x_%d',index);

if isfield(sysStruct,'StateName')
    if iscell(sysStruct.StateName) & length(sysStruct.StateName)>=index,
        string = sysStruct.StateName{index};
    elseif ~iscell(sysStruct.StateName) & ischar(sysStruct.StateName),
        string = sysStruct.StateName;
    end
end

set(handles.StateName, 'String', string);

% --- Executes during object creation, after setting all properties.
function StateName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StateName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

string = '';
if isfield(sysStruct,'StateName')
    if iscell(sysStruct.StateName),
        string = sysStruct.StateName{1};
    else
        string = sysStruct.StateName;
    end
else
    string = 'x_1';
end

set(hObject, 'String', string);


function StateName_Callback(hObject, eventdata, handles)
% hObject    handle to StateName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StateName as text
%        str2double(get(hObject,'String')) returns contents of StateName as a double

index = get(handles.SelectState, 'Value');
mpt_setSysStruct(index, 'StateName', get(hObject, 'String'));

% --- Executes during object creation, after setting all properties.
function SelectInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

if ~isstruct(sysStruct),
    return
end

nu = mpt_setSysStruct('nu');
xstrings = '';
for ii=1:nu
    xstrings = strvcat(xstrings, sprintf('Input %d',ii));
end
set(hObject, 'String', xstrings);
set(hObject, 'Value', 1);


% --- Executes on selection change in SelectInput.
function SelectInput_Callback(hObject, eventdata, handles)
% hObject    handle to SelectInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SelectInput contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectInput

sysStruct = mpt_setSysStruct;
index = get(hObject, 'Value');

string = sprintf('u_%d',index);

if isfield(sysStruct,'InputName')
    if iscell(sysStruct.InputName) & length(sysStruct.InputName)>=index,
        string = sysStruct.InputName{index};
    elseif ~iscell(sysStruct.InputName) & ischar(sysStruct.InputName),
        string = sysStruct.InputName;
    end
end

set(handles.InputName, 'String', string);


% --- Executes during object creation, after setting all properties.
function InputName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

string = '';
if isfield(sysStruct,'InputName')
    if iscell(sysStruct.InputName),
        string = sysStruct.InputName{1};
    else
        string = sysStruct.InputName;
    end
else
    string = 'u_1';
end

set(hObject, 'String', string);



function InputName_Callback(hObject, eventdata, handles)
% hObject    handle to InputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputName as text
%        str2double(get(hObject,'String')) returns contents of InputName as a double

whichU = get(handles.SelectInput,'Value');
mpt_setSysStruct(whichU,'InputName',get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function SelectOutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

if ~isstruct(sysStruct),
    return
end

ny = mpt_setSysStruct('ny');

xstrings = '';
for ii=1:ny
    xstrings = strvcat(xstrings, sprintf('Output %d',ii));
end
set(hObject, 'String', xstrings);
set(hObject, 'Value', 1);


% --- Executes on selection change in SelectOutput.
function SelectOutput_Callback(hObject, eventdata, handles)
% hObject    handle to SelectOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SelectOutput contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectOutput

sysStruct = mpt_setSysStruct;
index = get(hObject, 'Value');

string = sprintf('y_%d',index);

if isfield(sysStruct,'OutputName')
    if iscell(sysStruct.OutputName) & length(sysStruct.OutputName)>=index,
        string = sysStruct.OutputName{index};
    elseif ~iscell(sysStruct.OutputName) & ischar(sysStruct.OutputName),
        string = sysStruct.OutputName;
    end
end

set(handles.OutputName, 'String', string);


% --- Executes during object creation, after setting all properties.
function OutputName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

sysStruct = mpt_setSysStruct;

string = '';
if isfield(sysStruct,'OutputName')
    if iscell(sysStruct.OutputName),
        string = sysStruct.OutputName{1};
    else
        string = sysStruct.OutputName;
    end
else
    string = 'y_1';
end

set(hObject, 'String', string);


function OutputName_Callback(hObject, eventdata, handles)
% hObject    handle to OutputName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputName as text
%        str2double(get(hObject,'String')) returns contents of OutputName as a double

whichU = get(handles.SelectOutput,'Value');
mpt_setSysStruct(whichU,'OutputName',get(hObject,'String'));

% --- Executes on button press in Return.
function Return_Callback(hObject, eventdata, handles)
% hObject    handle to Return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.TextLabels);

% --- Executes on button press in button_help.
function button_help_Callback(hObject, eventdata, handles)
% hObject    handle to button_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(['Here you can assign text labels to every state, input and output variable. These '...
        'labels will be used to denote axis in simulation plots.']);


%============================================================
function sub_helpdlg(msg)

msgbox(msg, 'Help', 'modal');