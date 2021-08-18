function varargout = mpt_guiInputs(varargin)
% MPT_GUIINPUTS M-file for mpt_guiInputs.fig
%      MPT_GUIINPUTS, by itself, creates a new MPT_GUIINPUTS or raises the existing
%      singleton*.
%
%      H = MPT_GUIINPUTS returns the handle to a new MPT_GUIINPUTS or the handle to
%      the existing singleton*.
%
%      MPT_GUIINPUTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIINPUTS.M with the given input arguments.
%
%      MPT_GUIINPUTS('Property','Value',...) creates a new MPT_GUIINPUTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiInputs_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiInputs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiInputs

% Last Modified by GUIDE v2.5 09-Feb-2005 15:03:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiInputs_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiInputs_OutputFcn, ...
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


% --- Executes just before mpt_guiInputs is made visible.
function mpt_guiInputs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiInputs (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiInputs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

list_inputs_Callback(handles.list_inputs, eventdata, handles);

% UIWAIT makes mpt_guiInputs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiInputs_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function list_inputs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

global mpt___sysStruct

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

nu = mpt_setSysStruct('nu');
if isempty(nu),
    nu = 2;
end
if isempty(nu)
    errordlg('Please defined state update equation first!', 'Error', 'modal');
    return
end
inp_str = '';
for ii = 1:nu,
    inp_str = strvcat(inp_str, sprintf('Input %d', ii));
end
set(hObject, 'String', inp_str);

if ~isfield(mpt___sysStruct, 'Uset')
    mpt___sysStruct.Uset = cell(1, nu);
    for iu = 1:nu,
        mpt___sysStruct.Uset{iu} = [-Inf Inf];
    end
end
           
            
% --- Executes on selection change in list_inputs.
function list_inputs_Callback(hObject, eventdata, handles)
% hObject    handle to list_inputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_inputs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_inputs

global mpt___sysStruct

curinput = get(hObject, 'Value');

set(handles.radio_continuous, 'Value', 1);
set(handles.radio_binary, 'Value', 0);
set(handles.radio_alphabet, 'Value', 0);

inptype = 1;
if isstruct(mpt___sysStruct)
    if isfield(mpt___sysStruct, 'Uset'),
        Uset = mpt___sysStruct.Uset;
        if ~iscell(Uset),
            Uset = {Uset};
        end
        if length(Uset) >= curinput,
            if all(isinf(Uset{curinput})),
                % first input is continuous
                inptype = 1;
            elseif any(Uset{curinput}~=ceil(Uset{curinput})) | any(isinf(Uset{curinput}))
                % finite alphabet input
                inptype = 3;
            else
                minU = min(Uset{curinput});
                maxU = max(Uset{curinput});
                if ~isempty(setdiff(minU:maxU,Uset{curinput}))
                    % rule out "gaps" in integers, e.g. [-2 -1 1 3] -> error
                    inptype = 3;
                else
                    inptype = 2;
                end
            end
        end
    end
end

if inptype == 1,
    % continuous input
    set(handles.radio_continuous, 'Value', 1);
    set(handles.radio_binary, 'Value', 0);
    set(handles.radio_alphabet, 'Value', 0);
    set(handles.textfield_alphabet, 'String', '[]');
    set(handles.textfield_alphabet, 'Enable', 'Off');
elseif inptype == 2,
    % boolean input
    set(handles.radio_continuous, 'Value', 0);
    set(handles.radio_binary, 'Value', 1);
    set(handles.radio_alphabet, 'Value', 0);
    set(handles.textfield_alphabet, 'String', '[]');
    set(handles.textfield_alphabet, 'Enable', 'Off');
else
    % finite alphabet input
    set(handles.radio_continuous, 'Value', 0);
    set(handles.radio_binary, 'Value', 0);
    set(handles.radio_alphabet, 'Value', 1);
    set(handles.textfield_alphabet, 'String', mat2str(Uset{curinput}));
    set(handles.textfield_alphabet, 'Enable', 'On');
end


% --- Executes on button press in radio_continuous.
function radio_continuous_Callback(hObject, eventdata, handles)
% hObject    handle to radio_continuous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_continuous

val = get(hObject, 'Value');
if val == 1,
    set(handles.radio_binary, 'Value', 0);
    set(handles.radio_alphabet, 'Value', 0);
    set(handles.textfield_alphabet, 'Enable', 'Off');
    curinput = get(handles.list_inputs, 'Value');
    mpt_setSysStruct(curinput, 'Uset', [-Inf Inf]);
end


% --- Executes on button press in radio_binary.
function radio_binary_Callback(hObject, eventdata, handles)
% hObject    handle to radio_binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_binary

val = get(hObject, 'Value');
if val == 1,
    set(handles.radio_continuous, 'Value', 0);
    set(handles.radio_alphabet, 'Value', 0);
    set(handles.textfield_alphabet, 'Enable', 'Off');
    curinput = get(handles.list_inputs, 'Value');
    mpt_setSysStruct(curinput, 'Uset', [0 1]);
end


% --- Executes on button press in radio_alphabet.
function radio_alphabet_Callback(hObject, eventdata, handles)
% hObject    handle to radio_alphabet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_alphabet


val = get(hObject, 'Value');
if val == 1,
    set(handles.radio_continuous, 'Value', 0);
    set(handles.radio_binary, 'Value', 0);
    set(handles.textfield_alphabet, 'Enable', 'On');
end

% --- Executes during object creation, after setting all properties.
function textfield_alphabet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_alphabet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_alphabet_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_alphabet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_alphabet as text
%        str2double(get(hObject,'String')) returns contents of textfield_alphabet as a double

curinput = get(handles.list_inputs, 'Value');
mpt_setSysStruct(curinput, 'Uset', get(hObject, 'String'));


% --- Executes on button press in button_return.
function button_return_Callback(hObject, eventdata, handles)
% hObject    handle to button_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);
