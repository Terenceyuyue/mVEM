function varargout = mpt_guiWeights(varargin)
% MPT_GUIWEIGHTS M-file for mpt_guiWeights.fig
%      MPT_GUIWEIGHTS, by itself, creates a new MPT_GUIWEIGHTS or raises the existing
%      singleton*.
%
%      H = MPT_GUIWEIGHTS returns the handle to a new MPT_GUIWEIGHTS or the handle to
%      the existing singleton*.
%
%      MPT_GUIWEIGHTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIWEIGHTS.M with the given input arguments.
%
%      MPT_GUIWEIGHTS('Property','Value',...) creates a new MPT_GUIWEIGHTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiWeights_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiWeights_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiWeights

% Last Modified by GUIDE v2.5 22-Feb-2005 12:55:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiWeights_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiWeights_OutputFcn, ...
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


% --- Executes just before mpt_guiWeights is made visible.
function mpt_guiWeights_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiWeights (see VARARGIN)

global mpt__probStruct

mpt_guiCenter(hObject);

if ~isfield(mpt__probStruct, 'Q'),
    mpt__probStruct.Q = [];
end
if ~isfield(mpt__probStruct, 'R'),
    mpt__probStruct.R = [];
end

set(handles.textfield_Q, 'String', mat2str(mpt__probStruct.Q));
set(handles.textfield_R, 'String', mat2str(mpt__probStruct.R));
if isfield(mpt__probStruct, 'Qy'),
    set(handles.textfield_Qy, 'String', mat2str(mpt__probStruct.Qy));
end
if isfield(mpt__probStruct, 'Rdu'),
    set(handles.textfield_Rdu, 'String', mat2str(mpt__probStruct.Rdu));
end


% Choose default command line output for mpt_guiWeights
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiWeights wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiWeights_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function textfield_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_R_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_R as text
%        str2double(get(hObject,'String')) returns contents of textfield_R as a double

global mpt___sysStruct mpt__probStruct

[nx, nu, ny] = mpt_sysStructInfo(mpt___sysStruct);
value = get(hObject, 'String');
try,
    R = evalin('base', value);
    if ~isa(R, 'double'),
        R = str2num(R);
    end
    if mpt__probStruct.norm == 2,
        if size(R, 1) ~= size(R, 2) | size(R, 1) ~= nu,
            errordlg(sprintf('Penalty on inputs must be a %dx%d matrix!', nu, nu), 'Error', 'modal');
            return
        end
    else
        if size(R, 2) ~= nu,
            errordlg(sprintf('Penalty on outputs must have %d columns!', nu), 'Error', 'modal');
            return
        end
    end
catch
    sub_errordlg(lasterr);
    return
end
mpt__probStruct.R = R;

% --- Executes during object creation, after setting all properties.
function textfield_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_Q_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_Q as text
%        str2double(get(hObject,'String')) returns contents of textfield_Q as a double

global mpt___sysStruct mpt__probStruct

[nx, nu, ny] = mpt_sysStructInfo(mpt___sysStruct);
value = get(hObject, 'String');
try,
    Q = evalin('base', value);
    if ~isa(Q, 'double'),
        Q = str2num(Q);
    end
    if mpt__probStruct.norm == 2,
        if size(Q, 1) ~= size(Q, 2) | size(Q, 1) ~= nx,
            errordlg(sprintf('Penalty on states must be a %dx%d matrix!', nx, nx), 'Error', 'modal');
            return
        end
    else
        if size(Q, 2) ~= nx,
            errordlg(sprintf('Penalty on states must have %d columns!', nx), 'Error', 'modal');
            return
        end
    end
catch
    sub_errordlg(lasterr);
    return
end
mpt__probStruct.Q = Q;

% --- Executes during object creation, after setting all properties.
function textfield_Qy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_Qy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_Qy_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_Qy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_Qy as text
%        str2double(get(hObject,'String')) returns contents of textfield_Qy as a double

global mpt___sysStruct mpt__probStruct

[nx, nu, ny] = mpt_sysStructInfo(mpt___sysStruct);
value = get(hObject, 'String');
try,
    Qy = evalin('base', value);
    if ~isa(Qy, 'double'),
        Qy = str2num(Qy);
    end
    if mpt__probStruct.norm == 2,
        if size(Qy, 1) ~= size(Qy, 2) | size(Qy, 1) ~= ny,
            errordlg(sprintf('Penalty on outputs must be a %dx%d matrix!', ny, ny), 'Error', 'modal');
            return
        end
    else
        if size(Qy, 2) ~= ny,
            errordlg(sprintf('Penalty on outputs must have %d columns!', ny), 'Error', 'modal');
            return
        end
    end
catch
    sub_errordlg(lasterr);
    return
end
mpt__probStruct.Qy = Qy;

% --- Executes during object creation, after setting all properties.
function textfield_Rdu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_Rdu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_Rdu_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_Rdu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_Rdu as text
%        str2double(get(hObject,'String')) returns contents of textfield_Rdu as a double

global mpt___sysStruct mpt__probStruct

[nx, nu, ny] = mpt_sysStructInfo(mpt___sysStruct);
value = get(hObject, 'String');
try,
    R = evalin('base', value);
    if ~isa(R, 'double'),
        R = str2num(R);
    end
    if mpt__probStruct.norm == 2,
        if size(R, 1) ~= size(R, 2) | size(R, 1) ~= nu,
            errordlg(sprintf('Penalty on delta U must be a %dx%d matrix!', nu, nu), 'Error', 'modal');
            return
        end
    else
        if size(R, 2) ~= nu,
            errordlg(sprintf('Penalty on delta U must have %d columns!', nu), 'Error', 'modal');
            return
        end
    end
catch
    sub_errordlg(lasterr);
    return
end
mpt__probStruct.Rdu = R;


% --- Executes on button press in button_return.
function button_return_Callback(hObject, eventdata, handles)
% hObject    handle to button_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(handles.figure1);


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


% --- Executes on button press in button_help_Q.
function button_help_Q_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['Defined penalty on states Q in the cost function:\n\n' ...
        'J = sum( x(k)''*Q*x(k) + u(k)''*R*u(k))\n\n or\n\nJ = sum(||Q*x(k)|| + ||R*u(k)||)']));


% --- Executes on button press in button_help_U.
function button_help_U_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_U (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['Defined penalty on inputs R in the cost function:\n\n' ...
        'J = sum( x(k)''*Q*x(k) + u(k)''*R*u(k))\n\n or\n\nJ = sum(||Q*x(k)|| + ||R*u(k)||)']));

% --- Executes on button press in button_help_Qy.
function button_help_Qy_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_Qy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['Defined penalty on outputs Qy in the cost function:\n\n' ...
        'J = sum( x(k)''*Q*x(k) + u(k)''*R*u(k) + y(k)''*Qy*y(k))\n\n or\n\nJ = sum(||Q*x(k)|| + ||R*u(k)|| + ||Qy*y(k)||)']));


% --- Executes on button press in button_help_Rdu.
function button_help_Rdu_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_Rdu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


sub_helpdlg(sprintf(['Defined penalty on slew rate of U in the cost function:\n\n' ...
        'J = sum( x(k)''*Q*x(k) + (u(k)-u(k+1))''*Rdu*(u(k)-u(k+1)) )\n\n or \n\nJ = sum(||Q*x(k)|| + ||R*(u(k)-u(k+1))||)']));

%============================================================
function sub_helpdlg(msg)

msgbox(msg, 'Help', 'modal');