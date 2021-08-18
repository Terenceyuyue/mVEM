function varargout = mpt_guiWeights_export(varargin)
% MPT_GUIWEIGHTS_EXPORT M-file for mpt_guiWeights_export.fig
%      MPT_GUIWEIGHTS_EXPORT, by itself, creates a new MPT_GUIWEIGHTS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIWEIGHTS_EXPORT returns the handle to a new MPT_GUIWEIGHTS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIWEIGHTS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIWEIGHTS_EXPORT.M with the given input arguments.
%
%      MPT_GUIWEIGHTS_EXPORT('Property','Value',...) creates a new MPT_GUIWEIGHTS_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiWeights_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiWeights_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiWeights_export

% Last Modified by GUIDE v2.5 23-Feb-2005 21:54:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiWeights_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiWeights_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiWeights_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiWeights_export is made visible.
function mpt_guiWeights_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiWeights_export (see VARARGIN)

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


% Choose default command line output for mpt_guiWeights_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiWeights_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiWeights_export_OutputFcn(hObject, eventdata, handles)
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

% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiWeights_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

h1 = figure(...
'Units','characters',...
'Color',[0.847058823529412 0.847058823529412 0.847058823529412],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','Define Penalties',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[103.8 47.4615384615385 69 14],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', [], ...
'taginfo', struct(...
'figure', 2, ...
'text', 6, ...
'edit', 6, ...
'pushbutton', 8), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiWeights.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[1.8 9.15384615384616 18 1.15384615384615],...
'String','Penalty on inputs',...
'Style','text',...
'Tag','text1');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiWeights_export(''textfield_Q_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[22.4 11.3846153846154 38.4 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiWeights_export(''textfield_Q_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_Q');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiWeights_export(''textfield_R_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[22.4 8.92307692307692 38.4 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiWeights_export(''textfield_R_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_R');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiWeights_export(''textfield_Qy_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[22.4 6.53846153846154 38.4 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiWeights_export(''textfield_Qy_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_Qy');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[2.2 11.6153846153846 17.6 1.15384615384615],...
'String','Penalty on states',...
'Style','text',...
'Tag','text2');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[2 6.76923076923077 17.8 1.15384615384615],...
'String','Penalty on outputs',...
'Style','text',...
'Tag','text4');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[0.8 4.15384615384615 19 1.15384615384615],...
'String','Penalty on delta U',...
'Style','text',...
'Tag','text5');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiWeights_export(''textfield_Rdu_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[22.4 3.84615384615385 38.4 1.76923076923077],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiWeights_export(''textfield_Rdu_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_Rdu');


h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiWeights_export(''button_return_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[25.4 1 18.6 1.76923076923077],...
'String','Return',...
'Tag','button_return');


h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiWeights_export(''button_help_Q_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 11.3076923076923 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_Q');


h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiWeights_export(''button_help_U_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 8.84615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_U');


h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiWeights_export(''button_help_Qy_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 6.46153846153846 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_Qy');


h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiWeights_export(''button_help_Rdu_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 3.84615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_Rdu');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUIWEIGHTS_EXPORT, by itself, creates a new MPT_GUIWEIGHTS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIWEIGHTS_EXPORT returns the handle to a new MPT_GUIWEIGHTS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIWEIGHTS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIWEIGHTS_EXPORT.M with the given input arguments.
%
%      MPT_GUIWEIGHTS_EXPORT('Property','Value',...) creates a new MPT_GUIWEIGHTS_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2005/02/23 21:13:14 $

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [getfield(gui_State, gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % MPT_GUIWEIGHTS_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUIWEIGHTS_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUIWEIGHTS_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) & ischar(varargin{ind+1}) & ...
                strncmpi(varargin{ind},'visible',len1) & len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try, set(gui_hFigure, varargin{index}, varargin{index+1}), catch, break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)
try
    gui_hFigure = openfig(name, singleton, 'auto');
catch
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end

