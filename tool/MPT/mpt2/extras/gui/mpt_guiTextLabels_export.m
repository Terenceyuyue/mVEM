function varargout = mpt_guiTextLabels_export(varargin)
% MPT_GUITEXTLABELS_EXPORT M-file for mpt_guiTextLabels_export.fig
%      MPT_GUITEXTLABELS_EXPORT, by itself, creates a new MPT_GUITEXTLABELS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUITEXTLABELS_EXPORT returns the handle to a new MPT_GUITEXTLABELS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUITEXTLABELS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUITEXTLABELS_EXPORT.M with the given input arguments.
%
%      MPT_GUITEXTLABELS_EXPORT('Property','Value',...) creates a new MPT_GUITEXTLABELS_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiTextLabels_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiTextLabels_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to button_help mpt_guiTextLabels_export

% Last Modified by GUIDE v2.5 28-Feb-2005 10:56:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiTextLabels_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiTextLabels_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiTextLabels_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiTextLabels_export is made visible.
function mpt_guiTextLabels_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiTextLabels_export (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiTextLabels_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiTextLabels_export wait for user response (see UIRESUME)


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
function varargout = mpt_guiTextLabels_export_OutputFcn(hObject, eventdata, handles)
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

% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiTextLabels_export_LayoutFcn(policy)
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
'Name','Assign Text Labels',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[59.8 46.0769230769231 74.2 12.4615384615385],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','TextLabels',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', [], ...
'taginfo', struct(...
'figure', 2, ...
'popupmenu', 4, ...
'edit', 4, ...
'pushbutton', 6), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiTextLabels.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''SelectState_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.8 9.76923076923077 33.8 1.69230769230769],...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiTextLabels_export(''SelectState_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','SelectState');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''StateName_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[37.4 9.76923076923077 34 1.61538461538462],...
'String','x_1',...
'Style','edit',...
'CreateFcn','mpt_guiTextLabels_export(''StateName_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','StateName');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''SelectInput_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.8 7.53846153846154 33.8 1.69230769230769],...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiTextLabels_export(''SelectInput_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','SelectInput');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''InputName_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[37.4 7.53846153846154 34 1.61538461538462],...
'String','u_1',...
'Style','edit',...
'CreateFcn','mpt_guiTextLabels_export(''InputName_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','InputName');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''SelectOutput_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.8 5.30769230769231 33.8 1.69230769230769],...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiTextLabels_export(''SelectOutput_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','SelectOutput');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiTextLabels_export(''OutputName_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[37.4 5.30769230769231 34 1.61538461538462],...
'String','y_1',...
'Style','edit',...
'CreateFcn','mpt_guiTextLabels_export(''OutputName_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','OutputName');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiTextLabels_export(''Return_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[10.2 1.61538461538462 18.6 2.38461538461538],...
'String','Return',...
'Tag','Return');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiTextLabels_export(''button_help_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[45 1.61538461538462 18.6 2.38461538461538],...
'String','Help',...
'Tag','button_help');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUITEXTLABELS_EXPORT, by itself, creates a new MPT_GUITEXTLABELS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUITEXTLABELS_EXPORT returns the handle to a new MPT_GUITEXTLABELS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUITEXTLABELS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUITEXTLABELS_EXPORT.M with the given input arguments.
%
%      MPT_GUITEXTLABELS_EXPORT('Property','Value',...) creates a new MPT_GUITEXTLABELS_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.5 $ $Date: 2005/02/28 09:57:43 $

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
    % MPT_GUITEXTLABELS_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUITEXTLABELS_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUITEXTLABELS_EXPORT(...)
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

