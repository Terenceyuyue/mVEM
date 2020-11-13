function varargout = mpt_guiSimulate_export(varargin)
% MPT_GUISIMULATE_EXPORT M-file for mpt_guiSimulate_export.fig
%      MPT_GUISIMULATE_EXPORT, by itself, creates a new MPT_GUISIMULATE_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUISIMULATE_EXPORT returns the handle to a new MPT_GUISIMULATE_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUISIMULATE_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUISIMULATE_EXPORT.M with the given input arguments.
%
%      MPT_GUISIMULATE_EXPORT('Property','Value',...) creates a new MPT_GUISIMULATE_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiSimulate_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiSimulate_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiSimulate_export

% Last Modified by GUIDE v2.5 06-Apr-2005 11:54:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiSimulate_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiSimulate_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiSimulate_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiSimulate_export is made visible.
function mpt_guiSimulate_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiSimulate_export (see VARARGIN)

clear global mpt__sim_ctrl

global mpt__sim_ctrl

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiSimulate_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiSimulate_export wait for user response (see UIRESUME)
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
function varargout = mpt_guiSimulate_export_OutputFcn(hObject, eventdata, handles)
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

% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiSimulate_export_LayoutFcn(policy)
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
'Name','Simulator',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[-0.2 47.2307692307692 54.2 14.2307692307692],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'WindowStyle','modal',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', [], ...
'taginfo', struct(...
'figure', 2, ...
'popupmenu', 2, ...
'edit', 8, ...
'checkbox', 2, ...
'pushbutton', 2, ...
'radiobutton', 3, ...
'text', 4), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiSimulate.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiSimulate_export(''radio_closedloop_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[18.6 11.9230769230769 32.2 1.53846153846154],...
'String','Closed-loop trajectory',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_closedloop');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiSimulate_export(''radio_openloop_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[18.6 10.3846153846154 32.2 1.53846153846154],...
'String','Open-loop trajectory',...
'Style','radiobutton',...
'Tag','radio_openloop');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiSimulate_export(''textfield_x0_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[18.6 8.61538461538462 32.4 1.46153846153846],...
'String','[0;0]',...
'Style','edit',...
'CreateFcn','mpt_guiSimulate_export(''textfield_x0_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_x0');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiSimulate_export(''textfield_reference_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[18.6 6.92307692307692 32.4 1.46153846153846],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiSimulate_export(''textfield_reference_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_reference');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiSimulate_export(''textfield_horizon_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[18.6 5.07692307692308 32.4 1.46153846153846],...
'String','[]',...
'Style','edit',...
'CreateFcn','mpt_guiSimulate_export(''textfield_horizon_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_horizon');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiSimulate_export(''button_simulate_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[18 1.30769230769231 20.4 1.76923076923077],...
'String','Simulate',...
'Tag','button_simulate');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[1.8 8.69230769230769 15.4 1.30769230769231],...
'String','Initial state',...
'Style','text',...
'Tag','text1');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[1.8 7 15.4 1.30769230769231],...
'String','Reference point',...
'Style','text',...
'Tag','text_reference');


h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[1 5.23076923076923 16.2 1.15384615384615],...
'String','Simulation steps',...
'Style','text',...
'Tag','text3');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUISIMULATE_EXPORT, by itself, creates a new MPT_GUISIMULATE_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUISIMULATE_EXPORT returns the handle to a new MPT_GUISIMULATE_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUISIMULATE_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUISIMULATE_EXPORT.M with the given input arguments.
%
%      MPT_GUISIMULATE_EXPORT('Property','Value',...) creates a new MPT_GUISIMULATE_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.5 $ $Date: 2005/04/06 13:50:36 $

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
    % MPT_GUISIMULATE_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUISIMULATE_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUISIMULATE_EXPORT(...)
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

