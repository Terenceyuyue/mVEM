function varargout = mpt_guiLoadWorkspace_export(varargin)
% MPT_GUILOADWORKSPACE_EXPORT M-file for mpt_guiLoadWorkspace_export.fig
%      MPT_GUILOADWORKSPACE_EXPORT, by itself, creates a new MPT_GUILOADWORKSPACE_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUILOADWORKSPACE_EXPORT returns the handle to a new MPT_GUILOADWORKSPACE_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUILOADWORKSPACE_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUILOADWORKSPACE_EXPORT.M with the given input arguments.
%
%      MPT_GUILOADWORKSPACE_EXPORT('Property','Value',...) creates a new MPT_GUILOADWORKSPACE_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiLoadWorkspace_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiLoadWorkspace_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiLoadWorkspace_export

% Last Modified by GUIDE v2.5 23-Feb-2005 21:52:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiLoadWorkspace_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiLoadWorkspace_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiLoadWorkspace_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiLoadWorkspace_export is made visible.
function mpt_guiLoadWorkspace_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiLoadWorkspace_export (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiLoadWorkspace_export
handles.output = [];

if length(varargin)==1,
    handles.structtype = varargin{1};
end

%handles.output = [];

% Update handles structure
guidata(hObject, handles);



[varstr, varnames] = listWorkspace(handles);

set(handles.listbox_variables, 'String', varstr);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes mpt_guiLoadWorkspace_export wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiLoadWorkspace_export_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if isempty(handles)
    varargout{1} = [];
    return
end
varargout{1} = handles.output;
delete(handles.figure1);

% --- Executes during object creation, after setting all properties.
function listbox_variables_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in listbox_variables.
function listbox_variables_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_variables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_variables contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_variables


% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
% hObject    handle to button_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

variables = get(handles.listbox_variables, 'String');
selected = get(handles.listbox_variables, 'Value');

if isempty(variables),
    return
end

variable = variables(selected, :);
if isfield(handles, 'structtype'),
    try
        varvalue = evalin('base', variable);
        handles.output = varvalue;
        guidata(hObject, handles);
        uiresume(handles.figure1);
        return
    catch
        return
    end
end

brackpos = findstr(variable, '(');
if ~isempty(brackpos),
    varname = variable(1:brackpos-1);
    varvalue = evalin('base', varname);
    if ~isstruct(varvalue) & ~isa(varvalue, 'ss') & ~isa(varvalue, 'tf') & ...
            ~isa(varvalue, 'idss') & ~isa(varvalue, 'mpc'),
        error(sprintf('Cannot load variable of class "%s"!', class(varvalue)));
    end
    try
        if isstruct(varvalue) & isfield(varvalue, 'A'),
        else
            T = evalc('varvalue = mpt_sys(varvalue);');
        end
    catch
        errordlg(sprintf('An error occured during import:\n%s', lasterr), 'Error', 'modal');
        return
    end
    mpt_setSysStruct(varvalue);
    handles.output = mpt_setSysStruct;
end
    
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


%=================================================
function [varstr, varnames] = listWorkspace(handles)

if nargin==0,
    handles = [];
end
structtype = 'sysstruct';
    
argset = 0;
if isstruct(handles),
    if isfield(handles, 'structtype'),
        structtype = handles.structtype;
        argset = 1;
    end
end

% obtain list of variables
S = evalin('base', 'whos');

varstr = '';
varnames = {};
if ~isempty(S),
    for ii = 1:length(S),
        varclass = S(ii).class;
        varname = S(ii).name;
        if argset==0 & (strcmpi(varclass, 'ss') | strcmpi(varclass, 'tf') | strcmpi(varclass, 'idss') | strcmpi(varclass, 'mpc')),
            % an ss, tf, idss or an mpc toolbox object
            varstr = strvcat(varstr, sprintf('%s (%s)', varname, varclass));
            varnames{end+1} = varname;
            
        elseif strcmpi(varclass, 'mptctrl') & argset==1,
            if strcmpi(structtype, 'mptctrl'),
                varstr = strvcat(varstr, sprintf('%s', varname));
                varnames{end+1} = varname;
            end
    
        elseif strcmpi(varclass, 'struct'),
            if argset==1,
                if strcmpi(structtype, 'mptctrl'),
                    continue
                end
            end
            dummy = evalin('base', varname);
            if strcmpi(structtype, 'probstruct'),
                if isfield(dummy, 'Q') & isfield(dummy, 'R') & isfield(dummy, 'N'),
                    % probStruct structure
                    varstr = strvcat(varstr, sprintf('%s', varname));
                    varnames{end+1} = varname;
                end
            else
                if mpt_issysstruct(dummy),
                    % sysStruct structure
                    if argset
                        varstr = strvcat(varstr, sprintf('%s', varname));
                    else
                        varstr = strvcat(varstr, sprintf('%s (%s)', varname, varclass));
                    end
                    varnames{end+1} = varname;
                end
            end
            
        end
    end
end
        

% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiLoadWorkspace_export_LayoutFcn(policy)
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
'Name','Load from workspace',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[103.8 40.9230769230769 53.8 20.5384615384615],...
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
'listbox', 2, ...
'pushbutton', 2), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiLoadWorkspace.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiLoadWorkspace_export(''listbox_variables_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[1.6 3.53846153846154 50 16.3076923076923],...
'Style','listbox',...
'Value',1,...
'CreateFcn','mpt_guiLoadWorkspace_export(''listbox_variables_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','listbox_variables');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiLoadWorkspace_export(''button_load_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[18.4 0.846153846153846 16.4 1.76923076923077],...
'String','Load',...
'Tag','button_load');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUILOADWORKSPACE_EXPORT, by itself, creates a new MPT_GUILOADWORKSPACE_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUILOADWORKSPACE_EXPORT returns the handle to a new MPT_GUILOADWORKSPACE_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUILOADWORKSPACE_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUILOADWORKSPACE_EXPORT.M with the given input arguments.
%
%      MPT_GUILOADWORKSPACE_EXPORT('Property','Value',...) creates a new MPT_GUILOADWORKSPACE_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.3 $ $Date: 2005/02/23 21:13:14 $

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
    % MPT_GUILOADWORKSPACE_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUILOADWORKSPACE_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUILOADWORKSPACE_EXPORT(...)
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

