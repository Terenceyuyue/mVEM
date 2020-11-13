function varargout = mpt_guiInputs_export(varargin)
% MPT_GUIINPUTS_EXPORT M-file for mpt_guiInputs_export.fig
%      MPT_GUIINPUTS_EXPORT, by itself, creates a new MPT_GUIINPUTS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIINPUTS_EXPORT returns the handle to a new MPT_GUIINPUTS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIINPUTS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIINPUTS_EXPORT.M with the given input arguments.
%
%      MPT_GUIINPUTS_EXPORT('Property','Value',...) creates a new MPT_GUIINPUTS_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiInputs_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiInputs_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiInputs_export

% Last Modified by GUIDE v2.5 28-Feb-2005 10:54:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiInputs_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiInputs_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiInputs_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiInputs_export is made visible.
function mpt_guiInputs_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiInputs_export (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiInputs_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

list_inputs_Callback(handles.list_inputs, eventdata, handles);

% UIWAIT makes mpt_guiInputs_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiInputs_export_OutputFcn(hObject, eventdata, handles)
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


% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiInputs_export_LayoutFcn(policy)
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
'Name','Input Types',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[79.8 47.923076923077 33.6 13.5384615384615],...
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
'radiobutton', 4, ...
'edit', 2, ...
'text', 2, ...
'pushbutton', 5), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiInputs.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiInputs_export(''list_inputs_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.6 10.7692307692308 30.2 1.76923076923077],...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiInputs_export(''list_inputs_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','list_inputs');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiInputs_export(''radio_continuous_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.6 8.76923076923077 30.2 1.61538461538462],...
'String','Continuous',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_continuous');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiInputs_export(''radio_binary_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.6 7.30769230769231 30.2 1.61538461538462],...
'String','Binary (1 / 0)',...
'Style','radiobutton',...
'Tag','radio_binary');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiInputs_export(''radio_alphabet_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[1.6 5.84615384615385 30.2 1.61538461538462],...
'String','From finate alphabet:',...
'Style','radiobutton',...
'Tag','radio_alphabet');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiInputs_export(''textfield_alphabet_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[1.6 3.69230769230769 30.2 1.76923076923077],...
'String','[]',...
'Style','edit',...
'CreateFcn','mpt_guiInputs_export(''textfield_alphabet_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_alphabet');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiInputs_export(''button_return_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[8.8 0.769230769230769 16.4 1.76923076923077],...
'String','Return',...
'Tag','button_return');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUIINPUTS_EXPORT, by itself, creates a new MPT_GUIINPUTS_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIINPUTS_EXPORT returns the handle to a new MPT_GUIINPUTS_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIINPUTS_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIINPUTS_EXPORT.M with the given input arguments.
%
%      MPT_GUIINPUTS_EXPORT('Property','Value',...) creates a new MPT_GUIINPUTS_EXPORT or raises the
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
    % MPT_GUIINPUTS_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUIINPUTS_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUIINPUTS_EXPORT(...)
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

