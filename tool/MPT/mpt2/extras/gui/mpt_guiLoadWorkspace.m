function varargout = mpt_guiLoadWorkspace(varargin)
% MPT_GUILOADWORKSPACE M-file for mpt_guiLoadWorkspace.fig
%      MPT_GUILOADWORKSPACE, by itself, creates a new MPT_GUILOADWORKSPACE or raises the existing
%      singleton*.
%
%      H = MPT_GUILOADWORKSPACE returns the handle to a new MPT_GUILOADWORKSPACE or the handle to
%      the existing singleton*.
%
%      MPT_GUILOADWORKSPACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUILOADWORKSPACE.M with the given input arguments.
%
%      MPT_GUILOADWORKSPACE('Property','Value',...) creates a new MPT_GUILOADWORKSPACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiLoadWorkspace_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiLoadWorkspace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiLoadWorkspace

% Last Modified by GUIDE v2.5 10-Feb-2005 13:47:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiLoadWorkspace_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiLoadWorkspace_OutputFcn, ...
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


% --- Executes just before mpt_guiLoadWorkspace is made visible.
function mpt_guiLoadWorkspace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiLoadWorkspace (see VARARGIN)

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiLoadWorkspace
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

% UIWAIT makes mpt_guiLoadWorkspace wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiLoadWorkspace_OutputFcn(hObject, eventdata, handles)
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
        