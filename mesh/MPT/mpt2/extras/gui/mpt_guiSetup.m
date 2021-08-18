function varargout = mpt_guiSetup(varargin)
% MPT_GUISETUP M-file for mpt_guiSetup.fig
%      MPT_GUISETUP, by itself, creates a new MPT_GUISETUP or raises the existing
%      singleton*.
%
%      H = MPT_GUISETUP returns the handle to a new MPT_GUISETUP or the handle to
%      the existing singleton*.
%
%      MPT_GUISETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUISETUP.M with the given input arguments.
%
%      MPT_GUISETUP('Property','Value',...) creates a new MPT_GUISETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiSetup_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiSetup

% Last Modified by GUIDE v2.5 23-Feb-2005 20:26:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiSetup_OutputFcn, ...
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


% --- Executes just before mpt_guiSetup is made visible.
function mpt_guiSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiSetup (see VARARGIN)

global mptOptions

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if ~isstruct(mptOptions),
    try
        mpt_init;
    catch
        errordlg(lasterr, 'Error', 'modal');
        delete(handles.figure1);
        return
    end
end

opt = mpt_options;

a_lp = opt.solvers.lp_all;
a_qp = opt.solvers.qp;
a_milp = opt.solvers.milp;
a_miqp = opt.solvers.miqp;
a_extr = opt.solvers.extreme;

if isempty(a_lp),
    LPstr = {'Not installed'};
else
    LPstr = {'Fastest available'};
end
for ii = 1:length(a_lp),
    LPstr{end+1} = mpt_solverInfo('lp', a_lp(ii));
end

if isempty(a_qp),
    QPstr = {'Not installed'};
else
    QPstr = {'Fastest available'};
end
for ii = 1:length(a_qp),
    QPstr{end+1} = mpt_solverInfo('qp', a_qp(ii));
end

if isempty(a_milp),
    MILPstr = {'Not available'};
else
    MILPstr = {'Fastest available'};
end
for ii = 1:length(a_milp),
    MILPstr{end+1} = mpt_solverInfo('milp', a_milp(ii));
end

if isempty(a_miqp),
    MIQPstr = {'Not available'};
else
    MIQPstr = {'Fastest available'};
end
for ii = 1:length(a_miqp),
    MIQPstr{end+1} = mpt_solverInfo('miqp', a_miqp(ii));
end

if isempty(a_extr)
    EXTRstr = {'Not available'};
else
    EXTRstr = {'Fastest available'};
end
for ii = 1:length(a_extr),
    EXTRstr{end+1} = mpt_solverInfo('extreme', a_extr(ii));
end
set(handles.popup_lpsolver, 'String', LPstr);
set(handles.popup_qpsolver, 'String', QPstr);
set(handles.popup_milpsolver, 'String', MILPstr);
set(handles.popup_miqpsolver, 'String', MIQPstr);
set(handles.popup_extremesolver, 'String', EXTRstr);
set(handles.popup_lpsolver, 'Value', 1);
set(handles.popup_qpsolver, 'Value', 1);
set(handles.popup_milpsolver, 'Value', 1);
set(handles.popup_miqpsolver, 'Value', 1);
set(handles.popup_extremesolver, 'Value', 1);


set(handles.popup_debuglevel, 'Value', mptOptions.debug_level + 1);
set(handles.textfield_infbox, 'String', sprintf('%d', mptOptions.infbox));
set(handles.textfield_stepsize, 'String', sprintf('%.3g', mptOptions.step_size));
set(handles.textfield_abstol, 'String', sprintf('%.3g', mptOptions.abs_tol));
set(handles.textfield_reltol, 'String', sprintf('%.3g', mptOptions.rel_tol));

set(handles.checkbox_checkupdates, 'Value', mptOptions.checkupdates);

set(handles.textfield_hysdelpath, 'String', mptOptions.hysdelpath);

%handles.output = [];
% UIWAIT makes mpt_guiSetup wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiSetup_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isstruct(handles)
    %varargout{1} = handles.output;
else
    %varargout{1} = [];
end
varargout = [];

% --- Executes during object creation, after setting all properties.
function popup_lpsolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_lpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_lpsolver.
function popup_lpsolver_Callback(hObject, eventdata, handles)
% hObject    handle to popup_lpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_lpsolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_lpsolver


% --- Executes during object creation, after setting all properties.
function popup_qpsolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_qpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_qpsolver.
function popup_qpsolver_Callback(hObject, eventdata, handles)
% hObject    handle to popup_qpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_qpsolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_qpsolver


% --- Executes during object creation, after setting all properties.
function popup_extremesolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_extremesolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_extremesolver.
function popup_extremesolver_Callback(hObject, eventdata, handles)
% hObject    handle to popup_extremesolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_extremesolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_extremesolver


% --- Executes during object creation, after setting all properties.
function popup_debuglevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_debuglevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_debuglevel.
function popup_debuglevel_Callback(hObject, eventdata, handles)
% hObject    handle to popup_debuglevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_debuglevel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_debuglevel


% --- Executes during object creation, after setting all properties.
function textfield_stepsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_stepsize_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_stepsize as text
%        str2double(get(hObject,'String')) returns contents of textfield_stepsize as a double


% --- Executes on button press in button_help_lp.
function button_help_lp_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_lp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Specifies the default Linear Programming (LP) solver.');

% --- Executes on button press in button_help_qp.
function button_help_qp_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_qp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Specifies the default Quadratic Programming (QP) solver.');

% --- Executes on button press in button_help_extreme.
function button_help_extreme_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_extreme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Specifies the default Vertex enumeration method.');

% --- Executes on button press in button_help_debuglevel.
function button_help_debuglevel_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_debuglevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = ['Due to numerical problems, tiny regions are sometimes difficult to' ...
    'calculate, i.e. are not identified at all. This may create "gaps"'...
    'in the computed control law. For the exploration, these will be' ...
    'jumped over and the exploration in the state space will continue.'];


sub_helpdlg('Specifies to which extent results should be double-checked.');

% --- Executes on button press in button_help_stepsize.
function button_help_stepsize_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Sets default value of a "step size", i.e. a length of a step over a facet in multi-parametric solvers.');

% --- Executes during object creation, after setting all properties.
function textfield_abstol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_abstol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_abstol_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_abstol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_abstol as text
%        str2double(get(hObject,'String')) returns contents of textfield_abstol as a double


% --- Executes on button press in button_help_abstol.
function button_help_abstol_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_abstol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Sets default value of the absolute tolerance. This number should be at least of one order smaller than the relative tolerance.');

% --- Executes during object creation, after setting all properties.
function textfield_reltol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_reltol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_reltol_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_reltol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_reltol as text
%        str2double(get(hObject,'String')) returns contents of textfield_reltol as a double


% --- Executes on button press in button_help_reltol.
function button_help_reltol_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_reltol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Sets default value of the relative tolerance. This number should be at least of one order bigger than the absolute tolerance.');

% --- Executes during object creation, after setting all properties.
function popup_milpsolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_milpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_milpsolver.
function popup_milpsolver_Callback(hObject, eventdata, handles)
% hObject    handle to popup_milpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_milpsolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_milpsolver


% --- Executes during object creation, after setting all properties.
function popup_miqpsolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_miqpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_miqpsolver.
function popup_miqpsolver_Callback(hObject, eventdata, handles)
% hObject    handle to popup_miqpsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_miqpsolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_miqpsolver


% --- Executes on button press in button_help_milp.
function button_help_milp_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_milp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Specifies the default Mixed-Integer Linear Programming (MILP) solver.');

% --- Executes on button press in button_help_miqp.
function button_help_miqp_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_miqp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Specifies the default Mixed-Integer Quadratic Programming (MIQP) solver.');

% --- Executes during object creation, after setting all properties.
function textfield_infbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_infbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_infbox_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_infbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_infbox as text
%        str2double(get(hObject,'String')) returns contents of textfield_infbox as a double


% --- Executes on button press in button_help_infbox.
function button_help_infbox_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_infbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = ['Sets default value of the "infinite box". '...
        'The the R^n polyhedra is internally converted to a box of this size.'...
        sprintf('\n\n') ...
        'ATTENTION: we have observed certain numerical problems when value of the '...
        'infbox is too big, therefore we recommend not to change the default setting '...
        'unless necessary.'];

sub_helpdlg(msg);

% --- Executes on button press in button_return.
function button_return_Callback(hObject, eventdata, handles)
% hObject    handle to button_return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mptOptions

try
    abs_tol = evalin('base', get(handles.textfield_abstol, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

try
    rel_tol = evalin('base', get(handles.textfield_reltol, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

try
    step_size = evalin('base', get(handles.textfield_stepsize, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

try
    infbox = evalin('base', get(handles.textfield_infbox, 'String'));
catch
    errordlg(lasterr, 'Error', 'modal');
    return
end

if isempty(infbox),
    errordlg('Size of the infinite box must be a double!', 'Error', 'modal');
    return
elseif any(size(infbox)~=1) | ~isa(infbox, 'double'),
    errordlg('Size of the infinite box must be a double!', 'Error', 'modal');
    return
elseif infbox<=0,
    errordlg('Size of the infinite box must be positive!', 'Error', 'modal');
elseif infbox<=100,
    button = questdlg('Size of infite box appears to be small, you should usually use large values (e.g. 10000)', ...
        'Warning', 'Change value', 'Keep value', 'Change value');
    if ~isempty(findstr(button, 'Change')),
        return
    end
end

if isempty(step_size),
    errordlg('Step size must be a double!', 'Error', 'modal');
    return
elseif any(size(step_size)~=1) | ~isa(step_size, 'double'),
    errordlg('Step size must be a double!', 'Error', 'modal');
    return
end

if isempty(abs_tol),
    errordlg('Absolute tolerance must be a double!', 'Error', 'modal');
    return
elseif any(size(abs_tol)~=1) | ~isa(abs_tol, 'double'),
    errordlg('Absolute tolerance must be a double!', 'Error', 'modal');
    return
end

if isempty(rel_tol),
    errordlg('Relative tolerance must be a double!', 'Error', 'modal');
    return
elseif any(size(rel_tol)~=1) | ~isa(rel_tol, 'double'),
    errordlg('Relative tolerance must be a double!', 'Error', 'modal');
    return
end

debug_level = get(handles.popup_debuglevel, 'Value') - 1;

contents = get(handles.popup_lpsolver ,'String');
str = contents{get(handles.popup_lpsolver,'Value')};
try
    if strcmpi(str, 'Not available'),
        lpsol = [];
    else
        lpsol = mpt_solverInfo('lp', str);
    end
catch
    sub_errordlg(lasterr);
    return
end

contents = get(handles.popup_qpsolver ,'String');
str = contents{get(handles.popup_qpsolver,'Value')};
try
    if strcmpi(str, 'Not available'),
        qpsol = [];
    else
        qpsol = mpt_solverInfo('qp', str);
    end
catch
    sub_errordlg(lasterr);
    return
end

contents = get(handles.popup_milpsolver ,'String');
str = contents{get(handles.popup_milpsolver,'Value')};
try
    if strcmpi(str, 'Not available'),
        milpsol = [];
    else
        milpsol = mpt_solverInfo('milp', str);
    end
catch
    sub_errordlg(lasterr);
    return
end

contents = get(handles.popup_miqpsolver ,'String');
str = contents{get(handles.popup_miqpsolver,'Value')};
try
    if strcmpi(str, 'Not available'),
        miqpsol = [];
    else
        miqpsol = mpt_solverInfo('miqp', str);
    end
catch
    sub_errordlg(lasterr);
    return
end

contents = get(handles.popup_extremesolver ,'String');
str = contents{get(handles.popup_extremesolver,'Value')};
try
    if strcmpi(str, 'Not available'),
        extremesol = [];
    else
        extremesol = mpt_solverInfo('extreme', str);
    end
catch
    sub_errordlg(lasterr);
    return
end

hysdelpath = get(handles.textfield_hysdelpath, 'String');

if isempty(hysdelpath)
    button = questdlg('HYSDEL binary was not located, are you sure you want to disable HYSDEL support?', ...
        'HYSDEL binary not found', 'Locate HYSDEL', 'Disable HYSDEL', 'Locate HYSDEL');
    if strcmpi(button, 'locate hysdel'),
        button_browsehysdel_Callback(handles.button_browsehysdel, eventdata, handles)
        return
    elseif isempty(button),
        return
    end
end
    
checkupdates = get(handles.checkbox_checkupdates, 'Value');

button = questdlg('Make changes permanent or for current session only?', ...
    'Pemanent changes', 'Permanent', 'Current Session', 'Permanent');

if strcmpi(button, 'permanent'),
    p = path;
    mptpath = pwd;
    while ~isempty(p),
        [t, p] = strtok(p, pathsep);
        mptpos = findstr(t, 'mpt');
        if ~isempty(mptpos),
            mptpath = t(1:mptpos(1)-1);
            mptpath = [mptpath 'mpt'];
            break
        end
    end
    if ~exist(mptpath, 'dir'),
        mptpath = pwd;
    end
    [ffname, dirname] = uigetfile('mpt_init.m', 'Locate mpt_init.m');
    if ispc,
        initfile = [dirname '\mpt_init.m'];
    else
        initfile = [dirname '/mpt_init.m'];
    end
    if ~exist(initfile, 'file'),
        errordlg('This directory does not contain "mpt_init.m" file!', 'Error', 'modal');
        return
    end
    tmpfile = tempname;
    infid = fopen(initfile, 'r');
    if infid<0,
        errordlg(sprintf('Couldn''t open "%s" for reading!', initfile), 'Error', 'modal');
        return
    end
    outfid = fopen(tmpfile, 'w');
    if outfid<0,
        fclose(infid);
        errordlg(sprintf('Couldn''t open "%s" for writing!', tmpfile), 'Error', 'modal');
        return
    end
    
    replacepart = 1;
    while 1,
        tline = fgetl(infid);
        if ~ischar(tline),
            break
        end
        if ~isempty(findstr(tline, 'DO NOT EDIT BEYOND THIS LINE')),
            replacepart = 0;
        end
        if replacepart,
            if ~isempty(findstr(tline, 'mptOptions.lpsolver')),
                tline = sprintf('mptOptions.lpsolver = %s;', mat2str(lpsol));
            elseif ~isempty(findstr(tline, 'mptOptions.qpsolver')),
                tline = sprintf('mptOptions.qpsolver = %s;', mat2str(qpsol));
            elseif ~isempty(findstr(tline, 'mptOptions.milpsolver')),
                tline = sprintf('mptOptions.milpsolver = %s;', mat2str(milpsol));
            elseif ~isempty(findstr(tline, 'mptOptions.miqpsolver')),
                tline = sprintf('mptOptions.miqpsolver = %s;', mat2str(miqpsol));
            elseif ~isempty(findstr(tline, 'mptOptions.extreme_solver')),
                tline = sprintf('mptOptions.extreme_solver = %s;', mat2str(extremesol));
            elseif ~isempty(findstr(tline, 'mptOptions.abs_tol')),
                tline = sprintf('mptOptions.abs_tol = %.3g;', abs_tol);
            elseif ~isempty(findstr(tline, 'mptOptions.rel_tol')),
                tline = sprintf('mptOptions.rel_tol = %.3g;', rel_tol);
            elseif ~isempty(findstr(tline, 'mptOptions.step_size')),
                tline = sprintf('mptOptions.step_size = %.3g;', step_size);
            elseif ~isempty(findstr(tline, 'mptOptions.debug_level')),
                tline = sprintf('mptOptions.debug_level = %d;', debug_level);
            elseif ~isempty(findstr(tline, 'mptOptions.infbox')),
                tline = sprintf('mptOptions.infbox = %d;', infbox);
            elseif ~isempty(findstr(tline, 'mptOptions.hysdelpath')),
                tline = sprintf('mptOptions.hysdelpath = ''%s'';', hysdelpath);
            elseif ~isempty(findstr(tline, 'mptOptions.checkupdates')),
                tline = sprintf('mptOptions.checkupdates = %d;', checkupdates);
            end
        end
        fprintf(outfid, '%s\n', tline);
    end
    
    fclose(infid);
    fclose(outfid);
    [succ, msg] = copyfile(tmpfile, initfile);
    if succ==0,
        errordlg(msg);
        return
    end
    try
        delete(tmpfile);
    end
end

if strcmpi(button, 'current session') | strcmpi(button, 'permanent'),

    if isempty(lpsol),
        lpsol = mptOptions.solvers.lp(1);
    end
    if isempty(qpsol),
        qpsol = mptOptions.solvers.qp(1);
    end
    if isempty(milpsol),
        milpsol = mptOptions.solvers.milp(1);
    end
    if isempty(miqpsol),
        miqpsol = mptOptions.solvers.miqp(1);
    end
    if isempty(extremesol),
        extremesol = mptOptions.solvers.extreme(1);
    end

    mptOptions.lpsolver = lpsol;
    mptOptions.qpsolver = qpsol;
    mptOptions.milpsolver = milpsol;
    mptOptions.miqpsolver = miqpsol;
    mptOptions.extreme_solver = extremesol;
    mptOptions.abs_tol = abs_tol;
    mptOptions.rel_tol = rel_tol;
    mptOptions.debug_level = debug_level;
    mptOptions.step_size = step_size;
    mptOptions.infbox = infbox;
    mptOptions.hysdelpath = hysdelpath;
    mptOptions.checkupdates = checkupdates;
    
else
    return
end

uiresume(handles.figure1);
delete(handles.figure1);

% --- Executes on button press in button_news.
function button_news_Callback(hObject, eventdata, handles)
% hObject    handle to button_news (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = [sprintf('* Graphical user interface. Type ''mpt_studio'' to start it\n\n'), ...
        sprintf(['* Solution to on-line MPC problems where the optimization problem is solved ' ...
            'at every time step. Supports linear systems, linear systems with boolean/integer ' ...
            'inputs and hybrid systems modeled with HYSDEL. See ''help mpt_control'' for more details.\n\n']), ...
        sprintf('* Improved Simulink interface and ability to export explicit controllers to standalone C-code\n\n'), ...
        sprintf(['* Models of dynamical system can now be imported either from HYSDEL or from ' ...
            'SS (state-space), TF (transfer-function), IDSS (system identification toolbox) and ' ...
            'MPC (MPC toolbox) objects. See ''help mpt_sys'' for more details.\n\n']), ...
        sprintf(['* mpt_control now returns an instance of MPTCTRL object, which is a newly ' ...
            'introduced class to replace the ctrlStruct structures. all functions now accept ' ...
            'both the "mptctrl" object, as well as the ctrlStruct structures as inputs, such that '...
            'backwards compatibility is preserved. See ''help mptctrl'' and ''help mpt_control'' for more details\n\n']), ...
        sprintf(['* mpt_lyapunov.m - interface function for computation of Lyapunov functions. ' ...
            'See ''help mpt_lyapunov'' for more details.\n\n']), ...
        sprintf(['* mpt_invariantSet - calculates an invariant subset of an explicit controller. ' ...
            'See ''help mpt_invariantSet'' for more details.\n\n']), ...
        sprintf(['* Very fast conversion from an MLD (mixed logical-dynamical) representation to ' ...
            'an equivalent PWA (piecewise-affine) form. See ''help mpt_sys'' for more details.\n\n']), ...
        sprintf('* State constraints can now be used. Just define them in sysStruct.xmin and sysStruct.xmax')
];

msgbox(msg, 'What''s new in MPT 2.0', 'modal');

% --- Executes during object creation, after setting all properties.
function textfield_hysdelpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_hysdelpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_hysdelpath_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_hysdelpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_hysdelpath as text
%        str2double(get(hObject,'String')) returns contents of textfield_hysdelpath as a double


% --- Executes on button press in button_browsehysdel.
function button_browsehysdel_Callback(hObject, eventdata, handles)
% hObject    handle to button_browsehysdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mptOptions

if ispc,
    [fname, fpath] = uigetfile('*.exe', 'Locate HYSDEL binary');
else
    [fname, fpath] = uigetfile('*.*', 'Locate HYSDEL binary');
end

if isa(fname, 'double')
    % action was canceled
else
    fullpath = [fpath fname];
    [d1, d2, ext] = fileparts(fullpath);
    if ispc & ~strcmpi(ext, '.exe'),
        errordlg('This is not a valid executable file!', 'Error', 'modal');
        return
    elseif strcmpi(ext, '.m') | strcmpi(ext, '.mat') | strcmpi(ext, '.fig')
        errordlg('This is not a valid binary!', 'Error', 'modal');
        return
    end
    try
        if ispc,
            T = evalc(sprintf('system(''"%s" -V'');', fullpath));
        else
            T = evalc(sprintf('system(''%s -V'');', fullpath));
        end
    catch
        errordlg(sprintf('Couldn''t execute "%s"!', fname), 'Error', 'modal');
        return
    end
    if isempty(findstr(T, 'HYSDEL')),
        errordlg(sprintf('"%s" is not a valid HYSDEL executable!', fname), 'Error', 'modal');
        return
    end
    mptOptions.hysdelpath = [fpath fname];
    set(handles.textfield_hysdelpath, 'String', fullpath);
end




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


% --- Executes on button press in checkbox_checkupdates.
function checkbox_checkupdates_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_checkupdates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_checkupdates


% --- Executes on button press in button_help_updates.
function button_help_updates_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_updates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = ['MPT can automatically check the product website '...
        '(http://control.ee.ethz.ch/~mpt/) for updates. If new version is available for '...
        'download, the user will be notified. The check is run every time mpt_init is '...
        'called (typically at the start of new Matlab session or manually, by call to '...
        '''mpt_update'')'];

sub_helpdlg(msg);

% --- Executes on button press in button_help_hysdel.
function button_help_hysdel_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_hysdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = ['HYSDEL (HYbrid Systems DEscription Language) is a tool for rapid modelling of '...
        'hybrid systems. If you want to use this feature in MPT, please provide full '...
        'path to HYSDEL executable.'];

sub_helpdlg(msg);


% --- Executes on button press in button_checkupdates.
function button_checkupdates_Callback(hObject, eventdata, handles)
% hObject    handle to button_checkupdates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


T = evalc('mpt_update');

msgbox(T, 'Update information', 'modal', 'modal');



%============================================================
function sub_helpdlg(msg)

msgbox(msg, 'Help', 'modal');