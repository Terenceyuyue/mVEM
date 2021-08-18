function varargout = mpt_guiProbStruct_export(varargin)
% MPT_GUIPROBSTRUCT_EXPORT M-file for mpt_guiProbStruct_export.fig
%      MPT_GUIPROBSTRUCT_EXPORT, by itself, creates adyn new MPT_GUIPROBSTRUCT_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIPROBSTRUCT_EXPORT returns the handle to adyn new MPT_GUIPROBSTRUCT_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIPROBSTRUCT_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIPROBSTRUCT_EXPORT.M with the given input arguments.
%
%      MPT_GUIPROBSTRUCT_EXPORT('Property','Value',...) creates adyn new MPT_GUIPROBSTRUCT_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiProbStruct_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiProbStruct_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiProbStruct_export

% Last Modified by GUIDE v2.5 05-Sep-2005 18:23:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiProbStruct_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiProbStruct_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiProbStruct_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiProbStruct_export is made visible.
function mpt_guiProbStruct_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiProbStruct_export (see VARARGIN)

clear global mpt___sysStruct mpt__probStruct mpt__ctrl
global mpt__probStruct mpt___sysStruct mpt__ctrl mptOptions

% resolution = get(0, 'ScreenSize');
% res_width = resolution(3);
% res_height = resolution(4);
% set(hObject, 'Units', 'pixels');
% position = get(hObject, 'Position');
% height = position(3);
% width = position(4);
% xpos = round(res_width/2) - round(width/2);
% ypos = round(res_height/2) - round(width/2);
% set(hObject, 'Position', [xpos ypos position(3:4)]);

mpt_guiCenter(hObject);

if ~isstruct(mptOptions)
    try
        mpt_init;
    catch
        sub_errordlg(lasterr);
        delete(handles.figure1);
        return
    end
end
    

% Choose default command line output for mpt_guiProbStruct_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiProbStruct_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

mpt__probStruct.norm = 1;
Options.verbose = 0;
Options.guierrors = 1;
Options.forceverify = 1;

try
    sysStruct = evalin('base', 'sysStruct');
    T = evalc('sysStruct = mpt_verifySysStruct(sysStruct, Options);');
    button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles, sysStruct)
    mpt___sysStruct = sysStruct;
end
try
    probStruct = evalin('base', 'probStruct');
    T = evalc('probStruct = mpt_verifyProbStruct(probStruct, Options);');
    mpt__probStruct = probStruct;
    sub_setProbStructFields(handles, probStruct);
end

mpt__ctrl = [];


%=======================================================================
function res = CeqEye

global mpt___sysStruct

if isstruct(mpt___sysStruct)
    [nx,nu,ny,ndyn] = mpt_sysStructInfo(mpt___sysStruct);
    if ~iscell(mpt___sysStruct.C),
        res = isequal(mpt___sysStruct.C, eye(ny));
    else
        res = 1;
        for ii = 1:ndyn,
            % check if C matrix is square [issue128]
            res = isequal(mpt___sysStruct.C{ii}, eye(ny));
            if res == 0
                break
            end
        end
    end
else
    res = 1;
end


%=======================================================================
function sub_setProbStructFields(handles, probStruct);


sub_resetPopup('objective', handles);
if ~CeqEye,
    set(handles.popup_objective, 'String', {'Regulation towards origin', 'Regulation to output reference', ...
        'Free output trajectory tracking'});
end

if probStruct.subopt_lev==0 & isinf(probStruct.N),
    set(handles.popup_solutiontype, 'Value', 2);
    set(handles.textfield_N, 'Enable', 'off');
    set(handles.text_N, 'Enable', 'off');
    set(handles.popup_objective, 'String', {'Regulation towards origin'});
    set(handles.popup_objective, 'Value', 1);
elseif probStruct.subopt_lev==0 & ~isinf(probStruct.N),
    set(handles.popup_solutiontype, 'Value', 1);
elseif probStruct.subopt_lev==1,
    set(handles.popup_solutiontype, 'Value', 3);
    set(handles.textfield_N, 'Enable', 'off');
    set(handles.text_N, 'Enable', 'off');
    set(handles.popup_objective, 'String', {'Regulation towards origin'});
    set(handles.popup_objective, 'Value', 1);
    sub_resetPopup('norm', handles)
else
    set(handles.popup_solutiontype, 'Value', 4);
    set(handles.popup_objective, 'String', {'Regulation towards origin'});
    set(handles.popup_objective, 'Value', 1);
end
if probStruct.norm==2,
    sub_resetPopup('norm', handles);
    set(handles.popup_norm, 'Value', 3);
elseif probStruct.norm==1,
    set(handles.popup_norm, 'Value', 1);
else
    set(handles.popup_norm, 'Value', 2);
end
if length(get(handles.popup_objective, 'String')) > 1,
    if isfield(probStruct, 'Qy')
        if probStruct.tracking==0,
            set(handles.popup_objective, 'Value', 2);
            if isfield(probStruct, 'yref'),
                set(handles.textfield_yref, 'String', mat2str(probStruct.yref));
                set(handles.textfield_yref, 'Enable', 'on');
                set(handles.text_yref, 'Enable', 'on');
            end
        else
            set(handles.popup_objective, 'Value', length(get(handles.popup_objective, 'String')));
        end
    else
        if ~isfield(probStruct, 'tracking') | probStruct.tracking==0,
            set(handles.popup_objective, 'Value', 1);
        else
            set(handles.popup_objective, 'Value', 3);
        end
    end
end
set(handles.textfield_N, 'String', num2str(probStruct.N));
if ~isfield(probStruct, 'Qy'),
    set(handles.textfield_yref, 'String', '');
    set(handles.textfield_yref, 'Enable', 'off');
    set(handles.text_yref, 'Enable', 'off');
end

if ~CeqEye,
    set(handles.popup_objective, 'String', {'Regulation towards origin', 'Regulation to output reference', ...
        'Free output trajectory tracking'});
end


% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiProbStruct_export_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function popup_solutiontype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_solutiontype (see GCBO)
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have adyn white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_solutiontype.
function popup_solutiontype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_solutiontype (see GCBO)
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_solutiontype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_solutiontype

global mpt___sysStruct mpt__probStruct
value = get(hObject, 'Value');
% 1 - Finite horizon
% 2 - Infinite time
% 3 - Minimum time
% 4 - Low complexity

sysStruct = mpt___sysStruct;
if ~isstruct(sysStruct),
    sub_resetPopup('norm', handles);
    sub_resetPopup('objective', handles);
    sub_resetPopup('solutiontype', handles);
    return
end

ispwa = iscell(sysStruct.A);
isexplicit = get(handles.radio_explicit, 'Value');

switch value,
    case 1,
        % finite horizon solution
        if ispwa & isexplicit,
            % no 2-norm cost function
            set(handles.popup_norm, 'String', ...
              {'Linear cost function (1-norm)', 'Linear cost function (Inf-norm)'});
            set(handles.popup_norm, 'Value', 1);
        else
            sub_resetPopup('norm', handles);
        end
        sub_resetPopup('objective', handles);
        set(handles.textfield_N, 'Enable', 'on');
        set(handles.text_N, 'Enable', 'on');
        mpt__probStruct.subopt_lev = 0;
        mpt__probStruct.N = str2double(get(handles.textfield_N, 'String'));
        
    case 2,
        % infinite-time solution
        if ispwa,
            % no 2-norm cost function
            set(handles.popup_norm, 'String', ...
                {'Linear cost function (1-norm)', 'Linear cost function (Inf-norm)'});
            set(handles.popup_norm, 'Value', 1);
        else
            sub_resetPopup('norm', handles);
        end
        set(handles.popup_objective, 'String', {'Regulation towards origin'});
        set(handles.textfield_N, 'String', 'Inf');
        set(handles.textfield_N, 'Enable', 'off');
        set(handles.text_N, 'Enable', 'off');
        mpt__probStruct.subopt_lev = 0;
        mpt__probStruct.N = Inf;
        
    case 3,
        % minimum-time solution
        set(handles.popup_objective, 'String', {'Regulation towards origin'});
        set(handles.textfield_N, 'String', 'Inf');
        set(handles.textfield_N, 'Enable', 'off');
        set(handles.text_N, 'Enable', 'off');
        mpt__probStruct.subopt_lev = 1;
        mpt__probStruct.N = Inf;
        sub_resetPopup('norm', handles);
        
    case 4,
        % low-complexity solution
        if ispwa & isexplicit,
            % no 2-norm cost function
            set(handles.popup_norm, 'String', ...
              {'Linear cost function (1-norm)', 'Linear cost function (Inf-norm)'});
            set(handles.popup_norm, 'Value', 1);
        else
            sub_resetPopup('norm', handles);
        end
        set(handles.popup_objective, 'String', {'Regulation towards origin'});
        set(handles.textfield_N, 'String', '1');
        set(handles.textfield_N, 'Enable', 'on');
        set(handles.text_N, 'Enable', 'on');
        mpt__probStruct.subopt_lev = 2;
        mpt__probStruct.N = str2double(get(handles.textfield_N, 'String'));

end
    


% --- Executes during object creation, after setting all properties.
function popup_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_norm (see GCBO)
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have adyn white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_norm.
function popup_norm_Callback(hObject, eventdata, handles)
% hObject    handle to popup_norm (see GCBO)
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_norm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_norm

global mpt__probStruct
contents = get(hObject, 'String');
value = contents{get(hObject,'Value')};

if ~isempty(findstr(value, '1-norm')),
    mpt__probStruct.norm = 1;
elseif ~isempty(findstr(value, 'Inf-norm'))
    mpt__probStruct.norm = Inf;
else
    mpt__probStruct.norm = 2;
end

% --- Executes during object creation, after setting all properties.
function textfield_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_N_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_N as text
%        str2double(get(hObject,'String')) returns contents of textfield_N as a double

global mpt__probStruct

value = get(hObject, 'String');
try
    %value = evalin('base', value);
    if ~isa(value, 'double'),
        value = str2num(value);
    end
catch
    sub_errordlg(lasterr);
    return
end
if any(size(value)~=1),
    sub_errordlg('Prediction horizon must be a scalar!');
    return
end
mpt__probStruct.N = value;

% --- Executes during object creation, after setting all properties.
function popup_objective_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_objective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_objective.
function popup_objective_Callback(hObject, eventdata, handles)
% hObject    handle to popup_objective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_objective contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_objective

global mpt__probStruct

contents = get(hObject, 'String');
str = contents{get(hObject,'Value')};
if ~isempty(findstr(str, 'tracking')),
    if isfield(mpt__probStruct, 'yref'),
        mpt__probStruct = rmfield(mpt__probStruct, 'yref');
    end
    if ~isempty(findstr(str, 'state')),
        if isfield(mpt__probStruct, 'Qy'),
            mpt__probStruct = rmfield(mpt__probStruct, 'Qy');
        end
    end
end
if isempty(findstr(str, 'output reference')),
    set(handles.textfield_yref, 'Enable', 'off');
    set(handles.text_yref, 'Enable', 'off');
else
    set(handles.textfield_yref, 'Enable', 'on');
    set(handles.text_yref, 'Enable', 'on');
end

if ~isempty(findstr(str, 'tracking')),
    mpt__probStruct.tracking = 1;
else
    mpt__probStruct.tracking = 0;
end



% --- Executes during object creation, after setting all properties.
function textfield_yref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_yref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_yref_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_yref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_yref as text
%        str2double(get(hObject,'String')) returns contents of textfield_yref as a double

global mpt__probStruct mpt___sysStruct

[nx, nu, ny] = mpt_sysStructInfo(mpt___sysStruct);
value = get(hObject, 'String');
try
    %value = evalin('base', value);
    if ~isa(value, 'double'),
        value = str2num(value);
    end
catch
    sub_errordlg(lasterr);
    return
end
yref = value(:);
if length(yref) ~= ny,
    sub_errordlg(sprintf('Reference point must be a %dx1 vector!', ny));
    return
end
mpt__probStruct.yref = yref;


% --- Executes on button press in vanced.
function button_advanced_Callback(hObject, eventdata, handles)
% hObject    handle to button_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mpt_guiProbAdv_export;

% --- Executes on button press in HelpType.
function HelpType_Callback(hObject, eventdata, handles)
% hObject    handle to HelpType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['* Finite Horizon Solution\n\n* Infinite Time Solution\n\n* Minimum Time Solution\n\n* Low Complexity Solution']));

% --- Executes on button press in HelpNorm.
function HelpNorm_Callback(hObject, eventdata, handles)
% hObject    handle to HelpNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['* 1-norm cost function - linear cost function\n\n' ...
        '* Inf-norm cost function - linear cost function\n\n'...
        '* 2-norm cost function - quadratic cost function']));

% --- Executes on button press in HelpN.
function HelpN_Callback(hObject, eventdata, handles)
% hObject    handle to HelpN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf('Horizon over which predictions should be made'));

% --- Executes on button press in HelpObjective.
function HelpObjective_Callback(hObject, eventdata, handles)
% hObject    handle to HelpObjective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(sprintf(['* Regulation towards origin - controller will drive all systems towards zero\n\n' ...
        '* Regulation to output reference - controller will drive system outputs to a given (fixed) reference point\n\n' ...
        '* Free state trajectory tracking - such a controller will lead system STATES to track a certain free reference signal '...
        'which can be specified at the time the control law is about to be applied\n\n'...
        '* Free output trajectory tracking - instead of driving system states to certain reference signals, this controller '...
        'will force system OUTPUTS to track certain reference signals.']));

% --- Executes on button press in HelpYref.
function HelpYref_Callback(hObject, eventdata, handles)
% hObject    handle to HelpYref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('Value of fixed output reference you want to track.');


% --- Executes on button press in radio_explicit.
function radio_explicit_Callback(hObject, eventdata, handles)
% hObject    handle to radio_explicit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_explicit

set(handles.radio_explicit, 'Value', 1);
set(handles.radio_online, 'Value', 0);
sub_resetPopup('norm', handles);
sub_resetPopup('objective', handles);
sub_resetPopup('solutiontype', handles);

% --- Executes on button press in radio_online.
function radio_online_Callback(hObject, eventdata, handles)
% hObject    handle to radio_online (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_online

set(handles.radio_explicit, 'Value', 0);
set(handles.radio_online, 'Value', 1);
sub_resetPopup('norm', handles);
sub_resetPopup('objective', handles);
set(handles.popup_solutiontype, 'String', {'Finite Horizon Solution'});
set(handles.popup_solutiontype, 'Value', 1);

% --- Executes on button press in button_simplify.
function button_simplify_Callback(hObject, eventdata, handles)
% hObject    handle to button_simplify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

ctrlname = get(handles.textfield_ctrlname, 'String');
if ~isvarname(ctrlname)
    sub_errordlg(sprintf('"%s" is not a valid Matlab variable name!', ctrlname));
    return
end
if ~isa(mpt__ctrl, 'mptctrl')
    errodlg('Unknown controller type!');
    return
end
if mpt__ctrl.simplified,
    msgbox('The controller partition has already been simplified!', 'modal');
    return
end
nRbefore = length(mpt__ctrl);
set(handles.text_cdetails, 'String', 'Simplifying controller...');
try
    startt = clock;
    mpt__ctrl = mpt_simplify(mpt__ctrl, struct('statusbar', 1));
    runtime = etime(clock, startt);
catch
    mpt_statusbar;
    sub_errordlg(lasterr);
    set(handles.text_cdetails, 'String', strvcat(sprintf('Simplification error:\n', lasterr)));
    return
end
nRafter = length(mpt__ctrl);
msgbox(sprintf('Controller partition simplified from %d to %d regions in %d seconds.', nRbefore, nRafter, round(runtime)), 'modal');
set(handles.text_cdetails, 'String', char(mpt__ctrl));
set(handles.button_simplify, 'Enable', 'off');
set(handles.button_plotJ, 'Enable', 'off');
ctrlname = get(handles.textfield_ctrlname, 'String');
assignin('base', ctrlname, mpt__ctrl);


% --- Executes on button press in button_invariant.
function button_invariant_Callback(hObject, eventdata, handles)
% hObject    handle to button_invariant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

ctrlname = get(handles.textfield_ctrlname, 'String');
if ~isvarname(ctrlname)
    sub_errordlg(sprintf('"%s" is not a valid Matlab variable name!', ctrlname));
    return
end
if ~isa(mpt__ctrl, 'mptctrl')
    errodlg('Unknown controller type!');
    return
end
if isinvariant(mpt__ctrl),
    msgbox('The controller is already invariant!', 'modal');
    return
end
set(handles.text_cdetails, 'String', 'Computing invariant subset...');
try
    startt = clock;
    mpt__ctrl = mpt_invariantSet(mpt__ctrl, struct('statusbar', 1));
    runtime = etime(clock, startt);
catch
    mpt_statusbar;
    sub_errordlg(lasterr);
    %set(handles.text_cdetails, 'String', strvcat(sprintf('An error has occured:\n', lasterr)));
    return
end
msgbox('The invariant set has been successfully computed.', 'modal');
set(handles.text_cdetails, 'String', char(mpt__ctrl));
set(handles.button_invariant, 'Enable', 'off');
set(handles.button_simplify, 'Enable', 'on');
assignin('base', ctrlname, mpt__ctrl);


% --- Executes on button press in button_lyapunov.
function button_lyapunov_Callback(hObject, eventdata, handles)
% hObject    handle to button_lyapunov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

ctrlname = get(handles.textfield_ctrlname, 'String');
if ~isvarname(ctrlname)
    sub_errordlg(sprintf('"%s" is not a valid Matlab variable name!', ctrlname));
    return
end
if ~isa(mpt__ctrl, 'mptctrl')
    errodlg('Unknown controller type!');
    return
end
if isstabilizable(mpt__ctrl),
    msgbox('The controller stabilizes the given system, computation of a Lyapunov is not needed.', 'modal');
    return
end

lyaptype = mpt_guiLyapunov_export;

set(handles.text_cdetails, 'String', 'Computing lyapunov function...');
try
    startt = clock;
    mpt__ctrl = mpt_lyapunov(mpt__ctrl, lyaptype, struct('statusbar', 1));
    runtime = etime(clock, startt);
catch
    mpt_statusbar;
    sub_errordlg(lasterr);
    %set(handles.text_cdetails, 'String', strvcat(sprintf('An error has occured:\n', lasterr)));
    return
end
if isstabilizable(mpt__ctrl),
    msgbox('The Lyapunov function has been successfully computed.', 'modal');
else
    warndlg('Lypunov function not found, the closed-loop system may be unstable!');
    return
end
set(handles.text_cdetails, 'String', char(mpt__ctrl));
set(handles.button_lyapunov, 'Enable', 'off');
assignin('base', ctrlname, mpt__ctrl);

% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)
% hObject    handle to button_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

try
    Options.newfigure = 0;
    figure
    mpt_plotPartition(mpt__ctrl, Options);
catch
    sub_errordlg(lasterr);
end

% --- Executes on button press in button_plotU.
function button_plotU_Callback(hObject, eventdata, handles)
% hObject    handle to button_plotU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

try
    Options.newfigure = 0;
    figure
    mpt_plotU(mpt__ctrl, Options);
catch
    sub_errordlg(lasterr);
end

% --- Executes on button press in button_plotJ.
function button_plotJ_Callback(hObject, eventdata, handles)
% hObject    handle to button_plotJ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

try
    Options.newfigure = 0;
    figure
    mpt_plotJ(mpt__ctrl, Options);
catch
    sub_errordlg(lasterr);
end


% --- Executes on button press in button_generate.
function button_generate_Callback(hObject, eventdata, handles)
% hObject    handle to button_generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct mpt___sysStruct mpt__ctrl

ctrlname = get(handles.textfield_ctrlname, 'String');
if ~isvarname(ctrlname)
    sub_errordlg(sprintf('"%s" is not a valid Matlab variable name!', ctrlname));
    return
end

if isempty(mpt___sysStruct)
    errordlg('Define a model first!', 'Error', 'modal');
    return
end
if isempty(mpt__probStruct)
    errordlg('Define controller properties first!', 'Error', 'modal');
    return
end
try
    sysStruct = mpt___sysStruct;
    probStruct = mpt__probStruct;
    Options.guierrors = 1;
    Options.verbose = 0;
    Options.forceverify = 1;
    T = evalc('[sysStruct, probStruct] = mpt_verifySysProb(sysStruct, probStruct, Options);');
catch
    sub_errordlg(lasterr);
    return
end
    
try
    set(handles.text_cdetails, 'String', 'Computing controller...');
    sub_setPostProcessButtons('off', handles);
    online = get(handles.radio_online, 'Value');
    startt = clock;
    if online,
        ctrl = mpt_control(sysStruct, probStruct, 'online');
    else
        ctrl = mpt_control(sysStruct, probStruct, struct('statusbar', 1, 'guierrors', 1));
    end
    runtime = etime(clock, startt);
    assignin('base', ctrlname, ctrl);
    set(handles.text_cdetails, 'String', char(ctrl));
    sub_setPostProcessButtons('on', handles, ctrl);
    if online,
        msgbox('The controller has been successfully generated and saved to workspace.', 'modal');
    else
        msgbox(sprintf('The controller has been successfully generated in %d seconds and saved to workspace.', ceil(runtime)), 'modal');
    end
    mpt__ctrl = ctrl;
catch
    mpt_statusbar;
    set(handles.text_cdetails, 'String', strvcat(sprintf('Controller computation failed:\n'), lasterr));
    sub_errordlg(lasterr);
    return
end


% --- Executes on button press in button_weights.
function button_weights_Callback(hObject, eventdata, handles)
% hObject    handle to button_weights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mpt_guiWeights_export;

% --- Executes on button press in button_loadmodel.
function button_loadmodel_Callback(hObject, eventdata, handles, sysStruct)
% hObject    handle to button_invariant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct

if nargin<4,
    sysStruct = mpt_guiSysStruct_export;
end

if isstruct(sysStruct),
    [nx, nu, ny, ndyn] = mpt_sysStructInfo(sysStruct);
    ispwa = iscell(sysStruct.A);
    if ispwa,
        info_str = sprintf('PWA system (%d dynamics)\n', ndyn);
    else
        info_str = sprintf('LTI system\n');
    end
    plurx = ['state' repmat('s', 1, (nx>1))];
    pluru = ['input' repmat('s', 1, (nu>1))];
    plury = ['output' repmat('s', 1, (ny>1))];
    sysinfo_str = sprintf('%d %s\n%d %s\n%d %s', ...
    nx, plurx, nu, pluru, ny, plury);

    info_str = strvcat(info_str, sysinfo_str);
    set(handles.text_systeminfo, 'String', info_str);
else
    return
end
mpt___sysStruct = sysStruct;
sub_resetPopup('all', handles);
popup_solutiontype_Callback(handles.popup_solutiontype, eventdata, handles)
set(handles.button_weights, 'Enable', 'on');
set(handles.button_advanced, 'Enable', 'on');


%=============================================================
function sub_resetPopup(pType, handles)

if strcmp(pType, 'norm') | strcmp(pType, 'all'),
    set(handles.popup_norm, 'String', ...
        {'Linear cost function (1-norm)', 'Linear cost function (Inf-norm)', ...
        'Quadratic cost function (2-norm)'});
    set(handles.popup_norm, 'Value', 1);
end
if strcmp(pType, 'objective') | strcmp(pType, 'all'),
    if ~CeqEye,
        set(handles.popup_objective, 'String', {'Regulation towards origin', 'Regulation to output reference', ...
                'Free output trajectory tracking'});
        set(handles.popup_objective, 'Value', 1);
    else
        set(handles.popup_objective, 'String', ...
            {'Regulation towards origin', 'Regulation to output reference', ...
                'Free state trajectory tracking', 'Free output trajectory tracking'});
        set(handles.popup_objective, 'Value', 1);
    end
end
if strcmp(pType, 'solutiontype') | strcmp(pType, 'all'),
    set(handles.popup_solutiontype, 'String', ...
        {'Finite Horizon Solution', 'Infinite Time Solution', ...
        'Minimum Time Solution', 'Low Complexity Solution'});
    set(handles.popup_solutiontype, 'Value', 1);
    
end


%=============================================================
function sub_setPostProcessButtons(status, handles, ctrl)

if strcmp(status, 'off'),
    set(handles.button_simplify, 'Enable', 'off');
    set(handles.button_invariant, 'Enable', 'off');
    set(handles.button_lyapunov, 'Enable', 'off');
    set(handles.button_plot, 'Enable', 'off');
    set(handles.button_plotU, 'Enable', 'off');
    set(handles.button_plotJ, 'Enable', 'off');
    set(handles.button_simulate, 'Enable', 'off');
    return
end
if ~isexplicit(ctrl),
    set(handles.button_simplify, 'Enable', 'off');
    set(handles.button_invariant, 'Enable', 'off');
    set(handles.button_lyapunov, 'Enable', 'off');
    set(handles.button_plot, 'Enable', 'off');
    set(handles.button_plotU, 'Enable', 'off');
    set(handles.button_plotJ, 'Enable', 'off');
    set(handles.button_simulate, 'Enable', 'on');
    return
end

if ctrl.simplified,
    set(handles.button_simplify, 'Enable', 'off');
else
    set(handles.button_simplify, 'Enable', status);
end
if isinvariant(ctrl),
    set(handles.button_invariant, 'Enable', 'off');
else
    set(handles.button_invariant, 'Enable', status);
end
if isstabilizable(ctrl),
    set(handles.button_lyapunov, 'Enable', 'off');
else
    set(handles.button_lyapunov, 'Enable', status);
end
set(handles.button_plot, 'Enable', status);

if dimension(ctrl.Pn)<1 | dimension(ctrl.Pn)>2,
    set(handles.button_plotU, 'Enable', 'off');
else
    set(handles.button_plotU, 'Enable', status);
end

if ctrl.simplified | dimension(ctrl.Pn)<1 | dimension(ctrl.Pn)>2,
    set(handles.button_plotJ, 'Enable', 'off');
else
    set(handles.button_plotJ, 'Enable', status);
end
set(handles.button_simulate, 'Enable', status);

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


% --- Executes on button press in button_simulate.
function button_simulate_Callback(hObject, eventdata, handles)
% hObject    handle to button_simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl


if isempty(mpt__ctrl)
    sub_errordlg('No controller loaded!');
    return
end
if ~isa(mpt__ctrl, 'mptctrl')
    sub_errordlg('Unknown controller!');
end

clear global mpt__sim_options
mpt_guiSimulate_export(mpt__ctrl);


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_loadsession_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadsession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_newsession_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_newsession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct mpt___sysStruct mpt__ctrl

res = questdlg('Clear everything and start a new session?', 'Confirm action', ...
    'Yes', 'No', 'No');

if strcmpi(res, 'yes'),
    probStruct.norm = 1;
    probStruct.subopt_lev = 0;
    probStruct.N = 1;
    probStruct.Q = [];
    probStruct.R = [];
    probStruct.tracking = 0;
    mpt__probStruct = probStruct;
    sub_resetPopup('all', handles);
    sub_setProbStructFields(handles, probStruct);
    set(handles.radio_explicit, 'Value', 1);
    set(handles.radio_online, 'Value', 0);
    mpt___sysStruct = [];
    set(handles.text_systeminfo, 'String', 'No model defined');
    set(handles.text_cdetails, 'String', 'No controller available');
    mpt__ctrl = [];
    sub_setPostProcessButtons('off', handles)
end


% --------------------------------------------------------------------
function menu_file_savesession_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_savesession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct mpt__probStruct

res = inputdlg({'System structure:', 'Problem structure:'}, ...
    'Export system and problem setup to workspace', ...
    1, {'sysStruct', 'probStruct'});

if ~isempty(res),
    sysStructVar = res{1};
    probStructVar = res{2};
    if ~isvarname(sysStructVar),
        sub_errordlg(sprintf('"%s" is not a valid variable name!', sysStructVar));
        return
    end
    if ~isvarname(probStructVar),
        sub_errordlg(sprintf('"%s" is not a valid variable name!', probStructVar));
        return
    end
    assignin('base', sysStructVar, mpt___sysStruct);
    assignin('base', probStructVar, mpt__probStruct);
end

% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

res = questdlg('Do you really want to exit MPT Studio?', 'Confirm action', ...
    'Yes', 'No', 'No');

if strcmpi(res, 'yes'),
    delete(handles.figure1);
end

% --------------------------------------------------------------------
function menu_file_loadfile_Callback(hObject, eventdata, handles, fname)
% hObject    handle to menu_file_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct mpt__probStruct mpt__ctrl

if nargin<4,
    [fname, fpath] = uigetfile({'*.m', 'MATLAB files (*.m)'}, 'Pick an M-file');
    if fname==0,
        % the dialog box was closed manually
        return
    end
elseif nargin==4,
    fpath = '';
end

sysStruct = [];
probStruct = [];
if ~isempty(fname)
    sysStruct = execute_sysstructfile(fname, fpath);
    probStruct = execute_probstructfile(fname, fpath);
end
if ~isstruct(sysStruct) | isempty(sysStruct) | ~isstruct(probStruct) | isempty(probStruct),
    return
end

mpt__ctrl = [];
if isfield(probStruct, 'N'),
    set(handles.textfield_N, 'String', num2str(probStruct.N));
end
button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles, sysStruct);
set(handles.text_cdetails, 'String', 'No controller available');
sub_setPostProcessButtons('off', handles)
if ~isfield(probStruct, 'tracking'),
    probStruct.tracking = 0;
end
sub_resetPopup('all', handles);
sub_setProbStructFields(handles, probStruct);
mpt___sysStruct = sysStruct;
mpt__probStruct = probStruct;
mpt_setSysStruct(sysStruct);


% --------------------------------------------------------------------
function menu_file_savefile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct mpt__probStruct

if ~isstruct(mpt___sysStruct) | isempty(mpt___sysStruct)
    errordlg('Define a model first!', 'Error', 'modal');
    return
end
if ~isstruct(mpt__probStruct) | isempty(mpt__probStruct),
    errordlg('Define problem properties first!', 'Error', 'modal');
    return
end

[fname, fpath] = uiputfile('*.m', 'Save system as');

allowedfields = {'norm', 'N', 'subopt_lev', 'Q', 'R', 'Qy', 'Rdu', ...
        'y0bounds', 'yref', 'tracking', 'Tconstraint', 'useSymmetry', ...
        'feedback', 'FBgain', 'xref', 'uref', 'P_N'};

if isfield(mpt__probStruct, 'Tset')
    if isa(mpt__probStruct.Tset, 'polytope')
        if isfulldim(mpt__probStruct.Tset),
            allowedfields{end+1} = 'Tset';
        end
    end
end

if ischar(fname),
    try
        sub_saveSysStruct(mpt___sysStruct, fname, fpath);
        sub_saveSysStruct(mpt__probStruct, fname, fpath, allowedfields, 'a');
    catch
        sub_errordlg(lasterr);
    end
end


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg = {'Multi-Parametric Toolbox version 2.0', ...
        '', 'Copyright (C) 2003-2005 by Michal Kvasnica, Pascal Grieder, Mato Baotic', ...
        '', 'Send feedback, questions or bug reports to: mpt@control.ee.ethz.ch'};

msgbox(msg, 'About', 'modal');

% --------------------------------------------------------------------
function menu_help_gotompt_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_gotompt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if ispc | 1,
        web('http://control.ee.ethz.ch/~mpt/docs/');
    else
        web('http://control.ee.ethz.ch/~mpt/docs/', '-browser');
    end
catch
    sub_errordlg(lasterr);
end

% --------------------------------------------------------------------
function menu_help_version_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function textfield_ctrlname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_ctrlname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_ctrlname_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_ctrlname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_ctrlname as text
%        str2double(get(hObject,'String')) returns contents of textfield_ctrlname as a double


function menu_help_credits_Callback(hObject, eventdata, handles)

credits = {'List of contributors:', '', 'Baotic, Mato', ...
        'Baric, Miroslav', ...
        'Biswas, Pratik', ...
        'Cagienard, Raphael', ...
        'Christophersen, Frank J.', ...
        'Fukuda, Komei', ...
        'Geyer, Tobias', ...
        'Grieder, Pascal', ...
        'Jones, Colin N.', ...
        'Kvasnica, Michal', ...
        'Lagerberg, Adam', ...
        'Linder, Arne', ...
        'L?fberg, Johan', ...
        'Suard, Raphael', ...
        'Tetrev, Boris', ...
        'Torrisi, Fabio D.'};

msgbox(credits, 'Credits', 'modal');



% --------------------------------------------------------------------
function menu_file_loadw_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct mpt__probStruct mpt__ctrl

ctrl = mpt_guiLoadWorkspace_export('mptctrl');

if ~isempty(ctrl)
    ctrlStruct = struct(ctrl);
    if isfield(ctrlStruct.details, 'origSysStruct'),
        sysStruct = ctrlStruct.details.origSysStruct;
    else
        sysStruct = ctrlStruct.sysStruct;
    end
    if isfield(ctrlStruct.details, 'origProbStruct'),
        probStruct = ctrlStruct.details.origProbStruct;
    else
        probStruct = ctrlStruct.probStruct;
    end
    mpt__ctrl = ctrl;
    mpt___sysStruct = sysStruct;
    mpt__probStruct = probStruct;
    mpt_setSysStruct(sysStruct);

    button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles, sysStruct);
    
    if ~isfield(probStruct, 'tracking'),
        probStruct.tracking = 0;
    end
    sub_resetPopup('all', handles);
    sub_setProbStructFields(handles, probStruct);
    
    %set(handles.text_cdetails, 'String', 'No controller available');
    sub_setPostProcessButtons('on', handles, ctrl)
    set(handles.text_cdetails, 'String', char(mpt__ctrl));
end

% --------------------------------------------------------------------
function menu_file_savew_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_saveworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__ctrl

if isempty(mpt__ctrl)
    errordlg('Generate a controller first!', 'Error', 'modal');
    return
end

varname = inputdlg({'Export to workspace as'}, 'Export to workspace', 1, {'ctrl'});

if ~isempty(varname),
    try
        assignin('base', varname{1}, mpt__ctrl);
    catch
        sub_errordlg(lasterr);
    end
end

% --------------------------------------------------------------------
function menu_model_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_model_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles);

% --------------------------------------------------------------------
function menu_model_loadw_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model_loadworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct

sysStruct = mpt_guiLoadWorkspace_export('sysStruct');

if ~isempty(sysStruct)
    mpt___sysStruct = sysStruct;
    mpt_setSysStruct(sysStruct);
    button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles, sysStruct);
end

% --------------------------------------------------------------------
function menu_model_loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct

[fname, fpath] = uigetfile({'*.m', 'MATLAB files (*.m)'}, 'Pick an M-file');
if fname==0,
    % the dialog box was closed manually
    return
end
sysStruct = [];
if ~isempty(fname)
    sysStruct = execute_sysstructfile(fname, fpath);
end
if ~isstruct(sysStruct) | isempty(sysStruct),
    return
end

mpt___sysStruct = sysStruct;
mpt_setSysStruct(sysStruct);
button_loadmodel_Callback(handles.button_loadmodel, eventdata, handles, sysStruct);


% --------------------------------------------------------------------
function menu_model_savew_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model_saveworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct

if ~isstruct(mpt___sysStruct) | isempty(mpt___sysStruct)
    errordlg('Define a model first!', 'Error', 'modal');
    return
end

varname = inputdlg({'Export to workspace as'}, 'Export to workspace', 1, {'sysStruct'});

if ~isempty(varname),
    try
        assignin('base', varname{1}, mpt___sysStruct);
    catch
        sub_errordlg(lasterr);
    end
end

% --------------------------------------------------------------------
function menu_model_savefile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_model_savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt___sysStruct

if ~isstruct(mpt___sysStruct) | isempty(mpt___sysStruct)
    errordlg('Define a model first!', 'Error', 'modal');
    return
end

[fname, fpath] = uiputfile('*.m', 'Save system as');

if ischar(fname),
    try
        sub_saveSysStruct(mpt___sysStruct, fname, fpath);
    catch
        sub_errordlg(lasterr);
    end
end


% --------------------------------------------------------------------
function menu_problem_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_problem_loadw_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem_loadworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct

probStruct = mpt_guiLoadWorkspace_export('probstruct');

if ~isempty(probStruct)
    mpt__probStruct = probStruct;
    if ~isfield(probStruct, 'tracking'),
        probStruct.tracking = 0;
    end
    sub_resetPopup('all', handles);
    sub_setProbStructFields(handles, probStruct);
end

% --------------------------------------------------------------------
function menu_problem_loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct

[fname, fpath] = uigetfile({'*.m', 'MATLAB files (*.m)'}, 'Pick an M-file');
if fname==0,
    % the dialog box was closed manually
    return
end
probStruct = [];
if ~isempty(fname)
    probStruct = execute_probstructfile(fname, fpath);
end
if ~isstruct(probStruct) | isempty(probStruct),
    return
end

if ~isfield(probStruct, 'tracking'),
    probStruct.tracking = 0;
end
sub_resetPopup('all', handles);
sub_setProbStructFields(handles, probStruct);
mpt__probStruct = probStruct;

% --------------------------------------------------------------------
function menu_problem_savew_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem_saveworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct

if ~isstruct(mpt__probStruct) | isempty(mpt__probStruct)
    errordlg('Define problem properties first!', 'Error', 'modal');
    return
end

varname = inputdlg({'Export to workspace as'}, 'Export to workspace', 1, {'probStruct'});

if ~isempty(varname),
    try
        assignin('base', varname{1}, mpt__probStruct);
    catch
        sub_errordlg(lasterr);
    end
end


% --------------------------------------------------------------------
function menu_problem_savefile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem_savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct

if ~isstruct(mpt__probStruct) | isempty(mpt__probStruct)
    errordlg('Define problem properties first!', 'Error', 'modal');
    return
end

[fname, fpath] = uiputfile('*.m', 'Save problem as');

allowedfields = {'norm', 'N', 'subopt_lev', 'Q', 'R', 'Qy', 'Rdu', ...
        'y0bounds', 'yref', 'tracking', 'Tconstraint', 'useSymmetry', ...
        'feedback', 'FBgain', 'xref', 'uref', 'P_N'};

if isfield(mpt__probStruct, 'Tset')
    if isa(mpt__probStruct.Tset, 'polytope')
        if isfulldim(mpt__probStruct.Tset),
            allowedfields{end+1} = 'Tset';
        end
    end
end

if ischar(fname),
    try
        sub_saveSysStruct(mpt__probStruct, fname, fpath, allowedfields);
    catch
        sub_errordlg(lasterr);
    end
end


%============================================================
function sysStruct = execute_sysstructfile(fname, fpath)

origfname = fname;
dotpos = findstr(fname, '.');
if ~isempty(dotpos),
    try
        fname = fname(1:dotpos-1);
    catch
        sub_errordlg(sprintf('Corrupted file name "%s"', origfname));
    end
end
sysStruct = [];
curdir = pwd;
try
    if ~isempty(fpath),
        cd(fpath);
    end
    T = evalc(fname);
    if ~isempty(fpath),
        cd(curdir);
    end
    if ~mpt_issysstruct(sysStruct),
        sub_errordlg(sprintf('The file "%s" does not contain a valid system structure!', origfname));
    end
catch
    sub_errordlg(sprintf('The file "%s" does not contain a valid system structure!', origfname));
end
if ~isempty(fpath)
    cd(curdir);
end


%============================================================
function probStruct = execute_probstructfile(fname, fpath)

origfname = fname;
dotpos = findstr(fname, '.');
if ~isempty(dotpos),
    try
        fname = fname(1:dotpos-1);
    catch
        sub_errordlg(sprintf('Corrupted file name "%s"', origfname));
    end
end
probStruct = [];
curdir = pwd;
try
    if ~isempty(fpath),
        cd(fpath);
    end
    T = evalc(fname);
    if ~isempty(fpath),
        cd(curdir);
    end
    if ~sub_isprobstruct(probStruct),
        sub_errordlg(sprintf('The file "%s" does not contain a valid problem setup!', origfname));
    end
catch
    sub_errordlg(sprintf('The file "%s" does not contain a valid problem setup!', origfname));
end
if ~isempty(fpath),
    cd(curdir);
end


%============================================================
function sub_saveSysStruct(sysStruct, fname, fpath, allowedfields, filemode)

if isempty(findstr(fname, '.m')),
    fname = [fname '.m'];
end
fullname = [fpath fname];
if nargin<5,
    fid = fopen(fullname, 'w');
else
    fid = fopen(fullname, filemode);
end
if fid<0,
    sub_errordlg(sprintf('Error opening file "%s" for writing!', fullname));
    return
end

if nargin<4,
    allowedfields = {'A', 'B', 'C', 'D', 'f', 'g', 'guardX', ...
            'guardU', 'guardC', 'umax', 'umin', 'dumax', 'dumin', ...
            'ymax', 'ymin', 'dymax', 'dymin', 'xmax', 'xmin', ...
            'Uset', 'InputName', 'StateName', 'OutputName'};
    structname = 'sysStruct';
else
    structname = 'probStruct';
end
    
ssFields = fields(sysStruct);

if isfield(sysStruct, 'data')
    % this sysStruct came from a HYSDEL file
    hysdelfile = sysStruct.data.hysdel.fromfile;
    str = ['''Importing from "' hysdelfile '"...'''];
    fprintf(fid, 'disp(%s);\n', str);
    fprintf(fid, 'sysStruct = mpt_sys(''%s'');\n', hysdelfile);
    allowedfields = {'umax', 'umin', 'dumax', 'dumin', ...
        'ymax', 'ymin', 'dymax', 'dymin', 'xmax', 'xmin', ...
        'InputName', 'StateName', 'OutputName'};
end

for ii = 1:length(ssFields),
    fieldname = ssFields{ii};
    if ~any(ismember(allowedfields, fieldname)),
        continue
    end
    value = getfield(sysStruct, fieldname);
    if iscell(value),
        for jj = 1:length(value),
            cval = value{jj};
            if isa(cval, 'polytope'),
                [H, K] = double(cval);
                fprintf(fid, '%s.%s{%d} = polytope(%s, %s);\n', structname, fieldname, jj, mat2str(H), mat2str(K));
                
            elseif ischar(cval),
                cval = sprintf('''%s''', cval);

            elseif ~ischar(cval),
                cval = mat2str(cval);
            end
            
            if ~isa(cval, 'polytope'),
                fprintf(fid, '%s.%s{%d} = %s;\n', structname, fieldname, jj, cval);
            end
        end
    else
        if isa(value, 'polytope'),
            [H, K] = double(value);
            fprintf(fid, '%s.%s = polytope(%s, %s);\n', structname, fieldname, mat2str(H), mat2str(K));
        elseif ischar(value),
            value = sprintf('''%s''', value);
            
        elseif ~ischar(value),
            value = mat2str(value);
        end
        if ~isa(value, 'polytope') & ~isstruct(value),
            fprintf(fid, '%s.%s = %s;\n', structname, fieldname, value);
        end
    end
end

fclose(fid);


% --------------------------------------------------------------------
function menu_problem_reset_Callback(hObject, eventdata, handles)
% hObject    handle to menu_problem_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct

res = questdlg('Reset problem setup to default values?', 'Confirm action', ...
    'Yes', 'No', 'No');

if strcmpi(res, 'yes'),
    probStruct.norm = 1;
    probStruct.subopt_lev = 0;
    probStruct.N = 1;
    probStruct.Q = [];
    probStruct.R = [];
    probStruct.tracking = 0;
    mpt__probStruct = probStruct;
    sub_resetPopup('all', handles);
    sub_setProbStructFields(handles, probStruct);
    set(handles.radio_explicit, 'Value', 1);
    set(handles.radio_online, 'Value', 0);
end


%===============================================================
function res = sub_isprobstruct(str)

res = 0;
if isstruct(str),
   if isfield(str, 'N') & isfield(str, 'Q') & isfield(str, 'R'),
       res = 1;
   end
end


% --------------------------------------------------------------------
function menu_examples_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_ex_di_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples_di (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_file_loadfile_Callback(handles.menu_file_loadfile, eventdata, handles, 'Double_Integrator');

% --------------------------------------------------------------------
function menu_ex_3dlti_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples_3dlti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_file_loadfile_Callback(handles.menu_file_loadfile, eventdata, handles, 'ThirdOrder');

% --------------------------------------------------------------------
function menu_ex_pwasincos_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples_pwasincos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_file_loadfile_Callback(handles.menu_file_loadfile, eventdata, handles, 'pwa_sincos');

% --------------------------------------------------------------------
function menu_ex_pwacar_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples_pwacar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_file_loadfile_Callback(handles.menu_file_loadfile, eventdata, handles, 'pwa_car');

% --------------------------------------------------------------------
function menu_ex_twotanks_Callback(hObject, eventdata, handles)
% hObject    handle to menu_examples_twotanks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_file_loadfile_Callback(handles.menu_file_loadfile, eventdata, handles, 'two_tanks');


%============================================================
function sub_helpdlg(msg)

msgbox(msg, 'Help', 'modal');


% --- Executes during object creation, after setting all properties.
function text_cdetails_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_cdetails (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


% --- Executes on selection change in text_cdetails.
function text_cdetails_Callback(hObject, eventdata, handles)
% hObject    handle to text_cdetails (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns text_cdetails contents as cell array
%        contents{get(hObject,'Value')} returns selected item from text_cdetails




% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiProbStruct_export_LayoutFcn(policy)
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
'Name','MPT Studio',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[40 34.5384615384615 109.8 43.3846153846154],...
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
'popupmenu', 10, ...
'text', 62, ...
'edit', 39, ...
'frame', 20, ...
'checkbox', 9, ...
'togglebutton', 2, ...
'pushbutton', 60, ...
'radiobutton', 5, ...
'listbox', 2), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiProbStruct.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[1.4 20.1538461538462 106.4 3.23076923076923],...
'String',{  '' },...
'Style','frame',...
'Tag','frame19');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[1.4 10.0769230769231 106.4 9.15384615384616],...
'String',{  '' },...
'Style','frame',...
'Tag','frame18');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'CData',[],...
'ListboxTop',0,...
'Position',[41.8 24.5384615384615 66.2 18],...
'String',{  '' },...
'Style','frame',...
'Tag','frame13',...
'UserData',[]);


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'CData',[],...
'ListboxTop',0,...
'Position',[1.4 24.5384615384615 38 18],...
'String',{  '' },...
'Style','frame',...
'Tag','frame12',...
'UserData',[]);


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''popup_solutiontype_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[62.6 36.2307692307692 38.2 1.69230769230769],...
'String',{  'Finite Horizon Solution'; 'Infinite Time Solution'; 'Minimum Time Solution'; 'Low Complexity Solution' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiProbStruct_export(''popup_solutiontype_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','popup_solutiontype');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[43.8 36.4615384615385 17.2 1.23076923076923],...
'String','Type of solution',...
'Style','text',...
'Tag','text1');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[44.8 34.3076923076923 16 1.15384615384615],...
'String','Cost function',...
'Style','text',...
'Tag','text2');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''popup_norm_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[62.6 34 38.2 1.69230769230769],...
'String',{  'Linear cost function (1-norm)'; 'Linear cost function (Inf-norm)'; 'Quadratic cost function (2-norm)' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiProbStruct_export(''popup_norm_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','popup_norm');


h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[60.8 41.7692307692308 28.4 1.23076923076923],...
'String','Controller Setup',...
'Style','text',...
'Tag','text11');


h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_loadmodel_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[3.6 39.3076923076923 33.4 1.76923076923077],...
'String','Load / Edit model',...
'Tag','button_loadmodel');


h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[7.6 41.7692307692308 24 1.30769230769231],...
'String','Model Setup',...
'Style','text',...
'Tag','text39');


h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[43.6 32 17.2 1.15384615384615],...
'String','Control objective',...
'Style','text',...
'Tag','text41');


h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''popup_objective_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[62.6 31.7692307692308 38.2 1.69230769230769],...
'String',{  'Regulation towards origin'; 'Regulation to output reference'; 'Free state trajectory tracking'; 'Free output trajectory tracking' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiProbStruct_export(''popup_objective_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','popup_objective');


h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Enable','off',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[45.4 27.6923076923077 15.4 1.15384615384615],...
'String','Reference point',...
'Style','text',...
'Tag','text_yref');


h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''textfield_yref_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[62.6 27.4615384615385 38.2 1.61538461538462],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiProbStruct_export(''textfield_yref_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_yref');


h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_advanced_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[78.4 24.9230769230769 21.4 1.92307692307692],...
'String','Advanced Options',...
'Tag','button_advanced');


h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''HelpType_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[102.8 36.2307692307692 3.4 1.76923076923077],...
'String','?',...
'Tag','HelpType');


h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''HelpNorm_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[102.8 34.0769230769231 3.4 1.76923076923077],...
'String','?',...
'Tag','HelpNorm');


h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''HelpObjective_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[103 31.8461538461538 3.4 1.76923076923077],...
'String','?',...
'Tag','HelpObjective');


h21 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''HelpYref_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[103 27.4615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','HelpYref');


h22 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[43.6 30 17.4 1.15384615384615],...
'String','Prediction horizon',...
'Style','text',...
'Tag','text_N');


h23 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''textfield_N_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[62.6 29.7692307692308 38.2 1.61538461538462],...
'String','1',...
'Style','edit',...
'CreateFcn','mpt_guiProbStruct_export(''textfield_N_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_N');


h24 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''HelpN_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[103 29.6153846153846 3.4 1.76923076923077],...
'String','?',...
'Tag','HelpN');


h25 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''radio_explicit_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 40.0769230769231 38.2 1.38461538461538],...
'String','Explicit Controller',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_explicit');


h26 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''radio_online_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[63 38.5384615384615 38.2 1.38461538461538],...
'String','On-line Controller',...
'Style','radiobutton',...
'Tag','radio_online');


h27 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[44 40.0769230769231 17.2 1.23076923076923],...
'String','Type of controller',...
'Style','text',...
'Tag','text55');


h28 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[1.6 0.461538461538462 106.4 8.53846153846154],...
'String',{  '' },...
'Style','frame',...
'Tag','frame16');


h29 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[40 8.38461538461539 29.6 1.23076923076923],...
'String','Post-Processing',...
'Style','text',...
'Tag','text56');


h30 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_simplify_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[5 6.15384615384615 30.8 1.76923076923077],...
'String','Merge Controller Regions',...
'Tag','button_simplify');


h31 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_invariant_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[39.2 6.15384615384615 30.8 1.76923076923077],...
'String','Compute Invariant Subset',...
'Tag','button_invariant');


h32 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_lyapunov_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[73.6 6.15384615384615 30.8 1.76923076923077],...
'String','Compute Lyapunov Function',...
'Tag','button_lyapunov');


h33 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_plot_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[5 3.69230769230769 30.8 1.76923076923077],...
'String','Plot Controller Regions',...
'Tag','button_plot');


h34 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_plotU_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[39.2 3.69230769230769 30.8 1.76923076923077],...
'String','Plot Control Law',...
'Tag','button_plotU');


h35 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_plotJ_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[73.6 3.69230769230769 30.8 1.76923076923077],...
'String','Plot Value Function',...
'Tag','button_plotJ');


h36 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[4.2 25.4615384615385 33.4 13.0769230769231],...
'String','No model defined',...
'Style','text',...
'Tag','text_systeminfo');


h37 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_generate_Callback'',gcbo,[],guidata(gcbo))',...
'FontWeight','bold',...
'ForegroundColor',[1 0 0],...
'ListboxTop',0,...
'Position',[73.4 20.6153846153846 29 2.15384615384615],...
'String','Generate Controller',...
'Tag','button_generate');


h38 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_weights_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[52.2 25 21.4 1.92307692307692],...
'String','Define Penalties',...
'Tag','button_weights');


h39 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[39.8 18.6923076923077 29.6 1.23076923076923],...
'String','Controller Details',...
'Style','text',...
'Tag','text59');


h40 = uimenu(...
'Parent',h1,...
'Callback','mpt_guiProbStruct_export(''menu_file_Callback'',gcbo,[],guidata(gcbo))',...
'Label','File',...
'Tag','menu_file');

h41 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_newsession_Callback'',gcbo,[],guidata(gcbo))',...
'Label','New session...',...
'Tag','menu_file_newsession');

h42 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_loadfile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load setup from file...',...
'Separator','on',...
'Tag','menu_file_loadfile');

h43 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_savefile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save setup to file...',...
'Tag','menu_file_savefile');

h44 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_savesession_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save setup to workspace...',...
'Tag','menu_file_savesession');

h45 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_loadw_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load controller from workspace...',...
'Separator','on',...
'Tag','menu_file_loadw');

h46 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_savew_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save controller to workspace...',...
'Tag','menu_file_savew');

h47 = uimenu(...
'Parent',h40,...
'Callback','mpt_guiProbStruct_export(''menu_file_exit_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Exit',...
'Separator','on',...
'Tag','menu_file_exit');

h48 = uimenu(...
'Parent',h1,...
'Callback','mpt_guiProbStruct_export(''menu_model_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Model',...
'Tag','menu_model');

h49 = uimenu(...
'Parent',h48,...
'Callback','mpt_guiProbStruct_export(''menu_model_edit_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Edit...',...
'Tag','menu_model_edit');

h50 = uimenu(...
'Parent',h48,...
'Callback','mpt_guiProbStruct_export(''menu_model_loadw_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load from workspace...',...
'Separator','on',...
'Tag','menu_model_loadw');

h51 = uimenu(...
'Parent',h48,...
'Callback','mpt_guiProbStruct_export(''menu_model_loadfile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load from file...',...
'Tag','menu_model_loadfile');

h52 = uimenu(...
'Parent',h48,...
'Callback','mpt_guiProbStruct_export(''menu_model_savew_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save to workspace...',...
'Separator','on',...
'Tag','menu_model_savew');

h53 = uimenu(...
'Parent',h48,...
'Callback','mpt_guiProbStruct_export(''menu_model_savefile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save to file...',...
'Tag','menu_model_savefile');

h54 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''button_simulate_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[39.2 1.15384615384615 30.8 1.76923076923077],...
'String','Simulations',...
'Tag','button_simulate');


h55 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'FontWeight','bold',...
'ListboxTop',0,...
'Position',[36.8 22.6153846153846 35.6 1.46153846153846],...
'String','Controller Computation',...
'Style','text',...
'Tag','text60');


h56 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbStruct_export(''textfield_ctrlname_Callback'',gcbo,[],guidata(gcbo))',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[30.6 20.8461538461538 38.4 1.61538461538462],...
'String','ctrl',...
'Style','edit',...
'CreateFcn','mpt_guiProbStruct_export(''textfield_ctrlname_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_ctrlname');


h57 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[4 21.0769230769231 25.4 1.15384615384615],...
'String','Export to workspace as',...
'Style','text',...
'Tag','text61');


h58 = uimenu(...
'Parent',h1,...
'Callback','mpt_guiProbStruct_export(''menu_problem_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Problem',...
'Tag','menu_problem');

h59 = uimenu(...
'Parent',h58,...
'Callback','mpt_guiProbStruct_export(''menu_problem_reset_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Reset...',...
'Tag','menu_problem_reset');

h60 = uimenu(...
'Parent',h58,...
'Callback','mpt_guiProbStruct_export(''menu_problem_loadw_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load from workspace...',...
'Separator','on',...
'Tag','menu_problem_loadw');

h61 = uimenu(...
'Parent',h58,...
'Callback','mpt_guiProbStruct_export(''menu_problem_loadfile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Load from file...',...
'Tag','menu_problem_loadfile');

h62 = uimenu(...
'Parent',h58,...
'Callback','mpt_guiProbStruct_export(''menu_problem_savew_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save to workspace...',...
'Separator','on',...
'Tag','menu_problem_savew');

h63 = uimenu(...
'Parent',h58,...
'Callback','mpt_guiProbStruct_export(''menu_problem_savefile_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save to file...',...
'Tag','menu_problem_savefile');

h64 = uimenu(...
'Parent',h1,...
'Callback','mpt_guiProbStruct_export(''menu_examples_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Examples',...
'Tag','menu_examples');

h65 = uimenu(...
'Parent',h64,...
'Callback','mpt_guiProbStruct_export(''menu_ex_di_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Double Integrator',...
'Tag','menu_ex_di');

h66 = uimenu(...
'Parent',h64,...
'Callback','mpt_guiProbStruct_export(''menu_ex_3dlti_Callback'',gcbo,[],guidata(gcbo))',...
'Label','3rd order LTI system',...
'Tag','menu_ex_3dlti');

h67 = uimenu(...
'Parent',h64,...
'Callback','mpt_guiProbStruct_export(''menu_ex_pwasincos_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Periodical PWA system',...
'Tag','menu_exa_pwasincos');

h68 = uimenu(...
'Parent',h64,...
'Callback','mpt_guiProbStruct_export(''menu_ex_pwacar_Callback'',gcbo,[],guidata(gcbo))',...
'Label','PWA car',...
'Tag','menu_ex_pwacar');

h69 = uimenu(...
'Parent',h64,...
'Callback','mpt_guiProbStruct_export(''menu_ex_twotanks_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Two tanks system',...
'Tag','menu_ex_twotanks');

h70 = uimenu(...
'Parent',h1,...
'Callback','mpt_guiProbStruct_export(''menu_help_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Help',...
'Tag','menu_help');

h71 = uimenu(...
'Parent',h70,...
'Callback','mpt_guiProbStruct_export(''menu_help_gotompt_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Help on MPT homepage',...
'Tag','menu_help_gotompt');

h72 = uimenu(...
'Parent',h70,...
'Callback','mpt_guiProbStruct_export(''menu_help_credits_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Credits...',...
'Separator','on',...
'Tag','menu_help_credits');

h73 = uimenu(...
'Parent',h70,...
'Callback','mpt_guiProbStruct_export(''menu_help_about_Callback'',gcbo,[],guidata(gcbo))',...
'Label','About...',...
'Separator','on',...
'Tag','menu_help_about');

h74 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbStruct_export(''text_cdetails_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','inactive',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[2.4 10.3076923076923 105 8.30769230769231],...
'String','No controller loaded',...
'Style','text',...
'Value',1,...
'CreateFcn','mpt_guiProbStruct_export(''text_cdetails_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','text_cdetails');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUIPROBSTRUCT_EXPORT, by itself, creates a new MPT_GUIPROBSTRUCT_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIPROBSTRUCT_EXPORT returns the handle to a new MPT_GUIPROBSTRUCT_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIPROBSTRUCT_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIPROBSTRUCT_EXPORT.M with the given input arguments.
%
%      MPT_GUIPROBSTRUCT_EXPORT('Property','Value',...) creates a new MPT_GUIPROBSTRUCT_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.5 $ $Date: 2003/07/17 18:28:28 $

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
    % MPT_GUIPROBSTRUCT_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUIPROBSTRUCT_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUIPROBSTRUCT_EXPORT(...)
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

