function varargout = mpt_guiProbStruct(varargin)
% MPT_GUIPROBSTRUCT M-file for mpt_guiProbStruct.fig
%      MPT_GUIPROBSTRUCT, by itself, creates adyn new MPT_GUIPROBSTRUCT or raises the existing
%      singleton*.
%
%      H = MPT_GUIPROBSTRUCT returns the handle to adyn new MPT_GUIPROBSTRUCT or the handle to
%      the existing singleton*.
%
%      MPT_GUIPROBSTRUCT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIPROBSTRUCT.M with the given input arguments.
%
%      MPT_GUIPROBSTRUCT('Property','Value',...) creates adyn new MPT_GUIPROBSTRUCT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiProbStruct_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiProbStruct_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiProbStruct

% Last Modified by GUIDE v2.5 04-Mar-2005 10:49:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiProbStruct_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiProbStruct_OutputFcn, ...
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


% --- Executes just before mpt_guiProbStruct is made visible.
function mpt_guiProbStruct_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in adyn future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiProbStruct (see VARARGIN)

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
    

% Choose default command line output for mpt_guiProbStruct
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiProbStruct wait for user response (see UIRESUME)
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
            if res==0,
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
function varargout = mpt_guiProbStruct_OutputFcn(hObject, eventdata, handles)
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


