function varargout = mpt_guiProbAdv(varargin)
% mpt_guiProbAdv M-file for mpt_guiProbAdv.fig
%      mpt_guiProbAdv, by itself, creates a new mpt_guiProbAdv or raises the existing
%      singleton*.
%
%      H = mpt_guiProbAdv returns the handle to a new mpt_guiProbAdv or the handle to
%      the existing singleton*.
%
%      mpt_guiProbAdv('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in mpt_guiProbAdv.M with the given input arguments.
%
%      mpt_guiProbAdv('Property','Value',...) creates a new mpt_guiProbAdv or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiProbAdv_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiProbAdv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiProbAdv

% Last Modified by GUIDE v2.5 04-Apr-2005 11:04:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiProbAdv_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiProbAdv_OutputFcn, ...
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


% --- Executes just before mpt_guiProbAdv is made visible.
function mpt_guiProbAdv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiProbAdv (see VARARGIN)


global mpt__probStruct

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiProbAdv
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiProbAdv wait for user response (see UIRESUME)
% uiwait(handles.figure1);

probStruct = mpt__probStruct;
if isfield(probStruct, 'Tconstraint')
    switch probStruct.Tconstraint,
        case 0,
            set(handles.popup_Tconstraint, 'Value', 3);
            if isfield(probStruct, 'P_N'),
                set(handles.textfield_PN, 'String', mat2str(probStruct.P_N));
            end
            set(handles.textfield_PN, 'Enable', 'on');
            set(handles.text_PN, 'Enable', 'on');
            
        case 1,
            set(handles.popup_Tconstraint, 'Value', 1);
            
        case 2,
            set(handles.popup_Tconstraint, 'Value', 2);
            if isfield(probStruct, 'Tset'),
                Tset = probStruct.Tset;
                if isa(Tset, 'polytope'),
                    [H,K] = double(Tset);
                    set(handles.textfield_Tset, 'String', ...
                        sprintf('polytope(%s, %s)', mat2str(H), mat2str(K)));
                end
                set(handles.textfield_Tset, 'Enable', 'on');
                set(handles.textfield_PN, 'Enable', 'on');
                set(handles.text_PN, 'Enable', 'on');
                set(handles.text_Tset, 'Enable', 'on');
            end
    end
end
if isfield(probStruct, 'y0bounds'),
    if probStruct.y0bounds,
        set(handles.checkbox_y0bounds, 'Value', 1);
    else
        set(handles.checkbox_y0bounds, 'Value', 0);
    end
end

if isfield(probStruct, 'feedback')
    if probStruct.feedback,
        set(handles.checkbox_feedback, 'Value', 1);
    else
        set(handles.checkbox_feedback, 'Value', 0);
    end
end

if isfield(probStruct, 'FBgain'),
    set(handles.textfield_FBgain, 'String', mat2str(probStruct.FBgain));
end

set(handles.checkbox_deltaTracking, 'Enable', 'off');
if isfield(probStruct, 'tracking')
    if probStruct.tracking>0,
        set(handles.checkbox_deltaTracking, 'Enable', 'on');
    end
end

if isfield(probStruct, 'Nc'),
    set(handles.textfield_Nc, 'String', mat2str(probStruct.Nc));
end

% --- Outputs from this function are returned to the command line.
function varargout = mpt_guiProbAdv_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function textfield_PN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_PN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_PN_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_PN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_PN as text
%        str2double(get(hObject,'String')) returns contents of textfield_PN as a double

global mpt__probStruct
try
    mpt__probStruct.P_N = evalin('base', get(hObject, 'String'));
catch
    sub_errordlg(lasterr);
end

% --- Executes on button press in checkbox_y0bounds.
function checkbox_y0bounds_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_y0bounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_y0bounds

global mpt__probStruct
try
    mpt__probStruct.y0bounds = get(hObject, 'Value');
catch
    sub_errordlg(lasterr);
end


% --- Executes on button press in checkbox_feedback.
function checkbox_feedback_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_feedback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_feedback

global mpt__probStruct
value = get(hObject, 'Value');

if value==1,
    % checkbox_feedback gain requested, enable textfield_FBgain field
    set(handles.textfield_FBgain, 'Enable', 'on');
    set(handles.text_FBgain, 'Enable', 'on');
else
    % disable textfield_FBgain field
    set(handles.textfield_FBgain, 'Enable', 'off');
    set(handles.text_FBgain, 'Enable', 'off');
end

try
    mpt__probStruct.feedback = value;
catch
    sub_errordlg(lasterr);
end


% --- Executes during object creation, after setting all properties.
function popup_Tconstraint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Tconstraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popup_Tconstraint.
function popup_Tconstraint_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Tconstraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_Tconstraint contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Tconstraint

global mpt__probStruct
value = get(hObject,'Value');

switch value,
    case 1
        % stabilizing terminal set, disable textfield_PN and textfield_Tset fields
        set(handles.textfield_Tset, 'Enable', 'off');
        set(handles.textfield_PN, 'Enable', 'off');
        set(handles.text_Tset, 'Enable', 'off');
        set(handles.text_PN, 'Enable', 'off');
        mpt__probStruct.Tconstraint = 1;
    case 2
        % user defined set, enable textfield_PN and textfield_Tset fields
        set(handles.textfield_Tset, 'Enable', 'on');
        set(handles.textfield_PN, 'Enable', 'on');
        set(handles.text_Tset, 'Enable', 'on');
        set(handles.text_PN, 'Enable', 'on');
        mpt__probStruct.Tconstraint = 2;
    case 3
        % no terminal set, enable textfield_PN, disable textfield_Tset
        set(handles.textfield_Tset, 'Enable', 'off');
        set(handles.text_Tset, 'Enable', 'off');
        set(handles.textfield_PN, 'Enable', 'on');
        set(handles.text_PN, 'Enable', 'on');
        mpt__probStruct.Tconstraint = 0;
end

% --- Executes during object creation, after setting all properties.
function textfield_Tset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_Tset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_Tset_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_Tset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_Tset as text
%        str2double(get(hObject,'String')) returns contents of textfield_Tset as a double

global mpt__probStruct

value = get(hObject, 'String');

if isempty(value),
    if isfield(mpt__probStruct, 'Tset'),
        mpt__probStruct.Tset = polytope;
    end
    return
end

try
    Tset = evalin('base', value);
catch
    sub_errordlg(lasterr);
    return
end

if ~isa(Tset, 'polytope')
    errordlg('Target set must be a POLYTOPE object!', 'Error', 'modal');
    return
end
mpt__probStruct.Tset = Tset;

% --- Executes during object creation, after setting all properties.
function textfield_FBgain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_FBgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_FBgain_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_FBgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_FBgain as text
%        str2double(get(hObject,'String')) returns contents of textfield_FBgain as a double

global mpt__probStruct
value = get(hObject, 'String');
if strcmpi(value, 'lqr'),
    if isfield(mpt__probStruct, 'FBgain'),
        mpt__probStruct = rmfield(mpt__probStruct, 'FBgain');
    end
    return
end

try
    mpt__probStruct.FBgain = evalin('base', get(hObject, 'String'));
catch
    sub_errordlg(lasterr);
end


% --- Executes on button press in Return.
function Return_Callback(hObject, eventdata, handles)
% hObject    handle to Return (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mpt__probStruct mpt___sysStruct

try
    Options.verbose = 0;
    Options.guierrors = 1;
    Options.forceverify = 1;
    T = evalc('mpt_verifySysProb(mpt___sysStruct, mpt__probStruct, Options);');
catch
    sub_errordlg(lasterr);
    return
end

close(handles.figure1);

% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function popup_moveb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_moveb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function textfield_Nc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textfield_Nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function textfield_Nc_Callback(hObject, eventdata, handles)
% hObject    handle to textfield_Nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textfield_Nc as text
%        str2double(get(hObject,'String')) returns contents of textfield_Nc as a double

global mpt__probStruct

try
    valstr = get(hObject, 'String');
    if isempty(valstr),
        valstr = '[]';
    end
    value = evalin('base', valstr);
catch
    sub_errordlg(lasterr);
    return
end

if isempty(value),
    mpt__probStruct = rmfield(mpt__probStruct, 'Nc');
else
    mpt__probStruct.Nc = value;
end

% --- Executes on button press in button_help_Tconstraint.
function button_help_Tconstraint_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_Tconstraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg([sprintf('Which target set constraint to apply on the final state x(n)\n\n'), ...
        sprintf('* Stabilizing target set - enforces x(N) \\in LQRset (for 2-norm problems)\n\n'), ...
        sprintf('* User-defined target set - enforces x(N) \\in Tset, where Tset must be defined\n\n'), ...
        sprintf('* No target set - no target set constraint on x(N)')]);

% --- Executes on button press in button_help_Tset.
function button_help_Tset_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_Tset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg([sprintf('Defines the taget set which is to be applied as a set constraints on final state x(N), i.e.:\n\n'), ...
        sprintf('x(n) \\in Tset\n\n'), ...
        sprintf('The input must be a polytope object, e.g. ''unitbox(2, 1)''.')]);


% --- Executes on button press in button_help_PN.
function button_help_PN_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_PN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg([sprintf('Defines penalty P_N on final state x(N) in the cost function:\n\n'), ...
        sprintf('J = sum( x(k)''*Q*x(k) + u(k)''*R*u(k) + x(N)''*P_N*x(N)\n\n'), ...
        'Note: Not applicable if the "Terminal set" option is set to "Stabilizing target set".']);

% --- Executes on button press in button_help_y0bounds.
function button_help_y0bounds_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_y0bounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg('If disabled, cosntraints on outputs will be imposed only on y(1)...y(N) and not on y(0)');

% --- Executes on button press in button_help_feedback.
function button_help_feedback_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_feedback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(['If enabled, the control problem is augmnted such that "U = K*x + c", where "K" is '...
        'a state-feedback gain (typically an LQR controller) and the optimization is is performed '...
        sprintf('over the offset "c"\n\nNote: Only available for quadratic cost functions and LTI systems.')]);


% --- Executes on button press in button_help_FBgain.
function button_help_FBgain_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_FBgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(['If "Prestabilization with feedback" is enabled, a specific state-feedback gain '...
        'can be provided in this field.']);


% --- Executes on button press in help_moveb.
function help_moveb_Callback(hObject, eventdata, handles)
% hObject    handle to help_moveb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in help_Nc.
function help_Nc_Callback(hObject, eventdata, handles)
% hObject    handle to help_Nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(['Specifies number of free control moves, i.e. number of degrees of freedom in the optimization procedure.']);


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


% --- Executes on button press in checkbox_deltaTracking.
function checkbox_deltaTracking_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_deltaTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_deltaTracking

global mpt__probStruct

value = get(hObject, 'Value');

if isfield(mpt__probStruct, 'tracking')
    if mpt__probStruct.tracking>0,
        if value==1,
            mpt__probStruct.tracking = 1;
        else
            mpt__probStruct.tracking = 2;
        end
    end
end
            
        

% --- Executes on button press in button_help_tracking2.
function button_help_tracking2_Callback(hObject, eventdata, handles)
% hObject    handle to button_help_tracking2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sub_helpdlg(['By default, if a tracking controller is required, the state-space vector has to be '...
    'augmented to introduce past inputs. This can have a significant impact on run-time and '...
    'solution complexity. If you, however, know that zero input is an equilibrium for each possible '...
    'reference, set this option to False, thus greatly reducing computational complexity.' ...
    sprintf('\n\nNote: Offset-free tracking cannot be guaranteed if this option is set to FALSE!')]);



%============================================================
function sub_helpdlg(msg)

msgbox(msg, 'Help', 'modal');