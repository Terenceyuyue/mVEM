function varargout = mpt_guiProbAdv_export(varargin)
% mpt_guiProbAdv_export M-file for mpt_guiProbAdv_export.fig
%      mpt_guiProbAdv_export, by itself, creates a new mpt_guiProbAdv_export or raises the existing
%      singleton*.
%
%      H = mpt_guiProbAdv_export returns the handle to a new mpt_guiProbAdv_export or the handle to
%      the existing singleton*.
%
%      mpt_guiProbAdv_export('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in mpt_guiProbAdv_export.M with the given input arguments.
%
%      mpt_guiProbAdv_export('Property','Value',...) creates a new mpt_guiProbAdv_export or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_guiProbAdv_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_guiProbAdv_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt_guiProbAdv_export

% Last Modified by GUIDE v2.5 04-Apr-2005 11:19:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_guiProbAdv_export_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_guiProbAdv_export_OutputFcn, ...
                   'gui_LayoutFcn',  @mpt_guiProbAdv_export_LayoutFcn, ...
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


% --- Executes just before mpt_guiProbAdv_export is made visible.
function mpt_guiProbAdv_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mpt_guiProbAdv_export (see VARARGIN)


global mpt__probStruct

mpt_guiCenter(hObject);

% Choose default command line output for mpt_guiProbAdv_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mpt_guiProbAdv_export wait for user response (see UIRESUME)
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
function varargout = mpt_guiProbAdv_export_OutputFcn(hObject, eventdata, handles)
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

% --- Creates and returns a handle to the GUI figure. 
function h1 = mpt_guiProbAdv_export_LayoutFcn(policy)
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
'Name','Advanced Problem Options',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[50 48.3076923076923 72.6 24.7692307692308],...
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
'text', 11, ...
'edit', 7, ...
'checkbox', 6, ...
'popupmenu', 4, ...
'frame', 2, ...
'pushbutton', 12), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\MatlabFiles\mpt\mpt\extras\gui\mpt_guiProbAdv.m'));


h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Enable','off',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[2.6 17.1538461538462 20.8 1.15384615384615],...
'String','Penalty on final state',...
'Style','text',...
'Tag','text_PN');


h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbAdv_export(''popup_Tconstraint_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[25.2 21.9230769230769 38.2 1.69230769230769],...
'String',{  'Stabilizing terminal set'; 'User defined terminal set'; 'No terminal set' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn','mpt_guiProbAdv_export(''popup_Tconstraint_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','popup_Tconstraint');


h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbAdv_export(''textfield_Tset_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[25.2 19.4615384615385 38.4 1.84615384615385],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiProbAdv_export(''textfield_Tset_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_Tset');


h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbAdv_export(''textfield_PN_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[25.2 16.8461538461538 38.4 1.84615384615385],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiProbAdv_export(''textfield_PN_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_PN');


h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''checkbox_y0bounds_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[25.2 14.6923076923077 23 1.30769230769231],...
'String','Bounds on y(0)',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_y0bounds');


h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''checkbox_feedback_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[25.2 12.3846153846154 35.2 1.30769230769231],...
'String','Prestabilization with Feedback',...
'Style','checkbox',...
'Tag','checkbox_feedback');


h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[1.2 22.1538461538462 22.2 1.15384615384615],...
'String','Terminal set constraint',...
'Style','text',...
'Tag','text5');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Enable','off',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[13.4 19.7692307692308 10 1.15384615384615],...
'String','Target Set',...
'Style','text',...
'Tag','text_Tset');


h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Enable','off',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[2.6 9.92307692307692 20.8 1.15384615384615],...
'String','Feedback law',...
'Style','text',...
'Tag','text_FBgain');


h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbAdv_export(''textfield_FBgain_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[25.2 9.61538461538462 38.4 1.84615384615385],...
'String','LQR',...
'Style','edit',...
'CreateFcn','mpt_guiProbAdv_export(''textfield_FBgain_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_FBgain');


h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','mpt_guiProbAdv_export(''textfield_Nc_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[25.2 6.84615384615385 38.4 1.84615384615385],...
'String','',...
'Style','edit',...
'CreateFcn','mpt_guiProbAdv_export(''textfield_Nc_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','textfield_Nc',...
'UserData',[]);


h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''checkbox_deltaTracking_Callback'',gcbo,[],guidata(gcbo))',...
'Enable','off',...
'ListboxTop',0,...
'Position',[25.2 4.61538461538462 35 1.23076923076923],...
'String','Tracking with delta U formulation',...
'Style','checkbox',...
'Value',1,...
'Tag','checkbox_deltaTracking');


h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''Return_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[26 1.23076923076923 23.4 1.92307692307692],...
'String','Return',...
'Tag','Return');


h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'CData',[],...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[2.6 7.15384615384615 20.8 1.15384615384615],...
'String','Control horizon',...
'Style','text',...
'Tag','text_Nc',...
'UserData',[]);


h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_Tconstraint_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.4 21.8461538461538 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_Tconstraint');


h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_Tset_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.4 19.4615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_Tset');


h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_PN_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.4 16.8461538461538 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_PN');


h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_y0bounds_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.4 14.4615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_y0bounds');


h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_feedback_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.4 12.1538461538462 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_feedback');


h21 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_FBgain_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.2 9.61538461538462 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_FBgain');


h22 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''help_Nc_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.2 6.84615384615385 3.4 1.76923076923077],...
'String','?',...
'Tag','help_Nc');


h23 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','mpt_guiProbAdv_export(''button_help_tracking2_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[66.2 4.30769230769231 3.4 1.76923076923077],...
'String','?',...
'Tag','button_help_tracking2');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      MPT_GUIPROBADV_EXPORT, by itself, creates a new MPT_GUIPROBADV_EXPORT or raises the existing
%      singleton*.
%
%      H = MPT_GUIPROBADV_EXPORT returns the handle to a new MPT_GUIPROBADV_EXPORT or the handle to
%      the existing singleton*.
%
%      MPT_GUIPROBADV_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT_GUIPROBADV_EXPORT.M with the given input arguments.
%
%      MPT_GUIPROBADV_EXPORT('Property','Value',...) creates a new MPT_GUIPROBADV_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.3 $ $Date: 2005/04/04 09:24:51 $

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
    % MPT_GUIPROBADV_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % MPT_GUIPROBADV_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % MPT_GUIPROBADV_EXPORT(...)
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

