function [nph,nfh] = sbprogress(ph,x,msg,fh,title,modal)
%[nph,nfh] = function SBPROGRESS(ph,x,msg,fh,title,modal)
%   Creates a progress bar in the "status-bar section" of the figure.
%
%   - ph:   progress-bar axis handle; empty if one is to be created
%   - x:    can have a numerical value (-100..0..100) or any of the two
%           strings: 'close' or 'redraw'. Progress (x) with negative values means that
%           the percent number will not be displayed.
%           x given as the 'close' string closes an existing progress bar (ph) (and corresponding figure
%           if figure was created when sbprogress was called otherwise only  the progress bar will be removed from the figure).
%           x given as the 'redraw' string redraws an existing ph (progress
%           bar) - useful when figure changes its size for example.
%   - msg:  msg text (string)
%   - fh:   (optional) figure handle that progress bar lies on (if empty, a
%           new figure will be created)
%   - title:title when a new figure is raised (when fh is empty or
%           non-existing).
%   - modal:in case a new figure is raised, it can be raised as modal
%           window (default = 0).
%
%   - nph:  a handle to a newly created or existing progress-bar axis
%   - nfh:  a handle to a newly created or existing figure
%
%   EXAMPLES:   
%   sbprogress
%           Runs a demo
%   sbprogress([],10,'Progress:',[],'Test progress',0)
%           Creates a new figure with a progress bar...
%   ph = sbprogress([],-10,'Progress:',[],'Test progress',0)
%           Creates a new figure with a progress bar where the percent
%           number is not shown; a handle array to the progress-bar axis is returned
%   sbprogress(ph,'redraw')
%           Redraws the ph progress bar
%   sbprogress(ph,20)
%           Re-sets the progress-bar's progress
%   sbprogress([],20,'Test',fh)
%           Creates a progress bar on an existing figure fh
%   sbprogress(ph,'redraw')
%           Redraws an existing progress bar designated by ph handle-array
%
%   NOTES: when the progress-bar's percentage is changed, the bar will be redrawn too. 
%
%   Primoz Cermelj, 08.12.2002
%   (c) Primoz Cermelj, primoz.cermelj@email.si
%
%   Last revision: 28.01.2004
%--------------------------------------------------------------------------
global W H

W = 400;
H = 20;

if nargin == 0
    % Run the demo
    [nph,nfh] = sbprogress([],48.5,'Test:',[],'SBPROGRESS demo',0);
    return
elseif nargin < 2
    error('At least progress length should be given or ''close'' string');
end
if ~ishandle(ph)
    error('ph is not a handle');
end

if ischar(x)
    if strcmpi(x,'close')
        Data = get(ph,'UserData');
        if ~isempty(Data) && isfield(Data,'fhCreated') && isfield(Data,'fh') && isfield(Data,'frame')
            if Data.fhCreated
                delete(Data.fh);
            else
                delete(Data.frame);
                delete(ph);
            end
        else
            error('Wrong ph handle given');
        end
    elseif strcmpi(x,'redraw')
        Data = get(ph,'UserData');
        if ~isempty(Data) && isfield(Data,'fhCreated') && isfield(Data,'fh') && isfield(Data,'frame')
            updateprogress(ph,Data.x,Data.msg);
        else
            error('Wrong ph handle given');
        end
    else
        error('Use sbprogress(ph,''redraw'')');
    end
    return
end
phExists = 1;
fhExists = 1;
if nargin < 2 || isempty(ph) || ~ishandle(ph)
    phExists = 0;
else
    if ~ishandle(ph)
        phExists = 0;
    end        
end
if nargin < 4 || isempty(fh) || ~ishandle(fh)
    fhExists = 0;
else
    if ~ishandle(ph)
        phExists = 0;
    end     
end
if nargin < 5 || isempty(title)
    title = '';
end
if nargin < 6 || isempty(modal)
    modal = 0;
end

if ~phExists
    if nargin < 3 || isempty(msg)
        msg = '';
    end
    fhCreated = 0;
    if ~fhExists
       nfh = createnewfig(title,modal); 
       fh = nfh;
       fhCreated = 1;
    end
    nph = createnewprogress(msg,fh,fhCreated,x);
    ph = nph;
else
    if nargin < 3
        updateprogress(ph,x);
    else
        updateprogress(ph,x,msg);
    end    
end

if nargout > 0
    nph = ph;
    if nargout > 1
        nfh = fh;
    end
end



%-------------------------------------
function fh = createnewfig(title,modal)
global W H
set(0,'Units','Pixels');
scr = get(0,'ScreenSize');
fh = figure('Visible','off');
set(fh,'Units','Pixels',...
       'Tag','progressBarFigure',...
       'Position',[scr(3)/2-W/2 scr(4)/2-H/2 W H],...
       'MenuBar','none',...
       'Resize','off',...
       'NumberTitle','off',...
       'Name',title);
if modal
    set(fh,'WindowStyle','modal');
end
set(fh,'Visible','on');



%-------------------------------------
function [ph,hMsg] = createnewprogress(msg,fh,fhCreated,x)
global H

[fPos,pPos,framePos] = getcoords(fh);
w = pPos(3);
ph = axes;
set(ph,'Tag','progressBarAxes',...
       'Units','Pixels',...
       'Position',pPos,...
       'Box','off',...
       'Visible','off',...
       'XLim',[0 1],'YLim',[0 1]);
hMsg = text(1,0.45*H,msg,'Units','Pixels',...
                         'Tag','progressBarMsg');
tPos = get(hMsg,'Extent');
if (tPos(3)+tPos(1)) > 0.5*w
    pStart = 0.5;
else
    pStart = (tPos(3)+tPos(1)+10)/w;
end
Data.pStart = pStart;
Data.fh = fh;
Data.fhCreated = fhCreated;
if x < -100; x = -100; end;
if x > 100; x = 100; end;
Data.x = x;
Data.msg = msg;
set(ph,'UserData',Data);
pEnd = pStart + abs(x)/100*abs(1-pStart);
hRec = rectangle('Position',[pStart 0 1-pStart 1],'Tag','progressBarRec');
hPatch = patch([pStart; pEnd; pEnd; pStart],[0; 0; 1; 1],'b',...
                            'Tag','progressBarPatch');
hPercTxt = text( (pStart+1)/2,0.5,[num2str(abs(x)) '%'],'Units','Normalized','FontWeight','bold','Color','b',...
                            'EraseMode','xor','HorizontalAlignment','Center',...
                            'Tag','progressBarPercTxt','Visible','Off');
if x >= 0
    set(hPercTxt,'Visible','On');
end
frame = panel(fh,framePos);
Data = get(ph,'UserData');
Data.frame = frame;
set(ph,'UserData',Data);
drawnow;


%-------------------------------------
function updateprogress(ph,x,newMsg)
global H

Data = get(ph,'UserData');
fh = Data.fh;
hMsg = findobj(ph,'Tag','progressBarMsg');
if nargin > 2  & ~isempty(newMsg) 
    set(hMsg,'String',newMsg);
    Data.msg = newMsg;
end
[fPos,pPos,framePos] = getcoords(fh);
if x < -100; x = 0; end;
if x > 100; x = 100; end;

set(ph,'Position',pPos);
set(hMsg,'Position',[1 0.45*H]);
tPos = get(hMsg,'Extent');
w = pPos(3);
if (tPos(3)+tPos(1)) > 0.5*w
    pStart = 0.5;
else
    pStart = (tPos(3)+tPos(1)+10)/w;
end

Data.pStart = pStart;
Data.x = x;
set(ph,'UserData',Data);
pEnd = pStart + abs(x)/100*abs(1-pStart);
hRec = findobj(ph,'Tag','progressBarRec');
hPatch = findobj(ph,'Tag','progressBarPatch');
hPercTxt = findobj(ph,'Tag','progressBarPercTxt');
%
xf = framePos(1);
yf = framePos(2);
wf = framePos(3);
hf = framePos(4);
set(Data.frame(1),'Position',[xf yf+hf wf 1]);
set(Data.frame(2),'Position',[xf yf 1 hf]);
set(Data.frame(3),'Position',[xf+1 yf wf-1 1]);
set(Data.frame(4),'Position',[xf+wf-1 yf 1 hf]);
%
set(hRec,'Position',[pStart 0 1-pStart 1]);
set(hPatch,'XData',[pStart; pEnd; pEnd; pStart],'YData',[0; 0; 1; 1]);
set(hPercTxt,'Position',[(pStart+1)/2 0.5],'String',[num2str(abs(x)) '%']);
if x >= 0
    set(hPercTxt,'Visible','On');
else
    set(hPercTxt,'Visible','Off');
end
drawnow;
set(ph,'UserData',Data);



%-------------------------------------
function frame = panel(fh,pos);
% Creates a virtual panel surrounding the progress bar.
% See also the self-contained function BEVEL.
x = pos(1);
y = pos(2);
w = pos(3);
h = pos(4);
col = rgb2hsv(get(fh,'Color'));
downColor = col;   downColor(2) = 0.5*downColor(2);  downColor(3) = 0.9; downColor = hsv2rgb(downColor);
upColor = col;    upColor(3) = 0.4;  upColor = hsv2rgb(upColor);
frame = zeros(4,1);
frame(1) = uicontrol('Units','Pixels','Style','Frame','Position',[x y+h w 1],'ForeGroundColor',upColor,'BackgroundColor',upColor);
frame(2) = uicontrol('Units','Pixels','Style','Frame','Position',[x y 1 h],'ForeGroundColor',upColor,'BackgroundColor',upColor);
frame(3) = uicontrol('Units','Pixels','Style','Frame','Position',[x+1 y w-1 1],'ForeGroundColor',downColor,'BackgroundColor',downColor);
frame(4) = uicontrol('Units','Pixels','Style','Frame','Position',[x+w-1 y 1 h],'ForeGroundColor',downColor,'BackgroundColor',downColor);


%-------------------------------------
function [fPos,pPos,framePos] = getcoords(fh)
% Returns coordinates for axis, frame, and msg for our sbprogress (for
% creation or updating needs).
global H

% Progress-bar axis
figure(fh);
fUnits = get(fh,'Units');
set(fh,'Units','Pixels');
fPos = get(fh,'Position');
set(fh,'Units',fUnits);
pPos = fPos;
pPos(1) = 5;
pPos(2) = 4;
pPos(3) = pPos(3)-8;
pPos(4) = H-6;
% Frame
framePos = fPos;
framePos(1) = 2;
framePos(2) = 2;
framePos(3) = fPos(3)-3;
framePos(4) = H-3;
