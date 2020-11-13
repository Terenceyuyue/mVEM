function mpt_plotLyapunov(ctrl, meshgridpoints)
%MPT_PLOTLYAPUNOV Plots lyapunov function stored in a controlle
%
% mpt_plotLyapunov(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the polyhedral parition obtained as a solution of the 
% optimal control problem obtained by mpt_control
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl             - Explicit controller (MPTCTRL object)
% meshgridpoints   - number of grid points in one axis,
%                    (default: 30)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% see also MPT_PLOTPWA, MPT_PLOTPWQ
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

error(nargchk(1,2,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin < 2
    meshgridpoints = 30;
end

if ~isa(ctrl, 'mptctrl'),
    error('Input must be an MPTCTRL object.');
end
if ~isexplicit(ctrl),
    error('On-line controllers not supported by this function.');
end

if ~isfield(ctrl.details, 'lyapunov'),
    ctrlname = inputname(1);
    if ~isempty(ctrlname),
        fprintf('No Lyapunov function data stored in the object, use mpt_lyapunov(%s) first.\n', ctrlname);
    else
        fprintf('No Lyapunov function data stored in the object, use mpt_lyapunov first.\n');
    end
    return
end

Pn = ctrl.Pn;
if dimension(Pn)>2,
    error('Only Lyapunov functions defined over 2D partition can be plotted.');
end
lyap = ctrl.details.lyapunov;

switch lyap.type
    case 'quadratic',
        q = lyap.P;
        l = zeros(dimension(Pn), 1);
        c = 0;
        lenPn = length(Pn);
        Q = cell(lenPn, 1);
        L = cell(lenPn, 1);
        C = cell(lenPn, 1);
        for ii = 1:length(Pn),
            Q{ii} = q;
            L{ii} = l;
            C{ii} = c;
        end
        mpt_plotPWQ(Pn, Q, L, C, meshgridpoints);
        title('Lyapunov function');
        
    case 'pwq',
        Q = lyap.Q;
        L = lyap.L;
        C = lyap.C;
        mpt_plotPWQ(Pn, Q, L, C, meshgridpoints);
        title('Lyapunov function');
        
    case 'pwa',
        L = lyap.L;
        C = lyap.C;
        mpt_plotPWA(Pn, L, C, struct('showPn', 1));
        title('Lyapunov function');
        
    otherwise
        % cannot handle PWP and SOS lyapunov functions yet
        error('Unsopported type of Lyapunov function.');
end
