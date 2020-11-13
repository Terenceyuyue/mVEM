function result = isstabilizable(ctrl)
%STABILIZES Decides if a given controller is stabilizable
%
%   result = isstabilizable(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns 1 if a given explicit controller is stabilizable, i.e. it drives ALL
% system states/outputs to a given references. Returns 0 if stability of the
% closed-loop system cannot be determined (this means that the closed-loop
% system may, or may not be stable).
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPT controller (MPTCTRL object)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% result - 1 if the controller is stabilizable, 0 if no stability guarantees
%          can be given
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

probStruct = ctrl.probStruct;
sysStruct = ctrl.sysStruct;
result = 0;

if ~isexplicit(ctrl) & iscell(sysStruct.A)
    % no stability guarantees for on-line controllers for PWA systems
    result = 0;
    return
end

if probStruct.tracking>0
    % no stability guarantees for tracking problems
    return
end
if probStruct.norm==2 & probStruct.Tconstraint==1 & probStruct.subopt_lev==0
    % stability guaranteed for 2-norm problems with LQR target set
    result = 1;
    
elseif probStruct.subopt_lev == 1
    % stability guaranteed for minimum-time solutions
    result = 1;
    
elseif probStruct.subopt_lev==0 & isinf(probStruct.N)
    % stability guaranteed for infinite-time problems
    result = 1;
    
elseif isfield(ctrl.details, 'lyapunov')
    % stability is also guaranteed if a lyapunov function has been found
    result = ctrl.details.lyapunov.feasible;
    
end
if isfield(ctrl.probStruct, 'Nc'),
    if ctrl.probStruct.Nc==ctrl.probStruct.N & result==1,
        return
    else
        result = 0;
    end
end
if isfield(ctrl.probStruct, 'inputblocking') | isfield(ctrl.probStruct, 'deltablocking')
    % no stability guarantees for move blocking
    result = 0;
end
