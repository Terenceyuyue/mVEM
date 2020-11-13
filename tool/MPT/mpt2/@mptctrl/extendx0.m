function [x0, dumode] = extendx0(ctrl, x0, uprev, reference)
%EXTENDX0 Extends the initial state if needed
%
% [x0, dumode] = extendx0(ctrl, x0, uprev, reference)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% If tracking and/or deltaU constraints have been used to design a particular
% controller, the state vector needs to be augmented in order to be able to deal
% with these goals in closed-loop. This implies that the value of the initial
% state must also be augmented accordingly when a given controller is evaluated
% by
%
%   u = ctrl(x0)
%
% In order to automate this task, one can call this helper as follows:
%
%   [x0, dumode] = extendx0(ctrl, x0, uprev, reference)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl       - MPT controller object
% x0         - Current measured state
% uprev      - Previous value of the control action u(k-1)
% reference  - Reference trajectory (if tracking was used)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% x0         - The augmented state vector
% dumode     - If true, the controller returns deltaU=u(k)-u(k-1), i.e. the
%              "true" control action must be additionaly obtained by
%                utrue = u + uprev
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(2,4,nargin));
if nargin < 3,
    uprev = [];
    reference = [];
elseif nargin < 4,
    reference = [];
end

x0 = x0(:);
uprev = uprev(:);
reference = reference(:);

x0format = ctrl.details.x0format;
dumode = x0format.uprev > 0;

%=====================================================================
% exit quickly if the state already has correct length
if length(x0) == x0format.required,
    return
end

%=====================================================================
% augment the state vector if need
if x0format.uprev > 0,
    if length(uprev) ~= x0format.uprev,
        error('Wrong dimension of u(k-1).');
    end
    x0 = [x0; uprev];
end
if x0format.reference > 0,
    if length(reference) ~= x0format.reference,
        error('Wrong dimension of the reference.');
    end
    x0 = [x0; reference(:)];
end
if length(x0) ~= x0format.required,
    error('Wrong dimension of x0.');
end
