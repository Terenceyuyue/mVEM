function analyze(ctrl)
%ANALYZE Analyzes a given explicit controller and suggests further actions
%
%   analyze(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Analyzes a given explicit controller and suggests further actions
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl - explicit controller (an MPTCTRL object)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
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

if ~isa(ctrl, 'mptctrl')
    error('MPT_ANALYZE: First input must be a valid MPTCTRL object!');
end

if ~isexplicit(ctrl)
    fprintf('\nThis controller is an on-line MPC controller, no further simplifications can be made.\n\n');
    return
end

invariant = isinvariant(ctrl);
stable = isstabilizable(ctrl);
merged = ctrl.simplified;

inpname = inputname(1);

if ~merged & ctrl.probStruct.subopt_lev==0
    fprintf('\nThe controller can be simplified:\n');
    fprintf('   %s = mpt_simplify(%s)  will reduce number of regions.\n', ...
        inpname, inpname);
end

if ~invariant
    fprintf('\nThe closed-loop system may not be invariant:\n');
    fprintf('   %s = mpt_invariantSet(%s)  will identify the invariant subset.\n', ...
        inpname, inpname);
end
if ~stable
    fprintf('\nThe controller may not be stabilizable:\n');
    fprintf('   %s = mpt_lyapunov(%s,''pwq'')  will verify stability using PWQ Lyapunov function.\n',...
        inpname, inpname);
end
    
if (merged | ctrl.probStruct.subopt_lev~=0) & invariant & stable
    fprintf('\nYour controller is as good as it can be. Simple, invariant and stabilizable.\n\n');
else
    fprintf('\n');
end

