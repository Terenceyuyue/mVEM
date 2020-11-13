function fval = subsref(ctrl, X)
%SUBSREF Indexed referencing for MPTCTRL object
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Indexing of controller objects can perform one opf the following things:
%
% 1. If "u=ctrl(x0)" is used, the controller is evaluated for a given state x0,
%    i.e. it will return "u" such that u=mpt_getInput(ctrl, x0). 
%    If "ctrl(x0, Options)" is used, it corresponds to 
%    u=mpt_getInput(ctrl, x0, Options)
%
% 2. If "ctrl.fieldname" is used, it will return the value of field "fieldname",
%    e.g. "ctrl.Pn" returns the polyhedral partition of the controller. Type
%    'struct(ctrl)' to see all accessible fields.
%
% If a vector of NaNs is returned, the problem is infeasible for a given value
% of the initial state x0.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% see also SUBSASGN, MPT_GETINPUT
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

if numel(X)>1,
    XX = X(1);
else
    XX = X;
end
dotindex = strcmp(XX.type, '.');
evalindex = strcmp(XX.type, '()');
if ~(dotindex | evalindex),
    error(['Indexing with ''' X.type ''' not supported!']);
end

if evalindex,
    % evaluate the controller for given state
    x0 = X.subs{1};
    if length(X.subs) > 1,
        Options = X.subs{2};
    else
        Options = [];
    end
    try
        [fval, feasible] = mpt_getInput(ctrl, x0, Options);
        if ~feasible,
            fval = NaN*fval;
        end
    catch
        rethrow(lasterror);
    end
    
else
    % extract value of given field of the controller structure
    ctrl_str = struct(ctrl);
    if length(X)==1,
        fname = X.subs;
    else
        fname = X(1).subs;
    end
    if isfield(ctrl_str, fname)
        fval = getfield(ctrl_str, fname);
    else
        error(['??? Reference to non-existent field ''' fname '''.']);
    end
    if length(X)>1,
        try
            fval = subsref(fval, X(2:end));
        catch
            rethrow(lasterror);
        end
    end
end
