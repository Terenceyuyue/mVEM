function newvar = subsref(var, X)
%SUBSREF Indexed referencing for MPTVAR objects

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
    error('??? Attempt to reference field of non-structure array.');
else
    if ~strcmp(X.type,'()'),
        error(['Indexing with ''' X.type ''' not supported!']);
    end
end

indices = X.subs{1};
if isempty(indices),
    newvar = [];
    return
end

if (indices(1) > var.varlength) | (length(indices) > var.varlength),
    error('??? Index exceeds array dimension');
end

if length(indices)>1
    error('Indexing with more than one index not allowed!');
end

newvarname = sprintf('%s(%d)', var.varname, indices);
newvar = mptvar(var.vartype, 1, newvarname, indices, var);
newvar.multiplier = var.multiplier(indices);
