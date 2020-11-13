function obj = mptvar(vartype, varlength, varname, varindex, fromvar)
%MPTVAR Constructor of MPTVAR objects

% Copyright is with the following author(s):
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

error(nargchk(1,5,nargin));

if ~(strcmpi(vartype, 'state') | strcmpi(vartype, 'input')),
    error(sprintf('MPTVAR: unknown variable type "%s"!', vartype));
end
obj.vartype = vartype;

if nargin==1,
    varlength = 1;
    if strcmpi(vartype, 'state')
        varname = 'X';
    elseif strcmpi(vartype, 'input')
        varname = 'U';
    end
    varindex = 0;
    fromvar = [];
elseif nargin==2
    if strcmpi(vartype, 'state')
        varname = 'X';
    elseif strcmpi(vartype, 'input')
        varname = 'U';
    end
    varindex = 0;
    fromvar = [];
elseif nargin==3
    obj.varindex = 0;
    obj.fromvar = [];
elseif nargin==4
    fromvar = [];
end

obj.varlength = varlength;
obj.varname = varname;
obj.varindex = varindex;
obj.fromvar = fromvar;

obj.multiplier = ones(varlength, 1);
obj.Array = {};


obj = class(obj, 'mptvar');
