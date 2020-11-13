function newCS = mpt_mergeCS(CScell)
%MPT_GLUECS Merges a cell array of ctrlStruct structures
%
% newCS = mpt_glueCS(CScell)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Merges a cell array of ctrlStruct structures
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% CScell   - cell array of ctrlStruct structures
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% newCS    - merged ctrlStruct structure
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(1,1,nargin));

if ~iscell(CScell)
    newCS = CScell;
elseif length(CScell)==1,
    newCS = CScell;
else
    newCS = CScell{1};
    for ctr=2:length(CScell),
        if isempty(CScell{ctr}),
            continue
        end
        oneCS = CScell{ctr};
        newCS.Pn = [newCS.Pn oneCS.Pn];
        newCS.Pfinal = [newCS.Pfinal oneCS.Pfinal];
        newCS.Fi = {newCS.Fi{:}, oneCS.Fi{:}};
        newCS.Gi = {newCS.Gi{:}, oneCS.Gi{:}};
        newCS.Ai = {newCS.Ai{:}, oneCS.Ai{:}};
        newCS.Bi = {newCS.Bi{:}, oneCS.Bi{:}};
        newCS.Ci = {newCS.Ci{:}, oneCS.Ci{:}};
        newCS.dynamics = [newCS.dynamics oneCS.dynamics];
    end
end
