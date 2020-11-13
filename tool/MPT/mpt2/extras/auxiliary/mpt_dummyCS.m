function ctrlStruct=mpt_dummyCS(Pn,Pfinal)
%MPT_DUMMYCS Returns a dummy controller structure
%
% ctrlStruct=mpt_dummyCS(Pn,Pfinal)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns a dummy controller structure
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn     - Polytope array
% Pfinal - Polytope / Polytope array
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct - valid controller structure with ctrlStruct.Pn = Pn
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

if nargin<2
    Pfinal = Pn;
end

if ~isa(Pn, 'polytope')
    error('mpt_dummyCS: ''Pn'' must be a polytope object!');
end

if ~isa(Pfinal, 'polytope')
    error('mpt_dummyCS: ''Pfinal'' must be a polytope object!');
end

np = length(Pn);
nx = dimension(Pn);
dummyStruct.subopt_lev = 0;

ctrlStruct.Pn = Pn;
ctrlStruct.Pfinal = Pfinal;
ctrlStruct.sysStruct = dummyStruct;
ctrlStruct.probStruct = dummyStruct;

S = cell(1,np);
[S{:}] = deal(zeros(nx));
F = cell(1,np);
[F{:}] = deal(zeros(1,nx));
C = cell(1,np);
[C{:}] = deal(0);

ctrlStruct.Fi = F;
ctrlStruct.Gi = C;
ctrlStruct.Ai = S;
ctrlStruct.Bi = F;
ctrlStruct.Ci = C;
ctrlStruct.dynamics = repmat(0,np,1);
ctrlStruct.details = dummyStruct;
ctrlStruct.overlaps = 1;
