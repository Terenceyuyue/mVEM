function status=mpt_isValidCS(ctrlStruct,Options)
%MPT_ISVALIDCS Checks if input argument is a valid controller structure
%
% status=mpt_isValidCS(ctrlStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Checks if input argument is a valid controller structure
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrlStruct  - Controller structure as generated with mpt_control
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status      - 1 if input is a valid controller structure, 0 otherwise
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
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

error(nargchk(1,2,nargin));

if nargin<2
    Options = [];
end

if ~isfield(Options,'nowarnings'),
    Options.nowarnings = 0;
end

if isa(ctrlStruct, 'mptctrl') 
  ctrlStruct = struct(ctrlStruct);
end

status = 0;
if ~isfield(ctrlStruct,'Pn') | ...
        ~isfield(ctrlStruct,'sysStruct') | ...
        ~isfield(ctrlStruct,'probStruct') | ...
        ~isfield(ctrlStruct,'Fi') | ...
        ~isfield(ctrlStruct,'Gi') | ...
        ~isfield(ctrlStruct,'Ai') | ...
        ~isfield(ctrlStruct,'Bi') | ...
        ~isfield(ctrlStruct,'Ci') | ...
        ~isfield(ctrlStruct,'dynamics') | ...
        ~isfield(ctrlStruct,'details') | ...
        ~isfield(ctrlStruct,'overlaps') | ...
        ~isfield(ctrlStruct,'Pfinal')
    if Options.nowarnings == 0,
        warning('mpt_isValidCS: one of the fields missing! Check if ctrlStruct contains [sysStruct, probStruct, Fi, Gi, Ai, Bi, Ci, dynamics, details, overlaps, Pfinal]');
    end
    return
end
if ~isa(ctrlStruct.Pn,'polytope') | ...
        ~isa(ctrlStruct.Pfinal,'polytope') | ...
        ~isstruct(ctrlStruct.sysStruct) | ...
        ~isstruct(ctrlStruct.probStruct) | ...
        ~iscell(ctrlStruct.Fi) | ...
        ~iscell(ctrlStruct.Gi) | ...
        ~iscell(ctrlStruct.Ai) | ...
        ~iscell(ctrlStruct.Bi) | ...
        ~iscell(ctrlStruct.Ci)
    if Options.nowarnings == 0,
        warning('mpt_isValidCS: [Pn, Pfinal] must be polytopes, [sysStruct, probStruct] structures, [Fi,Gi,Ai,Bi,Ci] cells!');
    end
    return
end
if length(ctrlStruct.Pn)~=length(ctrlStruct.Fi) | ...
        length(ctrlStruct.Pn)~=length(ctrlStruct.Gi) | ...
        length(ctrlStruct.Pn)~=length(ctrlStruct.Ai) | ...
        length(ctrlStruct.Pn)~=length(ctrlStruct.Bi) | ...
        length(ctrlStruct.Pn)~=length(ctrlStruct.Ci) | ...
        length(ctrlStruct.Pn)~=length(ctrlStruct.dynamics),
    if Options.nowarnings == 0,
        warning('mpt_isValidCS: Pn, Fi, Gi, Ai, Bi, Ci, dynamics must be of the same lenght!');
    end
    status = 0;
    return
end
status = 1;
