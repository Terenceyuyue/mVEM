function ctrl = set(varargin)
%SET Set a field of MPTCTRL objects
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% ctrl = set(ctrl,'PropertyName', PropertyValue....) sets properties (fields) of
% the controller 'ctrl'
%
% USAGE:
%   ctrl = set(P,'PropertyName',PropertyValue)
%   e.g. ctrl = set(ctrl, 'Pfinal', P)
%
% List of properties that can be accessed:
%   Pfinal, Pn, Fi, Gi, Ai, Bi, Ci, dynamics, details, overlaps, simplified,
%   type, sysStruct, probStruct
%
% WARNING:
%   Using set() can be rather dangerous if you don't know what you are doing!
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl  - MPTCTRL object
% name  - String giving name of the property
% value - Value of that property
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrl  - Updated controller object
%
% see also GET
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

if ~isa(varargin{1}, 'mptctrl')
  error('SET: First argument must be an MPTCTRL object,');
end

ni = nargin;
if ~rem(ni, 2)
   error(['SET: Input arguments following the object name must be pairs', ...
          ' of the form PropertyName, PropertyValue']);
end

ctrl = struct(varargin{1});

for ii=2:2:ni
    if ~isfield(ctrl, varargin{ii})
        error('SET: Non-existing property (field)');
    end
    ctrl = setfield(ctrl, varargin{ii}, varargin{ii+1});
end
ctrl = class(ctrl, 'mptctrl');
