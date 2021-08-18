function P = set(varargin)
%SET Used to modify internal properties of a given polytope object
%
% P = set(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% P = SET(P,'PropertyName', PropertyValue....) sets properties (fields) of a polytope P
%
% USAGE:
%   VALUE = SET(P,'PropertyName',PropertyValue)
%   i.e. value = get(P,'minrep',0)
%
% List of properties that can be accessed:
%   H, K, minrep, normal, xCheb, RCheb, vertices
%
% WARNING:
%   Be sure that you set all fields to correct values! Many internal routines rely on
%   informations already stored in a structure, hence modifying them without a good 
%   knowledge about the internal structure may lead to serious consequences!
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P     - Polytopes
% name  - String giving name of the property
% value - Value of that property
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P     - Updated polytope P
%
% see also GET
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

ni = nargin;

if ~isa(varargin{1}, 'polytope')
  error('SET: First argument MUST be a polytope object');
end
P=struct(varargin{1});

if ~rem(ni, 2)
   error(['SET: Input arguments following the object name must be pairs', ...
          ' of the form PropertyName, PropertyValue']);
end

for ii=2:2:ni
    if ~isfield(P,varargin{ii})
        error('SET: Non-existing property (field)');
    end
    P=setfield(P,varargin{ii},varargin{ii+1});
end
P=class(P,'polytope');