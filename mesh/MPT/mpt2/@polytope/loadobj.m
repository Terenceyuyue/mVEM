function out = loadobj(obj)
%LOADOBJ load filter for POLYTOPE objects.
%
% [out] = loadobj(obj)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Invoked on any load command which loads polytope objects. Ensures backward
% compatibility by setting up default values for certain parameters (e.g.
% P.bbox)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%

% ---------------------------------------------------------------------------
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
% --------------------------------------------------------------------------

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

% convert the object to a structure if necessary
if ~isstruct(obj),
    obj = struct(obj);
end
if ~isfield(obj, 'bbox'),
    % versions prior to 1.3 did not store bounding box info in the polytope object
    obj.bbox = [];     
end
if isfield(obj, 'keptrows'),
    % versions after 1.3.1 do not store kept rows in the polytope object
    obj = rmfield(obj, 'keptrows');
end

% now convert the structure to an object if necessary
if ~isa(obj, 'polytope'),
    out = class(obj, 'polytope');
else
    out = obj;
end
