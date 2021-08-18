function rel = mpt_matlabrelease
%MPT_MATLABRELEASE Detects the release number of the running MATLAB instance
%
% release = mpt_matlabrelease
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns the release number of the running MATLAB instance. MATLAB
% R2006a and beyond is identified as release 15.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
%(C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

try
    rel = str2num(version('-release'));
catch
    % matlab versions prior to R12 do not allow to call version with the
    % '-release' argument, therefore if an error occures, it indicates an
    % old version of  matlab
    rel = 11;
    return
end
if isempty(rel)
    % Matlab R2006a, R2006b, ... return release as '2006a', '2006b', ...
    % hence str2num() will return an empty string
    rel = 15;
end
