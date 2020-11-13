function sysStruct = mpt_lti2pwa(sysStruct)
%MPT_LTI2PWA Converts an LTI system to a PWA system
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Converts an LTI-type sysStruct structure to a PWA system with one dynamics
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct  - system structure describing an LTI system
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sysStruct  - PWA system
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

infbox = mptOptions.infbox;

if iscell(sysStruct.A),
    % system is already in PWA format
    return
else
    % convert an LTI system to PWA description
    A{1} = sysStruct.A;
    B{1} = sysStruct.B;
    C{1} = sysStruct.C;
    D{1} = sysStruct.D;
    if isfield(sysStruct, 'Cy'),
        Cy{1} = sysStruct.Cy;
        sysStruct.Cy = Cy;
        Dy{1} = sysStruct.Dy;
        sysStruct.Dy = Dy;
    end
    if isfield(sysStruct,'f')
        f{1} = sysStruct.f;
        sysStruct.f = f;
    end
    if isfield(sysStruct,'g')
        g{1} = sysStruct.g;
        sysStruct.g = g;
    end
    nx = size(sysStruct.A,2);
    sysStruct.A = A;
    sysStruct.B = B;
    sysStruct.C = C;
    sysStruct.D = D;
    sysStruct.guardX{1} = [eye(nx); -eye(nx)];
    sysStruct.guardC{1} = [ones(2*nx,1)*infbox];
    verOpt.verbose=0;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end
