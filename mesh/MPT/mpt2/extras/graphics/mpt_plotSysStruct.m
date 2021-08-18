function mpt_plotSysStruct(sysStruct,Options)
%MPT_PLOTSYSSTRUCT Plots the system partitions of a PWA system
%
% mpt_plotSysStruct(sysStruct)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct    - System structure in the sysStruct format (help mpt_sysStruct)                  
%                        
%
% see also PLOT, MPT_PLOTPARTITION, MPT_PLOTJ, MPT_PLOTPWA, MPT_PLOTPWQ

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

if nargin<2,
    Options = [];
end
if ~isfield(Options,'shade'),
    Options.shade = 0.8;
end

if ~isfield(sysStruct,'verified'),
    sysStruct = mpt_verifySysStruct(sysStruct);
end

if(~iscell(sysStruct.A))
    %error('mpt_plotSysStruct function only works for PWA systems, i.e. dynamical matrices must be cells')
    sysStruct = mpt_lti2pwa(sysStruct);
end

[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

PP = polytope(intInfo.Pdyn);
if dimension(PP)==2,
    for ii=1:length(PP),
        if isfulldim(PP(ii))
            V = extreme(PP(ii),struct('extreme_solver',0));
            patch(V(:,1),V(:,2),ii*ones(size(V,1)),ii,'FaceAlpha',Options.shade);
        end
    end
    grid on
    view(3);
else
    plot(PP);
end
return
