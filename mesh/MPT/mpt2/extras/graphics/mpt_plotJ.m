function mpt_plotJ(ctrl,gridpoints,Options)
%MPT_PLOTJ Plots value function associated to a given controller
%
% mpt_plotJ(ctrl)
% mpt_plotJ(ctrl,gridpoints)
% mpt_plotJ(ctrl,gridpoints,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the value function defined over regions of a polyhedral partition
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl          - Explicit controller (MPTCTRL object)
% gridpoints    - Number of grid points in one axis (for PWQ case)
% Options       - Auxiliary options
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% see also MPT_PLOTU, MPT_PLOTPWA, MPT_PLOTPWQ
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

error(nargchk(1,3,nargin));

if nargin<3,
    Options = [];
end

if nargin<2,
    gridpoints = 30;
end
if nargin==2,
    if isstruct(gridpoints),
        Options=gridpoints;
        gridpoints=30;
    end
end
if ~isfield(Options,'showPn'),
    Options.showPn=1;
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl),
        error('This function supports only explicit controllers!');
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end

if ~mpt_isValidCS(ctrlStruct)
    error('mpt_plotJ: First argument has to be a valid controller structure! See mpt_control for details.');
end

if ~isa(gridpoints,'double'),
    error('mpt_plotJ: Second argument must be a number of grid points!');
end



if all(all([ctrlStruct.Ai{:}]==0)),
    % no quadratic term in the cost function, use PWA plot
    mpt_plotPWA(ctrlStruct.Pn, ctrlStruct.Bi, ctrlStruct.Ci, Options);
else
    mpt_plotPWQ(ctrlStruct.Pn, ctrlStruct.Ai, ctrlStruct.Bi, ctrlStruct.Ci, gridpoints, Options);
end
if dimension(ctrlStruct.Pn)==1,
    ylabel('J^{*}(x)');
else
    zlabel('J^{*}(x)');
end
title(sprintf('Value function over %d regions',length(ctrlStruct.Pn)),'FontSize',18);
% if isfield(ctrlStruct.sysStruct,'StateName'),
%     xlabels = ctrlStruct.sysStruct.StateName;
%     xlabel(xlabels{1}','Fontsize',16); % LaTeX math symbols are directly supported!
%     ylabel(xlabels{2},'Fontsize',16);
% end
