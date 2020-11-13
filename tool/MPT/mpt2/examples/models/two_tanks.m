%TWO_TANKS 2nd order PWA example
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Model of a two-tanks system with discrete and continuous inputs
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% none
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sysStruct, probStruct - system and problem definition structures stores
%                         in the workspace
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

clear sysStruct probStruct

dT = 10;  % sampling time
q_in = 0.1e-3;

%% physical parameters of the system
hv = 0.3; % height of valve 1
s = 10e-6; % scross-section of valves
hmax = 0.62; % maximum height of liquid in both tanks
g = 9.81;
h1s = 0.25; % steady-state level in first tank
h2s = 0.2; % steady-state level in second tank
k1 = s*sqrt(2*g/(hmax-hv));
k2 = s*sqrt(2*g/hmax);
F = 0.0143; % area of each tank
e=0;

F = F/dT;

%% common constraint on states:
% x2<=x1
Xcom = [-1 1];
Ucom = [0 0];
Ccom = [0];

% Linearization around h1s, h2s:

B = [q_in/F 0; 0 0];
C = [0 1];
D = [0 0];
f = [0; 0];

% V2 open, h1>=hv
sysStruct.A{1} = [1-k1/F -k2/F; (k1+k2)/F 1-k2/F];
sysStruct.B{1} = B;
sysStruct.C{1} = C;
sysStruct.D{1} = D;
sysStruct.f{1} = [hv/F*k1; -hv/F*k1];
sysStruct.guardX{1} = [-1 0; 0 0];
sysStruct.guardU{1} = [0 0; 0 -1];
sysStruct.guardC{1} = [-hv; -0.5];


% V2 open, h1<=hv
sysStruct.A{2} = [1-k2/F 0; k2/F 1-k2/F];
sysStruct.B{2} = B;
sysStruct.C{2} = C;
sysStruct.D{2} = D;
sysStruct.f{2} = f;
sysStruct.guardX{2} = [1 0; 0 0];
sysStruct.guardU{2} = [0 0; 0 -1];
sysStruct.guardC{2} = [hv; -0.5];

% V2 closed, h1>=hv
sysStruct.A{3} = [1-k1/F 0; k1/F 1-k2/F];
sysStruct.B{3} = B;
sysStruct.C{3} = C;
sysStruct.D{3} = D;
sysStruct.f{3} = [k1/F*hv; -k1/F*hv];
sysStruct.guardX{3} = [-1 0; 0 0];
sysStruct.guardU{3} = [0 0; 0 1];
sysStruct.guardC{3} = [-hv; 0.5];

% V2 closed, h1<=hv
sysStruct.A{4} = [1 0; 0 1-k2/F];
sysStruct.B{4} = B;
sysStruct.C{4} = C;
sysStruct.D{4} = D;
sysStruct.f{4} = f;
sysStruct.guardX{4} = [1 0; 0 0];
sysStruct.guardU{4} = [0 0; 0 1];
sysStruct.guardC{4} = [hv; 0.5];


for ii=1:4,
    sysStruct.guardX{ii} = [Xcom; sysStruct.guardX{ii}];
    sysStruct.guardU{ii} = [Ucom; sysStruct.guardU{ii}];
    sysStruct.guardC{ii} = [Ccom; sysStruct.guardC{ii}];
end

sysStruct.xmax = [hmax; hmax];
sysStruct.xmin = [0; 0];
sysStruct.ymax = [hmax];
sysStruct.ymin = [0];
sysStruct.umax = [1; 1+e];
sysStruct.umin = [0; 0];
sysStruct.Pbnd = unitbox(2,1);

%% 1st input is continuous, 2nd is from set {0,1}
sysStruct.Uset = {[-Inf Inf],[0 1]};

%% 1st input is from set {0, 0.5}, 2nd from set {0,1}
%sysStruct.Uset = {[0 0.5],[0 1]};

% string labels for states, inputs and outputs
sysStruct.StateName = { 'Level 1', 'Level 2' }; 
sysStruct.InputName = { 'Liquid input', 'Valve 2' };
sysStruct.OutputName = { 'Level 2' };


probStruct.norm = 1;
probStruct.N = 3;
probStruct.subopt_lev = 0;
probStruct.Q = eye(2);
probStruct.R = diag([1e-5,1e-5]);
probStruct.Qy = 200;

probStruct.y0bounds = 0;

% reference for output tracking
probStruct.yref = 0.2;


clear a11 a12 a21 a22 A B C D Ad Bd Cd Dd dT hv s hmax g h1s h2s k1 k2 F ii e q_in x1ref x2ref Ucom Xcom Ccom
return