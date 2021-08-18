function  [sys,x0,str,ts] = two_tanks_sim(t,x,u,flag,xinit)
%TWO_TANKS_SIM S-function of the two tank system
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Simulation of the two tank system
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


switch flag
    case 0                                                % Initialization
        sys = [2,      % number of continuous states
            0,      % number of discrete states
            2,      % number of outputs
            2,      % number of inputs
            0,      % reserved must be zero
            0,      % direct feedthrough flag
            1];     % number of sample times
        %x0  = [0.5; 0.2];
        x0 = xinit;
        str = [];
        ts  = [0 0];   % sample time: [period, offset]
        
    case 1                                                % Derivatives
        sys = state_update(x,u);
        
    case 2                                                % Discrete state update
        sys = []; % do nothing
        
    case 3
        if x(2)>0.3,
            x(2) = 0.3;
        end
        sys = x;
        
    case 9                                                % Terminate
        sys = []; % do nothing
        
    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end
return

function xnext = state_update(x0,u)

x1ref = 0.25;
x2ref = 0.2;

dT = 10;  % sampling time
q_in = 0.1e-3;

hv = 0.3; % height of valve 1
s = 10e-6; % scross-section of valves
hmax = 0.62; % maximum height of liquid in both tanks
g = 9.81;
h1s = 0.35; % steady-state level in first tank
h2s = 0.1; % steady-state level in second tank
k1 = s*sqrt(2*g/(hmax-hv));
k2 = s*sqrt(2*g/hmax);
F = 0.0143; % area of each tank

h1 = x0(1);
h2 = x0(2);
if h2<0,
    h2 = 0;
end
if h2>0.3,
    h2 = 0.3;
end
if h1<0
    h1 = 0;
end
if u(2)>0.5
    % valve 2 open
    if h1>=hv
        mode = 1;
        dh1 = 1/F * (q_in*u(1) - k1*sqrt(h1-hv) - k2*sqrt(h1-h2));
        dh2 = 1/F * (k1*sqrt(h1-hv) + k2*sqrt(h1-h2) - k2*sqrt(h2));
    else
        mode = 2;
        dh1 = 1/F * (q_in*u(1) - k2*sqrt(h1-h2));
        dh2 = 1/F * (k2*sqrt(h1-h2) - k2*sqrt(h2));
    end
else
    if h1>=hv
        mode = 3;
        dh1 = 1/F * (q_in*u(1) - k1*sqrt(h1-hv));
        dh2 = 1/F * (k1*sqrt(h1-hv) - k2*sqrt(h2));
    else
        mode = 4;
        dh1 = 1/F * q_in *u(1);
        dh2 = -1/F * k2*sqrt(h2);
    end
end
if ~isreal(dh1),
    dh1=real(dh1);
end
if ~isreal(dh2),
    dh2=real(dh2);
end

xnext = [dh1; dh2];