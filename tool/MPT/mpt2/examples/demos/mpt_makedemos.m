function mpt_makedemos
%MPT_MAKEDEMOS creates MAT files with data needed to run demos
%

% Copyright is with the following author(s):
%
%(C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

% ltiDIfh
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator
    probStruct.N=5; probStruct.norm=2;
    disp('making Double_Integrator; probStruct.N=5; probStruct.norm=2;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIfh ctrlStruct
end

% ltiDIih
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator
    probStruct.N=Inf; probStruct.norm=2;
    disp('making Double_Integrator; probStruct.N=Inf; probStruct.norm=2;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIih ctrlStruct
end

% ltiDIft
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator
    probStruct.N=5; probStruct.norm=2; probStruct.tracking=1;
    disp('making Double_Integrator; probStruct.N=Inf; probStruct.norm=2; probStruct.tracking=1;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIft ctrlStruct
end

% ltiDIlc
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator
    probStruct.N=1; probStruct.norm=2; probStruct.subopt_lev=2;
    disp('making Double_Integrator; probStruct.N=1; probStruct.norm=2; probStruct.subopt_lev=2;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIlc ctrlStruct
end

% ltiDIadd
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator_addU
    disp('making Double_Integrator_addU');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIadd ctrlStruct
end

% ltiDIrt
try
    clear sysStruct probStruct ctrlStruct
    Double_Integrator
    probStruct.yref = [2;0]; probStruct.Qy = 100*eye(2);
    disp('making Double_Integrator; probStruct.yref = [2;0]; probStruct.Qy = 100*eye(2);');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save ltiDIrt ctrlStruct
end

% pwa_sincos_n5
try
    clear sysStruct probStruct ctrlStruct
    pwa_sincos
    probStruct.norm = Inf; probStruct.subopt_lev = 0; probStruct.N = 5;
    disp('making pwa_sincos; probStruct.norm = Inf; probStruct.subopt_lev = 0; probStruct.N = 5;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save pwa_sincos_n5 ctrlStruct
end

% pwa_sincos_it
try
    clear sysStruct probStruct ctrlStruct
    pwa_sincos
    probStruct.norm = Inf; probStruct.subopt_lev = 0; probStruct.N = Inf;
    disp('making pwa_sincos; probStruct.norm = Inf; probStruct.subopt_lev = 0; probStruct.N = Inf;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save pwa_sincos_it ctrlStruct
end
    
% pwa_sincos_mt
try
    clear sysStruct probStruct ctrlStruct
    pwa_sincos
    disp('making pwa_sincos');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save pwa_sincos_mt ctrlStruct
end

% pwa_sincos_lc
try
    clear sysStruct probStruct ctrlStruct
    pwa_sincos
    probStruct.norm = 1; probStruct.N = 1; probStruct.subopt_lev=2;
    disp('making pwa_sincos; probStruct.norm = 1; probStruct.N = 1; probStruct.subopt_lev=2;');
    ctrlStruct = mpt_control(sysStruct, probStruct);
    save pwa_sincos_lc ctrlStruct
end

% pa4d
try
    clear sysStruct probStruct ctrlStruct
    FourthOrder
    probStruct.y0bounds = 0;
    ec=mpt_control(sysStruct, probStruct);
    pa4d = ec.Pn;
    save pa4d pa4d
end