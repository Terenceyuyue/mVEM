function [sys,x0,str,ts] = hylink(t,x,u,flag,filename,x_init)
%HYLINK: HYSDEL models in simulink via S-function
%
% use: drag a simulink block of type 'S-function' in your simulink model
%      set the function name to hylink and the filename as parameter.
%      The optional parameter x0 sets the initial state.
%      Use multiplex/demultiplex to access all the inputs/outputs
%   
% (C) 2000--2002 F.D. Torrisi,
% Automatic Control Laboratory, ETH Zentrum, CH-8092 Zurich, Switzerland
% torrisi@aut.ee.ethz.ch
%
% see license.txt for the terms and conditions.


switch flag,
  % Initialization
  case 0,                                                
    [sys,str,ts] = mdlInitializeSizes(filename);    
    if nargin > 5,
        x0 = x_init;
    else
        x0 = 0;
    end
  % Update 
  case 2,                                               
    sys = mdlUpdate(t,x,u,[filename '_sim']);
  % Output
  case 3,                                               
    sys = mdlOutputs(t,x,u,[filename '_sim']);    
  % Terminate 
  case 9,                                               
    sys = [];

  otherwise
    error(['unhandled flag = ',num2str(flag)]);
end

%end hylink

% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
function [sys,str,ts]=mdlInitializeSizes(filename)

hysdel(filename,[filename '_sim']);
eval(char(filename));

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = S.nx;
sizes.NumOutputs     = S.ny;
sizes.NumInputs      = S.nu;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

str = [];
ts  = [1 0]; % Sample period of 1 seconds (1Hz)

% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
function sys = mdlUpdate(t,x,u,simname)
eval(['sys = ' simname '(x,u);']);    

% mdlOutputs
% Return the output vector for the S-function
function sys = mdlOutputs(t,x,u,simname)
eval(['[dummy1, dummy2, dummy4, sys] = ' simname '(x,u);']);    

%end mdlOutputs

