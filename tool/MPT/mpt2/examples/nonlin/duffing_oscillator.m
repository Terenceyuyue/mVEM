function sys = duffing_oscillator(flag, x, u, mode)
% Duffing Oscillator example (nonlinear system)
%
% Example taken from:
% @InProceedings { FotEtal:2006:IFA_2224,
%     author={I.A. Fotiou and Ph. Rostalski and B. Sturmfels and M. Morari},
%     title={{An algebraic geometry approach to nonlinear parametric
% 	  optimization in control}},
%     booktitle={American Control Conference},
%     pages={},
%     year={2006},
%     address={},
%     month=jun,
%     url={http://control.ee.ethz.ch/index.cgi?page=publications;action=details;id=2224}
% }

ispiecewise = (nargin > 3);
Ts = 0.2;

switch flag,
    case 'info',
        % number of system states
        sys.nx = 2;
        
        % number of system inputs
        sys.nu = 1;
        
        % number of system outputs
        sys.ny = 1;
        
        % is this system piecewise non-linear? if yes, this parameter specifies
        % number of piecewise segments. set this value to zero if the system is
        % not piecewise non-linear. 
        sys.piecewise = 0;  
        
        % state constraints
        sys.xmin = [-5; -5];
        sys.xmax = [5; 5];
        
        % input constraints
        sys.umin = -10;
        sys.umax = 10;
        
        % output constraints
        sys.ymin = -5;
        sys.ymax = 5;
        
        % deltaU constraints
        sys.dumin = [];        
        sys.dumax = [];

        % are there any boolean or integer inputs? see "help mpt_sysStruct" for
        % details
        sys.Uset = [];
        
        % sampling time
        sys.Ts = Ts;
        
    case 'state',
        % define the state update equation. if the system is piecewise
        % non-linear, use the "mode" input to switch between different dynamics
        
        % pure non-linear system; state update equation
        h = Ts; ksi = 0.3;
        A = [1 h;-h 1-2*ksi*h];
        B = [0;h];
        sys = A*x + B*u + [0; -h*x(1)^3];
        
    case 'output',
        % define the output equation. if the system is piecewise non-linear, use
        % the "mode" input to switch between different dynamics 
        
        % output is the first state
        sys = x(1);
        
    case 'guards',
        % we have no guards in this example, because the system is not piecewise
        % nonlinear
        
end
