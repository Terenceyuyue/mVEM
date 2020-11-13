function sys = mpt_nonlinfcn(flag, x, u, mode)
% MPT_NONLINFCN Template function for (piecewise) nonlinear systems
%
% This is a template function for defining nonlinear or piecewise nonlinear
% systems. The function is divided into four main parts which are executed
% depending on the value of the "flag" input.
%
% In the "info" part, the user provides parameters which define system
% dimensions and constraints. If certain constraints are not known, or the
% user does not want to impose them, respective fields can be defined as empty
% matrices.
%
% In the "state" part, the state update equations are specified as (possibly
% nonlinear) functions of the states (denoted by "x") and inputs (denoted by
% "u"). If the system is piecewise nonlinear, one state update equation per each
% segment must be specified. Index of a respective segment is provided in the
% "mode" variable.
%
% The purpose of the "output" section is to define the output equation as a
% function of system states and system inputs. Again, one output equation per
% segment must be given if the system is piecewise nonlinear.
%
% Finally, the "guards" section specifies portions of the state-input space
% where a given segment (one piece of the piecewise nonlinear dynamics) is
% active. The equations can again be specified using the state and input
% variables.
%
% To create a system structure based on this template, simply call following
% command:
%
%   sysStruct = mpt_sys(@function_name)
%

ispiecewise = (nargin > 3);

switch flag,
    case 'info',
        % number of system states
        sys.nx = 0;
        
        % number of system inputs
        sys.nu = 0;
        
        % number of system outputs
        sys.ny = 0;
        
        % is this system piecewise non-linear? if yes, this parameter specifies
        % number of piecewise segments. set this value to zero if the system is
        % not piecewise non-linear. 
        sys.piecewise = 0;  
        
        % state constraints
        sys.xmin = [];
        sys.xmax = [];
        
        % input constraints
        sys.umin = [];
        sys.umax = [];
        
        % output constraints
        sys.ymin = [];
        sys.ymax = [];
        
        % deltaU constraints
        sys.dumin = [];        
        sys.dumax = [];

        % are there any boolean or integer inputs? see "help mpt_sysStruct" for
        % details
        sys.Uset = [];
        
        % sampling time
        sys.Ts = [];
        
    case 'state',
        % define the state update equation. if the system is piecewise
        % non-linear, use the "mode" input to switch between different dynamics
        
        
        if ispiecewise,
            % piecewise non-linear system
            switch mode,
                case 1,
                    % state update equation for mode 1
                    sys = x;
                    
                case 2,
                    % state update equation for mode 2
                    sys = x;
            end
            
        else
            % pure non-linear system; state update equation
            sys = x;
            
        end
        
    case 'output',
        % define the output equation. if the system is piecewise non-linear, use
        % the "mode" input to switch between different dynamics 
        
        if ispiecewise,
            % piecewise non-linear system
            switch mode,
                case 1,
                    % output equation for mode 1
                    sys = x;
                    
                case 2,
                    % output equation for mode 1
                    sys = x;
            end
            
        else
            % pure non-linear system; output equation
            sys = x;
            
        end
        
    case 'guards',
        % define the guardlines (only for piecewise non-linear systems)
        
        switch mode,
            case 1,
                sys = x(1) >= 1;
                
            case 2,
                sys = x(1) <= 1;
                
        end
end
