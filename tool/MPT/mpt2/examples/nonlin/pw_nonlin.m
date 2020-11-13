function sys = pw_nonlin(flag, x, u, mode)
% Example of an piecewise nonlinear model:
%
% x(k+1) = A1*x + B1*u    if x'*x <= 1
% x(k+1) = A2*x + B2*u    if x'*x >= 1
%   y(k) = C*x

ispiecewise = (nargin > 3);
Ts = 1;

switch flag,
    case 'info',
        % number of system states
        sys.nx = 2;
        
        % number of system inputs
        sys.nu = 1;
        
        % number of system outputs
        sys.ny = 1;

        % we have 2 different dynamics
        sys.piecewise = 2;
        
        % state constraints
        sys.xmin = [-5; -5];
        sys.xmax = [5; 5];
        
        % input constraints
        sys.umin = -1;
        sys.umax = 1;
        
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
        
        % one state-update equation per mode
        alpha = -pi/3;
        A1 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
        B1 = [0; 1];
        
        alpha = pi/3;
        A2 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
        B2 = [0; 1];
        
        switch mode,
            case 1,
                sys = A1*x + B1*u;
                
            case 2,
                sys = A2*x + B2*u;
                
        end
        
    case 'output',
        % define the output equation. if the system is piecewise non-linear, use
        % the "mode" input to switch between different dynamics 
        
        % one output equatiuon per mode
        switch mode,
            case 1,
                sys = x(1);
                
            case 2,
                sys = x(1);
                
        end
        
    case 'guards',
        % guards specify a region of the state-input space where a given "mode"
        % is active
        
        switch mode,
            case 1,
                sys = x'*x <= 1;
                
            case 2,
                sys = x'*x >= 1;
                
        end
        
end
