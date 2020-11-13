% "Ball and Plate" demo

% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

clear sysx sysy probStruct

% continuous-time model for the X-axis
Ax = [  0    1.0000         0         0
    0         0 -700.0000         0
    0         0         0    1.0000
    0         0         0  -33.1800  ];
Bx = [0; 0; 0; 3.7921 ];
Cx = [1 0 0 0]; 
Dx = 0;     

% continuous-time model for the Y-axis
Ay = [  0    1.0000         0         0
    0         0  700.0000         0
    0         0         0    1.0000
    0         0         0  -34.6900  ];
By = [0; 0; 0; 3.1119 ];
Cy = [1 0 0 0];
Dy = 0;

% sampling time
Ts = 0.03;

% define state-space objects
css_x = ss(Ax, Bx, Cx, Dx);
css_y = ss(Ay, By, Cy, Dy);

% discretize the systems
dss_x = c2d(css_x, Ts);
dss_y = c2d(css_y, Ts);

% convert state-space objects to system structures
sysx = mpt_sys(dss_x, Ts);
% input constraints
sysx.umax = 10;
sysx.umin = -10;
% output constraints
sysx.ymax = 30;
sysx.ymin = -30;
% state constraints
sysx.xmax = [30; 15; 15*pi/180; 1];
sysx.xmin = -sysx.xmax;

sysy = mpt_sys(dss_y, Ts);
% input constraints
sysy.umax = 10;
sysy.umin = -10;
% output constraints
sysy.ymax = 30;
sysy.ymin = -30;
% state constraints
sysy.xmax = [30; 15; 15*pi/180; 1];
sysy.xmin = -sysy.xmax;

% optimal solution
probStruct.subopt_lev = 0;

% penalty on states (will not be used since we have penalty on outputs)
probStruct.Q = eye(4);

% penalty on outputs
probStruct.Qy = 100;

% penalty on inputs
probStruct.R = 0.1;

% use tracking without deltaU formulation - saves one additional state
probStruct.tracking = 2;

% quadratic performance index
probStruct.norm = 2;

% prediction horizon
probStruct.N = 10;

% control horizon
probStruct.Nc = 2;

% penalty on final output
probStruct.P_N = 1;


% load explicit controllers
load ballandplateexp

% if you change any parameters of the dynamicals systems or the problem setup,
% you will need to recompute the explicit controllers by running
%   ecx = mpt_control(sysx, probStruct);
%   ecy = mpt_control(sysy, probStruct);


% you can also compute on-line controllers:
%onx = mpt_control(sysx, probStruct, 'online');
%ony = mpt_control(sysy, probStruct, 'online');


% load trajectory
load circle_10

% in the sequel we will use MEX implementation of mpt_getInput, therefore we
% need to extract information about the two explicit controllers:
ecx_data = mpt_mexData(ecx);
ecy_data = mpt_mexData(ecy);


fprintf('\nSimulating...\n');

% initial conditions
x0 = [0;0;0;0];
y0 = [0;0;0;0];

% initial values of input signals
vx = 0;
vy = 0;

% initial values of u(k-1)
vxprev = 0;
vyprev = 0;

X = [];
Y = [];
U = [];

for ii = 1:size(yref, 1),
    if mod(ii,10)==0,
        fprintf('.');
    end
    if mod(ii,300)==0,
        fprintf('\n');
    end
    
    XREF = xref(ii, 2);
    YREF = yref(ii, 2);
    
    if probStruct.tracking==1,
        % tracking=1 implies that the state vector is augmented to
        % xaug = [x; u(k-1); reference]
        x0_c = [x0; vxprev; XREF(:)];
        y0_c = [y0; vyprev; YREF(:)];
  
    elseif probStruct.tracking==2,
        % tracking=2 implies the state vector is augmented to
        % xaug = [x; reference]
        x0_c = [x0; XREF(:)];
        y0_c = [y0; YREF(:)];
    end
    
    % use mex_getInput - C-code implementation of the region identification
    % algorithm. this is much faster than calling mpt_getInput, but has limited
    % functionality (see mex_getInput.c for more details).
    [vx, regx] = mex_getInput(ecx_data, x0_c);
    [vy, regy] = mex_getInput(ecy_data, y0_c);

    if  regx==0 | regy==0,
        % if the region number is 0 (zero), no feasible control law was found
        % for a given state
        error('infeasible problem!');
    end

    if probStruct.tracking==1,
        % IMPORTANT: controllers generated with probStruct.tracking=1 return
        % deltaU instead of U, hence we need to add previous input u(k-1) to
        % obtain the "true" system input
        vx = vx + vxprev;
        vy = vy + vyprev;
        
        vxprev = vx;
        vyprev = vy;
    end
    
    % compute state update for the X-axis model
    xn = dss_x.A*x0 + dss_x.B*vx;

    % compute state update for the Y-axis model
    yn = dss_y.A*y0 + dss_y.B*vy;
    
    % check if any constraints are violated
    xmaxviol = find(xn > sysx.xmax);
    xn(xmaxviol) = sysx.xmax(xmaxviol);
    xminviol = find(xn < sysx.xmin);
    xn(xminviol) = sysx.xmin(xminviol);

    % check if any constraints are violated
    xmaxviol = find(yn > sysy.xmax);
    yn(xmaxviol) = sysy.xmax(xmaxviol);
    xminviol = find(yn < sysy.xmin);
    yn(xminviol) = sysy.xmin(xminviol);

    % store data to matrices
    x0 = xn;
    X = [X; x0'];

    y0 = yn;
    Y = [Y; y0'];
   
    U = [U; vx vy];
end
fprintf('\n');


% plot the closed-loop trajectories
plot((1:size(Y,1))*Ts,Y(:,1),(1:size(Y,1))*Ts,yref(:,2)); 
xlabel('Time [s]'); ylabel('y [cm]'); legend('Y position', 'Reference');
figure; plot((1:size(X,1))*Ts,X(:,1),(1:size(X,1))*Ts,xref(:,2))
xlabel('Time [s]'); ylabel('x [cm]'); legend('X position', 'Reference');
figure; plot(X(:,1),Y(:,1),XREF,YREF)
xlabel('x [cm]'); ylabel('y [cm]'); legend('Ball Position', 'Reference');


% run simulations in Simulink
fprintf('\n\nNow running simulation is Simulink...\n');
open_system('ballandplatesim');
sim('ballandplatesim', stoptime);