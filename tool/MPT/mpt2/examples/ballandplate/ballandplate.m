% this function only defines parameters of the "ball and plate" system,
% for more complete description, execute "ballandplatedemo"

% "Ball and Plate" system definitions file

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

% convert state-space objects to system structures
sysx = mpt_sys(css_x, Ts);
% input constraints
sysx.umax = 10;
sysx.umin = -10;
% output constraints
sysx.ymax = 30;
sysx.ymin = -30;
% state constraints
sysx.xmax = [30; 15; 15*pi/180; 1];
sysx.xmin = -sysx.xmax;

sysy = mpt_sys(css_y, Ts);
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

% load trajectory data
load star_40
