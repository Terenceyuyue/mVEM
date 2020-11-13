clear sysStruct probStruct

% import system from HYSDEL source
sysStruct = mpt_sys(which('two_tanks.hys'));

% state constraints (optional)
sysStruct.xmax = [0.62; 0.3];
sysStruct.xmin = [0; 0];

% output constraints (mandatory)
sysStruct.ymax = [0.3];
sysStruct.ymin = [0];

% input constraints (mandatory)
sysStruct.umax = [1; 1];
sysStruct.umin = [0; 0];

% sampling time (if not given, assuming Ts = 1)
sysStruct.Ts = 10;


% 1-norm type cost function
probStruct.norm = 1;

% prediction horizon
probStruct.N = 3;

% optimal solution
probStruct.subopt_lev = 0;

% penalty on states
probStruct.Q = eye(2);

% penalty on inputs
probStruct.R = diag([1e-5,1e-5]);

% penalty on outputs
probStruct.Qy = 200;

% system output will be regulated to this value
probStruct.yref = 0.2;



% compute the explicit controller
fprintf('Computing an explicit controller...\n\n');
expc = mpt_control(sysStruct, probStruct);

% compute the on-line controller
fprintf('Computing an on-line controller...\n\n');
onlc = mpt_control(sysStruct, probStruct, 'online');

% plot regions of the explicit controller
plot(expc);

% plot control action defined over regions
figure;
mpt_plotU(expc);

% plot closed-loop trajectories starting from given initial conditions
figure;
x0 = [0.5; 0];
[Xexp, Uexp, Yexp] = mpt_plotTimeTrajectory(expc, x0);
title('Simulation of the explicit controller');
figure;
[Xonl, Uonl, Yonl] = mpt_plotTimeTrajectory(onlc, x0);
title('Simulation of the on-line controller');


% open simulink model and start simulation
open_system('twotankssim');
sim('twotankssim', 1000);
