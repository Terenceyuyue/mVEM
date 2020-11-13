clear sysStruct probStruct
Double_Integrator

echo on
% Demonstration of time-varying models. In this example we enforce time-varying
% constraints on system outputs.

S1 = sysStruct; S1.ymax = [5; 5]; S1.ymin = [-5; -5];
S2 = sysStruct; S2.ymax = [4; 4]; S2.ymin = [-4; -4];
S3 = sysStruct; S3.ymax = [3; 3]; S3.ymin = [-3; -3];
S4 = sysStruct; S4.ymax = [2; 2]; S4.ymin = [-2; -2];

% create the time-varying model
model = {S1, S2, S3, S4};

% prediction horizon must be set to the number of models, i.e. to 4 in this
% example:
probStruct.N = 4;

% calculate an explicit controller
ctrl = mpt_control(model, probStruct);

% calculate the open-loop trajectory
[X, U, Y] = sim(ctrl, [5; 0], struct('openloop', 1));
Y

% Notice that each element of the system outputs satisfies constraints given at
% each step of the prediction.

echo off
