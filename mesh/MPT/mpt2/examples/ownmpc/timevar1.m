clear sysStruct probStruct
opt_sincos

echo on
% Demonstration of time-varying models

% "PWA" here represents a PWA system with 2 dynamics. We create an LTI
% model of the system by averaging the two dynamics:
PWA = sysStruct;
LTI.A = (PWA.A{1} + PWA.A{2})/2;
LTI.B = (PWA.B{1} + PWA.B{2})/2;
LTI.C = (PWA.C{1} + PWA.C{2})/2;
LTI.D = (PWA.D{1} + PWA.D{2})/2;
LTI.ymax = PWA.ymax;
LTI.ymin = PWA.ymin;
LTI.xmax = PWA.xmax;
LTI.xmin = PWA.xmin;
LTI.umax = PWA.umax;
LTI.umin = PWA.umin;

% for the first prediction we use the PWA model, for the rest we assume an LTI
% dynamics:
model = {PWA, LTI, LTI};

% prediction horizon must be set to the number of models, i.e. to 3 in this
% example:
probStruct.N = 3;

% Use linear cost function
probStruct.norm = 1;

% calculate an explicit controller
ctrl = mpt_control(model, probStruct);

% plot the controller partition
plot(ctrl)

% plot value of control actions
figure; plotu(ctrl);

% simulate the closed-loop system
figure; simplot(ctrl, [5; 5], 10);

echo off
