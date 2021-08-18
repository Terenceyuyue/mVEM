clear sysStruct probStruct
Double_Integrator

echo on
% Stability analysis of autonomous systems. In this example we consider an LTI
% system and we would like to find a Piecewise Quadratic Lyapunov function which
% would show stability of the autonomous system
%
% First we make the system autonomous by stabilizing the dynamics by an LQR
% controller
K = mpt_dlqr(sysStruct.A, sysStruct.B, probStruct.Q, probStruct.R);
sysStruct.A = sysStruct.A - sysStruct.B*K;

% Then we set the B matrices to zero to indicate that this system is autonumous
sysStruct.B = zeros(2, 1);

% Try to calculate a Piecewise Quadratic Lyapunov function
lyap = mpt_lyapunov(sysStruct, 'pwq');

% PWQ Lyapunov function was found, the system is stable
lyap.details.lyapunov.feasible

echo off
