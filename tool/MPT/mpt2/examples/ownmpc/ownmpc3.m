Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 4
probStruct.N = 4;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Now we want to add following move-blocking constraints:
% 1. u_0==u_1
% 2. (u_1-u_2)==(u_2-u_3)
% 3. u_3==K*x_3
%
% where K is the solution to the Ricatti equation

K = mpt_dlqr(sysStruct.A, sysStruct.B, probStruct.Q, probStruct.R);

% Add constraint u(0)==u(1)
C = C + set(V.u{1} == V.u{2});

% Add constraint (u(1)-u(2))==(u(2)-u(3))
C = C + set((V.u{2} - V.u{3}) == (V.u{3} - V.u{4}));

% Add constraint u(3)==K*x(3)
C = C + set(V.u{4} == K * V.x{4});

% Compute an on-line controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the open-loop system
x0 = [3; 0];
[X, U] = sim(ctrl, x0, struct('openloop', 1))

% Verify that u_0==u_1
abs(U(1) - U(2)) < 1e-14

% Verify that (u_1-u_2)==(u_2-u_3)
abs((U(2)-U(3)) - (U(3)-U(4))) < 1e-14

% Verify that u_3==K*x_3
abs(U(4) - K*X(4, :)') < 1e-14


echo off
