clear sysStruct probStruct
close all
oscilator

% In this demo we add time-varying constraints on the system output. Such
% constraints can be changed on-the-fly, since they are part of the state
% vector.

sysStruct.A = [sysStruct.A zeros(2, 1); zeros(1, 2) 0];
sysStruct.B = [sysStruct.B; 0];
sysStruct.C = [1 0 0];
sysStruct.D = 1;

sysStruct.xmax = [5; 5; 5];
sysStruct.xmin = [-5; -5; -5];
sysStruct.ymax = 5;
sysStruct.ymin = -5;

probStruct.N = 3;
probStruct.Q = eye(3);
probStruct.Tconstraint = 0;

[C, O, V] = mpt_ownmpc(sysStruct, probStruct);

for k = 1:length(V.y),
    C = C + set(-V.x{1}(3) <= V.y{k} <= V.x{1}(3));
end

ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

x0 = [3; 0];

ymax = 2
[X, U, Y] = sim(ctrl, [x0; ymax], struct('openloop', 1)); Y

ymax = 2.1
[X, U, Y] = sim(ctrl, [x0; ymax], struct('openloop', 1)); Y

ymax = 2.2
[X, U, Y] = sim(ctrl, [x0; ymax], struct('openloop', 1)); Y
