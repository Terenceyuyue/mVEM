% example with mpqp degeneracy

sysStruct.A = {[0.7 0.74; -0.44 -0.89], [-0.78 0.74; 0.81 -0.89]};
sysStruct.B = {[-0.64; -0.51], [-0.64; -0.51]};
sysStruct.C = {eye(2), eye(2)};
sysStruct.D = {[0;0], [0;0]};
sysStruct.guardX = { [1 0], [-1 0] };
sysStruct.guardC = { 0, 0};
sysStruct.umin = -1;
sysStruct.umax = 1;
sysStruct.ymax = 5*[1; 1];
sysStruct.ymin = 5*[-1; -1];

probStruct.subopt_lev = 1;
probStruct.N = Inf;
probStruct.Q = eye(2);
probStruct.R = 1;
probStruct.norm=2;
