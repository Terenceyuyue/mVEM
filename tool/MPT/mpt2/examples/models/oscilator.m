%OSCILATOR 2nd order LTI system

clear sysStruct probStruct

%x(k+1)=Ax(k)+Bu(k)
sysStruct.A= [   0.5403   -0.8415; 0.8415    0.5403];
sysStruct.B=[ -0.4597; 0.8415];

%y(k)=Cx(k)+Du(k)
sysStruct.C=[1 0];
sysStruct.D=0;

%set constraints on output
sysStruct.ymin = -10;
sysStruct.ymax = 10;

%set constraints on input
sysStruct.umin    =   -1;
sysStruct.umax    =   1;

probStruct.norm=2;
probStruct.Q=eye(2);
probStruct.R=1;
probStruct.N=5;
probStruct.subopt_lev=0;