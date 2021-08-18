clc
echo on
%*********************************************************
%
% Moment relaxations
%
%*********************************************************
%
% This example shows how to solve some polynomial problems
% using the theory of moment-relxations.
pause
clc

% Define variables
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
pause

% Define a polynomial to be analyzed...
p = -2*x1+x2-x3;
pause

% and a set of polynomial inequalities (concave constraint)
F = set(x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>0);
F = F + set(4-(x1+x2+x3)>0);
F = F + set(6-(3*x2+x3)>0);
F = F + set(x1>0);
F = F + set(2-x1>0);
F = F + set(x2>0);
F = F + set(x3>0);
F = F + set(3-x3>0);
pause

% A lower bound {min p(x) | F(x)>0} can be found using solvemoment 
solvemoment(F,p);
pause
clc
% The lower bound can be be found by checking the 
% value of the RELAXED variables. When the moment relaxation
% is calculated, x1 and x1^2 etc. are treated as independant
% variables. To extract the values of the relaxed variables,
% use the command relaxdouble
relaxdouble(p)
[double([x1;x2;x3;x1^2;x2^2;x3^2]) relaxdouble([x1;x2;x3;x1^2;x2^2;x3^2])]

pause
clc

% Better lower bound can be obtained by using higher order 
% relaxations
pause
solvemoment(F,p,[],2);
relaxdouble(p)

pause
clc

% Another way to obtain better bounds is to add more valid constraints
% One simple trick is to use the current lower bound by adding
% the constraint p>double(p)
pause
solvemoment(F+set(p>double(p)),p,[],2);
relaxdouble(p)

pause
clc
% ...or square some linear constraints
pause
solvemoment(F+set(9>x3^2)+set(4>x1^2),p,[],2);
relaxdouble(p)

pause
clc
% or use an even higher relaxation
pause
solvemoment(F,p,[],4);
relaxdouble(p)

pause

% The value -4 is known to be the global optimum, but our relaxed solution
% is still not a valid solution.
checkset(F)
pause

% YALMIP can try to recover globally optimal solutions, but to 
% do this, three output arguments are needed. The global solutions,
% if sucessfully extracted, are returned in the second output.
pause
[solution,xoptimal] = solvemoment(F,p,[],4);

% Use the first extracted global solution.
setsdpvar([x1;x2;x3],xoptimal{1});
double(p)
checkset(F)

pause

% In this example, there are only three variables and they
% correpond to the variables in the output x.
%
% In principle, the code setsdpvar(recover([depends(p) depends(F)]),x{1})
% should work, but to be on the safe side in a more complex scenario, 
% a third output argument should be used to keep track of the variables.
pause
[solution,xoptimal,momentdata] = solvemoment(F,p,[],4);

% The variables used in the relaxation are returned in 
% the fourth output, in the field 'x'
momentdata
pause

setsdpvar(momentdata.x,xoptimal{1});
pause

double(p)
checkset(F)
pause

clc

% Polynomial semidefinite constraints can be addressed also
sdpvar x1 x2
p = -x1^2-x2^2;
F = set([1-4*x1*x2 x1;x1 4-x1^2-x2^2]);
pause
[sol,xoptimal] = solvemoment(F,p,[],2);
setsdpvar([x1;x2],xoptimal{1});
double(p)
checkset(F)

pause
clc

% Note that the moment relaxation solver can be called in 
% a fashion that is more consistent with the rest of the YALMIP framework.
%
% To do this, just use the solver tag 'moment'.
pause
solution = solvesdp(F,p,sdpsettings('solver','moment','moment.order',2));
setsdpvar(solution.momentdata.x,solution.xoptimal{1});
double(p)
checkset(F)
pause

clc

% The extraction of the global solution is numerically sensitive and
% can easily lead to low accuracy, or even completly wrong solutions.
%
% Some tweaking can be done to improve upon this.
%
% Consider the following problem with known optima (0,2) and (0,-2)
clear all
yalmip('clear')
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
obj = -x1^2-x2^2;
g1 = 5-4*x1*x2-x1^2-x2^2;
g2 = 4-16*x1*x2-2*x1^2-x2^2+4*x1^3*x2+4*x1*x2^3;
F = set([g1 g2]);
pause

% We solve and extract solutions (requires a rather high order to find a solution)
% Fo reasons that will become clear later, we also specify the tolerance
% for the Gaussian elimination, used during the extraction process.
pause
[sol,xe] = solvemoment(F,obj,sdpsettings('moment.rceftol',1e-6),7);

% The accuracy can easily become rather poor for this problem 
% NOTE : There is a randomization step in the extraction algorithm, 
% so these numbers may differ.
[xe{1} xe{2}]
pause

% To improve the solutions, YALMIP can perform a number of Newton
% iterations on a set of polynomial equalities that are solved using 
% numerical linear algebra techniques during the extraction algorithm. 
% (Note, these are not the original polynomials, but polynomials that 
% define the global solutions, given the moment matrices)
pause
[sol,xe] = solvemoment(F,obj,sdpsettings('moment.refine',5),7);
[xe{1} xe{2}]
pause

% The main reason for the poor accuracy in this problem is actually
% the bad conditioning during a Gaussian elimination. To improve 
% this situation, it is possible to let YALMIP select a tolerance using
% some heuristics.
%
% This can be done by specifying the tolerance rceftol to -1. 
% In fact, this is the default setting, so we just run standard settings!
pause
[sol,xe] = solvemoment(F,obj,[],7);
[xe{1} xe{2}]
pause

% For this problem we needed a rather high relaxation order to be able to 
% extract the solution.
%
% It is possible to tell YALMIP to try to extract a solution, even though
% YALMIP belives this is impossible.
%
% To force YALMIP to try to extract global solutions, set the option 
% moment.extractrank to the desired number of solutions.
%
% Let us try to extract a solution from a 6th order relaxation.
pause
[sol,xe] = solvemoment(F,obj,sdpsettings('moment.extractrank',2),6);
[xe{1} xe{2}]
pause

% Hmm, that's poor. No wonder, YALMIP did not expect it to be possible to
% extract these solutions, so numerical problems are very likely during the
% extraction process. Let us try to refine it.
pause
[sol,xe] = solvemoment(F,obj,sdpsettings('moment.extractrank',2,'moment.refine',5),6);
[xe{1} xe{2}]
pause

% Pretty annoying to resolve the problem when you just want to change a
% setting in the extraction algorithm, right?
%
% To avoid this, use 4 output arguments
[sol,xe,momentdata] = solvemoment(F,obj,sdpsettings('moment.extractrank',2,'moment.refine',5),6);
pause

% The variable momentdata contains all information needed to extract
% solutions.
%
% Hence, we can compare solutions for different settings easily.
xe0 = extractsolution(momentdata,sdpsettings('moment.extractrank',2,'moment.refine',0));
xe1 = extractsolution(momentdata,sdpsettings('moment.extractrank',2,'moment.refine',2));
xe2 = extractsolution(momentdata,sdpsettings('moment.extractrank',2,'moment.refine',5));
[xe0{1} xe0{2} xe1{1} xe1{2} xe2{1} xe2{2}]
pause

echo off

