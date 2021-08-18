function [xopt,lambda,flag,fval]=esp_LP(f,H,He,x0)
%
% [xopt,lambda,flag,fval]=esp_LP(f,H,He,x0)
% 
% Solve the LP min f'x s.t. H*[x;-1]<0, He*[x;-1]=0
%
% Currently a wrapper for the MPT toolbox
%
% xopt      - The optimizer
% fval      - Value of the objective
% lambda    - Vector of Lagrangian multipliers
% exitflag  - An integer specifying result of the optimization
%          1 = optimal
%         -1 = infeasible
%         -2 = unbounded
% how       - States the result of optimization ('ok', 'unbounded', 'infeasible')
%
%
%  NOTE: Defaulting to linprog!
%

global mptOptions

if(nargin < 4) x0 = []; end;
if(nargin < 3) He = []; end;

[A,B]     = a2s(H);
[Aeq,Beq] = a2s(He);
% use the "rescue" mode which guarantees that if one solver fails, other LP
% solver will be used automatically. greatly improves reliability!
[xopt,fval,lambda,exitflag,how]=mpt_solveLPs(f,A,B,Aeq,Beq,x0,mptOptions.lpsolver);
lambda = lambda(1:length(B));

% Ensure we get the flags that esp is expecting independent of the lp solver used
if(strcmp(how,'ok') == 1)
  flag = 1;
elseif(strcmp(how,'unbounded') == 1)
  flag = -2;
elseif(strcmp(how,'infeasible') == 1)
  flag = -1;
else
  flag = -1;
end;
