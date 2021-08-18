function [upper,x_min] = rounder(p,relaxedsolution)

% Extremely simple heuristic for finding integer
% solutions.
%
% Rounds up and down, fixes etc.

% This was the relaxed solution
x = relaxedsolution.Primal;

% Assume we fail
upper = inf;
x_min = x;

% These should be integer
intvars = [p.integer_variables(:);p.binary_variables(:)];
p.options.rounding = [1 1 1 1];

% b=p.F_struc(:,1);
% A=-p.F_struc(:,2:end);
% c=p.c;
% 
% %x=rand(length(p.c),1)
% xtilde = round(x(intvars));
% xplus=sdpvar(length(xtilde),1);
% xminus=sdpvar(length(xtilde),1);
% xj=sdpvar(length(xtilde),1);
% t=sdpvar(1,1);
% 
% while 1
%     goneup   = find(xtilde==p.ub(intvars));
%     gonedown = find(xtilde==p.lb(intvars));
%     F=set([]);%xplus>0)+set(xminus>0)+set(xj==xtilde+xplus-xminus);
%     F=F + set(b(1:p.K.l) -A(1:p.K.l,:)*[t;xj]>0);
%     F = F + set(reshape(b(p.K.l+1:end) -A(p.K.l+1:end,:)*[t;xj],14,14) > 0)
%     solvesdp(F+set(0<xj<1),-sum(xj(goneup))+sum(xj(gonedown)));%+sum(xplus)+sum(xminus));
% 
%     xtilde=round(double(xj));
% end

if p.options.rounding(1)
    % Round, update nonlinear terms, and compute feasibility
    xtemp = x;xtemp(intvars) = round(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = setnonlinearvariables(p,xtemp);  
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = computecost(p.f,p.corig,p.Q,x_min,p);%p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return
    end
end

if p.options.rounding(2)
    % Do same using fix instead
    xtemp = x;xtemp(intvars) = fix(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));    
    xtemp = setnonlinearvariables(p,xtemp);  
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
         upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return
    end
end

if p.options.rounding(3)
    % ...or ceil
    xtemp = x;xtemp(intvars) = ceil(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));    
    xtemp = setnonlinearvariables(p,xtemp);  
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
         upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return
    end
end

if p.options.rounding(4)
    % or floor
    xtemp = x;xtemp(intvars) = floor(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));    
    xtemp = setnonlinearvariables(p,xtemp);  
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
         upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return
    end
end