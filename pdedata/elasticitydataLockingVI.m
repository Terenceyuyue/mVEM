function pde = elasticitydataLockingVI(para)

% example for no slipping condition

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

% --------- given by the symbolic computation ------
[u1,u2,u1x,u1y,u2x,u2y,sig11,sig12,sig22,f11,f12,g] = compute_rhs(lambda,mu);

para = struct('lambda',lambda, 'mu',mu);
    function val = f1(p)
        x = p(:,1); y = p(:,2);
        val = [f11(x,y)+0*x,f12(x,y)+0*x];
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y)+0*x, u2(x,y)+0*x];
    end
    function val = g_D(p)
        val = uexact(p);
    end
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        Du1 = [u1x(x,y)+0*x, u1y(x,y)+0*x];
        Du2 = [u2x(x,y)+0*x, u2y(x,y)+0*x];
        val = [Du1, Du2];
    end
    function val = g_N(p)
        x = p(:,1); y = p(:,2);
        val = [sig11(x,y)+0*x,sig22(x,y)+0*x,sig12(x,y)+0*x];
    end
    function val = g_C(p)
        x = p(:,1); y = p(:,2);
        val = g(x,y)+0*x;
    end
pde = struct('para',para, 'f1', @f1, 'uexact',@uexact,...
    'g_D',@g_D, 'g_N',@g_N, 'g_C',@g_C, 'Du',@Du);
end

function [u1,u2,u1x,u1y,u2x,u2y,sig11,sig12,sig22,f11,f12,g] ...
    = compute_rhs(lambda,mu)
    syms x y;
    % exact solution
    u1 = (-1+cos(2*pi*x))*sin(2*pi*y) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);
    u2 = -(-1+cos(2*pi*y))*sin(2*pi*x) + 1/(1+lambda)*sin(pi*x)*sin(pi*y);
    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);
    u2x = diff(u2,x);      u2y = diff(u2,y);
    % Eij
    E11 = u1x;
    E12 = 0.5*(u1y+u2x);
    E22 = u2y;
    % sigma(u)
    sig11 = (lambda+2*mu)*E11 + lambda*E22;
    sig12 = 2*mu*E12; % sig21 = sig12
    sig22 = (lambda+2*mu)*E22 + lambda*E11;
    % f1 = [f11,f12]
    f11 = -(diff(sig11,x) + diff(sig12,y));
    f12 = -(diff(sig12,x) + diff(sig22,y));
    % g
    g = 4*sqrt(sig11.^2 + sig12.^2 + sig12.^2 + sig22.^2);
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y});        u2 = matlabFunction(u2,'Vars',{x,y});
    u1x = matlabFunction(u1x,'Vars',{x,y});      u1y = matlabFunction(u1y,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});      u2y = matlabFunction(u2y,'Vars',{x,y});
    sig11 = matlabFunction(sig11,'Vars',{x,y});
    sig12 = matlabFunction(sig12,'Vars',{x,y});
    sig22 = matlabFunction(sig22,'Vars',{x,y});
    f11 = matlabFunction(f11,'Vars',{x,y});     f12 = matlabFunction(f12,'Vars',{x,y});
    g = matlabFunction(g,'Vars',{x,y});
end