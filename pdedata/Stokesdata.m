function pde = Stokesdata
% PDE: u = [u1,u2]
%     -nu*\Delta(u) - \nabla(p) = f,  in \Omega
%     div u = 0 ,  in \Omega
%     u = g_D,   on \partial \Omega

nu = 1;

% --------- given by the symbolic computation ------
[u1,u2,pe,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(nu);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y)+0*x,u2(x,y)+0*x];
    end
    function val = pexact(p)
        x = p(:,1); y = p(:,2);
        val = pe(x,y)+0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [f1(x,y)+0*x,f2(x,y)+0*x];
    end

% Dirichlet
    function val = g_D(p)
        val = uexact(p);
    end

% Du
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        Du1 = [u1x(x,y)+0*x, u1y(x,y)+0*x];
        Du2 = [u2x(x,y)+0*x, u2y(x,y)+0*x];
        val = [Du1, Du2];
    end

pde = struct('uexact',@uexact, 'pexact',@pexact, 'f',@f, ...
    'g_D',@g_D, 'Du',@Du, 'nu',nu);
end

function [u1,u2,pe,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(nu)
    syms x y;
    % exact solution: \Omega = [0,1]^2
    u1 = -0.5*cos(x)^2*cos(y)*sin(y);
    u2 = 0.5*cos(y)^2*cos(x)*sin(x);
    pe = sin(x)-sin(y);  % \int_Omega pe = 0
%     u1 = y^4 + 1;
%     u2 = x^4 + 2;
%     pe = x^3 - y^3;
    % derivative
    u1x = diff(u1,x);   u1y = diff(u1,y);
    u2x = diff(u2,x);   u2y = diff(u2,y);
    u1xx = diff(u1x,x); u1yy = diff(u1y,y);
    u2xx = diff(u2x,x); u2yy = diff(u2y,y);
    px = diff(pe,x);    py = diff(pe,y);
    % Lap(u)
    Lapu1 = u1xx + u1yy;
    Lapu2 = u2xx + u2yy;
    % f = -nu*Delta(u) - nabla(p)
    f1 = -nu*Lapu1 - px;
    f2 = -nu*Lapu2 - py;
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y});
    u2 = matlabFunction(u2,'Vars',{x,y});
    pe = matlabFunction(pe,'Vars',{x,y});
    f1 = matlabFunction(f1,'Vars',{x,y});
    f2 = matlabFunction(f2,'Vars',{x,y});
    u1x = matlabFunction(u1x,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});
    u1y = matlabFunction(u1y,'Vars',{x,y});
    u2y = matlabFunction(u2y,'Vars',{x,y});
end
