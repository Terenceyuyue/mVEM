function pde = NavierStokesdata_Kovasznay(Re)
% PDE: u = [u1,u2]
%     -nu*\Delta(u) + (u*\nabla)u - \nabla(p) = f,  in \Omega
%     div u = 0 ,  in \Omega
%     u = g_D,   on \partial \Omega

nu = 1/Re;
lam = Re/2 - sqrt(Re^2/4+4*pi^2);

% --------- given by the symbolic computation ------
[u1,u2,pe,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(nu,lam);

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

function [u1,u2,pe,f1,f2,u1x,u1y,u2x,u2y] = compute_rhs(nu,lam)
    syms x y;
    % exact solution: 
    % \Omega = [-0.5, 1] \times [-0.5, 1.5]
    u1 = 1-exp(lam*x)*cos(2*pi*y);
    u2 = lam/(2*pi)*exp(lam*x)*sin(2*pi*y);
    pe = -0.5*exp(2*lam*x); 
    p0 = -0.5*2*(exp(2*lam*1)-exp(-2*lam*0.5))/(2*lam);
    p0 = p0/(1.5*2);
    pe = pe - p0;
 
    % derivative
    u1x = diff(u1,x);   u1y = diff(u1,y);
    u2x = diff(u2,x);   u2y = diff(u2,y);
    u1xx = diff(u1x,x); u1yy = diff(u1y,y);
    u2xx = diff(u2x,x); u2yy = diff(u2y,y);
    px = diff(pe,x);    py = diff(pe,y);
    % Lap(u)
    Lapu1 = u1xx + u1yy;
    Lapu2 = u2xx + u2yy;
    % nonlinear term
    N1 = u1*u1x + u2*u1y;
    N2 = u1*u2x + u2*u2y;
    % f = -nu*Delta(u) - nabla(p)
    f1 = -nu*Lapu1 + N1 - px;
    f2 = -nu*Lapu2 + N2 - py;
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


