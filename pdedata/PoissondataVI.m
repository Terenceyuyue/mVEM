function pde = PoissondataVI(c)

if nargin == 0,  c = 0; end

% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs(c);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = u(x,y) + 0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y) + 0*x;
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p) + 0*p(:,1);
    end
% Neumann boundary conditions ( right side hand )
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y) + 0*x, uy(x,y) + 0*x];
    end

% Frictional condition
    function val = g_C(p)
        x = p(:,1); y = p(:,2);
        val = sqrt(ux(x,y).^2 + uy(x,y).^2) + 0*x;
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du,...
        'c', c, 'g_C',@g_C);
end

function [u,ux,uy,f] = compute_rhs(c)
    syms x y;
    % exact solution
    u = (sin(x)-x*sin(1))*sin(2*pi*y);
    %u = sin(pi*x)*sin(pi*y);
    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    % Lap(u)
    Lapu = diff(u,x,2)+diff(u,y,2);
    % f
    f = -Lapu + c*u;
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y});    
    ux = matlabFunction(ux,'Vars',{x,y});  uy = matlabFunction(uy,'Vars',{x,y}); 
end
