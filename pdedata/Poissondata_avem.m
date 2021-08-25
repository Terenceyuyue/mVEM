function pde = Poissondata_avem(c)

if nargin == 0,  c = 1; end

% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs(c);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = u(x,y);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% Neumann boundary conditions ( right side hand )
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y), uy(x,y)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du,'c', c);
end

function [u,ux,uy,f] = compute_rhs(c)
    syms x y;
    % exact solution
    u = x.*y.*(1-x).*(1-y).*exp(-1000*((x-0.5).^2+(y-0.117).^2));
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



