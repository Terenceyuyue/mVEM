function pde = Poisson3data(c)

if nargin == 0,  c = 0; end

% --------- given by the symbolic computation ------
[u,ux,uy,uz,rhs] = compute_rhs(c);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = u(x,y,z) + 0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = rhs(x,y,z) + 0*x;
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p) + 0*p(:,1);
    end
% Neumann boundary conditions
    function val = Du(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        val = [ux(x,y,z) + 0*x, uy(x,y,z) + 0*x, uz(x,y,z) + 0*x];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du,'c', c);
end

function [u,ux,uy,uz,f] = compute_rhs(c)
    syms x y z;
    % exact solution
    u = sin(2*x*y)*cos(z); % example 1
    %u = x^2+y^2+z^2; % example 2
    %u = sin(pi*x)*cos(pi*y)*cos(pi*z); % example 3   
    
    % derivative
    ux = diff(u,x);  uy = diff(u,y);  uz = diff(u,z);
    % Lap(u)
    Lapu = diff(u,x,2)+diff(u,y,2)+diff(u,z,2);
    % f
    f = -Lapu + c*u;
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y,z});    
    f = matlabFunction(f,'Vars',{x,y,z});    
    ux = matlabFunction(ux,'Vars',{x,y,z});  
    uy = matlabFunction(uy,'Vars',{x,y,z}); 
    uz = matlabFunction(uz,'Vars',{x,y,z}); 
end
