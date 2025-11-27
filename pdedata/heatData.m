function pde = heatData

% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs;

% exact solution
    function val = uexact(p,t)
        x = p(:,1); y = p(:,2); 
        val = u(x,y,t) + 0*x;
    end
% right side hand function
    function val = f(p,t)
        x = p(:,1); y = p(:,2); 
        val = rhs(x,y,t) + 0*x;
    end
% Dirichlet boundary conditions
    function val = g_D(p,t)
        val = uexact(p,t);
    end
% for Neumann boundary conditions
    function val = Du(p,t)
        x = p(:,1); y = p(:,2); 
        val = [ux(x,y,t) + 0*x, uy(x,y,t) + 0*x];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du',@Du);
end

function [u,ux,uy,f] = compute_rhs()
    syms x y t;
  
    % exact solution
    u = exp(t)*sin(pi*x)*sin(pi*y);  
    %u = exp(-x-y-t);    
    %u = exp(-t)*sin(2*x+0.5)*cos(y+0.3)+log(1+x*y);
    % derivative
    ut = diff(u,t);
    ux = diff(u,x);  uy = diff(u,y);
    Lapu = diff(ux,x) + diff(uy,y);
    % f
    f = ut-Lapu;

    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y,t});    
    f = matlabFunction(f,'Vars',{x,y,t});    
    ux = matlabFunction(ux,'Vars',{x,y,t});  
    uy = matlabFunction(uy,'Vars',{x,y,t}); 
end
