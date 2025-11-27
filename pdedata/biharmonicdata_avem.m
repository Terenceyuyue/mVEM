function pde = biharmonicdata_avem(id)

if nargin==0, id = 1; end

% --------- given by the symbolic computation ------
[u,ux,uy,uxx,uxy,uyy,rhs] = compute_rhs(id);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = u(x,y)+0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y)+0*x;
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% derivative
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y)+0*x, uy(x,y)+0*x];
    end
% 2-order derivative
    function val = DDu(p)
        x = p(:,1); y = p(:,2);
        val = [uxx(x,y)+0*x, uxy(x,y)+0*x, uxy(x,y)+0*x, uyy(x,y)+0*x];
    end


pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du, 'DDu', @DDu);
end

function [u,ux,uy,uxx,uxy,uyy,f] = compute_rhs(id)
    syms x y; 
    switch id
        case 1
            u = 10*x^2*y^2*(1-x)^2*(1-y)^2*sin(pi*x) + x^2+y^2;
        case 2
            u = x.*y.*(1-x).*(1-y).*exp(-1000*((x-0.5).^2+(y-0.117).^2));
        case 3
            eps = 1e-4;
            u = ((x-0.5)^2+(y-0.5)^2+eps)^(5/6);
    end   
    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    uxx = diff(ux,x);    uxy = diff(ux,y);    uyy = diff(uy,y);
    Lapu = uxx + uyy;
    LapLapu = diff(Lapu,x,2) + diff(Lapu,y,2);
    % f
    f = LapLapu;    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});
    ux = matlabFunction(ux,'Vars',{x,y});
    uy = matlabFunction(uy,'Vars',{x,y});
    uxx = matlabFunction(uxx,'Vars',{x,y});
    uxy = matlabFunction(uxy,'Vars',{x,y});   
    uyy = matlabFunction(uyy,'Vars',{x,y}); 
    f = matlabFunction(f,'Vars',{x,y});
end
