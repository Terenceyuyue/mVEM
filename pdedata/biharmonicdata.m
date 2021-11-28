function pde = biharmonicdata
%
%     f = 64*pi^4*sin(2*pi*x)*cos(2*pi*y);
%     u = sin(2*pi*x)*cos(2*pi*y);
%     Du = (2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y));
%     w = -Delta u = 8*pi^2*sin(2*pi*x)*cos(2*pi*y);
%     du/dn = g_N on [0,1]^2.
%
pde = struct('f',@f,'uexact',@uexact,'Du',@Du,'wexact',@wexact,'Dw',@Dw,...
    'g_D',@g_D);

% right hand side function
    function val =  f(p)
        x = p(:,1); y = p(:,2);
        val =  64*pi^4*sin(2*pi*x).*cos(2*pi*y);
    end

% exact solution
    function val =  uexact(p)
        x = p(:,1); y = p(:,2);
        val =  sin(2*pi*x).*cos(2*pi*y);
    end

% derivative of the exact solution
    function val =  Du(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = 2*pi*cos(2*pi*x).*cos(2*pi*y);
        val(:,2) = -2*pi*sin(2*pi*x).*sin(2*pi*y);
    end

% exact solution of w = -Delta u
    function val =  wexact(p)
        x = p(:,1); y = p(:,2);
        val = 8*pi^2*sin(2*pi*x).*cos(2*pi*y);
    end

% derivative of the exact solution w
    function val =  Dw(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = 16*pi^3*cos(2*pi*x).*cos(2*pi*y);
        val(:,2) = -16*pi^3*sin(2*pi*x).*sin(2*pi*y);
    end

% Dirichlet boundary condition
    function val = g_D(p)
        val = uexact(p);
    end

end