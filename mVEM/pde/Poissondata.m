function pde = Poissondata()

c = 1;
pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'g_N', @g_N, 'c', c);


% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = y.^2.*sin(pi*x);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = (pi^2*y.^2-2).*sin(pi*x) + c*uexact(p);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% Neumann boundary conditions ( right side hand )
    function val = g_N(p)
        x = p(:,1); y = p(:,2);
        val = [pi*y.^2.*cos(pi*x), 2*y.*sin(pi*x)];
        %val = pi*y.^2.*cos(pi*x);
    end
end
