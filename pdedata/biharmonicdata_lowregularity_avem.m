function pde = biharmonicdata_lowregularity_avem()

alpha = 4/3;  beta = 4/3;
eps = 0;
% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        alpha = 4/3; beta = 4/3;
        r = sqrt(x.^2+y.^2 + 1e-5);
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        val = r.^alpha.*sin(beta*theta);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        r = sqrt(x.^2+y.^2+eps);
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        val = (alpha^2-beta^2)*(alpha^2-beta^2-4*alpha+4)*r.^(alpha-4).*sin(beta*theta);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% derivative
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        r = sqrt(x.^2+y.^2 + 1e-5);
        theta = atan2(y, x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        ur = alpha * r.^(alpha-1) .* sin(beta * theta);
        utheta = r.^alpha * beta .* cos(beta * theta);
        ux = ur.* (x./r) + utheta.*(-y./(r.^2));
        uy = ur.* (y./r) + utheta.*(x./(r.^2));
        val = [ux,uy];
    end
pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du);
end


