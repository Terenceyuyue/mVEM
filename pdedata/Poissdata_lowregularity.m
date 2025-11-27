function pde = Poissdata_lowregularity()
%% Lshape Problem
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% Detailed explanation can be found in 
% https://www.math.uci.edu/~chenlong/ifemdoc/afem/afemdoc.html

c = 0;

% exact solution
    function val = uexact(p) % exact solution
    eps = 1e-5;
    r = sqrt(sum(p.^2,2) + eps);
    theta = atan2(p(:,2),p(:,1));
    theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
    val = r.^(2/3).*sin(2*theta/3);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); 
        val = 0*x;
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% Neumann boundary conditions ( right side hand )
    function val = Du(p) % exact solution
    x = p(:,1); y = p(:,2);
    eps = 1e-5;
    r = sqrt(sum(p.^2,2) + eps);
    theta = atan2(y,x);
    theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
    val(:,1) = 2/3*r.^(-1/3).*sin(2*theta/3).*x./r ...
                - 2/3*r.^(2/3).*cos(2*theta/3).*y./r.^2;
    val(:,2) = 2/3*r.^(-1/3).*sin(2*theta/3).*y./r ...
                + 2/3*r.^(2/3).*cos(2*theta/3).*x./r.^2;
    end
pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du,'c', c);
end