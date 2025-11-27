function pde = NavierStokesdata_cavity(Re)
% PDE: u = [u1,u2]
%     -nu*\Delta(u) + (u*\nabla)u - \nabla(p) = f,  in \Omega
%     div u = 0 ,  in \Omega
%     u = g_D,   on \partial \Omega

nu = 1/Re;

% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [0*x,0*y];
    end

% Dirichlet
    function val = g_D(p)
        x = p(:,1); y = p(:,2);
        u1 = 1*(abs(y-1)<1e-4);  u2 = 0*(abs(y-1)<1e-4);
        val = [u1, u2];
    end

rho = 1.2;  alpha = 7*Re/10;
pde = struct('f',@f, 'g_D',@g_D, 'nu',nu, 'rho', rho, 'alpha',alpha);
end


