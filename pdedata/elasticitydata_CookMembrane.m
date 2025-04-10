function pde = elasticitydata_CookMembrane(varargin)

% ------ Lame constants ------
if nargin==0
    E = 250;  % Young's modulus
    nu = 0.3;  % Poisson's ratio
else
    [E, nu] = deal(varargin{:});
end
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu)); 


para = struct('lambda',lambda, 'mu',mu);
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [0*x, 0*y];
    end
    function val = g_D(p)
        x = p(:,1); y = p(:,2);
        val = [0*x, 0*y];
    end
    function val = g_N(p)
        x = p(:,1); 
        g = 6.25;  % vertical traction force or shear load per length 
        % The shear load is imposed on the right side of the membrane
        % sigma*n = [0; g]   gives  sig11 = 0, sig12 = g  since n = [1 0].
        val = [0*x, 0*x, g+0*x]; % [sig11,sig22,sig12]
    end
pde = struct('para',para, 'f', @f, 'g_D',@g_D, 'g_N',@g_N);
end