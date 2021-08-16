function pde = PlateBendingData4(varargin)

% data for plate
if nargin == 0
    t = 0.1; % plate thickness (m)
    E = 10920; % Young's modulus (Nm/kg^2)
    nu = 0.3;  % Poisson's ratio
    D = E*t^3/12/(1-nu^2); % flexural rigidity
    c = 0; % 弹性耦合常数
    para = struct('t',t, 'E',E, 'nu',nu, 'D',D, 'c',c);    
end

% exact solution
    function val =  uexact(p)
        x = p(:,1); y = p(:,2);
        val = (1-x-y).^2.*x.^2.*y.^2;
    end

% load data (right hand side function)
    function val =  f(p)
        x = p(:,1); y = p(:,2);
        val = 64*x.*y + 8*(x + y - 1).^2 + 32*x.^2 + 32*y.^2 + 16*x.*(2*x + 2*y - 2) ...
        + 16*y.*(2*x + 2*y - 2) + c*uexact(p);
    end

% derivative of the exact solution
    function val =  Du(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = x.^2.*y.^2.*(2*x + 2*y - 2) + 2*x.*y.^2.*(x + y - 1).^2;
        val(:,2) = x.^2.*y.^2.*(2*x + 2*y - 2) + 2*x.^2.*y.*(x + y - 1).^2;
    end

% Dirichlet boundary condition
    function val = g_D(p)
        val = uexact(p);
    end

pde = struct('para', para, 'f',@f, 'uexact',@uexact, 'g_D',@g_D, 'Du',@Du);
end