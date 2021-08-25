function pde = PlateBendingData3(varargin)

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
        val = (sin(pi*x).*sin(pi*y)).^2;
    end

% load data (right hand side function)
    function val =  f(p)
        x = p(:,1); y = p(:,2);
        val = 24*pi^4*sin(pi*x).^2.*sin(pi*y).^2 + 8*pi^4*cos(pi*x).^2.*cos(pi*y).^2 ...
            - 16*pi^4*cos(pi*x).^2.*sin(pi*y).^2 - 16*pi^4*cos(pi*y).^2.*sin(pi*x).^2 ...
            + c*uexact(p);
    end

% derivative of the exact solution
    function val =  Du(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = 2*pi*cos(pi*x).*sin(pi*x).*sin(pi*y).^2;
        val(:,2) = 2*pi*cos(pi*y).*sin(pi*x).^2.*sin(pi*y);
    end

% Dirichlet boundary condition
    function val = g_D(p)
        val = uexact(p);
    end

pde = struct('para', para, 'f',@f, 'uexact',@uexact, 'g_D',@g_D, 'Du',@Du);
end