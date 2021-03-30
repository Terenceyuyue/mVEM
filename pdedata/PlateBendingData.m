function pde = PlateBendingData(varargin)

% data for plate
if nargin == 0
    t = 0.1; % plate thickness (m)
    E = 10920; % Young's modulus (Nm/kg^2)
    nu = 0.3;  % Poisson's ratio (nu \neq 1)
    D = E*t^3/12/(1-nu^2); % flexural rigidity
    c = 0; % 弹性耦合常数
    para = struct('t',t, 'E',E, 'nu',nu, 'D',D, 'c',c);    
end

% --------- given by the symbolic computation ------
[u,rhs,ux,uy,uxx,uxy,uyy] = compute_rhs(para);

% exact solution
    function val =  uexact(p)
        x = p(:,1); y = p(:,2);
        %val =  sin(2*pi*x).*cos(2*pi*y);
        val = u(x,y);
    end

% load data (right hand side function)
    function val =  f(p)
        x = p(:,1); y = p(:,2);
        %val =  D*64*pi^4*sin(2*pi*x).*cos(2*pi*y) + c*uexact(p);
        val = rhs(x,y);
    end

% derivative of the exact solution
    function val =  Du(p)
        x = p(:,1); y = p(:,2);
        %val(:,1) = 2*pi*cos(2*pi*x).*cos(2*pi*y);
        %val(:,2) = -2*pi*sin(2*pi*x).*sin(2*pi*y);
        val = [ux(x,y), uy(x,y)];
    end

% Dirichlet boundary condition
    function val = g_D(p)
        val = uexact(p);
    end

% second order derivative of the exact solution
    function val =  DDu(p) 
        x = p(:,1); y = p(:,2);
        val(:,1) = uxx(x,y);          val(:,2) = uxy(x,y);  
        val(:,3) = uxy(x,y);          val(:,4) = uyy(x,y); 
    end

pde = struct('para', para, 'f',@f, 'uexact',@uexact, 'g_D',@g_D, 'Du',@Du, 'DDu',@DDu);
end

function [u,f,ux,uy,uxx,uxy,uyy] = compute_rhs(para)    
    syms x y;    
    
    % exact solution
    u = sin(2*pi*x).*cos(2*pi*y);

    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    
    % second derivative
    uxx = diff(ux,x); uxy = diff(ux,y); uyy = diff(uy,y);

    % f
    Lapu = uxx + uyy;
    LapLapu = diff(Lapu,x,2) + diff(Lapu,y,2);
    f = para.D*LapLapu + para.c*u;
    
    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y});   
    ux = matlabFunction(ux,'Vars',{x,y});  uy = matlabFunction(uy,'Vars',{x,y});
    
    uxx = matlabFunction(uxx,'Vars',{x,y});  
    uxy = matlabFunction(uxy,'Vars',{x,y});  
    uyy = matlabFunction(uyy,'Vars',{x,y});
end