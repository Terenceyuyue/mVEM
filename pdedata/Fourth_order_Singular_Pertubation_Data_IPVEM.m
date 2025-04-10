function pde = Fourth_order_Singular_Pertubation_Data_IPVEM(id,varargin)

% data for plate
if nargin == 0
    epsilon = 1;    
else 
    epsilon = varargin{1};
end
para = struct('epsilon',epsilon);  

[u,rhs,ux,uy,uxx,uxy,uyy] = compute_rhs(para,id);

% exact solution
    function val =  uexact(p)
        x = p(:,1); y = p(:,2);
        val =  u(x,y);
    end

% load data (right hand side function)
    function val =  f(p)
        x = p(:,1); y = p(:,2);
        val =  rhs(x,y);
    end

% derivative of the exact solution
    function val =  Du(p)
        x = p(:,1); y = p(:,2);
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
function [u,f,ux,uy,uxx,uxy,uyy] = compute_rhs(para,id)    
    eps = para.epsilon;
    syms x y;    
    
    % exact solution    
    switch id
        case 1
            u = 10*x^2*y^2*(1-x)^2*(1-y)^2*sin(pi*x);
        case 2
            u = sin(pi*x)^2*sin(pi*y)^2;
        case 3
            u = eps*(exp(-x/eps) + exp(-y/eps)) - x^2*y^3;
        case 4
            iota = eps;
            u = (exp(sin(pi*x))-1-pi*iota*(cosh(1/(2*iota))-cosh((2*x-1)/(2*iota)))/sinh(1/(2*iota)))...
        *(exp(sin(pi*y))-1-pi*iota*(cosh(1/(2*iota))-cosh((2*y-1)/(2*iota)))/sinh(1/(2*iota)));
           % u = (sin(pi*x)-pi*iota*(cosh(1/(2*iota))-cosh((2*x-1)/(2*iota)))/sinh(1/(2*iota)))...
        %*(sin(pi*y)-pi*iota*(cosh(1/(2*iota))-cosh((2*y-1)/(2*iota)))/sinh(1/(2*iota)));
        case 5
            u = sin(pi*x)*sin(pi*y); % Poisson
            eps = 0;
    end   

    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    
     % second derivative
    uxx = diff(ux,x); uxy = diff(ux,y); uyy = diff(uy,y);
    
    % Laplace u
    Lap = uxx+uyy;

    % Laplace^2 u
    LapLap = diff(Lap,x,2) + diff(Lap,y,2);

    % f
    f = eps^2*LapLap - Lap; %f = simplify(f);
    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y}); 
    ux = matlabFunction(ux,'Vars',{x,y});  uy = matlabFunction(uy,'Vars',{x,y});
    
    uxx = matlabFunction(uxx,'Vars',{x,y});  
    uxy = matlabFunction(uxy,'Vars',{x,y});  
    uyy = matlabFunction(uyy,'Vars',{x,y});
end