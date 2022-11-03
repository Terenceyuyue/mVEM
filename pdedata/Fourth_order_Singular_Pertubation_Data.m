function pde = Fourth_order_Singular_Pertubation_Data(varargin)

% data for plate
if nargin == 0
    epsilon = 1;    
else 
    epsilon = varargin{1};
end
para = struct('epsilon',epsilon);    

[u,rhs,ux,uy,uxx,uxy,uyy] = compute_rhs(para);

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
function [u,f,ux,uy,uxx,uxy,uyy] = compute_rhs(para)    
    epsilon = para.epsilon;
    syms x y;    
    
    % exact solution
    u = (sin(pi*x)*sin(pi*y))^2;
    %u = (exp(cos(2*pi*x))-exp(1)).*(exp(cos(2*pi*y))-exp(1));

    % derivative
    ux = diff(u,x);      uy = diff(u,y);
    
     % second derivative
    uxx = diff(ux,x); uxy = diff(ux,y); uyy = diff(uy,y);
    
    % Laplace u
    Lap = uxx+uyy;

    % Laplace^2 u
    LapLap = diff(Lap,x,2) + diff(Lap,y,2);

    % f
    f = epsilon^2*LapLap - Lap; %f = simplify(f);
    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y}); 
    ux = matlabFunction(ux,'Vars',{x,y});  uy = matlabFunction(uy,'Vars',{x,y});
    
    uxx = matlabFunction(uxx,'Vars',{x,y});  
    uxy = matlabFunction(uxy,'Vars',{x,y});  
    uyy = matlabFunction(uyy,'Vars',{x,y});
end