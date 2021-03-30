function pde = Strain_gradient_elasticity_data2(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;  iota = 1;
else
    lambda = para.lambda; mu = para.mu; iota = para.iota;
end

% --------- given by the symbolic computation ------
para = struct('lambda',lambda, 'mu',mu, 'iota',iota);
[u1,u2,f1,f2,u1x,u1y,u2x,u2y,u1xx,u1xy,u1yy,u2xx,u2xy,u2yy] = compute_rhs(para);

% ---------- functions in PDE ---------------
% rhs
    function val = f(p)  
        x = p(:,1); y = p(:,2);
        val = [f1(x,y), f2(x,y)];
    end
% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = u1(x,y);
        val(:,2) = u2(x,y);
    end
% Dirichlet
    function val = g_D(p)
        val = uexact(p);
    end
% derivative of the exact solution
    function val =  Du(p)  % w = [u1,u2]
        x = p(:,1); y = p(:,2);
        val(:,1) = u1x(x,y);  val(:,2) = u1y(x,y); 
        val(:,3) = u2x(x,y);  val(:,4) = u2y(x,y);
    end
% second order derivative of the exact solution
    function val =  DDu(p)  % w = [u1,u2]
        x = p(:,1); y = p(:,2);
        DDu1(:,1) = u1xx(x,y);  DDu1(:,2) = u1xy(x,y);  DDu1(:,3) = u1xy(x,y); DDu1(:,4) = u1yy(x,y); 
        DDu2(:,1) = u2xx(x,y);  DDu2(:,2) = u2xy(x,y);  DDu2(:,3) = u2xy(x,y); DDu2(:,4) = u2yy(x,y); 
        val = [DDu1, DDu2];
    end

pde = struct('para', para, 'f', @f, 'uexact',@uexact,'g_D',@g_D, 'Du',@Du, 'DDu',@DDu);
end


function [u1,u2,f1,f2,u1x,u1y,u2x,u2y,u1xx,u1xy,u1yy,u2xx,u2xy,u2yy] = compute_rhs(para)    
    lambda = para.lambda; mu = para.mu; iota = para.iota;
    syms x y;    
    
    % exact solution
    u1 = (exp(sin(pi*x))-1-pi*iota*(cosh(1/(2*iota))-cosh((2*x-1)/(2*iota)))/sinh(1/(2*iota)))...
        *(exp(sin(pi*y))-1-pi*iota*(cosh(1/(2*iota))-cosh((2*y-1)/(2*iota)))/sinh(1/(2*iota)));
    u2 = (sin(pi*x)-pi*iota*(cosh(1/(2*iota))-cosh((2*x-1)/(2*iota)))/sinh(1/(2*iota)))...
        *(sin(pi*y)-pi*iota*(cosh(1/(2*iota))-cosh((2*y-1)/(2*iota)))/sinh(1/(2*iota)));
    
    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);
    u2x = diff(u2,x);      u2y = diff(u2,y);

    % Laplace u
    Lap1 = diff(u1,x,2)+diff(u1,y,2);
    Lap2 = diff(u2,x,2)+diff(u2,y,2);

    % grad (grad .u)
    divu = u1x+u2y;
    graddiv1 = diff(divu,x);
    graddiv2 = diff(divu,y);

    % Linear elasticity
    Ls1 = mu*Lap1 + (lambda+mu)*graddiv1;
    Ls2 = mu*Lap2 + (lambda+mu)*graddiv2;

    % Laplace (Ls)
    LapLs1 = diff(Ls1,x,2) + diff(Ls1,y,2);
    LapLs2 = diff(Ls2,x,2) + diff(Ls2,y,2);

    % f
    f1 = iota^2*LapLs1 - Ls1; %f1 = simplify(f1);
    f2 = iota^2*LapLs2 - Ls2; %f2 = simplify(f2);
    
    % second order derivatives
    u1xx = diff(u1,x,2);  u1xy = diff(diff(u1,x),y);  u1yy = diff(u1,y,2);
    u2xx = diff(u2,x,2);  u2xy = diff(diff(u2,x),y);  u2yy = diff(u2,y,2);
    
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y});    u2 = matlabFunction(u2,'Vars',{x,y});
    f1 = matlabFunction(f1,'Vars',{x,y});    f2 = matlabFunction(f2,'Vars',{x,y}); 
    u1x = matlabFunction(u1x,'Vars',{x,y});  u1y = matlabFunction(u1y,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});  u2y = matlabFunction(u2y,'Vars',{x,y}); 
    
    u1xx = matlabFunction(u1xx,'Vars',{x,y});  u1xy = matlabFunction(u1xy,'Vars',{x,y});  u1yy = matlabFunction(u1yy,'Vars',{x,y});
    u2xx = matlabFunction(u2xx,'Vars',{x,y});  u2xy = matlabFunction(u2xy,'Vars',{x,y});  u2yy = matlabFunction(u2yy,'Vars',{x,y});
end