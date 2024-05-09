function pde = elasticitydata(para)
% Omega = [0,1] * [0,1]

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

% --------- given by the symbolic computation ------
[u1,u2,u1x,u1y,u2x,u2y,sig11,sig12,sig22,f1,f2] = compute_rhs(lambda,mu);

para = struct('lambda',lambda, 'mu',mu);
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [f1(x,y)+0*x,f2(x,y)+0*x];
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y)+0*x, u2(x,y)+0*x];
    end
    function val = g_D(p)
        val = uexact(p);
    end
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        Du1 = [u1x(x,y)+0*x, u1y(x,y)+0*x];
        Du2 = [u2x(x,y)+0*x, u2y(x,y)+0*x];
        val = [Du1, Du2];
    end
    function val = g_N(p)
        x = p(:,1); y = p(:,2);
        val = [sig11(x,y)+0*x,sig22(x,y)+0*x,sig12(x,y)+0*x];
    end
pde = struct('para',para, 'f', @f, 'uexact',@uexact,...
    'g_D',@g_D, 'g_N',@g_N, 'Du',@Du);
end

function [u1,u2,u1x,u1y,u2x,u2y,sig11,sig12,sig22,f1,f2] = compute_rhs(lambda,mu)
    syms x y;
    % exact solution
    %     u1 = cos(pi*x)*cos(pi*y);
    %     u2 = sin(pi*x)*sin(pi*y);
    u1 = pi*cos(pi*y)*sin(pi*x)^2*sin(pi*y);
    u2 = -pi*cos(pi*x)*sin(pi*x)*sin(pi*y)^2;

    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);
    u2x = diff(u2,x);      u2y = diff(u2,y);

    % modify the solution such that
    %    int_\DOmega u dx = 0,   int_\Omega rot u dx = 0
    %    DOmega = \parital\Omega
    % that is, we substract a function p\in RM
    % Let p = c0*[1;0] + c1*[0;1] + c2*[-y;x].
    % Then
    %     c2 = 0.5*int_\Omega rot u dx
    %     c0 = ( int_\DOmega u1 dx + c2*\int_\DOmega y dx ) / |DOmega|
    %     c1 = ( int_\DOmega u2 dx - c2*\int_\DOmega x dx ) / |DOmega|
    %     Omega = [0,1] * [0,1]
    rotu = u2x - u1y;
    DOmega = 4;
    Iu1 = double(int(subs(u1,y,0),x,0,1) + int(subs(u1,y,1),x,0,1) ...
        + int(subs(u1,x,0),y,0,1) + int(subs(u1,x,1),y,0,1));
    Iu2 = double(int(subs(u2,y,0),x,0,1) + int(subs(u2,y,1),x,0,1) ...
        + int(subs(u2,x,0),y,0,1) + int(subs(u2,x,1),y,0,1));
    Ix = 0.5+1+0.5+0;
    Iy = 0+0.5+1+0;
    c2 = 0.5*double(int(int(rotu,x,0,1),y,0,1));
    c0 = (Iu1+c2*Iy)/DOmega;
    c1 = (Iu2-c2*Ix)/DOmega;
    u1 = u1-c0+c2*y;
    u2 = u2-c1-c2*x;

    % Eij
    E11 = u1x;
    E12 = 0.5*(u1y+u2x);
    E22 = u2y;
    % sigma(u)
    sig11 = (lambda+2*mu)*E11 + lambda*E22;
    sig12 = 2*mu*E12; % sig21 = sig12
    sig22 = (lambda+2*mu)*E22 + lambda*E11;
    % f1 = [f11,f12]
    f1 = -(diff(sig11,x) + diff(sig12,y));
    f2 = -(diff(sig12,x) + diff(sig22,y));
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y});     u2 = matlabFunction(u2,'Vars',{x,y});
    u1x = matlabFunction(u1x,'Vars',{x,y});   u1y = matlabFunction(u1y,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});   u2y = matlabFunction(u2y,'Vars',{x,y});
    sig11 = matlabFunction(sig11,'Vars',{x,y});
    sig12 = matlabFunction(sig12,'Vars',{x,y});
    sig22 = matlabFunction(sig22,'Vars',{x,y});
    f1 = matlabFunction(f1,'Vars',{x,y});    f2 = matlabFunction(f2,'Vars',{x,y});
end