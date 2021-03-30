function pde = Darcydata1
% u = K*grad(p)
% Boundary condition: u*n = g

K = [2 1; 1 2];

% --------- given by the symbolic computation ------
[pe,Kpx,Kpy,rhs] = compute_rhs(K);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [Kpx(x,y)+0*x,Kpy(x,y)+0*y];  % u = K*grad(p)
    end
    function val = pexact(p)
        x = p(:,1); y = p(:,2);
        val = pe(x,y)+0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y)+0*x;
    end

pde = struct('uexact',@uexact, 'pexact',@pexact, 'f',@f, 'K',K);
end

function [pe,Kpx,Kpy,f] = compute_rhs(K)
    syms x y;
    % exact solution
    pe = (x+1)*(y+1)-1/(x+y+3)-1/3*(7/4-16*log(2)+9*log(3));
    % derivative
    px = diff(pe,x);      py = diff(pe,y);
    % K*grad(p)
    Kpx = K(1,1)*px + K(1,2)*py;
    Kpy = K(2,1)*px + K(2,2)*py;
    % f = -div(K*grad(p))
    f = -(diff(Kpx,x) + diff(Kpy,y));
    % convert to anonymous functions
    pe = matlabFunction(pe,'Vars',{x,y});    
    f = matlabFunction(f,'Vars',{x,y});    
    Kpx = matlabFunction(Kpx,'Vars',{x,y});  
    Kpy = matlabFunction(Kpy,'Vars',{x,y}); 
end
