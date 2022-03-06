function x = cvx2cvx(A,b,w,mu)
% This function finds the minimizer of the sum of two convex functions:
%
%        min   f(x) + g(x),
%              f(x) = 0.5*x'*A*x - b'*x 
%              g(x) = w'*abs(x)
% 
% using the fixed point algorithm based on proximity operator.
% The proximal algorithm is:
%
%        x_{n+1} = prox_{lambda*g} o (I - lambda*grad(f)) (x_n)
%        grad(f)(x) = Ax - b
%        prox_{lambda*g}(x) = [prox_1, prox_2, ..., prox_n]'
%        prox_i = sgn(x_i)*max( |x_i| - w_i*lambda, 0 )
%
%   References
%   H. K. Xu, "Properties and iterative methods for  the Lasso and Its
%   Variants", Chin. Ann. Math., Vol 35B. No 3., pp. 501¨C518, 2014.
%   (See Eq. (3.4) there)
%
% Copyright (C)  Terence Yu.

if nargin==3, mu = 1e-2; end

x0 = zeros(size(b));
err = 1; tol = 1e-8;
while err > tol
    r = A*x0-b;    % grad f(x) 
    x = x0-mu*r; % (I - lambda*grad(f))(x)
    xw = abs(x)-mu*w;
    x = sign(x).*xw.*(xw>0); % prox_{lambda*g}(x)
    err = norm(x-x0);
    if err>1e5 || sum(isnan(x))  % restart the loop for nonconvergence
        x = zeros(size(b));
        mu = 1/max(eigs(A));
    end
    x0 = x; 
end