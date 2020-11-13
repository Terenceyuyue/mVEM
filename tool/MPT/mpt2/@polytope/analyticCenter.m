function x = analyticCenter(P,Options)
%ANALYTICCENTER Computes an analytic center of a polytope
%
% x = center(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes an analytic center of a polytope. The input polytope must be bounded
% and fully dimensional! Newton iteration based algorithm is used.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - Polytope
% Options.x0       - a point inside relative interior of the polytope P
%                    (default: Chebyshev center)
% Options.tol      - tolerance used as a stopping criterion, corresponds
%                    to the step of the Newton iteration
%                    (default: Chebyshev radius/10)   
% Options.maxIter  - maximal number of Newton step iterations 
%                    (default 10)    
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% x        - Analytic center of the polytope P
%
% see also CHEBYBALL
%

% Copyright is with the following author(s):
%
% (C) 2005 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          loefberg@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------


if ~isa(P,'polytope')
    error('ANALYTICCENTER: Input argument must be a polytope!');
end

if ~isempty(P.Array)
    error('ANALYTICCENTER: This function does not work with polytope arrays!');
end

if ~isfulldim(P)
    error('ANALYTICCENTER: Polytope must be fully dimensional!');
end

if ~isbounded(P)
    error('ANALYTICCENTER: Polytope must be bounded!');
end

if ~isminrep(P)
    P = reduce(P);
end

if ~isnormal(P)
    P = normalize(P);
end

xCheby = P.xCheb;
rCheby = P.RCheb;


if ( nargin == 1 )
    % using default settings
    nIter = 10;
    EPS = rCheby / 10;
    x0 = xCheby;
else
    if ~isfield(Options,'x0')        
        x0 = xCheby; % use Chebyshev center by default
    else
        x0 = Options.x0;
        if ~isinside(P,x0)
            warning(['ANALYTICCENTER: specified initial point is not inside the ' ...
                     'polytope, using Chebyshev center.']);
            x0 = xCheby;
        end
    end
    
    if ~isfield(Options,'tol')
        EPS = rCheby/10;
    else
        EPS = Options.tol;
    end
    
    if ~isfield(Options,'nIter')
        nIter = 10;
    else
        nIter = Options.maxIter;
    end
end

ALPHA = 0.2;
BETA  = 0.717;

[A,b]   = double(P); 
x       =  x0;
nablaf  = -A'*(1./(b-A*x));
Hessian =  A'*diag((b-A*x).^(-2))*A;
fk      = -sum(log(b-A*x));

while ( (nablaf'*(Hessian\nablaf) > 2*EPS) && (nIter > 0) )
    dx = Hessian\nablaf;
       
    % initial step length (maximum step which keeps us within the bounds of
    % the polytope)
    %
    inProd = A * dx;
    idxGTZero = find(inProd > 0);
    inProd = inProd(idxGTZero); % take only > 0
    slacks = b(idxGTZero) - A(idxGTZero,:) * x;
    t0 = min ( slacks * norm(dx) ./ inProd );
    t0 = min(0.99*t0,1);
    
    % do the backtracking linesearch
    %
    t = t0;
    while ( 1 )
        fkpp   = -sum(log(b-A*(x+t*dx)));
        if ( (fk - fkpp + ALPHA*t*nablaf'*dx) > 0 )
            break;
        end
        t = BETA * t;
    end
    x = x + t*dx;

    fk     = -sum(log(b-A*x));
    s      = 1./(b-A*x);
    nablaf = -A'*s;
    SA = diag(s)*A;    
    Hessian = SA'*SA;    
    nIter = nIter - 1;
end