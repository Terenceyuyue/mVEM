function [xc,R,lambda] = chebyball(P,Options,force)
%CHEBYBALL Computes center and radius of the largest ball inscribed in a polytope
%
% [xc,R,lambda] = chebyball(P,Options,force)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns center and radius of a largest ball that can be inscribed in a polytope P
%
% USAGE:
%   [xc, Rc] = chebyball(P)
%   [xc, Rc] = chebyball(P, Options)
%   [xc, Rc, lambda] = chebyball(P)
%
% NOTE: usually the 'POLYTOPE' constructor computes these values automatically.
%       If they are already stored in the internal structure, no computation will
%       be performed unless force=1
%
% For nice plots, call:
%   mpt_plotArrangement(P);
%   hold on
%   [xc,R] = chebyball(P);
%   plotellip(eye(2)/R^2,xc);
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                 - Polytope
% Options.lpsolver  - LP solver to use (see help mpt_solveLP)
% force             - force recomputing of the Chebyshev's ball
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% xc                - Center of Chebyshev's ball
% R                 - Radius of Chebyshev's ball
% lambda            - Set of lagrangian multipliers
%
% see also POLYTOPE

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

% if ~isa(P, 'polytope')
%   error('First argument MUST be a polytope object');
% end

global mptOptions;

if nargin<2,
    Options=[];
end
if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end

if nargin<3
    force=0;
end

lenP=length(P.Array);
if lenP>0,
    % if P is a polyarray, we recursively solve the problem by calling chebyball for each element of P
    R=zeros(lenP,1);                       % initialize the output
    xc=zeros(lenP,size(P.Array{1}.H,2));   % initialize the output
    lambda=[];                             % initialize the output
    for ii=1:lenP,
        [x,R(ii),l]=chebyball(P.Array{ii},Options,force);  % solve the problem for each element of P
        xc(ii,:)=x';                                      
        lambda=[lambda; l];
    end
    return
end
    

if ~isempty(P.RCheb) & force==0,
    % if the chebyshev's data are already stored in the structure, we return them without recomputing
    % if force=1, we skip this step
    if P.RCheb>-Inf
        xc=P.xCheb;
        R=P.RCheb;
        lambda=[];
        return
    end
end

H = P.H;
K = P.K;
[nc,nx]=size(H);

if nc==0 | nx==0
    error('CHEBYBALL: empty polytope!');
end

% If any boundary is -Inf polytope P is empty
%--------------------------------------------
if any(K==-Inf)
    %[nc,nx]=size(P.H);
    xc=zeros(nx,1);
    R=-Inf;
    lambda=zeros(nc,1);
    return;
end

% Remove rows with Inf boundaries
%--------------------------------
ii=(K==Inf);
H(ii,:)=[];
K(ii)=[];

[nc,nx]=size(H);

if nc==0,
    xc=zeros(nx,1);
    R=Inf;
    lambda=zeros(nc,1);
    return;
end

Aconstr = [H, -sqrt(sum(H.*H,2))];

% use 'rescue' function - resolve an LP automatically if it is infeasible
% with default solver
[xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],Aconstr,K,[],[],[],Options.lpsolver);

if ~strcmp(how,'ok'),
    % maybe there is a numerical problem, thus we normalize H and K
    
    if ~P.normal,
        P = normalize(P);
        H = P.H;
        K = P.K;
        Aconstr = [H, -sqrt(sum(H.*H,2))];
    end
    x0 = [zeros(nx,1); 1000];         % hard-coded initial conditions
    
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],Aconstr,K,[],[],x0,Options.lpsolver);
end

nag3 = (nargout>3);
if nag3,
    lambda = Aconstr*xopt - K;  % active constraints are given by  lambda<=eps
end

xc=xopt(1:nx); % Center of the ball
R=-xopt(nx+1); % Radius of the ball
