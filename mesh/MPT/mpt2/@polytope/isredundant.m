function [how] = isredundant(P,row,Options)
%ISREDUNDANT Check if a constraint is redundant
%
% [how] = isredundant(P,row,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Check if constraint defined by f*x<=g (f=P.H(row,:), g=P.K(row)) is redundant.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - Polytope
% row              - index of a row in P.H
% Options.abs_tol  - absolute tolerance
% Options.lpsolver - LP solver to be used (see help mpt_solveLP)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% how   = 1 (TRUE) redundant
%       = 0 (FALSE) non-redundant
%       = 2 empty polyhedron.
%
% see also REDUCE
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options=[];
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end

if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end

% if ~isa(P,'polytope'),
%     error('ISREDUNDANT: input argument MUST be a polytope object!');
% end

if length(P.Array)>0,
    error('ISREDUNDANT: this function does not support polytope arrays!');
end

[m,n]=size(P.H);

if m ~= size(P.K,1) | row > m | row < 1,
   error('ISREDUNDANT: Bad sizes in ISREDUNDANT.'); 
end

%LPi%f  = P.H(row,:)';
f  = P.H(row,:);

ii = [1:row-1, row+1:m];

if ~isempty(ii),
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(-f,P.H(ii,:),P.K(ii),[],[],[],Options.lpsolver);
   %LPi%obj=f'*xopt-P.K(row);
   obj=f*xopt-P.K(row);
else
   status='unbounded';
end

if strcmp(status,'unbounded') | (obj>Options.abs_tol & strcmp(status,'ok')),  
   how=0; % non-redundant
elseif(strcmp(status,'infeasible'))
   %double check to make sure that the problem is really infeasible
   [x,R]=chebyball(polytope(P.H(ii,:),P.K(ii),0,2),Options.lpsolver);
   if(max(P.H(ii,:)*x-P.K(ii))>0)
     how=2; % infeasible
   else
     how=1;
   end   
else 
   how=1; % redundant
end
