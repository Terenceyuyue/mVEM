function [F,P,feasible] = mpt_getStabFeedback(A,B,Q1,R,Options)
%MPT_GETSTABFEEDBACK Computes a stabilizing feedback law for a PWA system
%
% [F,P,feasible] = mpt_getStabFeedback(A,B,Q,R)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Function computes a stabilizing feedback law for a PWA system defined by the
% cells A{}, B{}. It is assumed that the dynamics are valid in the entire state-
% space. The obtained feedback law also minimizes the cost J=x'Q1x + u'Ru.
% Computation solves an SDP (Yalmip required). 
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% A{},B{}          - cell array defining the dynamics; x^+=A{i}x+B{i}u
% Q,R              - cost objective J=x'Qx + u'Ru.
% Options.verbose  - level of verbosity
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% F{}       - cell array containing the feedback law for each dynamic i
% P         - common quadratic Lyapunov function V(x)=x'Px for feedback laws F{}
% feasible  - 0/1 flag. Set to 1 if computation was successful
%
% see also MPT_GETCOMMONLYAPFCT, MPT_GETPWQLYAPFCT

% Copyright is with the following author(s):
%
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

error(nargchk(4,5,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<5,
    Options=[];
end

if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end

if(~iscell(A))
    tmp=A;
    clear A;
    A{1}=tmp;
end
if(~iscell(B))
    tmp=B;
    clear B;
    B{1}=tmp;
end

if(isempty(A{1}))
    error('A matrix is empty !')
end

nx=length(A{1});
nu=size(B{1},2);

myprog=lmi;
Q = sdpvar(nx,nx,'symmetric'); 
gam = 1;
for i=1:length(A)
    Y{i} = sdpvar(nu,nx,'full');  
    myprog=myprog+set([Q (A{i}*Q+B{i}*Y{i})' (Q1^0.5*Q)' (R^0.5*Y{i})';
        (A{i}*Q+B{i}*Y{i})  Q    zeros(nx,nx)  zeros(nx,nu);
        (Q1^0.5*Q)         zeros(nx,nx)        eye(nx)  zeros(nx,nu);
        (R^0.5*Y{i})       zeros(nu,nx)        zeros(nu,nx) eye(nu)]>0);
end
myprog = myprog + set(Q>0);  
if ~isempty(mptOptions.sdpsettings)
    options = mptOptions.sdpsettings;
else
    options=sdpsettings('Verbose',0);
end
solution = solvesdp(myprog,-trace(Q),options);
dQ=double(Q);
P=inv(dQ);
for i=1:length(A)
    F{i}=double(Y{i})*P;
end

if(solution.problem~=0)  
    [lP,feasible]=mpt_getCommonLyapFct(F,A,B,Options);
    if(~feasible)
        fprintf('\n\n');
        disp('mpt_getStabFeedback: NO SOLUTION FOUND')
        feasible=0;
    end
else
    feasible=1;
end

return