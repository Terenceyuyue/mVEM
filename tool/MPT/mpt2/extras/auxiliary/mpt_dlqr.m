function [K, P, E] = mpt_dlqr(A, B, Q, R)
%MPT_DLQR Linear-quadratic regulator design for discrete-time systems.
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the optimal gain matrix K such that the state-feedback law
% u = -K*x minimizes the cost function
%
%           J = Sum {x'Qx + u'Ru}
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% A, B     - Matrices defining the dynamics; x^+=A*x+B*u
% Q, R     - Weighting matrices
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% K        - Optimal state-feedback
% P        - Common quadratic Lyapunov function V(x)=x'Px for feedback law K
% E        - E = EIG(A - B*K)
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

error(nargchk(4,4,nargin));

if exist('dlqr', 'file'),
    % use the Control Toolbox if possible
    [K, P, E] = dlqr(A, B, Q, R);
    
else
    % use mpt_getStabFeedback
    [F, P] = mpt_getStabFeedback(A, B, Q, R);
    K = -F{1};
    if nargout>2,
        E = eig(A - B*K);
    end
end