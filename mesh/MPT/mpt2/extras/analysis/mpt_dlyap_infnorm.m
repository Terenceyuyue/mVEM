function [W,fail,Q,CPU] = mpt_dlyap_infnorm(A,Options)
% Lyapunov equation with inf-norm for discrete-time LTI systems
%
%  [W,fail,Q,cpu] = mpt_dlyap_infnorm(A,Options)
%
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% A                system matrix x(k+1) = A * x(k), with  |eig(A)| < 1
% Options.         optional arguments
%  .abs_tol        absolute tolerance (default mptOptions.abs_tol)
%  .lpsolver       equivalent to the MPT settings (default mptOptions.lpsolver)
%  .complexW       if 1, a complex W is allowed. (default 0)
%
% ---------------------------------------------------------------------------
% OUTPUT
% --------------------------------------------------------------------------- 
% O:   W       final penalty marix of J(x)= ||W x||_inf
%      fail    if 0, found a solution
%      Q       find W and Q such that, W*A - Q*W = 0 and ||Q||_inf <1
%      cpu     CPU time
%
%
% see also  MPT_NORM2PWA   MPT_LYAPUNOV


% Reference:
%   Polanski: On infinity norms as Lyapunov functions for linear systems.
%     IEEE Trans. Aut.Contr., Vol. 40(7), Jul. 1995
%   Christophersen, Morari: Further Results on 'Infinity Norms as Lyapunov
%     Functions for Linear Systems', Tech. Report, Automatic Control Lab, ETH
%     Zurich, 2005

% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch

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

% 2006-03-14, fjc

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

error(nargchk(1,2,nargin));

if (nargin<2),
    Options = [];
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end
abs_tol = Options.abs_tol;

if ~isfield(Options,'lpsolver')
    Options.lpsolver = mptOptions.lpsolver;
end
if ~isfield(Options,'method_complexpair')
    Options.method_complexpair = 'algebraic';
    % alternative: fjc     (simple/effective LP based iteration)
    %              yalmip  (bad for eigenvalues close to stabolity boundary!!)
end
if ~isfield(Options,'complexW')
    Options.complexW = 0;
end
if ~isfield(Options,'minimal')
    %  .minimal        if 1 (default), the minimal size real W is computed (iteration),
    Options.minimal = 1;
end


W = [];
Q = [];
fail  = 1;


n     = size(A,1);
[V,lamA] = eig(A);
eigA  = diag(lamA).';

t0 = cputime;

if any(abs(eigA) >= 1 )
    disp('Discrete-time system is unstable, i.e. V(x) does not exist.')
    fail = 1;
    return
end


% resort the eigenvals/eigenvects
if ~isreal(eigA)
    indR = [];
    indC = [];
    for ii =1:n
        if isreal(eigA(ii))
            indR = [indR ii];
        else
            indC = [indC ii];
        end
    end%FOR ii
    nR = length(indR);
    nC = length(indC);
    V = [V(:,indR) V(:,indC)];
    eigA = [eigA(indR) eigA(indC)];
    lamA = diag(eigA);
else
    nR = n;
    nC = 0;
end


% detect possible multible eigenvalues and if they have the same eigenvector
% proper jordan block decomposition
flag = 0;
for ii = 1:n
    rest = setdiff(1:n,ii);
    ind = find(abs(eigA(ii)-eigA(rest)) <=abs_tol);
        
    %check if the egenvectors are independent
    if ~isempty(ind)
        for ll = 1:length(ind)
            
            % scale the other (possibly dependent) eigenvectors
            iiIdx = find(abs(V(:,ii))>abs_tol,1);
            V_ = V(:,rest(ind(ll)))/V(iiIdx,rest(ind(ll)))*V(iiIdx,ii);
            
            if norm(V(:,ii)-V_,inf) <= abs_tol
                
                if ~isreal(eigA(ii))
                    error('Multiple complex eigenvalues with non-independent eigenvectors is not supported.')
                    fail = 1;
                    return
                end
                
                flag = 1;
                break
            end
        end
    end
    if flag
        break
    end
end%FOR ii
if flag
    [V,lamA] = jordan(A);
    % unfortunately there is no real better solution;
end


if ~Options.complexW & ~isreal(eigA)
    [V,lamA] = cdf2rdf(V,lamA);
end



Q = [];
U = [];

% real eigenvalues
if nR>0
    [W_real, Q_real] = realEV(lamA(1:nR,1:nR),Options);
    
    if n>200
        U = sparse(blkdiag(U,W_real));
        Q = sparse(blkdiag(Q,Q_real));
    else
        U = blkdiag(U,W_real);
        Q = blkdiag(Q,Q_real);
    end
end% IF nR


%complex eigenvalues
if nC>0
    for ii = nR+1:2:nR+nC-1
        if Options.complexW
            W_complex = inv(lamA(ii:ii+1,ii:ii+1));
            Q_complex = lamA(ii:ii+1,ii:ii+1);
        else
            [W_complex, Q_complex] = complexpair(lamA(ii:ii+1,ii:ii+1),Options);
        end
        
        U = blkdiag(U,W_complex);
        Q = blkdiag(Q,Q_complex);
    end%FOR ii
    
    
    if n>200
        U = sparse(U);
        Q = sparse(Q);
    end
end% IF nC

    
W = U/V;
    

if ~Options.complexW
    if norm(imag(W),inf) < abs_tol
        W = real(W);
    else
        warning('Numerical errors occured. W should be real but has a significant complex part.')
    end
    
    if norm(imag(Q),inf) < abs_tol
        Q = real(Q);
    else
        warning('Numerical errors occured. Q should be real but has a significant complex part.')
    end
end



% final test
if  (norm(W*A - Q*W,inf) < abs_tol)  & ( norm(Q,inf) < 1 + abs_tol ) &  ( rank(W)==n )
    fail = 0;
else
    fail = 1;
end

if n>200
    W = full(W);
    Q = full(Q);
end


CPU = cputime -t0;

return





% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SUBROUTINES
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [W_real, Q_real] = realEV(A_real,Options)

n = size(A_real,1);
abs_tol = Options.abs_tol;

% decomposition was done before
U = eye(n);
lamA = A_real;
lam = diag(lamA);

if n>1
    
    if 1  % scale ||Q||_inf < 1
    
        p = diag(lamA,1);
        Idxp = find(abs(p)<=abs_tol);
        Idxp = [Idxp(:)' n];
        
        
        if length(Idxp)==n
            % no multiple eigenvalues
            d = ones(n,1);
            
        else
            % multiple eigenvalues: scaling of Q
            
            H = [diag(abs(p)), zeros(n-1,1)] + [zeros(n-1,1), diag(abs(lam(1:end-1)) -1 + abs_tol)];
            H = [H;-eye(n)];
            K = zeros(n-1,1);
            K = [K; -abs_tol*ones(n,1)];
            
            EE = eye(n);
            Heq = EE(Idxp,:);
            Keq = ones(length(Idxp),1);
            
            [d,fval,lambda,exitflag,how] = mpt_solveLP(-ones(n,1),H,K,Heq,Keq,[],Options.lpsolver);
            
            % handle numerical problems, d>0
            idx_d = find(abs(d)<abs_tol);
            if ~isempty(idx_d)
                d(idx_d) = abs_tol*ones(length(idx_d),1);
            end
        end

    else  % minimize ||Q||_inf  and ||Q||_inf < 1
        
        % not to be used!!!: produces very ill-conditioned W_real, Q_real!!!
        
        r = diag(lamA,1);

        f = [1; zeros(n-1,1)];
        H = [f'; -f'; -f'; [zeros(n-1,1) -eye(n-1)]; [-ones(n-1,1) diag(r)]];
        K = [1-abs_tol, -abs_tol, -abs(lam(n)), -abs_tol*ones(1,n-1), -abs(lam(1:n-1))']';
        
        [z,q,lambda,exitflag,how] = mpt_solveLP(f,H,K,[],[],[],Options.lpsolver);
        
        kk = z(2:end);
        
        if strcmp(how,'ok')
            d(n) = 1e14;
            for nn = n-1:-1:1
                d(nn) = kk(nn)*d(nn+1);
            end
        else
           disp('Numerical problems.') 
        end
        
        
    end
    
    
elseif n==1
   % n==1
   d = 1;
   
else
   error('some problems in detecting the dimension')
end
        
W_real = diag(d)*U';
Q_real = W_real*A_real/W_real;
    
if norm(imag(W_real),inf) < abs_tol
    W_real = real(W_real);
else
    warning('Numerical errors occured. W should be real but has a significant complex part.')
end

if norm(imag(Q_real),inf) < abs_tol
    Q_real = real(Q_real);
else
    warning('Numerical errors occured. Q should be real but has a significant complex part.')
end

return


%++++++++++++++++++++++++++++

function [W_ab, Q_ab] = complexpair(A_ab,Options)
    


%find m_ab
% to solve: W_ab*A_ab = Q_ab*W_ab, ||Q_ab||_inf < 1

r = norm([A_ab(1,1) A_ab(1,2)]);
if A_ab(1,1) ~= 0
    alpha = atan(abs(A_ab(1,2)/A_ab(1,1)));
else
    alpha = pi/2;
end


m = 2;
m_max = ceil(pi/(acos(r)*2));

if ~Options.minimal % fastest way to compute W
    m = m_max;
end

    
if m_max<2
    error('numerical problems.')
end


if m_max==2
    W_ab = eye(2);
    Q_ab = A_ab;
    
else % iteration over m
    
    if strcmp(Options.method_complexpair, 'algebraic')    % method by FJC
        
        ok_flag = 0;
        while m <= m_max
                        
            k_low = floor(alpha*m/pi);
            % k_up  = k_low +1;
            
            % phi_low = k_low*pi/m;
            % phi_up  = k_up*pi/m;
            
            % rt = [cos(alpha) -cos(phi_low)+cos(phi_up); ...
            %      sin(alpha) -sin(phi_low)+sin(phi_up)]\[cos(phi_up);sin(phi_up)];
            
            rr = cos(pi/(2*m)) / cos(alpha - pi/(2*m)*(2*k_low+1));
            
            % if rt(1)<r | rt(2)<0 | rt(2)>1
            
            if rr < r
                m = m+1;
            else
                ok_flag = 1;
                break
            end
        end
        
        ww   = [0:m-1]*pi/m;
        W_ab = [cos(ww') sin(ww')];
        
        Aeq = [[W_ab' -W_ab']; ones(1,2*m)];
        Beq = [A_ab(1,1) A_ab(1,2) 1]';
                
        AA = [-eye(2*m)];
        BB = [zeros(2*m,1)-Options.abs_tol/(2*m)];
                
        [t,fval,lambda,exitflag,how]=mpt_solveLP(zeros(2*m,1),AA,BB,...
            Aeq,Beq,[],Options.lpsolver);
        
    else % solving via LP/optimization
        
        if strcmp(Options.method_complexpair, 'yalmip')
            sdp_opt = sdpsettings('verbose',0,'solver','nag');
            % sdp_opt = sdpsettings('verbose',1);
        end
        
        ok_flag = 0;
        while m <= m_max
            
            ww   = [0:m-1]*pi/m;
            W_ab = [cos(ww') sin(ww')];
            
            if strcmp(Options.method_complexpair, 'yalmip')
                % bad for eigenvalues close to stabolity boundary!!
                
                x = sdpvar(m,1);
                
                % construct the structure of Q_ab explicitly
                Q_ab = zeros(m);
                for kk = 1:m
                    Q_ab = Q_ab + diag(x(kk)*ones(1,m-kk+1),kk-1);
                end%kk
                
                for kk = 2:m
                    Q_ab = Q_ab - diag(x(kk)*ones(1,kk-1),-m+kk-1);
                end%kk
                
                % alternative
                %  Q_ab = sdpvar(m,m,'full');
                
                %            F = set(W_ab*A_ab==Q_ab*W_ab) + set(norm(Q_ab,inf)<1 - Options.abs_tol);
                F = set(W_ab*A_ab==Q_ab*W_ab) + set(norm(x,1)<1 - Options.abs_tol);
                sdp_info = solvesdp(F,[],sdp_opt);
                
                if sdp_info.problem==0
                    ok_flag = 1;
                    break
                else
                    m = m +1;
                end
                
            else
                % method by FJC
                
                Aeq = [[W_ab' -W_ab']; ones(1,2*m)];
                Beq = [A_ab(1,1) A_ab(1,2) 1]';
                
                AA = [-eye(2*m)];
                BB = [zeros(2*m,1)-Options.abs_tol/(2*m)];
                
                [t,fval,lambda,exitflag,how]=mpt_solveLP(zeros(2*m,1),AA,BB,...
                    Aeq,Beq,[],Options.lpsolver);
                
                    if strcmp(how,'ok')
                        ok_flag = 1;
                        break
                    else
                        m = m+1;
                    end
            end%IF: method
        end% WHILE: m_ab
    end%IF: method
    
    if m == m_max & ~ok_flag
        warning(['maximum nubmer of possible iterations reached without' ...
            ' getting a proper solution'])
    end
    
    
    if strcmp(Options.method_complexpair, 'yalmip')
        Q_ab = double(Q_ab);
        
    else % LP based iteration or algebraic
        x = [eye(m) -eye(m)]*t;
        
        Q_ab = zeros(m);
        for kk = 1:m
            Q_ab = Q_ab + diag(x(kk)*ones(1,m-kk+1),kk-1);
        end%kk
        
        for kk = 2:m
            Q_ab = Q_ab - diag(x(kk)*ones(1,kk-1),-m+kk-1);
        end%kk
        
    end%IF: method
    

end% IF m_max==2


if norm(imag(W_ab),inf) < Options.abs_tol
    W_ab = real(W_ab);
else
    warning('Numerical errors occured. W should be real but has a significant complex part.')
end

if norm(imag(Q_ab),inf) < Options.abs_tol
    Q_ab = real(Q_ab);
else
    warning('Numerical errors occured. Q should be real but has a significant complex part.')
end
    
return


%+++++++++++++++++++++++++++++++++++++++++++++++++

function D = scalingD(A,Options)

n = size(A,1);
abs_tol = Options.abs_tol;

% decomposition was done before
lamA = A;
lam = diag(lamA);

if n==1
    d = 1;
    
else%%IF n > 1
    p = diag(lamA,1);
    Idxp = find(abs(p)<=abs_tol);
    Idxp = [Idxp(:)' n];

    if length(Idxp)==n
        % no multiple eigenvalues
        d = ones(n,1);
        
    else
        % multiple eigenvalues: scaling of Q

    end%IF: length(Idxp)==n 
    
end%%IF n


if n>1
   
    
    
        
        H = [diag(abs(p)), zeros(n-1,1)] + [zeros(n-1,1), diag(abs(lam(1:end-1)) -1 + abs_tol)];
        H = [H;-eye(n)];
        K = zeros(n-1,1);
        K = [K; -abs_tol*ones(n,1)];
        
        EE = eye(n);
        Heq = EE(Idxp,:);
        Keq = ones(length(Idxp),1);
        
        [d,fval,lambda,exitflag,how] = mpt_solveLP(-ones(n,1),H,K,Heq,Keq,[],Options.lpsolver);
        
        % handle numerical problems, d>0
        idx_d = find(abs(d)<abs_tol);
        if ~isempty(idx_d)
            d(idx_d) = abs_tol*ones(length(idx_d),1);
        end
    end
    


D = diag(d);
return


