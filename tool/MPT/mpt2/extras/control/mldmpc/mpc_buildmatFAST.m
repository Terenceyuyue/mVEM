function [S1, S2, S3, F1, F2, F3, c1, c2, c3, IntIndex, Ext] = mpc_buildmatFAST(horizon, ...
              SYSTEM, WEIGHT, x1, u1, d1, z1, y1, eps2, xtt, Options) 

%===============================================================================
%
% Title:       mpc_buildmatFAST.m                                              
%
% Version:     2.0
%                                                                       
% Project:     Control of MLD systems
%                                                                       
% Sub-Project: Constrained Finite Time Optimal Control (CFTOC) of MLD systems
%                                                                       
% Purpose:     Formulate multi-parametric (mixed integer) program
%              equivalent to the CFTOC problem of MLD system
%                                                                        
% Author:      (C) Mato Baotic, Zurich, March 17, 2003
%              (C) Michal Kvasnica, Zurich, June 23, 2005
%
% History:     date        subject
%              2003.11.19  first public release                 
%              2003.07.29  fast implementation (static allocation of memory space)
%                          Tobias Geyer = Version 1.2
%              2003.03.18  Various bugs fixed
%                          = Version 1.1
%              2003.03.17  Initial Version
%                          = Version 1.0
%                      
% Description: Given:   horizon
%                       SYSTEM,
%                       WEIGHT, 
%                       x1, u1, d1, z1, y1,
%                       eps2, xtt, Options
%              Returns: S1, S2, S3, F1, F2, F3, c1, c2, c3, IntIndex
%
%              The optimization for the control problem is formulated as:
%                         
%                   min  v' S1 v + 2(S2 + x0' S3) v
%                   s.t. F1 v   <= F2  + F3 x0
%              Here we denote:
%
%                   horizon - vector of time segments (steps) through which
%                             specific SYSTEM is used in state update equation.
%                   NT      - prediction horizon [1, inf), NT=sum(horizon)
%                   SYSTEM  - cell structure with MLD system descriptions:
%                             A, B1, B2, B3, B5, C, D1, D2, D3, D5, E1, E2, E3, E4, E5
%                             Note: if input is not cell structure (but just a structure)
%                                   then we have a time invariant system
%                   WEIGHT  - cell structure with weights on u, d, z, x, y:
%                             Qu, Qd, Qz, Qx, Qy
%                   IntIndex - indices of integer variables in optimizer vector v
%
%                   v     = [ U',D',Z',EU',ED',EZ',EX',EY']'
%                   U     = [ u'(0), ... , u'(T-1)]'
%                   D     = [ d'(0), ... , d'(T-1)]'
%                   Z     = [ z'(0), ... , z'(T-1)]'
%                   EU    = [ epsu'(0), ... , epsu'(T-1)]'
%                   ED    = [ epsd'(0), ... , epsd'(T-1)]'
%                   EZ    = [ epsz'(0), ... , epsz'(T-1)]'
%                   EX    = [ epsx'(0), ... , epsx'(T-1)]'
%                   EY    = [ epsy'(0), ... , epsy'(T-1)]
%
%              The steady state values are: x1, d1, u1, z1, y1. These variables
%              can also be time-varying trajectories. 
%              If the steady state values have not the length n_*NT (n_ is the
%              corresponding number of components), then the last value of the
%              steady state vector is repeated.
%
%              The terminal state is given by xtt
%
%              c1,c2,c3: parameters for the constant term of the optimization
%                        the complete constant term is obtained as: 
%                        Jc = 0.5*x(0)'*c1*x(0) + c2*x(0) + c3
%
%              It is assumed, that the prediction horizon and the control
%              horizon are the same
%
% Contact:     Mato Baotic
%              Automatic Control Laboratory                
%              ETH Zentrum,
%              Zurich, Switzerland
%
%              baotic@control.ee.ethz.ch
%
%              Comments and bug reports are highly appreciated                 
%

%===============================================================================
%
% Legal note:   This program is free software; you can redistribute it and/or
%               modify it under the terms of the GNU General Public License as 
%               published by the Free Software Foundation; either version 2.1 of
%               the License, or (at your option) any later version. 
%
%               This library is distributed in the hope that it will be useful,
%               but WITHOUT ANY WARRANTY; without even the implied warranty of
%               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%               Lesser General Public License for more details.
% 
%               You should have received a copy of the GNU Lesser General Public
%               License along with this library; if not, write to the 
%               Free Software Foundation, Inc., 
%               59 Temple Place, Suite 330, 
%               Boston, MA  02111-1307  USA
%
%===============================================================================

if nargin==0
    disp('Version 1.1');
    return;
end

if nargin<11
    Options=[];
end
if ~isfield(Options,'norm')
    Options.norm=inf;
end
if ~isfield(Options,'TerminalConstraint')
    Options.TerminalConstraint=1;
end
if ~isfield(Options, 'reduceMLD')
    % if set to true, removes redundant constraints defining AD constraints from
    % the MLD representation. 
    % NOTE! this only works if AD items are in the form "d = f(x) >= 0" !!!
    Options.reduceMLD = 0;
end
if ~isfield(Options, 'useSparse')
    % if set to true, uses sparse matrices for constraints
    % NOTE! some solvers might not support the sparse format!
    Options.useSparse = 1;
end
if ~isfield(Options, 'reduceSlacks')
    % if set to true, does not introduce slacks for variables which are not
    % penalized
    Options.reduceSlacks = 1;
end
if ~isfield(Options, 'allXconstrained'),
    % if set to true, state constraints will be added on every state x(k). if
    % set to false, state constraints are imposed only on the final state x(N).
    Options.allXconstrained = 1;
end

% NT = sum of all horizons (if there are more than 1)
NT=sum(horizon);
if NT < 1
    error('Prediction horizon must be greater 0!');
end
if NT ~= round(NT),
    error('Prediction horizon must be an integer number!');
end

if Options.reduceMLD,
    % remove redundant rows. see help of the subfunction for more details.
    for ii = 1:length(SYSTEM),
        SYSTEM{ii} = sub_remove_redundant_ineq(SYSTEM{ii});
    end
end

% Dimensions
nu = size(SYSTEM{1}.B1,2);   % nu=dimension of u
nd = size(SYSTEM{1}.B2,2);   % nd=dimension of \delta
nz = size(SYSTEM{1}.B3,2);   % nz=dimension of z
nx = size(SYSTEM{1}.A,2);    % nx=dimension of x
ne = size(SYSTEM{1}.E5,1);   % ne=dimension of E_5
ny = size(SYSTEM{1}.C,1);    % ny=number of outputs

nub = SYSTEM{1}.nub;
nur = SYSTEM{1}.nur;
nxb = SYSTEM{1}.nxb;
nxr = SYSTEM{1}.nxr;
nyb = SYSTEM{1}.nyb;
nyr = SYSTEM{1}.nyr;

if nxb>0 | nyb>0
    %error('Tobias we have to think about this!!!!');
    % we do not need to take special care about the binary states and inputs. 
    % They do not show up in the optimizer.
end

if length(SYSTEM)~=length(horizon)
    error('SYSTEM should be a cell structure of length equal to the number of prediction horizons in horizon');
end
if ~(length(WEIGHT)==length(horizon) | length(WEIGHT)==1)
    error('WEIGHT should be a cell structure of length equal to 1 or the number of prediction horizons in horizon');
end

% relation time-to-system
t2s=zeros(1,NT);
beg=0;
for ii=1:length(horizon)
    t2s(beg+1:beg+horizon(ii))=ii;
    beg=beg+horizon(ii);
end

% relation time-to-weight
if length(WEIGHT)==1
    t2w=ones(1,NT);
else
    t2w=t2s;
end



% Define number of slack variables
if Options.norm==inf
    neu = min(1,nu);
    ned = min(1,nd);
    nez = min(1,nz);
    nex = min(1,nx);
    ney = min(1,ny);
elseif Options.norm==1
    neu = nu; 
    ned = nd;
    nez = nz;
    nex = nx;
    ney = ny;
elseif Options.norm==2
    neu = 0;
    ned = 0;
    nez = 0;
    nex = 0;
    ney = 0;
else
    error('Unknown norm has been specified!');
end

if Options.reduceSlacks,
    if Options.norm~=2,
        % set number of slacks to zero if all corresponding weights are zero
        allzero = 1;
        for ii = 1:length(WEIGHT),
            if ~isempty(WEIGHT{ii}.Qd),
                if any(any(WEIGHT{ii}.Qd~=0)),
                    allzero = 0;
                    break
                end
            end
        end
        if allzero,
            ned = 0;
        end
        allQdzero = allzero;
        
        allzero = 1;
        for ii = 1:length(WEIGHT),
            if ~isempty(WEIGHT{ii}.Qx),
                if any(any(WEIGHT{ii}.Qx~=0)),
                    allzero = 0;
                    break
                end
            end
        end
        if allzero,
            nex = 0;
        end
        allQxzero = allzero;
        
        allzero = 1;
        for ii = 1:length(WEIGHT),
            if ~isempty(WEIGHT{ii}.Qu),
                if any(any(WEIGHT{ii}.Qu~=0)),
                    allzero = 0;
                    break
                end
            end
        end
        if allzero,
            neu = 0;
        end
        allQuzero = allzero;
        
        allzero = 1;
        for ii = 1:length(WEIGHT),
            if ~isempty(WEIGHT{ii}.Qz),
                if any(any(WEIGHT{ii}.Qz~=0)),
                    allzero = 0;
                    break
                end
            end
        end
        if allzero,
            nez = 0;
        end
        allQzzero = allzero;
        
        allzero = 1;
        for ii = 1:length(WEIGHT),
            if ~isempty(WEIGHT{ii}.Qy),
                if any(any(WEIGHT{ii}.Qy~=0)),
                    allzero = 0;
                    break
                end
            end
        end
        if allzero,
            ney = 0;
        end
        allQyzero = allzero;
    end
else
    allQdzero = 0;
    allQxzero = 0;
    allQuzero = 0;
    allQyzero = 0;
    allQzzero = 0;
end
    
% Define total number of variables nV
%------------------------------------
% Note: we collapse optimizer to the smallest dimension needed.
% For instance if there are no d variables do not introduce slacks for it.
nV  = NT*(nu+nd+nz+neu+ned+nez+nex+ney);


% Handling of the dimensions for the steady state values
%========================================================
u1 = u1(:); % make column
if length(u1) == NT*nu
    % the complete trajectory for u is passed to the file
    UU1 = u1;
elseif length(u1) < NT*nu
    nuref  = length(u1)/nu;  % number of time steps for which u1 was given
    nurefr = round(nuref);   % rounding necessary because of division
    if abs(nurefr-nuref) > 1e-6
        % in u1 there is a number of components that is not a multiple of the
        % number of inputs nu: the vector u1 is incomplete
        error('dimension of u1 is wrong')
    else
        % repetition of the last value of the steady state for u1
        UU1 = [u1; kron(ones(NT-nurefr,1), u1((nurefr-1)*nu+1:length(u1)))];
    end 
else
    warning('u1 too long')
    UU1 = u1(1:NT*nu);
end     

d1 = d1(:); % make column
if length(d1) == NT*nd
    % the complete trajectory for u is passed to the file
    DD1 = d1;
elseif length(d1) < NT*nd
    ndref  = length(d1)/nd;  % number of time steps for which d1 was given
    ndrefr = round(ndref);   % rounding necessary because of division
    if abs(ndrefr-ndref) > 1e-6
        % in d1 there is a number of components that is not a multiple of the
        % number of inputs nd: the vector d1 is incomplete
        error('dimension of d1 is wrong')
    else
        % repetition of the last value of the steady state for d1
        DD1 = [d1; kron(ones(NT-ndrefr,1), d1((ndrefr-1)*nd+1:length(d1)))];
    end   
else
    warning('d1 too long')
    DD1 = d1(1:NT*nd);
end     

z1 = z1(:); % make column
if length(z1) == NT*nz
    % the complete trajectory for u is passed to the file
    ZZ1 = z1;
elseif length(z1) < NT*nz
    nzref  = length(z1)/nz;  % number of time steps for which z1 was given
    nzrefr = round(nzref);   % rounding necessary because of division
    if abs(nzrefr-nzref) > 1e-6
        % in z1 there is a number of components that is not a multiple of the
        % number of inputs nz: the vector z1 is incomplete
        error('dimension of z1 is wrong')
    else
        % repetition of the last value of the steady state for z1
        ZZ1 = [z1; kron(ones(NT-nzrefr,1), z1((nzrefr-1)*nz+1:length(z1)))];
    end   
else
    warning('z1 too long')
    ZZ1 = z1(1:NT*nz);
end 

y1 = y1(:); % make column
if length(y1) == NT*ny
    % the complete trajectory for u is passed to the file
    YY1 = y1;
elseif length(y1) < NT*ny
    nyref  = length(y1)/ny;  % number of time steps for which y1 was given
    nyrefr = round(nyref);   % rounding necessary because of division
    if abs(nyrefr-nyref) > 1e-6
        % in y1 there is a number of components that is not a multiple of the
        % number of inputs ny: the vector y1 is incomplete
        error('dimension of y1 is wrong')
    else
        % repetition of the last value of the steady state for y1
        YY1 = [y1; kron(ones(NT-nyrefr,1), y1((nyrefr-1)*ny+1:length(y1)))];
    end   
else
    warning('y1 too long')
    YY1 = y1(1:NT*ny);
end     

x1 = x1(:); % make column
if length(x1) == NT*nx
    % the complete trajectory for u is passed to the file
    XX1 = x1;
elseif length(x1) < NT*nx
    nxref  = length(x1)/nx;  % number of time steps for which x1 was given
    nxrefr = round(nxref);   % rounding necessary because of division
    if abs(nxrefr-nxref) > 1e-6
        % in x1 there is a number of components that is not a multiple of the
        % number of inputs nx: the vector x1 is incomplete
        error('dimension of x1 is wrong')
    else
        % repetition of the last value of the steady state for x1
        XX1 = [x1; kron(ones(NT-nxrefr,1), x1((nxrefr-1)*nx+1:length(x1)))];
    end   
else
    warning('x1 too long')
    XX1 = x1(1:NT*nx);
end   


%===========================================================
%
% Build cost function and constraint matrices
%
%===========================================================
%===========================================================
% Cost function for the optimizer
%
%       min   1/2 v' S1 v + S2 v + x0' S3 v
%       s.t.  F1 v   <= F2  + F3 x0
%
% We have:
%  U(i)=[u1(i) u2(i) ... u_nu(i)]',
%  D(i)=[d1(i) d2(i) ... d_nd(i)]',
%  Z(i)=[z1(i) z2(i) ... z_nz(i)]',
% with i=0,...,NT-1
%
% Optimizer is:
%  v=[
%     U(0);    U(1); ...    U(NT-1);      %  nu*NT elements
%     D(0);    D(1); ...    D(NT-1);      %  nd*NT elements
%     Z(0);    Z(1); ...    Z(NT-1);      %  nz*NT elements
%     epsu(0); epsu(1); ... epsu(NT-1);   % neu*NT elements
%     epsd(0); epsd(1); ... epsd(NT-1);   % ned*NT elements
%     epsz(0); epsz(1); ... epsz(NT-1);   % nez*NT elements
%     epsx(0); epsx(1); ... epsx(NT-1);   % nex*NT elements
%     epsy(0); epsy(1); ... epsy(NT-1)    % ney*NT elements
%    ]
%
%===========================================================

%============================================
% Constraints matrices  F1 v <= F2 + F3 x(0)
%============================================
F1 = [];
F2 = [];
F3 = [];

% we introduce matrices Ax, Bx and fx for storing expression for x(k)
% x(k) = Ax{k+1} * x(0) + Bx{k+1} * v + fx{k+1}, k=0,...,NT
Ax = {};
Bx = {};
fx = {};
% we introduce matrices Ay, By and fy for storing expression for y(k)
% y(k) = Ay{k+1} * x(0) + By{k+1} * v + fy{k+1}, k=0,...,NT-1
Ay = {};
By = {};
fy = {};


% define offsets for various variables
uoff  = 0;
doff  = uoff  + nu*NT;
zoff  = doff  + nd*NT;
euoff = zoff  + nz*NT;
edoff = euoff + neu*NT;
ezoff = edoff + ned*NT;
exoff = ezoff + nez*NT;
eyoff = exoff + nex*NT;

% indices of integer variables
%IntIndex=doff+1:doff+nd*NT;

aux=[repmat([zeros(nur,1); ones(nub,1)],NT,1);
     repmat(ones(nd,1),NT,1);
     repmat(zeros(nz,1),NT,1)];

IntIndex=find(aux);

% it is helpful to define following Empty Rows
ER_G = zeros(1,nV);
ER_W = zeros(1,1);
ER_S = zeros(1,nx);
ER_A = zeros(1,nx);
ER_B = zeros(1,nV);
ER_f = zeros(1,1);


% initialize Ax{1}, Bx{1} and fx{1} i.e., for the time 0
Ax{1}=eye(nx);
Bx{1}=repmat(ER_B,nx,1);
fx{1}=repmat(ER_f,nx,1);
% well that was easy!


% allocate memory space for F1, F2 and F3 to speed buildmat up
if Options.norm==1 | Options.norm==inf
    %nR = NT*( ne + 2*(nu+nd+nz+nx+ny) );    % number of rows
    nR = NT*( ne + 2*(neu+ned+nez+nex+ney) );    % number of rows
elseif Options.norm==2
    nR = NT*ne;                             % number of rows
else
    error('unknown norm')
end;
if Options.TerminalConstraint
    nR = nR + 2*nx;
end;
if isfield(Options, 'Tset'),
    % increase space by number of constraints of the terminal set
    if isa(Options.Tset, 'polytope'),
        nR = nR + nconstr(Options.Tset);
    end
end

addXconstraints = 0;
if isfield(Options, 'xmax'),
    % increase space to include state constraints
    if all(isinf(Options.xmax)) & all(isinf(Options.xmin)),
        % do not include +/- Inf constraints
    else
        if Options.allXconstrained,
            % add constraints on x(k), k = 0..NT
            nR = nR + 2*nx*(NT+1);
        else
            % add constraints only on x(NT)
            nR = nR + 2*nx;
        end
        addXconstraints = 1;
    end
end

addYconstraints = 0;
if isfield(Options, 'ymax'),
    % increase space to include output constraints
    if all(isinf(Options.ymax)) & all(isinf(Options.ymin)),
        % do not include +/- Inf constraints
    else
        % add constraints on y(k), k = 0..NT-1
        nR = nR + 2*ny*NT;
        addYconstraints = 1;
    end
end

cR = 1;                                 % pointer to the current row
if Options.useSparse,
    F1 = sparse(nR, nV);
    F2 = sparse(nR, 1);
    F3 = sparse(nR, nx);
else
    F1 = NaN*ones(nR, nV);
    F2 = NaN*ones(nR, 1);
    F3 = NaN*ones(nR, nx);
end

% Now let's do a more tricky part: loop through all time steps
% Have in mind that MLD system is a time varying system....

for kk=1:NT    % Actual time is given as k = kk-1, so this corresponds to k = 0:NT-1

    % make shorcuts
    SYS=SYSTEM{t2s(kk)};
    WEI=WEIGHT{t2w(kk)};
    
    % Note that expression for x(k) already exists
    %---------------------------------------------
    
    % Express y(k) = SYS.C  * x(k) + SYS.D1 * u(k) + SYS.D2 * d(k) +...
    %              + SYS.D3 * z(k) + SYS.D5
    %-----------------------------------------------------------------------------------------
    Ay{kk} = SYS.C*Ax{kk};
    By{kk} = SYS.C*Bx{kk};
    By{kk}(:,uoff+(kk-1)*nu+1:uoff+kk*nu) = By{kk}(:,uoff+(kk-1)*nu+1:uoff+kk*nu) + SYS.D1;
    By{kk}(:,doff+(kk-1)*nd+1:doff+kk*nd) = By{kk}(:,doff+(kk-1)*nd+1:doff+kk*nd) + SYS.D2;
    By{kk}(:,zoff+(kk-1)*nz+1:zoff+kk*nz) = By{kk}(:,zoff+(kk-1)*nz+1:zoff+kk*nz) + SYS.D3;
    fy{kk} = SYS.C*fx{kk} + SYS.D5;
    
    % Express x(k+1) = SYS.A  * x(k) + SYS.B1 * u(k) + SYS.B2 * d(k) +...
    %                + SYS.B3 * z(k) + SYS.B5
    %-----------------------------------------------------------------------------------------
    Ax{kk+1} = SYS.A*Ax{kk};
    Bx{kk+1} = SYS.A*Bx{kk};
    Bx{kk+1}(:,uoff+(kk-1)*nu+1:uoff+kk*nu) = Bx{kk+1}(:,uoff+(kk-1)*nu+1:uoff+kk*nu) + SYS.B1;
    Bx{kk+1}(:,doff+(kk-1)*nd+1:doff+kk*nd) = Bx{kk+1}(:,doff+(kk-1)*nd+1:doff+kk*nd) + SYS.B2;
    Bx{kk+1}(:,zoff+(kk-1)*nz+1:zoff+kk*nz) = Bx{kk+1}(:,zoff+(kk-1)*nz+1:zoff+kk*nz) + SYS.B3;
    fx{kk+1} = SYS.A*fx{kk} + SYS.B5;
    
    %=========================================================================
    % Let's formulate inequalities!
    % For this purpose we will use temp variables: A x(0) + B * v + f <= 0
    % Note that the outside form of the constraints is: F1 v <= F2 + F3 * x(0)
    %=========================================================================
    
    % MLD inequalities
    %=================
    
    % SYS.E2 * d(k) + SYS.E3 * z(k) <=   SYS.E1 * u(k) +
    %                                              + SYS.E4 * x(k) + SYS.E5
    %--------------------------------------------------------------------------------------
%     A = -SYS.E4*Ax{kk};
%     B = -SYS.E4*Bx{kk};
%     B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) = B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) - SYS.E1;
%     B(:,doff+(kk-1)*nd+1:doff+kk*nd) = B(:,doff+(kk-1)*nd+1:doff+kk*nd) + SYS.E2;
%     B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) = B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) + SYS.E3;
%     f = -SYS.E4*fx{kk} - SYS.E5;
%     
%     F1(cR:cR+ne-1,:) = B;
%     F2(cR:cR+ne-1,:) = -f;
%     F3(cR:cR+ne-1,:) = -A;
%     cR = cR + ne;
    
    % faster implementation of the same:
    A = -SYS.E4*Ax{kk};
    F1(cR:cR+ne-1,:) = -SYS.E4*Bx{kk};
    F1(cR:cR+ne-1,uoff+(kk-1)*nu+1:uoff+kk*nu) = F1(cR:cR+ne-1,uoff+(kk-1)*nu+1:uoff+kk*nu) - SYS.E1;
    F1(cR:cR+ne-1,doff+(kk-1)*nd+1:doff+kk*nd) = F1(cR:cR+ne-1,doff+(kk-1)*nd+1:doff+kk*nd) + SYS.E2;
    F1(cR:cR+ne-1,zoff+(kk-1)*nz+1:zoff+kk*nz) = F1(cR:cR+ne-1,zoff+(kk-1)*nz+1:zoff+kk*nz) + SYS.E3;
    f = -SYS.E4*fx{kk} - SYS.E5;
    
    F2(cR:cR+ne-1,:) = -f;
    F3(cR:cR+ne-1,:) = -A;
    cR = cR + ne;
    
    % slack inequalities
    %===================
    
    % in the case of Options.norm==2 there are no slack variables
    if Options.norm==2
        continue;
    end
    
    if ~allQuzero,
        % - eu(k) <= + WEI.Qu * (u(k)-u_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nu,1);
        B = repmat(ER_B,nu,1);
        B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) = B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) - WEI.Qu;
        B(:,euoff+(kk-1)*neu+1:euoff+kk*neu) = B(:,euoff+(kk-1)*neu+1:euoff+kk*neu) - eye(neu);
        f = WEI.Qu * UU1((kk-1)*nu+1:kk*nu);
        
        F1(cR:cR+nu-1,:) = B;
        F2(cR:cR+nu-1,:) = -f;
        F3(cR:cR+nu-1,:) = -A;
        cR = cR + nu;
        
        % - eu(k) <= - WEI.Qu * (u(k)-u_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nu,1);
        B = repmat(ER_B,nu,1);
        B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) = B(:,uoff+(kk-1)*nu+1:uoff+kk*nu) + WEI.Qu;
        B(:,euoff+(kk-1)*neu+1:euoff+kk*neu) = B(:,euoff+(kk-1)*neu+1:euoff+kk*neu) - eye(neu);
        f = -WEI.Qu * UU1((kk-1)*nu+1:kk*nu);
        
        F1(cR:cR+nu-1,:) = B;
        F2(cR:cR+nu-1,:) = -f;
        F3(cR:cR+nu-1,:) = -A;
        cR = cR + nu;
    end

    if ~allQdzero,
        % - ed(k) <= + WEI.Qd * (d(k)-d_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nd,1);
        B = repmat(ER_B,nd,1);
        B(:,doff+(kk-1)*nd+1:doff+kk*nd) = B(:,doff+(kk-1)*nd+1:doff+kk*nd) - WEI.Qd;
        B(:,edoff+(kk-1)*ned+1:edoff+kk*ned) = B(:,edoff+(kk-1)*ned+1:edoff+kk*ned) - eye(ned);
        f = WEI.Qd * DD1((kk-1)*nd+1:kk*nd);
        
        F1(cR:cR+nd-1,:) = B;
        F2(cR:cR+nd-1,:) = -f;
        F3(cR:cR+nd-1,:) = -A;
        cR = cR + nd;
    
        % - ed(k) <= - WEI.Qd * (d(k)-d_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nd,1);
        B = repmat(ER_B,nd,1);
        B(:,doff+(kk-1)*nd+1:doff+kk*nd) = B(:,doff+(kk-1)*nd+1:doff+kk*nd) + WEI.Qd;
        B(:,edoff+(kk-1)*ned+1:edoff+kk*ned) = B(:,edoff+(kk-1)*ned+1:edoff+kk*ned) - eye(ned);
        f = -WEI.Qd * DD1((kk-1)*nd+1:kk*nd);
        
        F1(cR:cR+nd-1,:) = B;
        F2(cR:cR+nd-1,:) = -f;
        F3(cR:cR+nd-1,:) = -A;
        cR = cR + nd;
    end
    
    if ~allQzzero,
        % - ez(k) <= + WEI.Qz * (z(k)-z_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nz,1);
        B = repmat(ER_B,nz,1);
        B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) = B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) - WEI.Qz;
        B(:,ezoff+(kk-1)*nez+1:ezoff+kk*nez) = B(:,ezoff+(kk-1)*nez+1:ezoff+kk*nez) - eye(nez);
        f = WEI.Qz * ZZ1((kk-1)*nz+1:kk*nz);
        
        F1(cR:cR+nz-1,:) = B;
        F2(cR:cR+nz-1,:) = -f;
        F3(cR:cR+nz-1,:) = -A;
        cR = cR + nz;
        
        % - ez(k) <= - WEI.Qz * (z(k)-z_e(k))
        %-------------------------------------------
        A = repmat(ER_A,nz,1);
        B = repmat(ER_B,nz,1);
        B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) = B(:,zoff+(kk-1)*nz+1:zoff+kk*nz) + WEI.Qz;
        B(:,ezoff+(kk-1)*nez+1:ezoff+kk*nez) = B(:,ezoff+(kk-1)*nez+1:ezoff+kk*nez) - eye(nez);
        f = -WEI.Qz * ZZ1((kk-1)*nz+1:kk*nz);

        F1(cR:cR+nz-1,:) = B;
        F2(cR:cR+nz-1,:) = -f;
        F3(cR:cR+nz-1,:) = -A;
        cR = cR + nz;
    end
    
    if ~allQxzero,
        % - ex(k) <= + WEI.Qx * (x(k)-x_e(k))
        %-------------------------------------------
        A = -WEI.Qx * Ax{kk};
        B = -WEI.Qx * Bx{kk};
        B(:,exoff+(kk-1)*nex+1:exoff+kk*nex) = B(:,exoff+(kk-1)*nex+1:exoff+kk*nex) - eye(nex);
        f = -WEI.Qx * (fx{kk} - XX1((kk-1)*nx+1:kk*nx));
        F1(cR:cR+nx-1,:) = B;
        F2(cR:cR+nx-1,:) = -f;
        F3(cR:cR+nx-1,:) = -A;
        cR = cR + nx;
    
        % - ex(k) <= - WEI.Qx * (x(k)-x_e(k))
        %-------------------------------------------
        A = WEI.Qx * Ax{kk};
        B = WEI.Qx * Bx{kk};
        B(:,exoff+(kk-1)*nex+1:exoff+kk*nex) = B(:,exoff+(kk-1)*nex+1:exoff+kk*nex) - eye(nex);
        f = WEI.Qx * (fx{kk} - XX1((kk-1)*nx+1:kk*nx));
    
        F1(cR:cR+nx-1,:) = B;
        F2(cR:cR+nx-1,:) = -f;
        F3(cR:cR+nx-1,:) = -A;
        cR = cR + nx;
    end

    if ~allQyzero,
        % - ey(k) <= + WEI.Qy * (y(k)-y_e(k))
        %-------------------------------------------
        A = -WEI.Qy * Ay{kk};
        B = -WEI.Qy * By{kk};
        B(:,eyoff+(kk-1)*ney+1:eyoff+kk*ney) = B(:,eyoff+(kk-1)*ney+1:eyoff+kk*ney) - eye(ney);
        f = -WEI.Qy * (fy{kk} - YY1((kk-1)*ny+1:kk*ny));
        
        
        F1(cR:cR+ny-1,:) = B;
        F2(cR:cR+ny-1,:) = -f;
        F3(cR:cR+ny-1,:) = -A;
        cR = cR + ny;
    
        % - ey(k) <= - WEI.Qy * (y(k)-y_e(k))
        %-------------------------------------------
        A = WEI.Qy * Ay{kk};
        B = WEI.Qy * By{kk};
        B(:,eyoff+(kk-1)*ney+1:eyoff+kk*ney) = B(:,eyoff+(kk-1)*ney+1:eyoff+kk*ney) - eye(ney);
        f = WEI.Qy * (fy{kk} - YY1((kk-1)*ny+1:kk*ny));
        
        F1(cR:cR+ny-1,:) = B;
        F2(cR:cR+ny-1,:) = -f;
        F3(cR:cR+ny-1,:) = -A;
        cR = cR + ny;
    end
end

if Options.TerminalConstraint
    % Terminal state constraint
    %==========================
    
    % x(NT) <= xtt + eps2
    %--------------------
    A = Ax{NT+1};
    B = Bx{NT+1};
    f = fx{NT+1} - xtt - eps2;
    
    F1(cR:cR+nx-1,:) = B;
    F2(cR:cR+nx-1,:) = -f;
    F3(cR:cR+nx-1,:) = -A;
    cR = cR + nx;
    
    % x(NT) >= xtt - eps2
    %--------------------
    A = -Ax{NT+1};
    B = -Bx{NT+1};
    f = -fx{NT+1} + xtt - eps2;
    
    F1(cR:cR+nx-1,:) = B;
    F2(cR:cR+nx-1,:) = -f;
    F3(cR:cR+nx-1,:) = -A;
    cR = cR + nx;
end

if addXconstraints,
    % state constraints
    %==================
    
    if Options.allXconstrained,
        % add state constraints on every state
        C_horizon = 1:NT+1;
    else
        % otherwise just add constraints on final state
        C_horizon = NT+1;
    end
    for iN = C_horizon,
        % x(k) <= xmax
        %-------------
        A = Ax{iN};
        B = Bx{iN};
        f = fx{iN} - Options.xmax;
        
        F1(cR:cR+nx-1,:) = B;
        F2(cR:cR+nx-1,:) = -f;
        F3(cR:cR+nx-1,:) = -A;
        cR = cR + nx;
        
        % x(k) >= xmin
        %--------------
        A = -Ax{iN};
        B = -Bx{iN};
        f = -fx{iN} + Options.xmin;
        
        F1(cR:cR+nx-1,:) = B;
        F2(cR:cR+nx-1,:) = -f;
        F3(cR:cR+nx-1,:) = -A;
        cR = cR + nx;
    end
end

if addYconstraints,
    % output constraints
    %===================
    
    for iN = 1:NT,
        % y(k) <= ymax
        %-------------
        A = Ay{iN};
        B = By{iN};
        f = fy{iN} - Options.ymax;
        
        F1(cR:cR+ny-1,:) = B;
        F2(cR:cR+ny-1,:) = -f;
        F3(cR:cR+ny-1,:) = -A;
        cR = cR + ny;
        
        % y(k) >= ymin
        %--------------
        A = -Ay{iN};
        B = -By{iN};
        f = -fy{iN} + Options.ymin;
        
        F1(cR:cR+ny-1,:) = B;
        F2(cR:cR+ny-1,:) = -f;
        F3(cR:cR+ny-1,:) = -A;
        cR = cR + ny;
    end
end

if isfield(Options, 'Tset')
    % Terminal set constraint
    %==========================

    if ~isa(Options.Tset, 'polytope'),
        error('Terminal set must be a polytope object.');
    end
    
    [Hfinal, Kfinal] = double(Options.Tset);
    Tset_nr = nconstr(Options.Tset);
    
    % Hf * x(NT) <= Kf
    %-----------------
    A = Ax{NT+1};
    B = Bx{NT+1};
    f = fx{NT+1};
    
    F1(cR:cR+Tset_nr-1,:) = Hfinal * B;
    F2(cR:cR+Tset_nr-1,:) = Kfinal - Hfinal * f;
    F3(cR:cR+Tset_nr-1,:) = -Hfinal * A;
    cR = cR + Tset_nr;
end    
    

%===================================
% Let's formulate a cost function!
%===================================

%===========================================
% J (v,x0) = 1/2 v' S1 v + S2 v + x0' S3 v
%===========================================
if Options.useSparse,
    S1 = sparse(nV, nV);
    S2 = sparse(1, nV);
    S3 = sparse(nx, nV);
else
    S1=zeros(nV,nV);
    S2=zeros(1,nV);
    S3=zeros(nx,nV);
end

% c1,c2,c3: parameters for the constant term of the optimization
%           the complete constant term is obtained as: 
%           Jc = 1/2 x(0)'*c1*x(0) + c2'*x(0) + c3
%---------------------------------------------------------------
c1=zeros(nx,nx);
c2=zeros(1,nx);
c3=0;

if Options.norm==inf | Options.norm==1
    S2(:,(nu+nd+nz)*NT+1:end)=1;
    
elseif Options.norm==2
    
    % somebody should spend some time on this issue and check it
    for kk=1:NT
        
        % make a shorcut
        WEI=WEIGHT{t2w(kk)};

        
        % u(k)'*WEI.Qu*u(k) - 2*u_e(k)'*WEI.Qu*u(k) + u_e(k)'*WEI.Qu*u_e(k)
        %------------------------------------------------------------------------------------
        ii=uoff+(kk-1)*nu+1:uoff+kk*nu;
        S1(ii,ii) = S1(ii,ii) + 2 * WEI.Qu;
        S2(:,ii) = S2(:,ii) - 2 * UU1((kk-1)*nu+1:kk*nu)'*WEI.Qu;
        c3 = c3 + UU1((kk-1)*nu+1:kk*nu)'*WEI.Qu*UU1((kk-1)*nu+1:kk*nu);
        
        % d(k)'*WEI.Qd*d(k) - 2*d_e(k)'*WEI.Qd*d(k) + d_e(k)'*WEI.Qd*d_e(k)
        %------------------------------------------------------------------------------------
        ii=doff+(kk-1)*nd+1:doff+kk*nd;
        S1(ii,ii) = S1(ii,ii) + 2 * WEI.Qd;
        S2(:,ii) = S2(:,ii) - 2 * DD1((kk-1)*nd+1:kk*nd)'*WEI.Qd;
        c3 = c3 + DD1((kk-1)*nd+1:kk*nd)'*WEI.Qd*DD1((kk-1)*nd+1:kk*nd);
        
        % z(k)'*WEI.Qz*z(k) - 2*z_e(k)'*WEI.Qz*z(k) + z_e(k)'*WEI.Qz*z_e(k)
        %------------------------------------------------------------------------------------
        ii=zoff+(kk-1)*nz+1:zoff+kk*nz;
        S1(ii,ii) = S1(ii,ii) + 2 * WEI.Qz;
        S2(:,ii) = S2(:,ii) - 2 * ZZ1((kk-1)*nz+1:kk*nz)'*WEI.Qz;
        c3 = c3 + ZZ1((kk-1)*nz+1:kk*nz)'*WEI.Qz*ZZ1((kk-1)*nz+1:kk*nz);
        
        
        % For the following have in mind that
        % x(k) = Ax{k} x(0) + Bx{k} v + fx{k}
        % y(k) = Ay{k} x(0) + By{k} v + fy{k}
        
        % x(k)'*WEI.Qx*x(k) - 2*x_e(k)'*WEI.Qx*x(k) + x_e(k)'*WEI.Qx*x_e(k)
        %------------------------------------------------------------------------------------
        % after a lot of checking...

        aux = fx{kk} - XX1((kk-1)*nx+1:kk*nx);
        
        S1 = S1 + 2 * Bx{kk}' * WEI.Qx * Bx{kk};
        
        % small hack, without it (x-x_e) is not being penalized!
        % before it was:
        % S2 = S2 + 2 * fx{kk}' * WEI.Qx * Bx{kk};
        S2 = S2 + 2 * aux' * WEI.Qx * Bx{kk};
        
        S3 = S3 + 2 * Ax{kk}' * WEI.Qx * Bx{kk};
        

        c1 = c1 + 2 * Ax{kk}' * WEI.Qx * Ax{kk};
        c2 = c2 + 2 * aux' * WEI.Qx * Ax{kk};
        c3 = c3 + aux' * WEI.Qx * aux;
        
        % y(k)'*WEI.Qy*y(k) - 2*y_e(k)'*WEI.Qy*y(k) + y_e(k)'*WEI.Qy*y_e(k)
        %------------------------------------------------------------------------------------
        % after a lot of checking...
        
        aux = fy{kk} - YY1((kk-1)*ny+1:kk*ny);
        
        S1 = S1 + 2 * By{kk}' * WEI.Qy * By{kk};
        
        % small hack, without it (y-y_e) is not being penalized!
        % before:
        % S2 = S2 + 2 * fy{kk}' * WEI.Qy * By{kk};
        S2 = S2 + 2 * aux' * WEI.Qy * By{kk};
        
        S3 = S3 + 2 * Ay{kk}' * WEI.Qy * By{kk};
        
        c1 = c1 + 2 * Ay{kk}' * WEI.Qy * Ay{kk};
        c2 = c2 + 2 * aux' * WEI.Qy * Ay{kk};
        c3 = c3 + aux' * WEI.Qy * aux;
        
    end
    
else
    error('Unknown norm has been specified!');
end



% return additional information:
%------------------------------------------------------------------------------------

% x(k) = Ax{k+1} * x(0) + Bx{k+1} * v + fx{k+1}, k=0,...,NT
Ext.OL.A = Ax;
Ext.OL.B = Bx;
Ext.OL.f = fx;

% y(k) = Ay{k+1} * x(0) + By{k+1} * v + fy{k+1}, k=0,...,NT-1
Ext.OL.C = Ay;
Ext.OL.D = By;
Ext.OL.g = fy;


%------------------------------------------------------------------------------------
function S = sub_remove_redundant_ineq(S)
% removes redundant rows introduced in AD section

neb = S.ne;
removerows = [];
for ii = 1:length(S.rowinfo.ineq),
    ineq = S.rowinfo.ineq{ii};
    if strcmp(ineq.section, 'AD'),
        if ineq.subindex == 2,
            % remove rows with 'subindex=2' because they define
            %   [delta = 1] -> [f(x) <= 0]
            % while we only keep rows which say
            %  [f(x) <= 0] -> [delta = 1]
            removerows = [removerows; ii];
        end
    end
end

S.E1(removerows,:) = [];
S.E2(removerows,:) = [];
S.E3(removerows,:) = [];
S.E4(removerows,:) = [];
S.E5(removerows,:) = [];
S.ne = size(S.E1, 1);

%fprintf('%d rows originally, %d after reduction\n', neb, S.ne);
