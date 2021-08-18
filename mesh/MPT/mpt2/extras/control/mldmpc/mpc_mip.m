function [ut, dt, zt, Eflag] = mpc_mip( S, xt, r, Q, pr, co, Options )

%===============================================================================
%
% Title:       mpc_mip.m
%
% Version:     2.0
%                                                                       
% Project:     Control of MLD systems
%                                                                       
% Purpose:     MPC of MLD systems                  
%                                                                        
% Authors:     Tobias Geyer, Domenico Mignone, Mato Baotic
%              based on mipc1.m by Alberto Bemporad, 12/5/1998 
%              (C) 2005 by Michal Kvasnica
%
% History:     date:       subject:                             author:
%              -----       --------                             -------
%              2003.11.19  first public release                 Tobias Geyer
%              2003.03.17  time-varying models                  Tobias & George
%              2002.05.21  Added flag, removed output messages  Tobias
%              2002.05.17  Added constraints on z variables     Mato & Tobias
%              2002.04.19  Relaxation of MLD and terminal       Mato Baotic
%                          state constraints enabled via
%                          Options.relaxind
%              2002.04.13  Included infinity norm case          Mato Baotic
%                          (Options.norm = inf)
%              2000.11.16  Initial version adapted from         Domenico Mignone
%                          mipc_main.m
%
% Inputs:     S:  the MLD model in compact or valid format
%                 If S is a cell array, every cell corresponds to a different
%                 MLD model with corresponding prediction horizons given in the
%                 vectors pr and co.
%             x:  the current state
%             r:  the reference value (see below for the possibilities)
%             Q:  structure with the weight matrices (see below for details)
%             pr: vector of prediction horizons
%             co: vector of control horizons
%             Options: further inputs (see below for the possibilities)
%
%             All inputs, except Options are mandatory
%
%             Details on some inputs:
%
%             r: if r has dimension ny, it is understood as constant reference
%                     for y
%                if r has dimension pr*ny, it is understood as a time varying
%                     reference trajectory for y
%                if r is a structure, the following fields are accepted as 
%                     references for the involved variables:
%                     r.x reference for the state
%                     r.u reference for the input
%                     r.d reference for delta
%                     r.z reference for z
%                     r.y reference for the output
%                     each field is assumed to be a constant reference, if it is
%                     of dimension nx, nu, nd, nz, ny and is assumed to be a
%                     time varying reference if the dimension is pr*nx, pr*nu...
%
%             Q: the following fields can be specified
%                Q.Qx weights on x
%                Q.Qu weights on u
%                Q.Qd weights on d
%                Q.Qz weights on z
%                Q.Qy weights on y
%
%             Options accepts the following fields:
%                Options.bigB             infinity for bounds on variables
%                                         (default = 1e6)
%                Options.epsil            global softening of ax <= b
%                                         (default = 1e-7)
%                Options.solver           MIQP solver for panmiqp
%                                         (default = 'miqp')
%                Options.solver.options   specific solver options
%                                         (default = [])
%                Options.saveflag         flag denoting, whether workspace is
%                                         saved before/after optimization
%                                         (default = 0)
%                Options.savestr          name of files for saving
%                                         (default = 'data')
%                Options.eps2             tolerance on fulfillment of terminal 
%                                         state constraints
%                                         (default = 0)
%                Options.umin             lower bounds for inputs
%                                         (default = -bigB)
%                Options.umax             upper bounds for inputs
%                                         (default = bigB)
%                Options.zmin             lower bounds for aux variable
%                                         (default = -bigB)
%                Options.zmax             upper bounds for aux variable
%                                         (default = bigB)
%                Options.flagihard        if 1 hard constraints on Ax<=b
%                                         if 0 epsilon softening of Ax<=b
%                                         (default = 1)
%                Options.weps             weight in cost function if flagihard=0
%                                         (default = 1)
%                Options.relaxind         indices of constraints Ax<=b to be
%                                         relaxed if flagihard==0. Only MLD and
%                                         terminal state constraints can be relaxed.
%                                         (default = [1:ne*pr+2*nx], i.e., relax
%                                         *ALL* MLD and terminal state constraints
%                                         and use only *ONE* slack variable)
%                                         Note: if Options.relaxindex==[], relax
%                                               *ALL* MLD and terminal state
%                                               constraints and use *DIFFERENT*
%                                               slack variable for each constraint
%                Options.xopt             initial guess for the optimization
%                                         (default = 0)
%                Options.norm             norm used in MPC formulation
%                                         {2 | inf | 1}
%                                         (default = 2)
%                Options.verbose          display additional information
%                                         (default = 0)
%                Options.removeInfBounds  if true, removes constraints of the
%                                         form z <= Inf and replaces +/- Inf
%                                         bounds on variables by +/- Options.bigB.
%                                         (default = 0)
%                Options.bounds2ineq      if true, converts bounds of variables
%                                         into inequality constraints
%                                         (default = 0)
%
%
% Outputs:     ut    : current input to be applied to the system
%              dt    : current value of the variable delta
%              zt    : current value of the variable z
%              Eflag : structure for further output informations
%                      Eflag.slkeps: slack variables of constraints if flagihard==0
%
% Notes:       Soft constraints on the whole set of inequalities (flagihard=0)
%              The problem becomes:  min  J + ||weps, epsil||
%                                    s.t. Ax-epsil*ones(nlin,1) <= B
%                                                         epsil >= 0 
%              constraint relaxations on u (flaguhard = 0)
%              are only effective, if all constraints on u are expressed via
%              umin and umax. If they are already included in the MLD matrices
%              Ei, this relaxation has no effect.
%
%              This is an extension to implement multiple sampling times within
%              the prediction horizon. Therefore, each model has the same dimensions
%              but different associated prediction horizons (and sampling times).
%
% Requires:    mpc_checkSysWeight
%              mpc_buildmatFAST
%
% Contact:     Mato Baotic, baotic@control.ee.ethz.ch
%              Tobias Geyer, geyer@control.ee.ethz.ch
%              Michal Kvasnica, kvasnica@control.ee.ethz.ch
%
%              Automatic Control Laboratory
%              ETH Zentrum,
%              Zurich, Switzerland
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

global mptOptions

if ~isstruct(mptOptions)
    mpt_error;
end

error(nargchk(6,7,nargin));

if nargin < 7
   Options = [];
end      



% defaults for undefined variables, according to comments above
% ---------------------------------------------------------------------

if ~isfield(Options,'bigB')
   bigB = 1e6;
else
   bigB = Options.bigB;   
end  
if ~isfield(Options,'epsil')
   epsil = 1e-10;
    if ~isempty(mptOptions.milpsolver),
        if mptOptions.milpsolver == 2,
            % set "epsil" to zero for GLPK, noticed some numerical problems when a
            % non-zero value is used
            epsil = 0;
        end
    end
else
   epsil = Options.epsil;   
end
if ~isfield(Options,'solver')
   Options.solver = 'miqp'; 
end
if ~isfield(Options,'saveflag')
   saveflag = 0;
else
   saveflag = Options.saveflag;   
end   
if ~isfield(Options,'savestr')
   savestr = 'data';
else
   savestr = Options.savestr;   
end 
if ~isfield(Options,'eps2')
   eps2 = 0;
else
   eps2 = Options.eps2;   
end
if ~isfield(Options,'flagihard')
   flagihard = 1;
else
   flagihard = Options.flagihard;   
end
if ~isfield(Options,'weps')
   weps = 1;
else
   weps = Options.weps;   
end
if ~isfield(Options,'norm')
    Options.norm = 2;
end
if ~isfield(Options,'verbose')
    Options.verbose = 0;
end;
if ~isfield(Options, 'returnproblem'),
    % if true, mpc_mip will exit immediatelly after it constructs matrices of
    % the problem
    Options.returnproblem = 0;
end
if ~isfield(Options, 'problemmatrices'),
    % if set, it has to contain matrices of the problem constructed beforehand
    Options.problemmatrices = [];
end
if ~isfield(Options, 'convert2eq'),
    % if true, detects inequality constraints which form equality constraints
    Options.convert2eq = 0;
end
if ~isfield(Options, 'removeInfBounds'),
    % if true, removes constraints of the form z <= Inf and replaces +/- Inf
    % bounds on variables by +/- 1e9, respectively
    Options.removeInfBounds = 0;
end
if ~isfield(Options, 'bounds2ineq'),
    % if true, converts bounds of variables ( bl <= z <= bu) into inequality
    % constraints
    Options.bounds2ineq = 0;
end


% check dimensions and turn S and Q into cell structures
% ---------------------------------------------------------------------

% check dimensions                             
[S, Q] = mpc_checkSysWeight(S, Q, Options.verbose);
% S and Q are now cell structures for sure

% Dimensions
nu = size(S{1}.B1,2);   % nu=dimension of u
nd = size(S{1}.B2,2);   % nd=dimension of \delta
nz = size(S{1}.B3,2);   % nz=dimension of z
nx = size(S{1}.A,2);    % nx=dimension of x
ne = size(S{1}.E5,1);   % ne=dimension of E_5
ny = size(S{1}.C,1);    % ny=number of outputs

if pr == co
    horizon = pr;
else
    error('different horizons for control and prediction are not implemented yet')
end;

if length(S) ~= length(horizon)
    error('S should be a cell structure of length equal to the number of prediction horizons in pr or co');
end

if ~(length(Q)==length(horizon) | length(Q)==1)
    error('Q should be a cell structure of length equal to 1 or the number of prediction horizons in pr or co');
end

% total length of horizon
NT=sum(horizon);




% some more defaults for undefined variables, according to comments above
% ---------------------------------------------------------------------

if ~isfield(Options,'umin')
   umin = -bigB*ones(1,nu);
else
   umin = Options.umin;   
end
if ~isfield(Options,'umax')
   umax = bigB*ones(1,nu);
else
   umax = Options.umax;   
end
if ~isfield(Options,'zmin')
   zmin = -bigB*ones(1,nz);
else
   zmin = Options.zmin;   
end
if ~isfield(Options,'zmax')
   zmax = bigB*ones(1,nz);
else
   zmax = Options.zmax;   
end




% extract value for the reference
% ---------------------------------------------------------------------

if isstruct(r) % r has been recognized to be a structure
    if isfield(r, 'y'),
        y1 = r.y;
    else
        y1 = zeros(ny,1);
    end
    if isfield(r, 'x')
        x1 = r.x;
    else
        x1 = zeros(nx,1);
    end
    if isfield(r, 'u')
        u1 = r.u;
    else
        u1 = zeros(nu,1);
    end
    if isfield(r, 'd'),
        d1 = r.d;
    else
        d1 = zeros(nd,1);
    end
    if isfield(r, 'z'),
        z1 = r.z;
    else
        z1 = zeros(nz,1);
    end
else
    % do not use terminal state constraint if reference on outputs is provided.
    % otherwise there is a contradiction - we want to reach certain reference
    % but we enforce x_N = 0 by default in mpc_buildmatFAST
    Options.TerminalConstraint = 0;
    y1 = r;
    x1 = zeros(nx,1); 
    u1 = zeros(nu,1);
    d1 = zeros(nd,1);
    z1 = zeros(nz,1);
end



% get row vectors
% ---------------------------------------------------------------------

umin = umin(:);
umax = umax(:);
zmin = zmin(:);
zmax = zmax(:);

xt = xt(:);      % current state x(t) 



% build optimization problem
% ---------------------------------------------------------------------

% terminal state constraint
xtt = x1(:,end);

F1eq = [];
F2eq = [];
F3eq = [];
matrices = [];
if isempty(Options.problemmatrices),
    % Build optimization matrices
    [S1, S2, S3, F1, F2, F3, c1, c2, c3, IntIndex, Ext] = mpc_buildmatFAST(horizon, ...
        S, Q, x1, u1, d1, z1, y1, eps2, xtt, Options) ;
    
    if Options.returnproblem,
        % only return matrices of the constructed problem
        matrices.S1 = S1;
        matrices.S2 = S2;
        matrices.S3 = S3;
        matrices.F1 = F1;
        matrices.F2 = F2;
        matrices.F3 = F3;    
        matrices.c1 = c1;
        matrices.c2 = c2;
        matrices.c3 = c3;    
        matrices.IntIndex = IntIndex;
        matrices.Ext = Ext;
        ut = matrices;
        return
    end
else
    % use matrices constructed beforehand
    matrices = Options.problemmatrices;
    S1 = matrices.S1;
    S2 = matrices.S2;
    S3 = matrices.S3;
    F1 = matrices.F1;
    F2 = matrices.F2;
    F3 = matrices.F3;
    c1 = matrices.c1;
    c2 = matrices.c2;
    c3 = matrices.c3;
    if isfield(matrices, 'F1eq'),
        if ~isfield(matrices, 'F2eq') | ~isfield(matrices, 'F3eq'),
            error('F2eq or F3eq missing in Options.problemmatrices.');
        end
        F1eq = matrices.F1eq;
        F2eq = matrices.F2eq;
        F3eq = matrices.F3eq;
    end
    IntIndex = matrices.IntIndex;
    Ext = matrices.Ext;
end

% number of (linear) constraints
nlin = length(F2);

% number of variables (=length of optimizer)
nvar = length(S2);

% Since panmip solves the problem:
%	     
% minimize 0.5 x' G x  + C' x
% subject to  Ax <= B
%
% S1 would have to be multiplied by 2. By definition (see mipc_optmat.m)
% also CC must be multiplied by 2. Both factors are ommited, since the 
% optimizer xt is not affected by them.
G  = S1;
CC = (S2+xt'*S3)';
AA = F1;
B  = F2+F3*xt;
AAeq = [];
BBeq = [];
if ~isempty(F1eq),
    AAeq = F1eq;
    BBeq = F2eq + F3eq*xt;
else
    eqconverted = 0;
    if isfield(matrices, 'info'),
        if isfield(matrices.info, 'eqconverted'),
            eqconverted = matrices.info.eqconverted;
        end
    end
    if ~eqconverted,
        if Options.convert2eq,
            ne = size(AA, 1);
            [AA, B, AAeq, BBeq] = mpt_ineq2eq(AA, B);
            fprintf('%d inequalities originally / %d inequalities replaced by %d equalities\n', ne, 2*length(BBeq), length(BBeq));
            nlin = length(B);
        else
            AAeq = [];
            BBeq = [];
        end
    end
end


% initial solution of optimizer
% ---------------------------------------------------------------------
if ~isfield(Options,'xopt') % initial solution of optimizer
    Options.xopt = zeros(nvar,1); 
end

xopt = Options.xopt;   

if length(xopt)<nvar
   warning('Provided initial solution is to short! Populating it with zeros!');
   xopt = [xopt(:); zeros(nvar-length(xopt),1)];
end

% Initial guess for MIQP=previous optimal V, shifted
% (actually it should be not completed by the new u1,d1,z1, but by the
% last u1,d1,z1 which gave feasible solution. No difference if the
% disturbance is constant after t>=pr)
% check this!!!

Vt1 = [xopt(nu+1:nu*NT); u1(:,1);
       xopt(nu*NT+nd+1:nu*NT+nd*NT); d1(:,1);
       xopt(nu*NT+nd*NT+nz+1:nu*NT+nd*NT+nz*NT); z1(:,1)];

% x0=initial guess for optimization

x0 = [Vt1; zeros(nvar-length(Vt1),1)];   % we introduce zeros in the places of (potential) slack variables



% fix additional constraints on u and z (deltas are ALWAYS in [0 1])
% ---------------------------------------------------------------------

% constraints on u

bl = [kron(ones(NT,1),umin); zeros(nd*NT,1)];
bu = [kron(ones(NT,1),umax);  ones(nd*NT,1)];

% constraints on z    

bl = [bl; kron(ones(NT,1),zmin)];
bu = [bu; kron(ones(NT,1),zmax)];  

% extend number of constraints if necessary (for 1- and inf-norm)

bl = [bl;     zeros(nvar-length(bl),1)];
bu = [bu; bigB*ones(nvar-length(bu),1)];




% relaxation of some/all of the MLD and/or terminal state constraints
% -------------------------------------------------------------------
% min J+ ||weps, ESP1||
% Ax-EPS1*ones(nlin,nEPS1)<=B
% 0<=EPSl<bigB

if ~flagihard
    
    %error('this part of the code needs to be rewritten')
    
    if ~isfield(Options,'relaxind')
        relaxind = [1:ne*NT+2*nx];
        nEPS1=1;
    elseif isempty(Options.relaxind);
        relaxind = [1:ne*NT+2*nx];
        nEPS1=ne*NT+2*nx;
    else
        relaxind = Options.relaxind;
        nEPS1=length(relaxind);
        if any(relaxind>ne*NT+2*nx | relaxind<1)
            error('relaxation indices out of range');
        end
    end
    EPS1 = zeros(nlin,nEPS1); 
    for ii=1:nEPS1
        EPS1(relaxind(ii),ii)=1;
    end;    
    
    AA	 = [AA,-EPS1];
    bl	 = [bl; repmat(   0,nEPS1,1)];    
    bu	 = [bu; repmat(bigB,nEPS1,1)];
    
    if length(weps)==1,
        weps=repmat(weps,nEPS1,1);
    elseif length(weps)~= nEPS1,
        error('length of weps and relaxind do not match');
    end
    
    % add one component in each matrix
    x0	 = [x0; zeros(nEPS1,1)];
    if Options.norm==2
        G	 = [G,zeros(nvar,nEPS1); zeros(nEPS1,nvar),diag(weps)];
        CC	 = [CC; zeros(nEPS1,1)];
    else % infinity and 1 norm have the same expression
        G	 = [G,zeros(nvar,nEPS1); zeros(nEPS1,nvar),zeros(nEPS1,nEPS1)];
        CC   = [CC; weps];
    end
    
    % optimization variable x now contains former x and slacks
end



% solve optimization problem
% ---------------------------------------------------------------------

ctype             = char('L'*ones(size(AA,1),1));
vartype           = char('C'*ones(size(G,1),1));
vartype(IntIndex) = 'B';

if isfield(Options, 'dumax')
    % include deltaU constraints
    % WARNING: these constraints will only be satisfied in open-loop solution if
    %          Options.Uprev is not given!!!
    dumin = Options.dumin;
    dumax = Options.dumax;
    ndu = length(dumax);
    [nAr, nAc] = size(AA);
    for ipr = 1:sum(horizon)-1,
        aa = [zeros(2*ndu, ndu*(ipr-1)) [-eye(ndu) eye(ndu); eye(ndu) -eye(ndu)] ...
                zeros(2*ndu, ndu*(sum(horizon)-ipr-1))];
        bb = [dumax; -dumin];
        AA = [AA; aa zeros(2*ndu, nAc-size(aa,2))];
        B = [B; bb];
        nlin = nlin+2*ndu;
    end
    if isfield(Options, 'Uprev'),
        % include deltaU constraints using previous input, such that we
        % guarantee satisfaction of deltaU constraints in closed loop
        if length(Options.Uprev) ~= length(dumax),
            error('Wrong dimension of Options.Uprev');
        end
        aa = [eye(ndu) zeros(ndu, nAc - ndu)];
        AA = [AA; aa; -aa];
        B = [B; dumax+Options.Uprev; -(dumin + Options.Uprev)];
        nlin = nlin + 2*ndu;
    else
        fprintf('Warning: No "Options.Uprev" specified, deltaU constraints will only be satisfied in open-loop solution!\n');
    end
end

MIoptions = [];
% use initial guess of xmin if provided
if isfield(Options, 'usex0'),
    MIoptions.usex0 = Options.usex0;
end
% propagate Options.save_prob to MIoptions
if isfield(Options, 'save_prob'),
    MIoptions.save_prob = Options.save_prob;
end
if isfield(Options, 'save_only'),
    MIoptions.save_only = Options.save_only;
end
if isfield(Options, 'logfile'),
    MIoptions.logfile = Options.logfile;
end

if isfield(Options, 'nocost'),
    % set objective function to zero if this flag is set
    % this case is used in reachability analysis where we are only looking for a
    % feasible solution
    if Options.nocost==1,
        if Options.norm==2,
            G = zeros(size(G));
        end
        CC = zeros(size(CC));
    end
end

if Options.removeInfBounds,
    % remove constraints of the form z <= Inf
    infB = find(B==Inf);
    if ~isempty(infB),
        AA(infB, :) = [];
        B(infB) = [];
    end

    % replace all lower bounds equal to -Inf by -Options.bigB
    if ~isempty(bl),
        infbl = find(bl==-Inf);
        if ~isempty(infbl),
            bl(infbl) = -bigB;
        end
    end
    
    % replace all upper bounds equal to Inf by Options.bigB
    if ~isempty(bu),
        infbu = find(bu==Inf);
        if ~isempty(infbu),
            bu(infbu) = bigB;
        end
    end
end

if Options.bounds2ineq
    % convert bounds to constraints
    if ~isempty(bl) & ~isempty(bu),
        nvar = size(AA, 2);
        AA = [AA; eye(nvar); -eye(nvar)];
        B = [B; bu; -bl];
        bu = [];
        bl = [];
        nlin = nlin + 2*nvar;
    end
end

startt = clock;
if Options.norm==2,
    [xopt, fopt, Eflagm, flag] = mpt_solveMIQP(G, CC, AA, B+epsil*ones(nlin,1), AAeq, BBeq, ...
        bl, bu, vartype, [], MIoptions);
else
    [xopt, fopt, Eflagm, flag] = mpt_solveMILP(CC, AA, B+epsil*ones(nlin,1), AAeq, BBeq, ...
        bl, bu, vartype, [], MIoptions);
end
runtime = etime(clock, startt);

full_xopt = xopt;
if ~flagihard
    %error('this part of the code needs to be rewritten')
    slkeps = xopt(nvar+1:nvar+nEPS1);   % slack variable for soft-constraints
    xopt = xopt(1:nvar);                % cut slacks off
else
    slkeps = 0;            % no relaxation has been performed
end            

% fix binary variables to 1/0
for ii=1:length(IntIndex)
    if xopt(IntIndex(ii))<0.5
        xopt(IntIndex(ii))=0;
    else
        xopt(IntIndex(ii))=1;
    end
end

% display
% ---------------------------------------------------------------------

if Options.verbose
    if flag == 1,
        if Options.norm==2
            disp(['MIQP solved']);
        else
            disp(['MILP solved']);
        end
    else
        if Options.norm==2
            disp(['flag of last MIQP = ', num2str(flag)]);
            error('MIQP solution FAILED');
        else
            disp(['flag of last MILP = ', num2str(flag)]);
            error('MILP solution FAILED');
        end
    end
end


% build outputs
% ---------------------------------------------------------------------

% row i contains i-th vector at all times
% column j contains vector at time j

uaux   = xopt(1:nu*NT);   
u      = reshape(uaux,nu,NT);

daux   = xopt(nu*NT+1:(nu+nd)*NT);
d      = reshape(daux,nd,NT);

zaux   = xopt((nu+nd)*NT+1:(nu+nd+nz)*NT);
z      = reshape(zaux,nz,NT);

ut     = u(:,1);
dt     = d(:,1);
zt     = z(:,1); 

Eflag            = []; 
Eflag.slkeps     = slkeps;
Eflag.xopt		 = xopt(:);
Eflag.fopt       = fopt;
Eflag.solverflag = Eflagm;  
Eflag.fl         = flag;
Eflag.u          = u;
Eflag.d          = d;
Eflag.z          = z;
Eflag.OL         = Ext.OL;
Eflag.full_xopt  = full_xopt;
Eflag.runtime    = runtime;
