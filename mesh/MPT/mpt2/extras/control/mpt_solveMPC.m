function [U,feasible,fval,sysStruct,probStruct,Matrices,X]=mpt_solveMPC(x0,sysStruct,probStruct,Matrices,Options)
%MPT_SOLVEMPC Solves the on line optimization MPC problem
%
% [U,feasible,fval,sysStruct,probStruct,Matrices]=mpt_solveMPC(x0,sysStruct,probStruct,Matrices,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%    
% Solves the optimization problem associated to the initial state x0 via LP, QP or LMI
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% x0            -   current state / initial state x(0)
% sysStruct     -   system Structure (type "help mpt_sysStruct")
% probStruct    -   problem Structure (type "help mpt_probStruct")
% Matrices      -   Optional: constraint matrices, if already known
% Options       -   Options for mpt_solveLP / mpt_solveQP
%   .solveLMI   -   1/0 if set to "1", the SOCP solver is used to solve the 2-norm
%                   problem; 
%                   NOTE: first call to SOCP solver takes significantly longer than
%                         subsequent calls.
%
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% U             -   optimizer, i.e. optimal input sequence
% feasible      -   1/0 flag; set to 1 if problem is feasible
% fval          -   Value function associated to optimizer
% sysStruct     -   updated (verified) sysStruct 
% probStruct    -   updated (verified) probStruct
% Matrices      -   constraint and objective Matrices; can be used as input for the 
%                   next call to this function
%
% see also mpt_mplp, mpt_mpqp for explicit feedback solutions 

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
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

global mptOptions;

error(nargchk(3,5,nargin));

if ~isa(x0, 'double'),
    error('mpt_solveMPC: first argument must be a vector of initial states');
end

if ~isstruct(mptOptions),
    mpt_error;
end
if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end
if ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end
if(nargin<4 | isempty(Matrices))
    Matrices=[];
end
if(nargin<5)
    Options=[];
end
if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end
if(isempty(Matrices))
    Matrices=mpt_constructMatrices(sysStruct,probStruct,Options);
end
if(iscell(sysStruct.A))
    error('MPT_SOLVEMPC: This function only works for LTI systems!')
end

x0=x0(:);

G=Matrices.G;
W=Matrices.W;
E=Matrices.E;
H=Matrices.H;
F=Matrices.F;
Y=Matrices.Y;

nuH=(size(sysStruct.B,2)*probStruct.N);
nx=size(x0,1);

[valid,nu,nbool,ubool,uint] = validboolinput(sysStruct);

if nbool>0,
    if ~valid
        error('mpt_solveMPC: all inputs must be either boolean (0/1) or integer with no "gaps"');
    else
        % handle systems with boolean and/or integer inputs
        bvartype = repmat('C',nu,1);
        LB = -1e9*ones(size(G,2),1);
        UB = 1e9*ones(size(G,2),1);
        for ii=1:length(uint),
            if uint{ii}.min==0 & uint{ii}.max==1,
                % this input is boolean
                bvartype(uint{ii}.uindex) = 'B';
            else
                % otherwise the input is integer
                bvartype(uint{ii}.uindex) = 'I';
            end
            for in=1:probStruct.N,
                % add proper bounds for boolean/integer variables
                LB((in-1)*nu+uint{ii}.uindex) = uint{ii}.min;
                UB((in-1)*nu+uint{ii}.uindex) = uint{ii}.max;
            end
        end
        vartype = [repmat(bvartype, probStruct.N, 1); ...
                repmat('C', size(G,2)-probStruct.N*nu, 1)]; 
        
        if probStruct.norm==2,
            [U,fmin,how,exitflag]=mpt_solveMIQP(H,x0'*F,G,W+E*x0, [], [], LB, UB, vartype);
            fval=0.5*U'*H*U+x0'*F*U+x0'*Y*x0;
        else
            [U,fmin,how,exitflag]=mpt_solveMILP(H, G, W+E*x0, [], [], LB, UB, vartype);
            fval=H*U;
        end
        umin=U(1:nu*probStruct.N);
        umin=reshape(umin,nu,probStruct.N)';            
        U = umin;
        Umat = U;
    end
else
    
    if(probStruct.norm==2 & isfield(Options,'solveLMI') & Options.solveLMI==1)
        if ~isfield(Options,'sdpoptions'),
            Options.sdpoptions = mptOptions.sdpsettings;
            %Options.sdpoptions = sdpsettings('Verbose',Options.verbose,'cachesolvers',1,'solver','sedumi');
        end
        H=H/2; %adjust construct matrices to LMI setup
        F=F/2; %adjust construct matrices to LMI setup
        
        alpha=sdpvar(1,1);      %define variables
        U=sdpvar(nuH,1);
        if 0,  %use LMIs
            myprog=LMI;             %build matrices       
            myprog=myprog+set(W+E*x0-G*U>0);    %constraint matrices
            %objective matrices
            myprog=myprog+set([alpha        (H*U+F'*x0)'    x0';...
                    (H*U+F'*x0)        H          zeros(nuH,nx);     ...
                    x0           zeros(nx,nuH)   inv((Y-F*inv(H)*F'))]>0); 
                
            else  %use cones... that's faster              
                myprog=set;             %build matrices       
                myprog=myprog+set(W+E*x0-G*U>0);    %constraint matrices
                
                Qhalf  = chol(inv(H));
                Rhalf  = chol(Y-F*inv(H)*F');
                
                myprog = myprog + set(cone([Qhalf*(H*U+F'*x0);Rhalf*x0],alpha));
            end
            
            solution = solvesdp(myprog,alpha,Options.sdpoptions);   %find solution 
            
            exitflag=~solution.problem; %extract solution 
            U=double(U);
            fval=double(alpha);
            
        elseif(probStruct.norm==2)
            [U,lambda,how,exitflag,fval]=mpt_solveQP(H,x0'*F+Matrices.Cf,G,W+E*x0);
            fval=0.5*U'*H*U+x0'*F*U+x0'*Y*x0+Matrices.Cf*U+Matrices.Cc+Matrices.Cx*x0; %+Matrices.Cx*x0+Matrices.Cc;
        else
            [U,fval,lambda,exitflag,how]=mpt_solveLP(H,G,W+E*x0);
            fval=H*U;
            U=U(1:nuH);
            if strcmpi(how, 'ok'),
                exitflag = 1;
            else
                exitflag = 0;
            end
        end
        Umat = reshape(U, nu, size(U,1)/nu)';
    end
    
if(exitflag<=0)
    feasible=0;
    U=[];
    disp(sprintf('no feasible control law for state x = [%s]',num2str(x0')));
else
    feasible=1;
end
    
X = x0';
if nargout > 6,
    if ~iscell(sysStruct.A),
        for ii=1:size(Umat,1),
            x0 = mpt_simSys(sysStruct, x0, Umat(ii,:)');
            x0 = x0';
            X = [X; x0'];
        end
    end
end


%------------------------------------------------------------------------
function [valid,nu,nbool,ubool,uint] = validboolinput(sysStruct)
% checks if all inputs specified in sysStruct.Uset are purely boolean
% i.e. no finite alphabet involving doubles or integer inputs

[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);
if ~isfield(sysStruct,'Uset')
    valid = 0;
    nbool = [];
    ubool = [];
    uint = {};
    return
end

uint = {};
for ii=1:length(ubool)
    ib = ubool(ii);
    Uset = sysStruct.Uset{ib};
    if any(Uset~=ceil(Uset)) | any(isinf(Uset))
        % first rule out non-integer inputs
        valid = 0;
        return
    else
        minU = min(Uset);
        maxU = max(Uset);
        if ~isempty(setdiff(minU:maxU,Uset))
            % rule out "gaps" in integers, e.g. [-2 -1 1 3] -> error
            valid = 0;
            return
        else
            uints.min = minU;
            uints.max = maxU;
            uints.uindex = ib;
            uint{ii} = uints;
        end
    end
end
valid = 1;