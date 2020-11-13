function [solution]=mpt_getCommonSOSLyapFct(ctrl,ndeg,Options)
%MPT_GETCOMMONSOSLYAPFCT Calculates Common SOS Lyapunov function for system with additive disturbance 
%
% [solution]=mpt_getCommonSOSLyapFct(ctrlStruct,ndeg,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function attempts to compute a higher order Sum of Squares Lyapunov
% function
% V(x) which guarantees exponential stability. The following is satisfied
% alpha * x^2 <= V(x(k))<= beta * x^2  
% V(x(k+1)) - V(x(k))<= - rho * x^2      
% (alpha,beta,rho>=0)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% ctrl                  - Explicit controller object
% ndeg                  - Degree of Lyapunov function desired (must be even number >=2)
% Options.sossolver     - SOS solver to be used , 0 : Yalmip (default) or 1: SOSTools
% Options.abs_tol       - Absolute tolerance
% Options.epsilon       - This is a tolerance factor which is introduced to turn
%                         LMI inequalities into strict inequalities.
% Options.debug_level   - If this is set to 1, the solution provided by the LMI
%                         solver will be double-checked manually. We strongly
%                         advise to set this to 1, since we've experienced
%                         numerous numerical issues with certain LMI solvers.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% V          - SOS Lyapunov function: Overall polynomial on partition
%                
% feasible   - 1: stable or 0: no statement about stability possible
% rho        - the Lyapunov decay rate over the partition
%              (is this is less than zero, stability cannot be guaranteed)
%              The Lyapunov value decrease Delta V <= - rho * ||x||^2
%
% see also MPT_GETPWQLYAPFCT,MPT_GETPWQLYAPFCT, MPT_GETQUADLYAPFCT, MPT_GETPWSOSLYAPFCT,MPT_GETSTABFEEDBACK
%

% ---------------------------------------------------------------------------
%   LITERATURE:
% ---------------------------------------------------------------------------
%
% Proceedings of the American Control Conference (ACC), Denver, CO. 2003.
% "Analysis of Switched and Hybrid Systems - Beyond Piecewise Quadratic
% Methods",
% S. Prajna, A. Papachristodoulou. 

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%          baric@control.ee.ethz.ch

%
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

error(nargchk(2,3,nargin));

if(mod(ndeg,2)~=0 | ndeg<2)
    error('mpt_getCommonSOSLyapFct: Order must be even number greater than or equal 2');
end;

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if (nargin<3),
    Options = [];
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
    	error('This function supports only explicit controllers!');
	end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end

% TOLERANCES
%
MAX_ZERO_TOLERANCE = 1e-7;
MIN_DECAY_RATE     = 1e-5;

% Options
%-------------
%
% SOS solver
YALMIP   = 0;
SOSTOOLS = 1;

localOptions = Options;

if ~isfield(localOptions,'debug_level'),
    localOptions.debug_level = 0;
end

if ~isfield(Options,'epsilon'),
    localOptions.epsilon=mptOptions.abs_tol;   %epsilon is added to guarantee that a quadratic lower bound on PWQ exists
end

if ~isfield(localOptions,'sossolver'),
    % YALMIP used by default
    %
    localOptions.sossolver = YALMIP;
end

if ~isfield(localOptions,'abs_tol'),
    % zero tolerance
    %
    localOptions.abs_tol = min(mptOptions.abs_tol,MAX_ZERO_TOLERANCE);
else
    localOptions.abs_tol = min(localOptions.abs_tol,MAX_ZERO_TOLERANCE);
end

if ~isfield(localOptions,'checkmulticonvexity'),
    %
    % checking if the Lyapunov function is 'multiconvex', by default,
    % don't do the cheking
    %
    localOptions.checkmulticonvexity = 0;
end

if ~isfield(localOptions,'sdpsolver')
    % this field is ignored in case SOSTools is selected as sossolver
    %
    localOptions.sdpsolver = 'sedumi';
else
    lowSolver = lower(localOptions.sdpsolver);
    if ( strcmp(lowSolver,'sedumi') ~= 1 & ...
            strcmp(lowSolver,'sdpt3')  ~= 1 ),
        warning ('Unsupported solver, switching to SeDuMi!');
        localOptions.sdpsolver = 'sedumi';
    end
end

if ( isfield(localOptions,'SProcForNoise') ),
    % if set, this flag enforces using S-procedure for the
    % noise. Otherwise, stability conditions are checked for
    % extreme realisations of the noise, enforcing the convexity of
    % the Lyapunov function as a necessity.
    %
    SProcForNoise = localOptions.SProcForNoise;
else
    SProcForNoise = 1;
    localOptions.SProcForNoise = 1;
end

% initial output values
%
timing   = struct('setupTime',0,...
    'compTime',0);

solution = struct('found'       ,0,...
    'V'           ,[],...
    'details'     ,struct('message','','timing',timing));

% input data
%
sysStruct = ctrlStruct.sysStruct;
verOptions.verbose=0;
if ~isfield(sysStruct,'verified')
    sysStruct = mpt_verifySysStruct(sysStruct, verOptions);
end

Pn        = ctrlStruct.Pn;
Fi        = ctrlStruct.Fi;
Gi        = ctrlStruct.Gi;
Pfinal    = ctrlStruct.Pfinal;
[nu,nx]   = size(Fi{1});

setupStartTime = cputime;

% get system information and extract relevant data
%
sysInfo   = subfun_getSysInfo(ctrlStruct,localOptions);
isNoisy   = sysInfo.isNoisy;
idxActDyn = sysInfo.idxActDyn;
Hnoise    = sysInfo.Hnoise;
Knoise    = sysInfo.Knoise;
decayPn   = sysInfo.decayPn;

%%%%%%%%%%%%
% SOSTools %
%%%%%%%%%%%%
%
if ( localOptions.sossolver == SOSTOOLS ),
    % define variables
    %
    varStrx = [];
    varStrw = [];
    for ii = 1:nx,
        auxStr = sprintf('x_%d ',ii);
        varStrx = [varStrx auxStr];
        if ( sysInfo.isNoisy & SProcForNoise ),
            auxStr = sprintf('w_%d ',ii);
            varStrw = [varStrw auxStr]; 
        end
    end
    varStr = [varStrx varStrw];
    eval(['syms ' varStr ' real']);
    eval(['vars = [' varStr '];']);
    vars = vars';
    
    try
        lyap = sosprogram(vars);   % Setup the SOS Program
    catch
        error(['mpt_getCommonSOSLyapFct: mpt_getPWSOSLyapFun: You need to download and  install' ...
                ' SOSTools for this function to work!']);
    end
    
    % create vector of monomials, starting from 2 (no need for
    % the upper bound constraint for the Lyapunov function)
    %
    xMon = monomials (vars(1:nx),[2:ndeg]);
    
    % variables determining decay rate and lower bound
    %
    [lyap, ALPHA] = sossosvar(lyap,1);
    [lyap, RHO]   = sossosvar(lyap,1);
    
    % define polynomial Lyapunov function
    %
    [lyap,Vk] = sospolyvar (lyap, xMon);
    
    % set negative decay of the Lyapunov function on a subset of the
    % feasible region (decayPn)
    %--------------------------------------------------
    %
    for region = 1:length(decayPn),
        [HH,KK] = double(sysInfo.decayPn(region));
        rowRange = [1:nx] + (idxActDyn(region)-1) * nx;
        ABF = sysInfo.bigABF(rowRange,:);
        BG  = sysInfo.bigBG(rowRange,:);
        if ( isNoisy & SProcForNoise ),
            % extend the state vector with the disturbance
            %
            ABF = [ABF eye(nx)];
            HH  = [HH zeros(size(HH,1),nx); ...
                    zeros(size(Hnoise,1),nx) Hnoise];
            KK  = [KK;Knoise];
        end
        %
        if ( isNoisy & ~SProcForNoise ),
            %
            % TODO
            %
            error (['mpt_getCommonSOSLyapFct: Robust stability without SProcedure for noise not yet' ...
                    ' implemented!']);
        else
            xk   = num2cell(vars(1:nx));
            xkpp = num2cell(ABF * vars + BG);
            Vkpp = subs (Vk,xk,xkpp);
            deltaV{region} = Vk - Vkpp - RHO * vars' * vars;%JOHAN
            if ( any(KK<0) | isNoisy ),
                %
                % S-procedure for extended (x,w) polytopes
                %
                dimN   = size(HH,1);
                Gdecay = KK - ( HH * vars );            
                for ii = 1:dimN,
                    for jj = ii+1:dimN,
                        [lyap, N{region}(ii,jj)] = sossosvar(lyap,1);
                        deltaV{region} = deltaV{region} - N{region}(ii,jj) ...
                            * Gdecay(ii) * Gdecay(jj);
                    end
                end
            end
        end    
        lyap = sosineq(lyap,deltaV{region});
    end
    
    % positivity of the Lyapunov function
    %---------------------------------------
    %
    for finalReg = 1:length(Pfinal),
        nConstrPfinal = nconstr(Pfinal(finalReg));
        [Hfin, Kfin] = double (Pfinal(finalReg));
        Vpos{finalReg} = Vk - ALPHA * vars(1:nx)' * vars(1:nx);
        if ( any(Kfin<0) ),
            Gpos = Kfin - ( Hfin * vars(1:nx) );
            for ii = 1:nConstrPfinal,
                for jj = ii+1:nConstrPfinal,
                    [lyap, M{finalReg}(ii,jj)] = sossosvar (lyap,1);
                    Vpos{finalReg} = Vpos{finalReg} - M{finalReg}(ii,jj) ...
                        * Gpos(ii) * Gpos(jj);
                end
            end
        end
        lyap = sosineq(lyap,Vpos{finalReg});
    end
    sossetobj(lyap,-RHO);
    
    setupTime = cputime - setupStartTime;
    
    % solve SOS problem
    %--------------------
    %
    disp ('mpt_getCommonSOSLyapFct: Solving SOS problem ...');
    
    compStartTime = cputime;
    [lyap,info] = sossolve(lyap,1,1e-12); 
    compTime = cputime - compStartTime;
    
    if ( (info.dinf == 1) | (info.pinf == 1) ),
        %
        % infeasible problem
        %
        disp ('mpt_getCommonSOSLyapFct: No common Lyapunov function has been found!');
        return;
    end
    
    % get whatever we have here and check if it's SOS
    %
    lyapFun  =  sosgetsol(lyap,Vk); 
    dALPHA  =  double(double(sosgetsol(lyap,ALPHA)));
    dRHO    =  double(double(sosgetsol(lyap,RHO)));
    [chkQ,chkZ] = findsos(lyapFun);
    if ( dALPHA < localOptions.abs_tol | dRHO < localOptions.abs_tol | ...
            isempty(chkQ) ),
        % we're not happy with the solution we have
        %
        if ( localOptions.debug_level ),
            warning (['mpt_getCommonSOSLyapFct: The obtained solution may not be valid. Try again using ' ...
                    'different parameters and/or solver.']);
        end
        return;
    end
    
    % fill up the output structure
    %
    solution.found    = 1;
    solution.V        = lyapFun;     % Lyap. function
    solution.details.rho   = dRHO;   %decay rate
    solution.details.alpha = dALPHA; % lower bound
    solution.details.timing = timing;
    solution.timing.setupTime = setupTime;
    solution.timing.compTime  = compTime;
    
    disp('mpt_getCommonSOSLyapFct: Commmon lyapunov function found');
    
else
    %%%%%%%%%%
    % YALMIP %
    %%%%%%%%%%
    %
    yalmip('clear');
    constrSet = set([]);    
    parametricVars = [];
    %
    strVecX = 'vecX = [';
    strVecW = 'vecW = [';
    %
    % declare all variables
    %
    for ii = 1:nx,
        varStr = sprintf('x%d',ii);
        eval([varStr ' = sdpvar(1,1);']);
        strVecX = [strVecX varStr ' '];
        if ( isNoisy & SProcForNoise ),
            varStrw = sprintf('w%d',ii);
            eval([varStrw ' = sdpvar(1,1);']);
            strVecW = [strVecW varStrw ' '];
        end
    end
    eval([strVecX ']'';']);
    eval([strVecW ']'';']);
    
    % variables defining decay rate and lower bound
    %
    ALPHA = sdpvar(1,1);
    RHO   = sdpvar(1,1);
    parametricVars = [parametricVars; ALPHA; RHO];
    
    % monomials: start from quadratic exponents (this way we don't
    % need upper bound constraints of the Lyapunov function)
    %
    vecMonoX = monolist(vecX,ndeg);
    
    % Skip constant and linear terms
    %
    vecMonoX = vecMonoX(2+length(vecX):length(vecMonoX));
    
    % define Grammian matrix and SOS Lyapunov function candidate
    %
    QQ = sdpvar(length(vecMonoX),1);
    V  = vecMonoX' * QQ;
    parametricVars = [parametricVars; QQ(:)];
    
    % set of constraints
    %
    %constrSet = constrSet + set('10 > ALPHA > 0') + set('10 > RHO > 0');
    constrSet = constrSet + set(ALPHA > 0) + set(10> RHO > 0);
    %                             
    for region = 1:length(decayPn),
        [HH,KK] = double(decayPn(region));
        rowRange = [1:nx] + (idxActDyn(region)-1) * nx;
        ABF = sysInfo.bigABF(rowRange,:);
        BG  = sysInfo.bigBG(rowRange,:);
        if ( isNoisy )
            if ( SProcForNoise ),
                % extend the state vector with the disturbance
                %
                ABF = [ABF eye(nx)];
                HH  = [HH zeros(size(HH,1),nx); ...
                        zeros(size(Hnoise,1),nx) Hnoise];
                KK  = [KK;Knoise];
            else
                % we'll introduce convexity constraints here
                %
            end
        end
        varsXW = [vecX;vecW];
        vecXpp = ABF * varsXW + BG;
        vecMonoXpp = monolist(vecXpp,ndeg);
        vecMonoXpp = vecMonoXpp(2+length(vecXpp):length(vecMonoXpp));
        vecMonoXpp = clean(vecMonoXpp,localOptions.abs_tol);            
        dimN  = size(HH,1);
        if ( any(KK < 0) | isNoisy ),
            % the region doesn't containt the origin or we're doing
            % robust stability analysis, add S-procedure
            % term
            %
            auxS = sdpvar(dimN,dimN);
            auxS = auxS - diag(diag(auxS));
            theSTerm = (KK - HH*varsXW)' * auxS * (KK - HH*varsXW);
            parametricVars = [parametricVars; auxS(:)];
            %constrSet = constrSet + set('norm(KK) > recover(depends(auxS)) > 0');
            constrSet = constrSet + set(recover(depends(auxS)) > 0);
        else
            %
            % check if the offset term BG is non-zero - if so, the origin
            % is not stable
            %
            if ( norm(BG,'inf') > localOptions.abs_tol ),
                solution.details.message = sprintf(['The origin is not the ' ...
                        'critical point.']);
                return;
            end    
            % remove constant monomials
            %
            vecMonoXpp = vecMonoXpp - full(getbasematrix(vecMonoXpp,0));
            theSTerm = 0;
        end
        deltaX = clean(vecMonoX - vecMonoXpp,localOptions.abs_tol);
        deltaV = deltaX' * QQ - theSTerm - RHO*vecX'*vecX;
        constrSet = constrSet + set(sos(deltaV));
    end
    
    % set the positivity of the Lyapunov function over the feasible
    % region
    %
    for finalReg = 1:length(Pfinal),
        nConstrPfinal = nconstr(Pfinal(finalReg));
        [Hfin, Kfin] = double (Pfinal(finalReg));
        Vpos = V - ALPHA*vecX'*vecX;        
        if ( any(Kfin < 0) ),
            % the region Pfinal(finalReg) doesn't contain the origin and
            % we'll need the S-procedure
            %
            auxS = sdpvar(nConstrPfinal,nConstrPfinal);
            auxS = auxS - diag(diag(auxS));
            Gpos = Kfin - ( Hfin * vecX );
            Vpos = Vpos - Gpos' * auxS * Gpos;
            parametricVars = [parametricVars; auxS(:)];
            constrSet = constrSet + set(norm(Kfin) > recover(depends(auxS)) > 0); 
        end
        constrSet = constrSet + set(sos(Vpos));
    end
    setupTime = cputime - setupStartTime;
    
    % set options
    %
    sosOptions = sdpsettings;
    sosOptions.clean   = localOptions.abs_tol;
    sosOptions.solver  = localOptions.sdpsolver;
    sosOptions.verbose = localOptions.debug_level;
    sosOptions.cachesolvers   = 1;
    sosOptions.sedumi.eps     = 1e-9;
    sosOptions.sos.extlp      = 1;
    sosOptions.sos.congruence = 0;
    sosOptions.sos.postprocess = 0;
    
    
    % solve SOS problem
    %--------------------
    %
    if ( localOptions.debug_level ),
        disp (['mpt_getCommonSOSLyapFct: Solving SOS problem using ' sosOptions.solver]);
    end
    
    compStartTime = cputime;
    
    % we maximize parameter RHO to get Lyapunov function of 'maximal
    % decrease rate' (doesn't change anything, we could also define a
    % feasibility problem as well)
    %   
    [sol,outMonomials,Grammian,primalSlacks] = solvesos(constrSet,-RHO,sosOptions,parametricVars);
    solution.details.message = yalmiperror(sol.problem);
    
    compTime = cputime - compStartTime;
    
    postprocStartTime = cputime;
    %---------------------------------------------------
    keepgoing=1;
    if 1%(solution.problem==4) | (solution.problem==3) | (solution.problem==0)
        mineig = zeros(length(Grammian),1);
        for i = 1:length(Grammian)
            mineig(i,1) = min(eig(Grammian{i}))/length(Grammian{i});
        end             
        if all(mineig>primalSlacks) & all(checkset(constrSet(find(~is(constrSet,'sos'))))>=0)
            disp('mpt_getCommonSOSLyapFct: certifiable SOS decomposition => Lyapunov function found!');            
            solution.found = 1;
            theFinalSolution = vecMonoX' * double(QQ);   
            strSolution = sdisplay(theFinalSolution);
            keepgoing=0;
        end
    end
    
    
    % analyse the output
    %
    if (keepgoing & (sol.problem == 4 | sol.problem == 0)),
        % we may have the solution
        %
        % verify the solution
        %
        solution.found = 1; %set to zero in case of failure
        
        elmwiseIdx  = find(is(constrSet,'elementwise'));
        sosIdx      = find(is(constrSet,'sos'));
        maxDecayIdx = length(decayPn);
        [primalSlacks,dualSlacks] = checkset(constrSet);
        
        % 0) check values of ALPHA and RHO
        %
        dALPHA = double(ALPHA);
        dRHO   = double(RHO);
        if ( dALPHA < localOptions.abs_tol | ...
                dRHO   < MIN_DECAY_RATE ),
            solution.details.message = ...
                sprintf('mpt_getCommonSOSLyapFct: No valid solution has been found.');
             solution.found = 0;
             keepgoing=0;
        end
        
        % 1) check elementwise constraints
        %
        if (keepgoing & (any(primalSlacks(elmwiseIdx) < -localOptions.abs_tol))),
            % the solution is crap
            %
            minSlack = min(primalSlacks(elmwiseIdx));
            solution.details.message = sprintf(['mpt_getCommonSOSLyapFct: Violated positivity ' ...
                    'constraint, residual: %E.\n'],minSlack);
            solution.found = 0;
            keepgoing=0;
        end
        
        % 2) check the eigenvalues of the SOS decomposition
        %
        sosOptions.verbose = 0;
        for ii = 1:length(sosIdx),
            epsMax = primalSlacks(sosIdx(ii));
            minEigValue = min(eig(Grammian{ii}));
            lenMonomials = length(outMonomials{ii});
            if ( minEigValue < -localOptions.abs_tol ),
                solution.details.message = ...
                    sprintf('mpt_getCommonSOSLyapFct: No valid solution has been found.');
                solution.found = 0;
                keepgoing=0;
            end
            if (keepgoing & (minEigValue < (lenMonomials * epsMax)) ),
                % possibly a bad solution
                %
                %
                % formulate SOS problem for the particular constraint
                %
                v = sosd(constrSet(sosIdx(ii)));
                auxPoly = v'*v;
                if ( ii <= maxDecayIdx ),
                    % decay constraint
                    %
                    auxPoly = auxPoly + (dRHO-MIN_DECAY_RATE) * vecX'*vecX;
                else
                    % positivity constraint
                    %
                    auxPoly = auxPoly + (dALPHA-MAX_ZERO_TOLERANCE) * vecX'*vecX;
                end
                
                auxSOS = set([]);
                auxSOS = set(sos(auxPoly));              
                [auxsol,auxmono,auxgramm] = solvesos(auxSOS,[], ...
                    sosOptions,[]);
                lenMonomials = length(auxmono{1});
                [auxEpsMax,dummy] = checkset(auxSOS);
                if ( min(eig(auxgramm{1})) < (lenMonomials * auxEpsMax)),
                    solution.details.message = sprintf(['mpt_getCommonSOSLyapFct: Violated SOS constraint: ' ...
                            'residual=%E, min. Grammian eigenvalue=%E.\n'], ...
                        epsMax,minEigValue);
                    solution.found = 0;
                end
            end
            %
        end
        theFinalSolution = vecMonoX' * double(QQ);   
        strSolution = sdisplay(theFinalSolution);

        % check multiconvexity
        %
        if ( localOptions.checkmulticonvexity ),
            solution.multiconvex = subfun_isMultiConvex(theFinalSolution,localOptions);
        end
    elseif(keepgoing),
        % the problem is not solved due to some other error
        % (infeasibility, solver failure etc.)
        %
        if ( localOptions.debug_level ),
            fprintf (2,'YALMIP ERROR: %s\n',yalmiperror(sol.problem));
        end
        strSolution = 'error';
    end
    
 
    
    solution.V             = strSolution;%subfun_getSym(strSolution{1});
    solution.details.rho   = double(RHO);   % decay rate
    solution.details.alpha = double(ALPHA); % lower bound
    
    
    solution.details.timing.postprocTime = cputime-postprocStartTime;
    solution.details.timing.setupTime = setupTime;
    solution.details.timing.compTime  = compTime;
    
    
end







if(solution.found == 1)
    disp(['mpt_getCommonSOSyapFun: SUCCESS:  Found Common Polynomial Lyapunov function of degree ' num2str(ndeg)])
else
    disp(['mpt_getCommonSOSyapFun: FAILURE:  Failed to find Common Polynomial Lyapunov function of degree ' num2str(ndeg)])
end




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the main subroutine %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

























function [symbSol] = subfun_getSym(strPoly)
%
% converts YALMIP's sdpvar polynomial expression into symbolic expression:
% again, a hack. It'll be supported by YALMIP in future releases.
%
% trim all the spaces
%
spaceIdx = strfind(strPoly,' '); 
strPoly(spaceIdx) = [];

% get occurences of 'x' variables
%
[starts,fs,tokens] = regexp(strPoly,'(x\d+)');
symVars = [];
newStrPoly = strPoly;
endStrIdx = length(strPoly);

% collect all x-es
%
for ii = 1:length(tokens),
    symVars = [symVars strPoly(tokens{ii}) ' '];
end

% add '*' for multiplications
%
newStrPoly = strrep(strPoly,'x','*x');
[auxs,auxf,auxt] = regexp(newStrPoly,'([++--]*)');
if ( ~isempty(auxs) ), % remove '-*' or '+*'
    newStrPoly(auxs+1) = [];
end

% declare symbolic variables (some may repeat, but it doesn't matter)
%
eval(['syms ' symVars ' symbSol;']);
eval(['symbSol = ' newStrPoly ';']);
%
%%%%%%%%%%%%%%%%%%%%%%%%
% end of subfun_getSym %
%%%%%%%%%%%%%%%%%%%%%%%%


function [isMConvex] = subfun_isMultiConvex(funLyap,Options)
%
% checks if the Lyapunov function is multiconvex
%

isMConvex = 0;        

if ( Options.sossolver ),
    %
    % SOStools
    %
    
    % not yet implemented
    return;
    
else
    %
    % YALMIP
    %
    diagPolys = diag(hessian(funLyap));
    auxSOSOpts = sosOptions;
    auxSOSOpts.verbose = 0;
    for ii  = 1:length(diagPolys),
        thePoly = diagPolys(ii);
        if ( isa(thePoly,'double') )
            % second derivative is constant
            %
            if ( double(thePoly) < -Options.abs_tol ),
                return;
            end
        else
            % we'll do SOS decomposition of the diagonal polynomials
            % of the Hessian
            %
            sosProb = set([]);
            sosProb = set(sos(thePoly));
            [auxsol,auxmonos,auxGramm] = solvesos(sosProb,[], ...
                auxSOSOpts,[]);
            [primSlacks,dualSlack] = checkset(sosProb);
            maxSlack = max(abs(primSlacks));
            if ( any(eig(auxGramm) < size(Gramm,1)*maxSlack ) ),
                return;
            end
        end
    end
    isMConvex = 1;
end       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  end of subfun_isMultiConvex %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sysInfo] = subfun_getSysInfo(ctrlStruct, Options)
% 
% transforms the system in a form convenient for the stability
% analysis
%
noiseVertices = [];
Hnoise  = []; 
Knoise  = [];
isNoisy = 0;
sysInfo = struct('isNoisy',      isNoisy,   ...
    'idxActDyn',    [],  ...
    'decayPn',      [],  ...
    'noiseVertices',noiseVertices,  ...
    'Hnoise',       Hnoise,  ...
    'Knoise',       Knoise,  ...
    'bigABF',       [],  ...
    'bigBG',        []);

sysStruct = ctrlStruct.sysStruct;
Pn        = ctrlStruct.Pn;
Fi        = ctrlStruct.Fi;
Gi        = ctrlStruct.Gi;
Pfinal    = ctrlStruct.Pfinal;

if isfield(ctrlStruct.probStruct, 'FBgain'),
    % handle pre-stabilization with feedback
    FBgain = ctrlStruct.probStruct.FBgain;
else
    FBgain = zeros(size(ctrlStruct.Fi{1}));
end

isPWASystem = 0;
if ( iscell(sysStruct.A) ),
    %
    % we have PWA system
    %
    if( isfield(sysStruct,'Aunc') & ~isempty(sysStruct.Aunc) )
        error('mpt_getCommonSOSLyapFct: Cannot handle PWA systems with polytopic uncertainty');
    end
    Acell = sysStruct.A;
    Bcell = sysStruct.B;
    isPWASystem = 1;
    if ( length(Acell) ~= length(Pn) ),
        noU = size(Bcell{1},2);
        cc=0;
        Acell={};
        Bcell={};
        for ii = 1:length(Pn),
            [xC,rC] = chebyball(Pn(ii));            
            Acell{ii}=[];
            for jj = 1:length(sysStruct.A),          
                guardX = sysStruct.guardX{jj};
                guardU = sysStruct.guardU{jj};
                guardC = sysStruct.guardC{jj};
                uu = (FBgain(1:noU,:) + Fi{ii}(1:noU,:))*xC + Gi{ii}(1:noU,:);
                if ( max(guardX*xC+guardU*uu - guardC) < Options.abs_tol ),    
                    % check which dynamics is active in the region  
                    %
                    Acell{ii}=sysStruct.A{jj};
                    Bcell{ii}=sysStruct.B{jj};
                    fcell{ii}=sysStruct.f{jj};
                end
            end
            if(isempty(Acell{ii}))
                error('mpt_getCommonSOSLyapFct: Faulty partition: Region could not be linked to any dynamic !!')
            end
        end
    else
        Acell = sysStruct.A;
        Bcell = sysStruct.B;
        fcell = sysStruct.f;
    end
    
elseif ~isfield(sysStruct,'Aunc') | ~isfield(sysStruct,'Bunc')
    % no polytopic uncertainty 
    % 
    Acell{1}=sysStruct.A;
    Bcell{1}=sysStruct.B;
else
    % polytopic uncertainty - check stability for all extreme
    % realizations of A and B
    %
    Acell=sysStruct.Aunc;
    Bcell=sysStruct.Bunc;
end

% relevant sizes
%
if ( isPWASystem ),
    nx = size(sysStruct.A{1},1);
    nu = size(sysStruct.B{1},2);
else
    nx = size(sysStruct.A,1);
    nu = size(sysStruct.B,2);
end
%
if ( isfield(sysStruct,'noise') ),
    Pw = sysStruct.noise;
    if ~isa(Pw, 'polytope'),
        error('Only polytopic noise is supported by this function.');
    end
    if ( isfulldim(Pw) ),
        isNoisy = 1;
        if ( ~(Options.SProcForNoise) ),
            noiseVertices = extreme(Pw);
        elseif ( isfulldim(Pw) )
            % we need these for the S-procedure
            %
            [Hnoise,Knoise] = double(Pw);
        end
    end
end    
%
nRegions = length(Pn);
if ( isPWASystem ),
    bigA = cat(1,Acell{:});
else
    bigA = repmat(Acell{1},nRegions,1);
end
bigABF = bigA;
bigBG  = zeros(size(bigA,1),1);
auxB   = Bcell{1};
auxf   = zeros(size(auxB,1),1);

idxContainsOrigin = [];
localOptions.abs_tol = Options.epsilon;
localOptions.fastbreak = 0;
for ii = 1:nRegions,
    rowRange = [1:nx] + (ii-1) * nx;
    if ( ii > 1 & isPWASystem ),
        auxB = Bcell{ii};
        auxf = fcell{ii};
    end
    bigABF(rowRange,:) = bigABF(rowRange,:) +  auxB * (Fi{ii}(1:nu,:) + FBgain(1:nu,:));
    bigBG(rowRange,:)  = auxB * Gi{ii}(1:nu) + auxf;
    
    % check if the origin is inside the closure of the current
    % region
    [isIn, dummy1, dummy2] = isinside (Pn,zeros(nx,1), ...
        localOptions);
    if ( isIn ),
        % if the region contains the origin, it has to be the
        % equilibrium point
        %
        if ( all(abs(bigBG(rowRange)) < Options.epsilon )),
            idxContainsOrigin(end+1) = ii;
        end
    end
end
decayPn   = Pn;
idxActDyn = 1:nRegions;
if ( isNoisy ),
    % additive noise
    %
    if ( ~isempty(idxContainsOrigin) ), 
        for ii = 1:length(idxContainsOrigin),
            rowRange = [1:nx] + (idxContainsOrigin(ii)-1) * nx;
            ABF{ii} = bigABF(rowRange,:);
            BG{ii}  = bigBG(rowRange,:);
        end
        
        % negative decay of the Lyapunov function should be
        % enforced outside maximal robust invariant set
        %
        [Xinfrob,dyn] = mpt_infsetPWA(Pn(idxContainsOrigin),ABF, ...
            BG,Pw);
        if ( isempty(Xinfrob) ),
            error (['mpt_getCommonSOSLyapFct: Cannot compute' ...
                    ' maximal robust invariant set!']);
        end
        Pnorig  = Pn(idxContainsOrigin) \ Xinfrob;
        decayPn(idxContainsOrigin) = [];
        %
        % slice up Pnorig
        %
        auxDyn = [];
        newPnorig = polytope;
        for ii = 1:length(Pnorig),
            for jj = idxContainsOrigin,
                auxIntersect = intersect(Pnorig(ii),Pn(jj));
                if ( ~isempty(auxIntersect) ),
                    newPnorig = [newPnorig auxIntersect];
                    auxDyn(end+1) = jj;
                end
            end
        end
        idxActDyn(idxContainsOrigin) = [];
        decayPn = [decayPn newPnorig];
        idxActDyn = [idxActDyn auxDyn];
    else
        % origin is not in the partition or it is not an
        % equilibrium point, quit
        %
        error (['mpt_getCommonSOSLyapFct: origin should be inside ' ...
                'the partiton.']);
    end
end
%
% set the values of output parameters
%
sysInfo.isNoisy       = isNoisy;
sysInfo.idxActDyn     = idxActDyn;
sysInfo.decayPn       = decayPn;
sysInfo.noiseVertices = noiseVertices;
sysInfo.Hnoise        = Hnoise;
sysInfo.Knoise        = Knoise;
sysInfo.bigABF        = bigABF;
sysInfo.bigBG         = bigBG;
