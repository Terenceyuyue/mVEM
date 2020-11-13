function [PolyOut,details] = mpt_greedyMerging(PolyIn, Options)
% MPT_GREEDYMERGING Greedy merging of polyhedra
%===============================================================================
%
% Title:        greedyMerging.m                                             	
%                                                                       
% Project:      merging of polyhedra
%                                                             
% Inputs:       PolyIn: set of polytopes either 
%                 as an array of polytopes using the MPT toolbox, or
%                 or as a structure using the halfspaces description Hi*x <= Ki
%               Options (optional):
%                 multiple: 0: one loop
%                           1: do multiple merging loops until no improvement (default)
%                 trials:   0: one trial (default)
%                          >0: do multiple trials. If a better solution is found,
%                              reset the counter to zero
%                 verbose:  0: verbosity off 
%                           1: verbosity on (default)
%                 union.lpsolver: LP solver used when computing the union
%                           0: NAG Toolbox (default)
%                           1: Matlab's linprog
%                 union.tolerance: tolerance when computing the union
%                           1e-6 (default)
%                 union.reduce: reduce union to minimal representation
%                           0: no (default)
%                           1: yes
%
% Outputs:      PolyOut: set of merged polytopes either
%                 as an array of polytopes using the MPT toolbox, or
%                 or as a structure using the halfspaces description Hi*x <= Ki.
%               The tool choses the same representation as for PolyIn
%
% Comments:     The algorithm tries to merge as many of the given polyhedra as 
%               possible using a greedy approach. The algorithm cycles through 
%               the regions and checks if any two regions form a convex union.
%               If so, the algorithm combines them in one region, and continues
%               checking the remaining regions.
%               To improve the solution, multiple merging loops are enabled
%               by default.
%               To reduce the problem of getting stuck in local minima, several 
%               trials can be used until the solution is not improved.
%
% Author:       Tobias Geyer
%                                                                       
% History:      date        subject                                       
%               2004.01.16  initial version based on mb_reducePWA by Mato Baotic
%               2004.05.06  debugged, comments added
%
% Requires:     multi parametric toolbox MPT,
%               mb_buildBClist
%               cell2mpt, mpt2cell,
%               mb_polyunion,
%               facetinnercircle
%
% Contact:      Tobias Geyer
%               Automatic Control Laboratory
%               ETH Zentrum, 
%               Zurich, Switzerland
%
%               geyer@control.ee.ethz.ch
%
%               Comments and bug reports are highly appreciated
%
%===============================================================================



global mptOptions

starttime = cputime;

% complement the options
if nargin<2
    Options=[];
end
if ~isfield(Options,'multiple')
    Options.multiple=1;
end
if ~isfield(Options,'trials')
    Options.trials = 1;
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options,'rel_tol')
    Options.rel_tol = mptOptions.rel_tol;
end

if ~isfield(Options, 'statusbar')
    Options.statusbar = 0;
end
if ~isfield(Options, 'status_min')
    Options.status_min = 0;
end
if ~isfield(Options, 'status_max')
    Options.status_max = 1;
end
if ~isfield(Options, 'closestatbar'),
    Options.closestatbar = 1;
end
if Options.statusbar,
    Options.verbose = -1;
end

Options.reduce = 0;

Options.emptypoly = mptOptions.emptypoly;

[Hn, Kn] = double(PolyIn);

Ri.nR = length(PolyIn);

if Ri.nR < 2,
    PolyOut = PolyIn;
    if nargout > 1,
        details.before = 1;
        details.after = 1;
        details.time = cputime - starttime;
        details.alg = 'greedy';
    end
    return
end

Ri.Pn = PolyIn;

details.before = Ri.nR;

% compute list of neighbours
Ri.BC = sub_buildBClist(Hn,Kn,PolyIn,Options);

% initialize the current best solution
Ri_best.nR = inf;
Ri_min = inf; Ri_max = 0;

% loop of trials
trial = 1;
while trial <= Options.trials
    if Options.verbose==2,
        fprintf('trial %i/%i:\n', trial, Options.trials); 
    elseif Options.verbose==1,
        fprintf('trial %i/%i: %i', trial, Options.trials, Ri.nR);
    end
    Ri_mer = merge(Ri, Options);
    if Options.verbose==1,
        fprintf(' --> %i', Ri_mer.nR); 
    end
    
    % update minimum and maximum
    Ri_min = min(Ri_min, Ri_mer.nR);
    Ri_max = max(Ri_max, Ri_mer.nR);
    if Ri_mer.nR < Ri_best.nR
        Ri_best = Ri_mer;
        trial = trial + 1;
        if Options.verbose==1,
            fprintf(' (new minimum)\n'); 
        end
    else
        trial = trial + 1;
        if Options.verbose==1,
            fprintf('\n'); 
        end
    end;
    if Ri_min<=1,
        break
    end
end;

PolyOut = Ri_best.Pn;

if Options.verbose>0,
    fprintf('  ==> min: %i  max: %i\n', Ri_min, Ri_max');
end

details.after = Ri_min;
details.runTime = cputime - starttime;
details.alg = 'greedy';

return


%-----------------------------------------------------------------------
function Ri = merge(Ri, Options);

% start the trial with the permutation K
iterateAgain = 1;
iter = 0;
while iterateAgain
    
    if Options.statusbar
        if isempty(mpt_statusbar(Options.status_handle, mod(iter, 10) / 10, Options.status_min, Options.status_max)),
            mpt_statusbar;
            error('Break...');
        end
    end
    
    iterateAgain = 0;
    iter = iter+1;
    if Options.verbose>0, 
        fprintf('  iteration %i: merging %i --> ', iter, Ri.nR); 
    end;

    % indicators whether region k has been (1) included in some union or not (0)
    issorted=zeros(Ri.nR,1);
    newRi.nR=0;
    newRi.Pn = Options.emptypoly;
    
    % find a random permutation of the indices (of polyhedra)
    K = randperm(Ri.nR);

    old2new = K;
    k_cnt = 0;
    for k = K
        if Options.statusbar
            k_cnt = k_cnt + 1;
            if mod(k_cnt, 5) == 0,        
                if isempty(mpt_statusbar(Options.status_handle, mod(iter, 10) / 10, Options.status_min, Options.status_max)),
                    mpt_statusbar;
                    error('Break...');
                end
            end
        end
        
        Pc = Ri.Pn(k);
        BCc=setdiff(Ri.BC{k},0);
        changed=0;
        if issorted(k)
            continue;
        end
        firstloop=1;
        while firstloop | changed
            firstloop=0;
            changed=0;
            BCc_old=BCc;
            BCc_old=BCc_old(randperm(length(BCc_old))); % permute BBc_old
            for ind_l=1:length(BCc_old), %1:Ri.nR
                ind2=BCc_old(ind_l);

                if(ind2==k | issorted(ind2))
                    continue;
                end

                [Pu,how] = union([Pc Ri.Pn(ind2)],Options);
                if how
                    Pc = Pu;
                    BCc=union(BCc, setdiff(Ri.BC{ind2},0));
                    if Options.verbose==2, disp(['regions ' num2str([k ind2]) ' are joined']); end;
                    issorted(ind2)=1;
                    old2new(ind2)=newRi.nR+1;
                    changed=1;
                    iterateAgain = 1;
                end
            end
        end
        newRi.nR=newRi.nR+1;
        old2new(k)=newRi.nR;
        newRi.Pn = [newRi.Pn Pc];
        newRi.BC{newRi.nR}=BCc;
        issorted(k)=1;
    end
    for i=1:newRi.nR
        newRi.BC{i}=unique(old2new(newRi.BC{i}));
    end

    if Options.verbose>0,
        percent = (Ri.nR-newRi.nR)/Ri.nR * 100;
        fprintf('%i (%2.1f percent)\n', newRi.nR, percent);
    end;

    Ri = newRi;
    clear newRi;

    if Options.multiple == 0
        iterateAgain = 0;
    end;
end;

if Options.closestatbar,
    mpt_statusbar;
end

return


%-----------------------------------------------------------------------
function BC = sub_buildBClist(Hn,Kn,Pn,Options);
% Inputs:  Ri.Hn, Ri.Kn: cell structure holding polytopes
%          Ri.nR: number of polytopes
% Outputs: BC: list of neighbours

alpha=5*Options.rel_tol;
tolerance=Options.abs_tol;

nR = length(Pn);
nx = dimension(Pn);

for ii=1:nR
    
    % display current number of polyhedron
    BC{ii}=[];
    % generate boundary points and check where they are
    
    Pnii = Pn(ii);
    Hnii = Hn{ii};
    for kk=1:nconstr(Pnii)
        [xfacet,rfacet]=facetcircle(Pnii,kk);
        if rfacet<Options.abs_tol & nx > 1,
            % facetcircle() always returns the radius as zero for 
            % 1-D polytopes. for dimensions above 1, zero radius indicates
            % a problem.
            continue;
        end
        xbeyond=xfacet+alpha*Hnii(kk,:)';
        
        found=0;
        for jj=1:nR
           if jj==ii
               continue;
           end
           if all(Hn{jj}*xbeyond-Kn{jj}<=Options.abs_tol)
               found=1;
               BC{ii}(end+1)=jj;
           end
        end
        if ~found
            BC{ii}(end+1)=0;
        end
    end
end