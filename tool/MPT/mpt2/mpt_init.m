function out=mpt_init(varargin)
%MPT_INIT Initializes the MPT toolbox
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Any routine of the MPT toolbox can be called with user-specified values
% of different parameters. To make usage of MPT toolbox as user-friendly as
% possible, we provide the option to store default values of the parameters
% in variable mptOptions, which is stored in MATLAB's workspace as a global
% variable (i.e. it stays there unless one types 'clear all'). 
%
% Please have a look at the code bellow to see what are the default values
% of the fields of mptOptions.
%
% Fields of mptOptions:
%
% lpsolver       - variable which sets the default LP solver. 
%                  allowed values:
%                    'nag'       - NAG LP solver (e04naf)
%                    'e04mbf'    - NAG LP solver (e04mbf)
%                    'cdd'       - CDD Criss-Cross
%                    'cplex'     - CPLEX 9 LP solver
%                    'cplex8'    - CPLEX 8 LP solver
%                    'glpk'      - GLPK solver (interfaced with glpkmex)
%                    'glpkcc'    - GLPK solver (interfaced with glpkcc)
%                    'linprog'   - Matlab's linprog
%                    'sedumi'    - SeDuMi
%                    'qsopt'     - QSopt
%                    'xpress'    - XPRESS
%                    'mosek'     - MOSEK
%                    'ooqp'      - OOQP
%                    'clp'       - CLP
%                    'bpmpd'     - BPMPD
%                    'cplexmex'  - CPLEX (interfaced with CPLEXMEX)
%
% qpsolver       - variable which sets the default QP solver
%                  allowed values:
%                    'nag'       - NAG QP solver
%                    'cplex'     - CPLEX 9 QP solver
%                    'cplex8'    - CPLEX 8 QP solver
%                    'quadprog'  - Matlab's quadprog
%                    'xpress'    - XPRESS
%                    'mosek'     - MOSEK
%                    'ooqp'      - OOQP
%                    'clp'       - CLP
%                    'bpmpd'     - BPMPD
%                    'cplexmex'  - CPLEX (interfaced with CPLEXMEX)
%
% milpsolver     - variable which sets the default MILP solver. 
%                  allowed values:
%                    'cplex'     - CPLEX 9 (interfaced with cplexint)
%                    'yalmip'    - YALMIP Branch & Bound algorithm
%                    'glpk'      - GLPK solver (interfaced with glpkmex)
%                    'glpkcc'    - GLPK solver (interfaced with glpkcc)
%                    'xpress'    - XPRESS
%                    'mosek'     - MOSEK
%                    'bintprog'  - bintprog.m
%                    'cplexmex'  - CPLEX (interfaced with CPLEXMEX)
%
% miqpsolver     - variable which sets the default MIQP solver. 
%                  allowed values:
%                    'cplex'     - CPLEX 9 (interfaced with cplexint)
%                    'yalmip'    - YALMIP Branch & Bound algorithm
%                    'xpress'    - XPRESS
%                    'mosek'     - MOSEK
%                    'cplexmex'  - CPLEX (interfaced with CPLEXMEX)
%
% extreme_solver - method to use for extreme points and convex hull enumeration
%                  allowed values:
%                     'matlab'   - analytical enumeration
%                     'lrs'      - LRS algorithm (matlab implementation)
%                     'cdd'      - CDD
%
% details        - defines how many details about the solution should be stored
%                  in the resulting controller structure.
% rel_tol        - relative tolerance for computation
% abs_tol        - absolute tolerance for computation
% step_size      - step of a size over a facet in mpQP/mpLP
% debug_level    - level of debugging in mpQP/mpLP
% infbox         - the R^n polyhedra is internally converted to a box,
%                  this option defines bounds of that box
% newfigure      - if set to 1, each picture will be plotted in a separate
%                  figure window
% verbose        - level of verbosity of the algorithms
%
%
% USAGE:
%
%   mpt_init
%   mpt_init('Property',Value,'Property',Value,...)
%   e.g.
%   mpt_init('qpsolver','quadprog','lpsolver','cdd','debug_level',2)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pair(s) Property, Value
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% mptOptions structure
%

% Copyright is with the following author(s):
%
% (C) 2007 Michal Kvasnica, Slovak University of Technology in Bratislava
%          michal.kvasnica@stuba.sk
% (C) 2003--2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%                kvasnica@control.ee.ethz.ch

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

mpt_ver = '2.6.3';

if mpt_matlabrelease < 12,
    error('MPT only works with Matlab R12 or newer.');
end

% returns version of MPT
if nargin>0,
    if strcmp(varargin{1},'version'),
        out = mpt_ver;
        return
    end
end

%-------------------------------------------
% Enable/disable GUI setup availability
%-------------------------------------------

announcegui = 1;


%-------------------------------------------
% Enable/disable automatic check for updates
%-------------------------------------------

% MPT can automatically check the product website
% (http://control.ee.ethz.ch/~mpt/) for updates. If new version is available for
% download, the user will be notified. The check is run every time mpt_init is
% called (typically at the start of new Matlab session or manually, by call to
% mpt_update).

mptOptions.checkupdates = 0;

% set checkupdates to 0 if you DON'T WANT mpt_init to connect to MPT website and
% check for new versions. No data will be send to the server (see mpt_update.m
% for source-code)

% Note: this feature is only available under Matlab 6.5 and Matlab 7



%---------------
% Path to HYSDEL
%---------------

% HYSDEL (HYbrid Systems DEscription Language) is a tool for rapid modelling of
% hybrid systems. If you want to use this feature in MPT and you are not running
% Windows, Linux or Solaris operating system, please provide full path to HYSDEL
% executable in "mptOptions.hysdelpath" variable.

comp = computer;

if ispc,
    mptOptions.hysdelpath = which('hysdel.exe');
elseif strcmpi(comp, 'glnx86'),
    mptOptions.hysdelpath = which('hysdel.linux');
elseif strcmpi(comp, 'sol2'),
    mptOptions.hysdelpath = which('hysdel.sol');
else
    % enter path to HYSDEL binary if you are not running Windows, Linux or
    % Solaris OS
    mptOptions.hysdelpath = '';
end



%------------------
% Default LP Solver
%------------------

% which LP solver to use:
%  0 - 'nag'       - NAG LP solver (e04naf)
%  9 - 'e04mbf'    - NAG LP solver (e04mbf)
%  3 - 'cdd'       - CDD Criss-Cross
%  2 - 'cplex'     - CPLEX 9 LP solver
%  8 - 'cplex8'    - CPLEX 8 LP solver
%  4 - 'glpk'      - GLPK solver
%  1 - 'linprog'   - Matlab's linprog
%  6 - 'sedumi'    - SeDuMi
%  7 - 'qsopt'     - QSopt
% 10 - 'xpress'    - XPRESS
% 11 - 'mosek'     - MOSEK
% 12 - 'ooqp'      - OOQP
% 13 - 'clp'       - CLP
% 14 - 'bpmpd'     - BPMPD
% 15 - 'cplexmex'  - CPLEXMEX

% fastest solver will be chosen if this field is empty (recommended)
mptOptions.lpsolver = [];

% by default 'mpt_init' will try to execute every solver which is found on your
% matlab path, to see if it is really working properly. however, due to a change
% in Matlab R14 API, some solver interfaces can crash Matlab if run on Linux.
% In such case provide the list of solvers which should be excluded from being
% checked in following variable.
%
% e.g. to exclude all CPLEX interfaces from checking, define
%   dontcheck.lp = [2 8 15]
dontcheck.lp = [];

% if you set mptOptions.rescueLP to 1 and some LP is infeasible with
% default solver, the problem will be automatically re-solved using a
% different solver which is available on your machine (if any)
mptOptions.rescueLP = 0;


%------------------
% Default QP Solver
%------------------

% which QP solver to use
%   0 - 'nag'       - NAG QP solver
%   2 - 'cplex'     - CPLEX 9 QP solver
%   4 - 'cplex8'    - CPLEX 8 QP solver
%   1 - 'quadprog'  - Matlab's quadprog
%   5 - 'xpress'    - XPRESS
%   6 - 'mosek'     - MOSEK
%   7 - 'ooqp'      - OOQP
%   3 - 'sedumi'    - SeDuMi
%   8 - 'clp'       - CLP
%   9 - 'bpmpd'     - BPMPD
%  10 - 'cplexmex'  - CPLEXMEX

% fastest solver will be chosen if this field is empty (recommended)
mptOptions.qpsolver = [];

% by default 'mpt_init' will try to execute every solver which is found on your
% matlab path, to see if it is really working properly. however, due to a change
% in Matlab R14 API, some solver interfaces can crash Matlab if run on Linux.
% In such case provide the list of solvers which should be excluded from being
% checked in following variable.
%
% e.g. to exclude all CPLEX interfaces from checking, define
%   dontcheck.qp = [2 4 10]
dontcheck.qp = [];

% if you set mptOptions.rescueQP to 1 and some QP is infeasible with
% default solver, the problem will be automatically re-solved using a
% different solver which is available on your machine (if any)
mptOptions.rescueQP = 0;



%--------------------
% Default MILP Solver
%--------------------

% which MILP solver to use
%   0 - 'cplex'     - CPLEX 9 (interfaced with cplexint)
%   1 - 'yalmip'    - YALMIP Branch & Bound algorithm
%   2 - 'glpk'      - GLPK solver (interfaced with glpkmex)
%   3 - 'xpress'    - XPRESS
%   4 - 'mosek'     - MOSEK
%   5 - 'bintprog'  - bintprog.m
%   6 - 'cplex8'    - CPLEX 8 (interfaced with milp_cplex)
%   7 - 'cplexmex'  - CPLEXMEX

% fastest solver will be chosen if this field is empty (recommended)
mptOptions.milpsolver = [];

% by default 'mpt_init' will try to execute every solver which is found on your
% matlab path, to see if it is really working properly. however, due to a change
% in Matlab R14 API, some solver interfaces can crash Matlab if run on Linux.
% In such case provide the list of solvers which should be excluded from being
% checked in following variable.
%
% e.g. to exclude all CPLEX interfaces from checking, define
%   dontcheck.milp = [0 6 7]
dontcheck.milp = [];


%--------------------
% Default MIQP Solver
%--------------------

% which MIQP solver to use
%   0 - 'cplex'     - CPLEX 9 (interfaced with cplexint)
%   1 - 'yalmip'    - YALMIP Branch & Bound algorithm
%   2 - 'xpress'    - XPRESS
%   3 - 'mosek'     - MOSEK
%   4 - 'cplex8'    - CPLEX 8 (interfaced with miqp_cplex)
%   5 - 'cplexmex'  - CPLEXMEX

% fastest solver will be chosen if this field is empty (recommended)
mptOptions.miqpsolver = [];

% by default 'mpt_init' will try to execute every solver which is found on your
% matlab path, to see if it is really working properly. however, due to a change
% in Matlab R14 API, some solver interfaces can crash Matlab if run on Linux.
% In such case provide the list of solvers which should be excluded from being
% checked in following variable.
%
% e.g. to exclude all CPLEX interfaces from checking, define
%   dontcheck.miqp = [0 4 5]
dontcheck.miqp = [];


%---------------------------------------------------
% Default solver for extreme points and convex hulls
%---------------------------------------------------

% which method to use for extreme points and convex hull computation
%   'cdd'      - CDD
%   'matlab'   - analytical enumeration
%   'lrs'      - LRS algorithm (matlab implementation)

% fastest solver will be chosen if this field is empty (recommended)
mptOptions.extreme_solver = [];




%---------------------------
% Default absolute tolerance
%---------------------------

% absolute tolerance
mptOptions.abs_tol = 1e-7;



%---------------------------
% Default relative tolerance
%---------------------------

% relative tolerance
mptOptions.rel_tol = 1e-6;




%-------------------------------------------
% Default step size for mpQP/mpLP Algorithms
%-------------------------------------------

% step of a size over a facet in mpQP/mpLP
mptOptions.step_size = 1e-5;


%---------------------------------
% Default level of double-checking
%---------------------------------

% 0 - No additional checks of results
% 1 - Perform double-checks when necessary
% 2 - Always double-check results

mptOptions.debug_level = 1;



%------------------------------------
% Default bounds for the Infinity-box
%------------------------------------

% the R^n polyhedra is internally converted to a box, this option defines bounds
% of that box 
% ATTENTION: we have observed certain numerical problems when value of the
% infbox is too big, therefore we recommend not to change the default setting
% unless necessary 

mptOptions.infbox = 10000;



%--------------------------------------------
% Default value for the "newfigure" parameter
%--------------------------------------------

% if set to 1, each picture will be plotted in a separate figure window

mptOptions.newfigure = 0;



%---------------------------
% Default value of verbosity
%---------------------------

% verbosity level (0 - important results only, 1 - intermediate results
% displayed, 2 - all results displayed) 

mptOptions.verbose = 1;


%----------------------------------
% Details stored about the solution
%----------------------------------

% defines how many details about the solution should be stored in the resulting
% controller structure. This can have a significant impact on the size of the
% structure. If you want to evaluate open-loop solution for PWA systems, set
% this to 1. Otherwise leave the default value.

mptOptions.details = 0;



%-----------------------
% Default YALMIP options
%-----------------------

% there are two possible options here:
%    'fixer'    - each branch&bound iteration takes longer, but this
%                 algorithm finds an initial feasible point quicker than
%                 'rounder' does
%    'rounder'  - faster iterations, but it takes longer to find a feasible
%                 initial point
bnb_upper_solver = 'fixer';

yalmip_options = {'verbose', 0, ...
    'warning', 0, ...
    'cachesolvers', 0, ...
    'bnb.branchrule', 'weight', ...
    'bnb.upper', bnb_upper_solver };



%-------------------------------------------------------------------------------
% DO NOT EDIT BEYOND THIS LINE!!!
%-------------------------------------------------------------------------------

% if set to 1, mpt_init will try to solve a problem with each solver to see
% which solvers are actually available.
% if set to 0, existence of a particular solver will be deduced upon existence
% of the corresponding interface file.
executesolver = 1;

% the "dosave" flag decides if global settings will be stored using setpref()
% function permanently. The flag will become true if mpt_init is called as:
%
%   mpt_init('save')
%
% or
%
%   mpt_init('lpsolver', 'cdd', ..., 'save')
dosave = 0;

% how many "rescue" LP solvers should be considered? increasing the value too
% musch can also lead to slow-downs!
nrescue.lp = 2;
nrescue.qp = 0;
nrescue.milp = 0;
nrescue.miqp = 0;
nrescue.extreme = 5;

% set prefered solvers, list organized according to solver's speed on a given
% problem
prefered.lp = [0 9 3 15 2 8 17 4 14 7 1 13 5 10 11 12 16];
prefered.qp = [0 1 9 2 8 4 5 6 7 10];
prefered.milp = [0 7 6 3 8 2 4 5 1];
prefered.miqp = [0 5 4 2 3 1];
prefered.extreme = [3 4 0 2 1];

% extend list of solvers which should be excluded from run tests
dontcheck.lp = [dontcheck.lp(:); 6]; % do not check SeDuMi - leads to a crash with Matlab 6.5
dontcheck.qp = [dontcheck.qp(:); 3]; % do not check SeDuMi - leads to a crash with Matlab 6.5
dontcheck.milp = dontcheck.milp(:);
dontcheck.miqp = dontcheck.miqp(:);
dontcheck.extreme = [];

nargs = nargin;
vargs = varargin;

% check if path is set properly
p = path;
mptpath = {
    ['extras' filesep 'analysis'],  ...
        ['extras' filesep 'auxiliary'], ...
        ['extras' filesep 'control'], ...
        ['extras' filesep 'control' filesep 'hys2pwa'], ...
        ['extras' filesep 'control' filesep 'mldmpc'], ...
        ['extras' filesep 'control' filesep 'optmerge'], ...
        ['extras' filesep 'geometry'], ...
        ['extras' filesep 'graphics'], ...
        ['extras' filesep 'gui'], ...
        ['extras' filesep 'simulink'], ...
        'solvers', ...
        'examples', ...
        ['examples' filesep 'demos'], ...
        'extras', ...
        'operators', ...
        ['modules' filesep 'global'], ...
        ['modules' filesep 'moment'], ...
        ['modules' filesep 'parametric'], ...
        ['modules' filesep 'sos']
};

% check if mpt_solverInfo is on the path
w = which('mpt_solverInfo');
if isempty(w),
    fprintf('\n\nPath is not set correctly!\n');
    fprintf('Did you use "addpath(genpath(''/path/to/mpt/''))"?\n\n');
    fprintf('Read the installation notes for more details:\n');
    fprintf('http://control.ee.ethz.ch/~mpt/docs/install.php\n\n');
    error('Cannot continue, you must set path to all subdirectories of MPT first.');
end

% check if YALMIP is on path just once
yalmiplocation = which('yalmip.m', '-all');
if any(size(yalmiplocation)>1),
    disp('Warning: You have multiple occurences of YALMIP on your path! This could lead to serious consequences.');
    yalmiplocation
end

% check if MPT is on path just once
mptlocation = which('mpt_control.m', '-all');
if any(size(mptlocation)>1),
    disp('Warning: You have multiple occurences of MPT on your path! This could lead to serious consequences.');
    mptlocation
end

% check if hashtable is on path
havehashtable = 1;
hashtablelocation = which('hashtable.m', '-all');
if isempty(hashtablelocation),
    disp('Warning: "mpt/extras/auxiliary" directory not in path or "mpt/extras/auxiliary/@hashtable" missing, check your installation.');
    havehastable = 0;
elseif any(size(hashtablelocation)>1),
    disp('Warning: multiple occurences of the HASHTABLE object on your path, this could lead to errors.');
    hashtablelocation
    havehashtable = 0;
end
    

% list of solvers interfaced by YALMIP
yalmip_solver_strings = {'xpress', 'mosek', 'bintprog', ...
        'cplex-milp-cplexint', 'cplex-miqp-cplexint', ...
        'bnb', 'sedumi', 'ooqp'};

for ii = 1:length(mptpath),
    if isempty(findstr(p, mptpath{ii})),
        % one of the necessary directories is not included in matlab path. this
        % can lead to missing and/or corrupted functionality. be sure to add the
        % whole MPT directory along with all subfolders to your matlab path.
        fprintf('\n\nPath to subdirectory "%s" is not set!\n\n', mptpath{ii});
        fprintf('Did you follow the installation instructions?\n');
        error('Path is not set correctly.');
    end
end

% try to load settings using "getpref()" and exit quickly without for checking
% available solvers.
% however we only do so if mpt_init was called with no input arguments
try
    if nargs==0 & ispref('MPT_toolbox'),
        fprintf('Restoring saved settings...\n');
        mptOptions = getpref('MPT_toolbox', 'mptOptions');

        if ~isstruct(mptOptions),
            % loaded settings is not a structure, remove it and instruct the
            % user to run mpt_init once more
            rmpref('MPT_toolbox');
            disp('Corrupted settings found, re-initializing...');
            dosave = 1;
            % this error will be catched by the master try-block
            error('Corrupted settings found.');
        end
        
        if ~isfield(mptOptions, 'version'),
            mptOptions.version = 'unknown';
        end
        if ~strcmp(mptOptions.version, mpt_ver),
            % this settings was saved for some other version, re-initialize
            % the toolbox
            disp('Settings saved for older version, re-initializing...');
            dosave = 1;
            error('Settings saved for older version.');
        end
        
        try
            % update mptOptions.sdpsettings if somebody installed a new version
            % of YALMIP
            mptOptions.sdpsettings = sdpsettings(yalmip_options{:});
        catch
            warning('YALMIP not found, some functionality may not be accessible.');
            mptOptions.sdpsettings = [];
        end
        
        try
            % prepare dummy model for fast execution of YALMIP-interfaced solvers
            % although this data is already save in mptOptions, we repeat the
            % assignment to deal with cases when the user installs a newer
            % version of YALMIP manually.
            mptOptions.yalmipdata = sub_prepareyalmipdata(yalmip_solver_strings);
        catch
            % if anything fails, remove the "yalmipdata" field just to be sure
            if isfield(mptOptions, 'yalmipdata'),
                mptOptions = rmfield(mptOptions, 'yalmipdata');
            end
        end

        sub_printcopyright(mpt_ver);
        sub_printsolvers(mptOptions);
        
        % check the MPT web-page for updates if desired
        % only available in Matlab 6.5 and Matlab 7
        if mptOptions.checkupdates & exist('urlread','file'),
            disp('Checking for updates (can be disabled by modifying mpt_init.m)...');
            mpt_update;
        elseif mptOptions.checkupdates
            disp('Automatic update works only with Matlab 6.5 and newer.');
        end
        
        if announcegui,
            fprintf(' Run ''mpt_studio'' to start the GUI. Run ''mpt_setup'' to set global parameters.\n\n');
        end
        
        if exist('mpt_verifySolution', 'file'),
            % "mpt_verifySolution.m" was removed in MPT 2.0 final, if it still exists,
            % something is wrong. before installing any version of MPT, you must first
            % remove any previous installation from your disk.
            warning('It appears that previous version of MPT was not removed from your computer. This can lead to serious problems! Please remove all previous versions of MPT from your disk and install latest version.');
        end
        
        out=mptOptions;
        if nargout<1,
            clear out
        end
        
        return
    end
end

% check if input arguments consist of pairs PropertyName, PropertyValue
if nargs==1,
    % mpt_init('rehash') can be used to re-scan available solvers and update
    % preferences stored by setpref()
    if isa(vargs{1}, 'char'),
        if strcmpi(vargs{1}, 'rehash'),
            if ispref('MPT_toolbox');
                rmpref('MPT_toolbox');
            end
            nargs = 0;
            dosave = 1;
        elseif strcmpi(vargs{1}, 'save'),
            dosave = 1;
            nargs = 0;
        end
    end
else
    for ii = 1:nargs,
        if isa(vargs{ii}, 'char'),
            if strcmpi(vargs{ii}, 'save'),
                % if the flag is equal to 'save', set dosave flag
                dosave = 1;
                
                % remove the 'save' flag from list of input arguments
                vargs = {vargs{setdiff(1:nargs, ii)}};
                
                % decrease number of input arguments by 1
                nargs = nargs - 1;

                % and exit this for-loop
                break
            end
        end
    end
end

if rem(nargs, 2)~=0,
    clear global mptOptions
    error(['mpt_init: Input arguments following the object name must be pairs', ...
            ' of the form PropertyName, PropertyValue']);
end

try
    mptOptions.sdpsettings = sdpsettings(yalmip_options{:});
catch
    warning('YALMIP not found, some functionality may not be accessible.');
    mptOptions.sdpsettings = [];
end

if executesolver,
    fprintf('looking for available solvers...\n');
end

solvers.lp = test_solvers('lp', prefered.lp, sub_maxsolvers('lp'), dontcheck.lp, nrescue.lp, executesolver);
solvers.lp_all = solvers.lp;
solvers.qp = test_solvers('qp', prefered.qp, sub_maxsolvers('qp'), dontcheck.qp, nrescue.qp, executesolver);
solvers.milp = test_solvers('milp', prefered.milp, sub_maxsolvers('milp'), dontcheck.milp, nrescue.milp, executesolver);
solvers.miqp = test_solvers('miqp', prefered.miqp, sub_maxsolvers('miqp'), dontcheck.miqp, nrescue.miqp, executesolver);

if isempty(solvers.lp),
    clear global mptOptions
    error('mpt_init: No supported LP solver available on your system, cannot proceed!');
end
if isempty(solvers.qp)
    solvers.qp = -1;
    mptOptions.qpsolver = -1; % NONE
    disp('There is no QP solver installed on your system.');
    disp('You will not be able to solve problems with quadratic cost function.');
    disp('Please be sure that you set ''probStruct.norm = 1'' before calling any control routine!');
end
if isempty(solvers.milp),
    warning('mpt_init: No supported MILP solver available on your system.');
end
if isempty(solvers.miqp),
    warning('mpt_init: No supported MIQP solver available on your system.');
end

if isempty(mptOptions.lpsolver),
    mptOptions.lpsolver = solvers.lp(1);
end
if isempty(mptOptions.qpsolver),
    mptOptions.qpsolver = solvers.qp(1);
end
if ~isempty(solvers.milp) & isempty(mptOptions.milpsolver),
    mptOptions.milpsolver = solvers.milp(1);
elseif isempty(mptOptions.milpsolver),
    mptOptions.milpsolver = -1;
end
if ~isempty(solvers.miqp) & isempty(mptOptions.miqpsolver),
    mptOptions.miqpsolver = solvers.miqp(1);
elseif isempty(mptOptions.miqpsolver),
    mptOptions.miqpsolver = -1;
end

% set the appropriate fields based on input arguments
for ii=1:2:nargs
    if strcmpi(vargs{ii}, 'save'),
        dosave = 1;
    elseif ~isfield(mptOptions,vargs{ii})
        clear global mptOptions
        error(['mpt_init: Non-existing property (' vargs{ii} ')']);
    end
    mptOptions=setfield(mptOptions,vargs{ii},vargs{ii+1});
end

if ischar(mptOptions.lpsolver),
    str_solver = mptOptions.lpsolver;
    [mptOptions.lpsolver, err] = mpt_solverInfo('lp', mptOptions.lpsolver);
    if err,
        clear global mptOptions
        error(['mpt_init: Unknown LP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.qpsolver),
    str_solver = mptOptions.qpsolver;
    [mptOptions.qpsolver, err] = mpt_solverInfo('qp', mptOptions.qpsolver);
    if err,
        clear global mptOptions
        error(['mpt_init: Unknown QP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.milpsolver),
    str_solver = mptOptions.milpsolver;
    [mptOptions.milpsolver, err] = mpt_solverInfo('milp', mptOptions.milpsolver);
    if err,
        clear global mptOptions
        error(['mpt_init: Unknown MILP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.miqpsolver),
    str_solver = mptOptions.miqpsolver;
    [mptOptions.miqpsolver, err] = mpt_solverInfo('miqp', mptOptions.miqpsolver);
    if err,
        clear global mptOptions
        error(['mpt_init: Unknown MIQP solver ''' str_solver ''' !']);
    end
end

if ~any(mptOptions.debug_level==[0 1 2]),
    clear global mptOptions
    error('mpt_init: debug_level can only have values 0, 1, or 2!');
end

if ~any(mptOptions.newfigure==[0 1]),
    clear global mptOptions
    error('mpt_init: newfigure can only have values 0 or 1!');
end

if ~any(mptOptions.verbose==[0 1 2]),
    clear global mptOptions
    error('mpt_init: verbose can only have values 0, 1, or 2!');
end


if ~test_lp(mptOptions.lpsolver);
    disp('The currently selected LP solver is not installed on your system.')
    disp('Please modify the entry "lpsolver" in "mpt_init.m"')
    if mptOptions.lpsolver==7,
        fprintf('\nDid you properly install QSopt?\n');
    end
    clear global mptOptions
    error('Error: cannot proceed with mpt_init');
    return
end

if mptOptions.qpsolver ~= -1,
    if ~test_qp(mptOptions.qpsolver)
        clear global mptOptions
        disp('The currently selected QP solver is not installed on your system.')
        disp('Please modify the entry "qpsolver" in "mpt_init.m"')
        error('Error: cannot proceed with mpt_init');
    end
else
    disp('There is no QP solver installed on your system.');
    disp('You will not be able to solve problems with quadratic cost function.');
    disp('Please be sure that you set ''probStruct.norm = 1'' before calling any control routine!');
end

mptOptions.emptypoly = polytope;

solvers.extreme = test_solvers('extreme', prefered.extreme, sub_maxsolvers('extreme'), dontcheck.extreme, nrescue.extreme, executesolver);

if ischar(mptOptions.extreme_solver),
    str_solver = mptOptions.extreme_solver;
    [mptOptions.extreme_solver, err] = mpt_solverInfo('extreme', mptOptions.extreme_solver);
    if err,
        clear global mptOptions
        error(['mpt_init: Unknown extreme point method ''' str_solver ''' !']);
    end
end

if isempty(mptOptions.extreme_solver)
    if ~isempty(solvers.extreme),
        mptOptions.extreme_solver = solvers.extreme(1);
    else
        error('mpt_init: No supported extreme point enumeration method found!');
        mptOptions.extreme_solver = -1;
    end
end

if ~test_extreme(mptOptions.extreme_solver)
    clear global mptOptions
    disp('The currently selected method for vertex enumeration failed.')
    disp('Please modify the entry "extreme_solver" in "mpt_init.m"')
    error('Error: cannot proceed with mpt_init');
end

mptOptions.solvers = solvers;

sub_printcopyright(mpt_ver);

if mptOptions.lpsolver==1,
    disp('WARNING: you have chosen linprog as a default LP solver.');
    disp('This solver is very slow and numerically not robust!');
    disp('We strongly advice you to use some other alternative if possible.');
    fprintf('\n');
end

out=mptOptions;

sub_printsolvers(mptOptions);

% check the MPT web-page for updates if desired
% only available in Matlab 6.5 and Matlab 7
if mptOptions.checkupdates & exist('urlread','file'),
    disp('Checking for updates (can be disabled by modifying mpt_init.m)...');
    mpt_update;
elseif mptOptions.checkupdates
    disp('Automatic update works only with Matlab 6.5 and newer.');
end

if announcegui,
    fprintf(' Run ''mpt_studio'' to start the GUI. Run ''mpt_setup'' to set global parameters.\n\n');
end

if exist('mpt_verifySolution', 'file'),
    % "mpt_verifySolution.m" was removed in MPT 2.0 final, if it still exists,
    % something is wrong. before installing any version of MPT, you must first
    % remove any previous installation from your disk.
    warning('It appears that previous version of MPT was not removed from your computer. This can lead to serious problems! Please remove all previous versions of MPT from your disk and install latest version.');
end

if nargout<1,
    clear out
end

if havehashtable,
    % prepare dummy model for fast execution of YALMIP-interfaced solvers
    mptOptions.yalmipdata = sub_prepareyalmipdata(yalmip_solver_strings);
end

% save version string to mptOptions
mptOptions.version = mpt_ver;

if dosave,
    % save settings such that they can be later restored by getpref()
    %
    % NOTE! once you save the options permanently, MPT will not look for
    % available solvers at every. This will lead to a faster initialization, but
    % will not reflect any new solvers added AFTER the settings are saved. In
    % order to rehash the actual state, call:
    %
    % mpt_init('rehash')
    %
    setpref('MPT_toolbox', 'mptOptions', mptOptions);
end

%------------------------------------------------------------------------
function yalmipdata = sub_prepareyalmipdata(yalmip_solver_strings)

% we use the new object - hashtable
yalmipdata = hashtable;
for ii = 1:length(yalmip_solver_strings),
    solstr = yalmip_solver_strings{ii};
    yalmipdata(['milp:' solstr]) = sub_getyalmipdata(solstr, 'milp');
    yalmipdata(['miqp:' solstr]) = sub_getyalmipdata(solstr, 'miqp');
end


%------------------------------------------------------------------------
function success=test_lp(solver, execute)

if nargin<2,
    % execute the solver by default
    execute = 1;
end

switch solver
    case 0, fname = 'e04naf';
    case 1, fname = 'linprog';
    case 2, fname = 'cplexint';
    case 3, fname = 'cddmex';
    case 4, fname = 'glpkmex';
    case 5, fname = 'cddmex';
    case 6, fname = 'solvesdp';
    case 7, fname = 'mexqsopt';
    case 8, fname = 'lp_cplex';
    case 9, fname = 'e04mbf';
    case 10, fname = 'mexpress';
    case 11, fname = 'mosekopt';
    case 12, fname = 'ooqp.m';
    case 13, fname = 'mexclp';
    case 14, fname = 'bp';
    case 15, fname = 'cplexmex';
    case 16, fname = 'pdco';
    case 17, fname = 'glpkcc';
    otherwise, fname = '';
end

if ~exist(fname,'file'),
    success = 0;
    return
end
if ~execute,
    % don't run the solver
    success = 1;
    return
end

success = 0;
try
    cmd='[xopt,fval,lambda,exitflag,how]=mpt_solveLP(1,[1;-1],[1;1],[],[],[],solver);';
    T = evalc(cmd);
    if isnan(xopt),
        return
    end
    if strcmpi(how,'ok')
        success = 1;
    end
end


%------------------------------------------------------------------------
function success=test_qp(solver, execute)

if nargin<2,
    % execute the solver by default
    execute = 1;
end

switch solver
    case 0, fname = 'e04naf';
    case 1, fname = 'quadprog';
    case 2, fname = 'cplexint';
    case 3, fname = 'solvesdp';
    case 4, fname = 'qp_cplex';
    case 5, fname = 'mexpress';
    case 6, fname = 'mosekopt';
    case 7, fname = 'ooqp.m';
    case 8, fname = 'mexclp';
    case 9, fname = 'bp';
    case 10, fname = 'cplexmex';
    otherwise, fname = '';
end

if ~exist(fname,'file'),
    success = 0;
    return
end
if ~execute,
    % don't run the solver
    success = 1;
    return
end

success = 0;
try
    cmd='[xopt,lambda,how,exitflag,objqp]=mpt_solveQP(1,1,[1;-1],[1;1],[],[],[],solver);';
    T = evalc(cmd);
    if isnan(xopt),
        return
    end
    if exitflag==1,
        success = 1;
    end
end


%------------------------------------------------------------------------
function success=test_milp(solver, execute)

if nargin<2,
    % execute the solver by default
    execute = 1;
end

switch solver
    case 0, fname = 'cplexint';
    case 1, fname = 'solvesdp';
    case 2, fname = 'glpkmex';
    case 3, fname = 'mexpress';
    case 4, fname = 'mosekopt';
    case 5, fname = 'bintprog';
    case 6, fname = 'milp_cplex';
    case 7, fname = 'cplexmex';
    case 8, fname = 'glpkcc';
    otherwise, fname = '';
end

if ~exist(fname,'file'),
    success = 0;
    return
end
if ~execute,
    % don't run the solver
    success = 1;
    return
end

success = 0;
try
    cmd='[xmin,fmin,how,exitflag]=mpt_solveMILP(1,[1;-1],[1;1],[],[],0,1,''B'',[],[],solver);';
    T = evalc(cmd);
    if isnan(xmin),
        return
    end
    if exitflag==1 | exitflag==-1,
        success = 1;
    end
end


%------------------------------------------------------------------------
function success=test_miqp(solver, execute)

if nargin<2,
    % execute the solver by default
    execute = 1;
end

switch solver
    case 0, fname = 'cplexint';
    case 1, fname = 'solvesdp';
    case 2, fname = 'mexpress';
    case 3, fname = 'mosekopt';
    case 4, fname = 'miqp_cplex';
    case 5, fname = 'cplexmex';
    otherwise, fname = '';
end

if ~exist(fname,'file'),
    success = 0;
    return
end
if ~execute,
    % don't run the solver
    success = 1;
    return
end

success = 0;
try
    cmd='[xmin,fmin,how,exitflag]=mpt_solveMIQP(1,1,[1;-1],[1;1],[],[],0,1,''B'',[],[],solver);';
    T = evalc(cmd);
    if isnan(xmin),
        return
    end
    if exitflag==1 | exitflag==-1,
        success = 1;
    end
end


%------------------------------------------------------------------------
function success=test_extreme(solver, executesolver)

options.extreme_solver = solver;
if ~exist('cddmex', 'file') & solver==3,
    success = 0;
    return
end
success = 0;
try
    P = hull([-1 -1; -1 1; 1 -1; 1 1],options);
    success = 1;
end


%------------------------------------------------------------------------
function solvers = fixpreferred(solvers,mptOptions)

if ~isempty(solvers.lp),
    solvers.lp(find(solvers.lp==mptOptions.lpsolver)) = [];
end
solvers.lp = [mptOptions.lpsolver solvers.lp];

if ~isempty(solvers.qp),
    solvers.qp(find(solvers.qp==mptOptions.qpsolver)) = [];
end
solvers.qp = [mptOptions.qpsolver solvers.qp];

if ~isempty(solvers.milp),
    solvers.milp(find(solvers.milp==mptOptions.milpsolver)) = [];
end
solvers.milp = [mptOptions.milpsolver solvers.milp];

if ~isempty(solvers.miqp),
    solvers.miqp(find(solvers.miqp==mptOptions.miqpsolver)) = [];
end
solvers.miqp = [mptOptions.miqpsolver solvers.miqp];

if ~isempty(solvers.extreme),
    solvers.extreme(find(solvers.extreme==mptOptions.extreme_solver)) = [];
end
solvers.extreme = [mptOptions.extreme_solver solvers.extreme];


%------------------------------------------------------------------------
function maxs = sub_maxsolvers(solverclass)
% checks how many solvers of a given class are interfaced with MPT

MAX_SOLVERS = 20;  % check this many solvers

maxs = MAX_SOLVERS;
for index = 1:MAX_SOLVERS,
    solvername = mpt_solverInfo(solverclass, index);
    if strcmpi(solvername, 'unknown'),
        maxs = index-1;
        return
    end
end


%------------------------------------------------------------------------
function sub_printcopyright(mpt_ver)
% prints the copyright notice

fprintf('\nMPT toolbox %s initialized...\n', mpt_ver);
fprintf('Copyright (C) 2003-2006 by M. Kvasnica, P. Grieder and M. Baotic\n');
fprintf('\nSend bug reports, questions or comments to mpt@control.ee.ethz.ch\n');
fprintf('For news, visit the MPT web page at http://control.ee.ethz.ch/~mpt/\n');


%------------------------------------------------------------------------
function sub_printsolvers(mptOptions)

s_lpsolver = mpt_solverInfo('lp', mptOptions.lpsolver);
s_qpsolver = mpt_solverInfo('qp', mptOptions.qpsolver);
s_milpsolver = mpt_solverInfo('milp', mptOptions.milpsolver);
s_miqpsolver = mpt_solverInfo('miqp', mptOptions.miqpsolver);
s_exsolver = mpt_solverInfo('extreme', mptOptions.extreme_solver);

disp(['         LP solver: ' s_lpsolver]);
disp(['         QP solver: ' s_qpsolver]);
disp(['       MILP solver: ' s_milpsolver]);
disp(['       MIQP solver: ' s_miqpsolver]);
disp(['Vertex enumeration: ', s_exsolver]);
fprintf('\n');


%------------------------------------------------------------------------
function S = test_solvers(type, prefered, maxsol, dontcheck, nrescue, executesolver)
% tests each solver

% connect list of prefered solvers and include all solvers not listed among
% 'prefered' solvers:
allsolvers = [prefered setdiff(0:maxsol, prefered)];

% now remove solvers marked for not checking:
toremove = [];
for ii = 1:length(allsolvers),
    if any(allsolvers(ii) == dontcheck),
        toremove = [toremove ii];
    end
end

% sort the set difference back
allsolvers(toremove) = [];

S = [];
for solver = allsolvers,
    switch lower(type),
        case 'lp',
            t = test_lp(solver, executesolver);
        case 'qp',
            t = test_qp(solver, executesolver);
        case 'milp',
            t = test_milp(solver, executesolver);
        case 'miqp',
            t = test_miqp(solver, executesolver);
        case 'extreme',
            t = test_extreme(solver, executesolver);
    end
    if t,
        % solver exists and works correctly
        S = [S solver];
        if length(S) > nrescue,
            % enough rescue solvers found, break
            break
        end
    end
end
