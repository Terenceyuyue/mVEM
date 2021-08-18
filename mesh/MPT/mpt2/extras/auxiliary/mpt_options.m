function out=mpt_options(varargin)
%MPT_OPTIONS Read / modify solver settings for MPT
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Prints information about surrent solver setting if called without any input
% and output arguments
%
% Returns mptOptions structure if called as out = mpt_options
%
% Sets field of mptOptions structure by providing a pair of 'property' and value.
%
% Fields of mptOptions:
%
% see help mpt_init for more information
%
% USAGE:
%
%   mpt_options
%   options = mpt_options;
%   mpt_options('Property',Value,'Property',Value,...)
%   e.g.
%   mpt_options('qpsolver',0,'lpsolver',3,'debug_level',2)
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
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_init;
end

if nargin==1 & isstruct(varargin{1}),
    % special case - input is a structure
    opt = varargin{1};
    f = fields(opt);
    for ii=1:length(f),
        mptOptions=setfield(mptOptions,f{ii},getfield(opt,f{ii}));
    end
elseif nargin==1 & isa(varargin{1}, 'char'),
    % return respective field of mptOptions
    try
        out = getfield(mptOptions, varargin{1});
    catch
        error(sprintf('No such field "%s" in mptOptions.', varargin{1}));
    end
    return
else
    % check if input arguments consist of pairs PropertyName, PropertyValue
    if rem(nargin, 2)~=0,
        error(['mpt_options: Input arguments following the object name must be pairs', ...
                ' of the form PropertyName, PropertyValue']);
    end
    
    % set the appropriate fields based on input arguments
    for ii=1:2:nargin
        if ~isfield(mptOptions,varargin{ii})
            error(['mpt_options: Non-existing property (' varargin{ii} ')']);
        end
        mptOptions=setfield(mptOptions,varargin{ii},varargin{ii+1});
    end
end

if ischar(mptOptions.lpsolver),
    str_solver = mptOptions.lpsolver;
    [mptOptions.lpsolver, err] = mpt_solverInfo('lp', mptOptions.lpsolver);
    if err,
        error(['mpt_options: Unknown LP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.qpsolver),
    str_solver = mptOptions.qpsolver;
    [mptOptions.qpsolver, err] = mpt_solverInfo('qp', mptOptions.qpsolver);
    if err,
        error(['mpt_options: Unknown QP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.milpsolver),
    str_solver = mptOptions.milpsolver;
    [mptOptions.milpsolver, err] = mpt_solverInfo('milp', mptOptions.milpsolver);
    if err,
        error(['mpt_options: Unknown MILP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.miqpsolver),
    str_solver = mptOptions.miqpsolver;
    [mptOptions.miqpsolver, err] = mpt_solverInfo('miqp', mptOptions.miqpsolver);
    if err,
        error(['mpt_options: Unknown MIQP solver ''' str_solver ''' !']);
    end
end

if ischar(mptOptions.extreme_solver),
    str_solver = mptOptions.extreme_solver;
    [mptOptions.extreme_solver, err] = mpt_solverInfo('extreme', mptOptions.extreme_solver);
    if err,
        error(['mpt_options: Unknown extreme point method ''' str_solver ''' !']);
    end
end

s_lpsolver = mpt_solverInfo('lp', mptOptions.lpsolver);
s_qpsolver = mpt_solverInfo('qp', mptOptions.qpsolver);
s_milpsolver = mpt_solverInfo('milp', mptOptions.milpsolver);
s_miqpsolver = mpt_solverInfo('miqp', mptOptions.miqpsolver);
s_exsolver = mpt_solverInfo('extreme', mptOptions.extreme_solver);

if mptOptions.newfigure==1,
    s_newfigure = 'yes';
else
    s_newfigure = 'no';
end

mptOptions.emptypoly = polytope;

if nargout>0,
    out=mptOptions;
elseif nargin==0
    disp(['         LP solver: ' s_lpsolver]);
    disp(['         QP solver: ' s_qpsolver]);
    disp(['       MILP solver: ' s_milpsolver]);
    disp(['       MIQP solver: ' s_miqpsolver]);
    disp(['Vertex enumeration: ' s_exsolver]);

    fprintf('Absolute tolerance: %.1e\n', mptOptions.abs_tol);
    fprintf('Relative tolerance: %.1e\n', mptOptions.rel_tol);
    fprintf(' Infinity box size: %d\n', mptOptions.infbox);
    fprintf('         Step size: %.1e\n', mptOptions.step_size);
    fprintf('       Debug level: %d (0/1/2)\n', mptOptions.debug_level);
    fprintf('Level of verbosity: %d (0/1/2)\n', mptOptions.verbose);
    fprintf('  Level of details: %d (0/1)\n', mptOptions.details);
    fprintf('   Open new figure: %s\n', s_newfigure);
    fprintf('\n');
end
