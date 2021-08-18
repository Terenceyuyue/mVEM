function info(ctrl)
%INFO Displays more detailed information a given MPT controller
%
% info(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Displays detailed information about a given controller, such as control
% objectives, memory requirements, etc.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl - MPT controller object
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

sys = char(ctrl);

inpname = inputname(1);
if isempty(inpname),
    inpname = 'ctrl';
end

fprintf('\n');
disp([inpname '=']);
fprintf('\n');

disp(sys);
fprintf('\n');

bytes_per_float = 8;   % each floating point number is represented by 8 bytes
bytes_per_integer = 4; % each integer is represented by 4 bytes

if isexplicit(ctrl),
    haveSearchTree = isfield(ctrl.details, 'searchTree');
    [nx, nu] = mpt_sysStructInfo(ctrl.sysStruct);
    
    % =========================================================
    % information about the runtime
    fprintf('    Computation time: %.0f secs\n', ceil(runtime(ctrl)));


    % =========================================================
    % information about performed simplification
    if isfield(ctrl.details, 'alg'),
        switch lower(ctrl.details.alg),
            case 'greedy',
                algtype = 'greedy';
            case 'optimal',
                algtype = 'optimal';
            otherwise,
                algtype = 'unknown';
        end
        fprintf(' Simplification done: %s merging (from %d to %d regions in %.0f secs)\n', ...
            algtype, ctrl.details.before_simpl, ctrl.details.after_simpl, ...
            ceil(ctrl.details.simpl_time));
    else
        fprintf(' Simplification done: none (run "%s=mpt_simplify(%s)" to perform simplification)\n', ...
            inpname, inpname);
    end
    
    fprintf('\n');
    
    % =========================================================
    % compute memory requirements
    
    % * control law
    nfloats = length(ctrl.Pn)*(nu*nx + nu);
    memfloats = nfloats * bytes_per_float;
    if memfloats > 4096,
        memsize = sprintf('%.2f kB', memfloats/1024);
    else
        memsize = sprintf('%d B', ceil(memfloats));
    end
    fprintf('         Control law: %s (%d floats)\n', ...
        memsize, nfloats);

    % regions / search tree
    [H, K] = pelemfun(@double, ctrl.Pn);
    if haveSearchTree,
        % * region storage
        nhyp = size(ctrl.details.searchTree, 1);
        nfloats = nhyp*(size(ctrl.details.searchTree, 2)-2);
        nint = nhyp*2;
        flag = '(if implemented using a search tree)';
    else
        % * region storage
        nfloats = sum(cellfun('prodofsize', H)) + sum(cellfun('prodofsize', K));
        % * facets per regions vector
        nint = length(ctrl.Pn); 
        flag = '(if implemented in a brute-force way)';
    end
    % * control law
    memfloats = nfloats * bytes_per_float;
    memint = nint * bytes_per_integer;
    if memfloats+memint > 4096,
        memsize = sprintf('%.2f kB', (memfloats+memint)/1024);
    else
        memsize = sprintf('%d B', ceil(memfloats+memint));
    end
    fprintf('      Memory storage: %s (%d floats, %d integers) %s\n', ...
        memsize, nfloats, nint, flag);
    
    
    % =========================================================
    % information about the search tree
    if haveSearchTree,
        st_nhp = size(ctrl.details.searchTree, 1);
        % maximum number of hyperplanes that will be checked before an active
        % region is found
        st_maxchecks = ceil(log2(max(max(abs(ctrl.details.searchTree(:, end-1:end))))))+1;
        
        fprintf('         Search tree: %d separating hyperplanes (depth %d)\n', ...
            st_nhp, st_maxchecks);
        
        % compute maximum number of floating point operations needed to identify
        % and compute the control action
        maxmult = st_maxchecks*nx + nu*nx;
        maxsums = st_maxchecks*(nx + 1) + nu*nx + nu;
        maxcomp = st_maxchecks;
        
        fprintf('Number of operations: %d mutliplications, %d additions, %d comparisons\n', ...
            maxmult, maxsums, maxcomp);
        
    else
        fprintf('         Search tree: none (run "%s=mpt_searchTree(%s)" to compute it)\n', ...
            inpname, inpname);
    end
    
end

fprintf('\n');
