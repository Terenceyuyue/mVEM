function ctrl = mpt_lyapunov(ctrl, funtype, ndeg, Options)
%MPT_LYAPUNOV Computes a Lapunov-type function for a given explicit controller
%
%   ctrl = mpt_lyapunov(ctrl, functiontype)
%   ctrl = mpt_lyapunov(ctrl, functiontype, ndeg)
%   ctrl = mpt_lyapunov(ctrl, functiontype, Options)
%   ctrl = mpt_lyapunov(ctrl, functiontype, ndeg, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes a Lyapunov function for a given explicit controller (if one exists).
%
% Three types of Lyapunov functions are supported:
%   * Common Quadratic Lyapunov Function       (funtype = 'quad')
%   * Common Sum of Squares Lyapunov Function  (funtype = 'sos')
%   * Piecewise-Affine Lyapunov Function       (funtype = 'pwa')
%   * Piecewise-Quadratic Lyapunov Function    (funtype = 'pwq')
%   * Piecewise-Polynomial Lyapunov Function   (funtype = 'pwp')
%
% NOTE: This function is NOT available for on-line controllers.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl     - explicit controller (an MPTCTRL object)
% funtype  - type of lyapunov function to compute. valid options are:
%               'quad' - quadratic lyapunov function
%               'sos'  - sum of squares lyapunov function
%               'pwa'  - piecewise-affine lyapunov function
%               'pwq'  - piecewise-quadratic lyapunov function
%               'pwp'  - piecewise-polynomial lyapunov function
% ndeg     - degree of piecewise-polynomial lyapunov function(should be even)
% Options  - additional options (see help mpt_getQuadLyapFct, mpt_getPWALyapFct,
%            mpt_getCommonSOSLyapFct and mpt_getPWQLyapFct for more details)
%
% Options.useTmap      - If set to true (default), transition map will
%                        be computed to rule out certain transitions
% Options.sphratio     - Gives factor which governs maximum number of separating
%                        hyperplanes computed in transition maps. Number of
%                        separating  hyperplnaes computed at each step is given
%                        by length(Pn)^2 / Options.ratio
%                        Default value is 20.
%                        Set this option to 0 if you don't want to impose any
%                        limit on number of separating hyperplanes. The higher
%                        the value of this parameter is, the less separating
%                        hyperplanes will be computed, resulting in (possibly)
%                        less efficiency.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrl     - updated controller with Lyapunov function stored in
%              ctrl.details.lyapunov
%
% see also MPT_GETQUADLYAPFCT, MPT_GETPWQLYAPFCT, MPT_GETPWALYAPFCT
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

error(nargchk(2,4,nargin));

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

if isstruct(ctrl) & mpt_issysstruct(ctrl),
    % convert the system structure into a dummy controller
    try
        ctrl = mptctrl(ctrl);
    catch
        rethrow(lasterror);
    end
        
elseif ~isa(ctrl, 'mptctrl')
    error('MPT_LYAPUNOV: First input argument must be a valid MPTCTRL object!');
end

if nargin==3,
    if isstruct(ndeg),
        Options = ndeg;
        ndeg = [];
    end
end
    
if nargin<3,
    Options = [];
    ndeg = [];
end
if ~ischar(funtype)
    error('MPT_LYAPUNOV: First input argument must be a string!');
end
if ~isexplicit(ctrl)
    error('MPT_LYAPUNOV: Cannot compute a Lyapunov function for on-line controllers!');
end
if ctrl.overlaps,
    error('MPT_LYAPUNOV: Controller regions overlap, use mpt_removeOverlaps first!');
end

if isstruct(ndeg),
    Options = ndeg;
    ndeg = [];
end

if ~exist('Options', 'var'),
    Options = [];
elseif ~isempty(Options),
    if ~isstruct(Options)
        error('MPT_LYAPUNOV: Third input argument must be a structure!');
    end
end

if ~isfield(Options, 'verbose')
    Options.verbose = mptOptions.verbose;
end

% transform the MPTCTRL object to structure
ctrl = struct(ctrl);

lyapunov.feasible = 0;

if strcmpi(funtype, 'quad')
    % common quadratic lyapunov function
    [lyapP,decay,feasible]=mpt_getQuadLyapFct(ctrl,Options);
    if feasible,
        lyapunov.type = 'quadratic';
        lyapunov.P = lyapP;
        lyapunov.decay = decay;
        lyapunov.feasible = 1;
    end

elseif strcmpi(funtype, 'pwq')
    % PWQ lyapunov function
    [lyapQ,lyapL,lyapC,feasible,decay]=mpt_getPWQLyapFct(ctrl,Options);
    if feasible,
        lyapunov.type = 'pwq';
        lyapunov.Q = lyapQ;
        lyapunov.L = lyapL;
        lyapunov.C = lyapC;
        lyapunov.decay = decay;
        lyapunov.feasible = 1;
    end
    
elseif strcmpi(funtype, 'pwa')
    % PWA lyapunov function
    % how many regions contain the origin?
    dim = dimension(ctrl.Pn(1));
    [isin, contain_origin] = isinside(ctrl.Pn, zeros(dim, 1));
    if (ctrl.probStruct.norm==2) && length(contain_origin)~=2^dim
        % no chance of finding a PWA lyap function for 2-norm problems
        fprintf('\nThere is no chance of finding a PWA Lyapunov function for a controller computed in 2-norm.\n');
        fprintf('Use mpt_lyapunov(%s, ''pwq'') or mpt_lyapunov(%s, ''quad'') instead.\n\n', inputname(1), inputname(1));
        error('MPT_LYAPUNOV: cannot continue, see message above.');
    end

    [lyapL,lyapC,feasible,decay]=mpt_getPWALyapFct(ctrl,Options);
    if feasible,
        lyapunov.type = 'pwa';
        lyapunov.L = lyapL;
        lyapunov.C = lyapC;
        lyapunov.decay = decay;
        lyapunov.feasible = 1;
    end

elseif strcmpi(funtype, 'pwp') | strcmpi(funtype, 'sos'),
    if isempty(ndeg),
        if strcmpi(funtype, 'pwp'),
            fprintf('No degree of piecewise-polynomail function specified, assuming ndeg=4\n\n');
        else
            fprintf('No degree of sum of squares function specified, assuming ndeg=4\n\n');
        end
        ndeg = 4;
    end
    if ndeg<2,
        if strcmpi(funtype, 'pwp'),
            error('"MPT_LYAPUNOV: Degree of piecewise-polynomial function must be greater 2!');
        else
            error('"MPT_LYAPUNOV: Degree of sum of squares function must be greater 2!');
        end
    end
    if mod(ndeg, 2)~=0,
        if strcmpi(funtype, 'pwp'),
            error('"MPT_LYAPUNOV: Degree of piecewise-polynomial function must be an even number!');
        else
            error('"MPT_LYAPUNOV: Degree of sum of squares function must be an even number!');
        end
    end
    
    if strcmpi(funtype, 'pwp'),
        [dQ,feasible,decay]=mpt_getPWPLyapFct(ctrl,ndeg,Options);
        if feasible,
            lyapunov.type = 'pwp';
            lyapunov.Q = dQ;
            lyapunov.decay = decay;
            lyapunov.feasible = 1;
        end
    else
        solution = mpt_getCommonSOSLyapFct(ctrl, ndeg, Options);
        if ~isempty(solution),
            if solution.found,
                lyapunov.type = 'sos';
                lyapunov.V = solution.V;
                lyapunov.feasible = 1;
                lyapunov.details = solution.details;
            end
        end
    end
    
else
    error(['MPT_LYAPUNOV: unknown Lyapunov type ''' funtype '''! Valid options are ''quad'', ''pwq'' and ''pwa''.']);
end

if isfield(ctrl.details, 'lyapunov')
    % remove any previously stored information
    % TODO: we might want to keep it and just copy it somewhere else instead
    ctrl.details = rmfield(ctrl.details, 'lyapunov');
end
ctrl.details.lyapunov = lyapunov;

%store information back to controller object
ctrl = mptctrl(ctrl);

% assign variable with new '.details.lyapunov' field in caller's workspace
if ~isempty(inputname(1)) & nargout==0,
    assignin('caller',inputname(1),ctrl);
end

lyap_type = upper(funtype);
if isequal(lyap_type, 'QUAD'),
    lyap_type = 'Quadratic';
end

if lyapunov.feasible
    fprintf('\n%s Lyapunov function found, closed-loop system is stable.\n', lyap_type);
    if isempty(inputname(1)) | nargout>0,
        fprintf('The Lyapunov function was stored to ctrl.details.lyapunov\n\n');
    else
        fprintf('The Lyapunov function was stored to %s.details.lyapunov\n\n', inputname(1));
    end
    
else
    fprintf('\n%s Lyapunov function NOT found, closed-loop system can be unstable.\n', lyap_type);
end
