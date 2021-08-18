function [pP,E]= esp(P,ax,testdim)
%
% [pP,E] = esp(P,ax,testdim)
%
% Compute the projection of the polytope P onto axes ax
%
% Inputs:
%   P           : Polytope P = [C D b]
%   ax [1 2]    : Axes to project to
%   testdim [0] : If 1 we test if the polytope is full-dimensional
%                 Can be very slow
%
% Outputs:
%   pP         : Projection of P (pP = [G g])
%   E          : Equality sets of each facet
%
%  2004/04/01
%     Colin Jones, Cambridge Control Laboratory, Cambridge UK
%     cnj22@cam.ac.uk

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if(nargin < 3) testdim = 0; end;

% Verbosity
global verbose
% verbose =  1  - list discovered/to search
% verbose >= 2  - include dual degenecy comments
% verbose >= 10 - include primal degeneracy comments
switch mptOptions.verbose,
    case 0, verbose=0;
    case 1, verbose=0;
    case 2, verbose=2;
end
    
%%verbose = 2;

% Numerical tolerance
global zerotol
zerotol = 1e-6;

% What to do when a dual-degenerate case is encountered
global DUAL_METHOD
DUAL_METHOD = 'max region'; % 'max region' or 'continuous law'

% Initialize hash table
esp_table('init');

% Convert from MPT format and call the main ESP function
if(isempty(P)) pP=[]; E=[]; return; end;
[pP,E] = esp_helper(P,ax,testdim);

if(verbose > 1)
    fprintf('\n');
end;
