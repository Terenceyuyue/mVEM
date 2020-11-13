function [isin, inwhich, closest] = isinside(P,x0,Options)
%ISINSIDE Checks if a given point lies in the interior of a given polytope
%
% [isin, inwhich, closest] = isinside(P,x0,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Checks if the given point x0 lies in the interior of the polytope P
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                 - Polytopes
% x0                - given point, if not specified: Options.clickx0=1
% Options.abs_tol   - absolute tolerance
% Options.fastbreak - if set to 1 and P is a polyarray, returns only the first
%                     element of P in which x0 lies in
% Options.clickx0   - find the region by clicking in the plot
% Options.newfigure - if a new figure for clicking is desired (default = 1)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% isin      - logical statement indicating if P.H*x0<=P.K was satisfied
% inwhich   - indices of polytopes in P (if P is an array of polytopes) in which
%             the points lies [i.e. P(inwhich) will return all elements of P which
%             contain the given point x0]
% closest   - if the point does not belong to P, index of the polytope which is
%             closest to the point x is returned
%
% see also POLYTOPE
%

% Copyright is with the following author(s):
%
% (C) 2004 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

if ~isa(P, 'polytope')
    error('ISINSIDE: input argument MUST be a polytope');
end

% if ~isstruct(mptOptions),
%     mpt_error;
% end

if nargin<3
    Options=[];
end

if ~isfield(Options,'abs_tol')
    % absolute tolerance
    abs_tol=mptOptions.abs_tol;
else
    abs_tol = Options.abs_tol;
end
if ~isfield(Options,'fastbreak') 
    % if set to 1 and P is a polyarray, returns only the first element of P
    % in which x0 lies in
    Options.fastbreak = 0;
end

if nargout < 2
    % we can always switch to fastbreak=1 if the user doesn't ask for
    % the "inwhich" output argument
    Options.fastbreak = 1;
end

if nargin<2
    if dimension(P)>2,
        error('polytope/isinside: Only 2D partitions are allowed if x0 is not provided!');
    end
    enough=0;
    cnt=0;
    
    X0      = [];
    isin    = [];
    inwhich = cell(0);
    closest = [];
    if ~isfield(Options,'newfigure'),
        Options.newfigure=1;
    end
    if Options.newfigure
      figure;
      plot(P);
    end
    
    title('Click on the figure to specify x_0');
    while ~enough
        [x1,x2,button]=ginput(1);
        if isempty(button) | button~=1,
            return;
        end
        x=[x1;x2];
        cnt = cnt + 1;
        X0  = [X0, x];
        opt = Options;
        opt.clickx0 = 0;
        [i1,i2,i3] = isinside(P,x,opt);
        isin(cnt)=i1;
        inwhich{cnt}=i2;
        closest{cnt}=i3;
        if isin(cnt)
            disp(['point x(' num2str(x') ') is in polytope [' num2str(inwhich{cnt}') '].']);
        else
            disp(['point x(' num2str(x') ') is no polytope. => closest polytope [' num2str(closest{cnt}') ']']);
        end
    end
end


lenP=length(P.Array);
if lenP==0
    H = P.H;
    K = P.K;
    if length(x0)~=size(H,2)
        if ~isfulldim(P)
            isin = 0;
            inwhich = [];
            closest = [];
            return
        end
        error('isinside: incorrect dimension of x0!');
    end

    closest = [];
    if all(H*x0-K<=abs_tol),   % check if point satisfies Hx<=K
        isin=1;
        inwhich=1;
    else
        isin=0;
        inwhich=[];
    end
else
    if length(x0(:))~=size(P.Array{1}.H,2),
        if ~isfulldim(P.Array{1}),
            isin = 0;
            inwhich = [];
            closest = [];
            return
        end
        error('isinside: incorrect dimension of x0!');
    end

    isin = 0;
    closest = [];
    PArray = P.Array;
    fastbreak = Options.fastbreak;
    inwhich = zeros(lenP, 1);
    for ii=1:lenP,                         % cycle through all elements of Pn
        PArray_ii = PArray{ii};
        if all(PArray_ii.H * x0 - PArray_ii.K <= abs_tol),    % for each P, check if Hx<=K holds
            isin=1;
            if fastbreak,
                % exit quickly if required
                inwhich = ii;
                return
            end
            inwhich(ii) = 1;            
        end
    end
    inwhich = find(inwhich);
end
if isin==0 && nargout>2
    if lenP==0
        closest = 1;
    else
        epsilons = [];
        for ii=1:lenP
            Q = P.Array{ii};
            if ~Q.normal
                Q = normalize(Q); 
            end
            epsilons = [epsilons; max(Q.H*x0-Q.K)];
        end
        [~,closest] = min(epsilons);
    end
end
