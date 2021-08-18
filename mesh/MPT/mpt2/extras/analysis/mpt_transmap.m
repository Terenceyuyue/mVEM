function [tmap,Pn,ex,explus] = mpt_transmap(Pn, Acell, fcell, Options)
%MPT_TRANSMAP Computes transition map
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Prunes infeasible transitions. Specifically:
%   domain(P2(j), A{i}, f{i}, P1(i)) is NOT feasible if tmap(i, j) is 0
%
% NOTE! If the function complains that extreme point enumeration failed, try to
% change the vertex enumeration method by:
%   mpt_options('extreme_solver', 'matlab')
% or
%   mpt_options('extreme_solver', 'cdd')
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn          - polytope array
% A, f        - cells containing affine map: x(k+1) = A{i}*x + f{i}
% Options.maxsph   - maximum number of separating hyperplnaes to compute
%                    (default is Inf, i.e. no limit)
% Options.targetPn - optional target partition if it is different from Pn
% Options.lpsolver - LP solver to use
% Options.abs_tol  - absolute tolerance
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% tmap        - sparse matrix with non-zero elements on indexes which correspond
%               to transitions whiche are POSSIBLY feasible.
% Pn          - input polyarray P1 with extreme points stored inside for each
%               polytope
% ex          - extreme points of Pn as a cell array
% explus      - extreme points of the affine transformation A*x+f
%

% Copyright is with the following author(s):
%
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          joloef@control.ee.ethz.ch
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

error(nargchk(3,4,nargin));

global mptOptions
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin < 4,
    Options = [];
end
if ~isfield(Options, 'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;
end
if ~isfield(Options, 'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options, 'maxsph'),
    Options.maxsph = Inf;
end

if isfield(Options, 'targetPn'),
    havetarget = 1;
    targetPn = Options.targetPn;
    if ~isa(targetPn, 'polytope'),
        error('Options.targetPn must be a polytope object.');
    end
    lenTarget = length(targetPn);
else
    havetarget = 0;
end

if ~isa(Pn, 'polytope'),
    error('First input must be a polytope object.');
end
if ~iscell(Acell),
    error('Second input must be a cell array.');
end
if ~iscell(fcell),
    error('Third input must be a cell array.');
end
if length(Acell) ~= length(Pn),
    error('Length of Acell must be equal to length of Pn.');
end
if length(fcell) ~= length(Pn),
    error('Length of fcell must be equal to length of Pn.');
end


lpsolver = Options.lpsolver;
abs_tol = Options.abs_tol;
large_tol = 1e3*abs_tol;
maxsph = Options.maxsph;

lenP = length(Pn);
if havetarget,
    tmap = ones(lenP, lenTarget);
else
    tmap = ones(lenP,lenP);
end

% Vertex enumeration
ext = {}; ex = {};
for is=1:lenP
    [E, R, Pn(is)] = extreme(Pn(is), Options);
    % round almost-zero elements
    E(find(abs(E) < 1e-12)) = 0;
    E = E';
    ex{is} = E;
    
    % affine map of extreme points
    Eplus = Acell{is}*E + repmat(fcell{is},1,size(E,2));
    % round almost-zero elements
    Eplus(find(abs(Eplus) < 1e-12)) = 0;
    explus{is} = Eplus;
end
if havetarget,
    for is=1:lenTarget
        E = extreme(targetPn(is), Options);
        % round almost-zero elements
        E(find(abs(E) < 1e-12)) = 0;
        ext{is} = E';
    end
end
if (havetarget & isempty(ext)) | isempty(ex)
    error('Extreme point enumeration failed, see "help mpt_transmap" for details.');
end

n = dimension(Pn);

% compute bounding box of every polytope
BoxMin = cell(1, lenP);
BoxMax = cell(1, lenP);
for i=1:lenP
    oneBoxMin = min(ex{i}')';
    oneBoxMax = max(ex{i}')';
    if ~havetarget,
        BoxMin{i} = oneBoxMin;
        BoxMax{i} = oneBoxMax;
    end
    %now extract all other extreme points of the bounding box
    for j=1:2^n
        index=dec2bin(j-1,n);
        for k=1:n
            if(index(k)=='1')
                boxPoint{i}{j}(k)=oneBoxMax(k);
            else
                boxPoint{i}{j}(k)=oneBoxMin(k);
            end
        end%for k=1:n
        if(size(boxPoint{i}{j},1)<size(boxPoint{i}{j},2))
            boxPoint{i}{j}= boxPoint{i}{j}';
        end
    end%for j=1:2^(n)
end%Pn
if havetarget,
    BoxMin = cell(1, lenTarget);
    BoxMax = cell(1, lenTarget);
    for i = 1:lenTarget,
        BoxMin{i} = min(ext{i}')';
        BoxMax{i} = max(ext{i}')';
    end
end


if havetarget,
    lenP2 = lenTarget;
    % save original ex points, replace them with ext
    ex_orig = ex;
    ex = ext;
else
    lenP2 = lenP;
end
% Initial prune, based on coordinate cuts
% Typically removes 10%
for r = 1:n
    a = eyev(n,r);
    b = 0;
    for prunei = 1:lenP
        iprunes(prunei) = all((a'*explus{prunei} < -large_tol));
    end
    for prunej = 1:lenP2
        jprunes(prunej) = all((a'*ex{prunej} > large_tol));
    end
    for prunei = 1:lenP
        for prunej = 1:lenP2
            if iprunes(prunei)
                if jprunes(prunej)
                    tmap(prunei,prunej) = 0;
                end
            end
        end
    end
end
if havetarget,
    % restore saved points
    ex = ex_orig;
end

if havetarget,
    ex = ext;
end
vecex = [ex{1:end}];
indiciesex = cumsum([1 cellfun('prodofsize',ex)/n]);
vecexplus = [explus{1:end}];
indiciesexplus = cumsum([1 cellfun('prodofsize',explus)/n]);

sphcount = 0;
for i=1:lenP

    %check if polytope is mapped onto a full dimensional polytope or onto
    %a facet (this has an impact on the subsequent reachability analysis)
    if(all(abs(eig(Acell{i}))>Options.abs_tol))
        fullMap=1;
    else
        fullMap=0;
    end
    
    if fullMap,
        lowerCur = min((Acell{i} * [boxPoint{i}{:}] + repmat(fcell{i},1,length(boxPoint{i})))')';
        upperCur = max((Acell{i} * [boxPoint{i}{:}] + repmat(fcell{i},1,length(boxPoint{i})))')';
    end
    
    for j = find(tmap(i,:))

        if tmap(i,j)
            if (fullMap)
                if any((upperCur + large_tol < BoxMin{j}) | (lowerCur - large_tol > BoxMax{j}))
                    tmap(i,j) = 0;
                end
            else
            end
        end

        if tmap(i,j)

            if sphcount > maxsph,
                % exit if maximum number of allowed separating haperplanes has
                % been reached
                return
            end
            [a,b,f] = separatinghp(ex{j},explus{i},lpsolver,n);
            sphcount = sphcount + 1;
            
            if (f == 1) % Found one!
                AA = a*vecexplus;
                BB = a*vecex;
                jprunes  = zeros(lenP,1);
                for prunej = 1:lenP2
                    jprunes(prunej) = all((BB(indiciesex(prunej):indiciesex(prunej+1)-1)-b)<=-large_tol);
                end
                jprunes = logical(jprunes);
                if any(jprunes)
                    iprunes  = zeros(lenP,1);
                    for prunei = i:lenP
                        iprunes(prunei) = all((AA(indiciesexplus(prunei):indiciesexplus(prunei+1)-1)-b)>=large_tol);
                    end

                    for prunei = find(iprunes)
                        tmap(prunei,jprunes) = 0;
                    end
                end
            end

        end
    end
end

%-------------------------------------------------------------------
function [a,b,how] = separatinghp(X,Y,lpsolver,nx)

mx = size(X,2);
my = size(Y,2);

A = [X' -ones(mx,1);-Y' ones(my,1)];
b = [-ones(mx,1);-ones(my,1)];
f = zeros(1,nx+1);
[xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,b,[],[],[],lpsolver);
a = xopt(1:nx)';
b = xopt(end);
how = ~isequal(how,'infeasible');
how = how & ~all(abs(xopt)<1e-8) & ~any(xopt==1e9);
