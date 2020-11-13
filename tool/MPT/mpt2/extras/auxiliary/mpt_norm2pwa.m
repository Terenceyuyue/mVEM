function pwafcn = mpt_norm2pwa(P,l,Options)
% MPT_NORM2PWA  transformes a linear norm into an equivalent PWA fcn representation
%
% pwafcn = mpt_norm2pwa(P,l,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Transformes a linear norm ||P*x||_l into an equivalent PWA function,
%   ||P*x||_l =  pwafcn.Bi{i}*x + pwafcn.Ci{i}  for x \in pwafcn.Pn(i)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P            - scaling matrix in ||P*x||_l
% l            - = 1 or Inf, standard vector norm (l=Inf will be assumed if not
%                set)
% Options      - optional arguments
%   .Pn        - one polytope over which the norm should be defined
%   .method    - if 1 (default) solution is found via enumeration 
%                if 2 solution will be found via YALMIP [should be used for larger problems!]
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% pwafcn       - descriptoin of the PWA function 
%    .Bi,.Ci    
%    .Pn       - polyhedral partition over which the PWA fcn is defined
%    .Pfinal   - domain of the PWA function 
%    .epi      - implicit epigraph description of ||P*x||_l using sdvar.
%                (more efficient than the PWA description for computation)
%                (only existent if Options.method=1)
%
% Note: for simple plotting use plot(pwafcn.epi)
%
%
%
% see also  MPT_DLYAP_INFNORM   MPT_LYAPUNOV


% Copyright is with the following author(s):
%
%(c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%         fjc@control.ee.ethz.ch
%(c) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%         loefberg@control.ee.ethz.ch

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


error(nargchk(1,3,nargin));

global mptOptions

if ~isstruct(mptOptions)
    mpt_error;
end

if nargin < 3
    Options = [];
end
if nargin < 2
    l = Inf;
end

if isempty(l) | (~isinf(l) & l~=1)
    error('only linear norms are supported, i.e. l = 1 or Inf');
end

[mP,nx] = size(P);  

if ~isfield(Options,'Pn')
    Pn = unitbox(nx,mptOptions.infbox);
else
    Pn = Options.Pn;
end
if length(Pn)>1
    error('it can only be defined over one polytope')
end

[H, K] = double(Pn);

if nx~= size(H,2)
    error('the dimension of P and Pn don''t match')
end
    
if ~isfield(Options,'method')
    Options.method = 1; %ENUMERATION
end


if Options.method == 2
    % use a more efficient epi graph way to describe the norm.
    % (more efficient than enumerating the possibilities)
    x = sdpvar(nx,1);
    [pwafcn.epi,pwafcn.Bi,pwafcn.Ci,pwafcn.Pn] = pwa(norm(P*x,l),set(H*x<=K));
    pwafcn.Pfinal = Pn;
    
elseif Options.method == 1
    % using 'enumeration' of the parts
    
        if l==1
            T=[];
            for ii=0:2^mP-1,
                a=dec2bin(ii, mP);
                t=[];
                for jj=1:length(a),
                    t=[t str2num(a(jj))];
                end,
                T=[T; t];
            end
            Tind = find(T(:)==0);
            T(Tind) = -1;
            
            cnt = 0;
            for cc= 1:2^mP
                pp = polytope([-repmat(T(cc,:)',1,nx).*P;H], [zeros(mP,1);K]);
                if isfulldim(pp)
                    cnt = cnt+1;
                    pwafcn.Bi{cnt} = T(cc,:)*P;
                    pwafcn.Ci{cnt} = 0;
                    pwafcn.Pn(cnt) = pp;
                end
            end%(FOR) cc
            pwafcn.Pfinal = Pn;
            
    else % l=Inf
        
        cnt = 0;
        for ii=1:mP
            ind = setdiff(1:mP,ii);
            
            Px = [P(ind,:); -P(ind,:)] - repmat(P(ii,:),(mP-1)*2,1);
            pp = polytope([Px;H],[zeros((mP-1)*2,1);K]);
            if isfulldim(pp)
                cnt = cnt+1;
                pwafcn.Bi{cnt} = P(ii,:);
                pwafcn.Ci{cnt} = 0;
                pwafcn.Pn(cnt) = pp;
            end
        end%(FOR) ii
        
        for ii=1:mP
            ind = setdiff(1:mP,ii);
            
            Px = [P(ind,:); -P(ind,:)] + repmat(P(ii,:),(mP-1)*2,1);
            pp = polytope([Px;H],[zeros((mP-1)*2,1);K]);
            if isfulldim(pp)
                cnt = cnt+1;
                pwafcn.Bi{cnt} = -P(ii,:);
                pwafcn.Ci{cnt} = 0;
                pwafcn.Pn(cnt) = pp;
            end
        end%(FOR) ii
        pwafcn.Pfinal = Pn;
                
    end%(IF) l
    
end% (IF) Options.method

return
