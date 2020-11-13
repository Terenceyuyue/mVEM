function status = eq(P,Q,Options)
%EQ Checks if two polytopes are equal
%
% status = eq(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = EQ(P,Q) returns TRUE (1) if P==Q
%
% USAGE:
%   P==Q
%   eq(P,Q)
%   eq(P,Q,Options)
%
% NOTE: If P and Q are empty polytopes, returns 1 (true)
%       If one of them is empty and second one is not, returns 0 (false)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q              - Polytopes or polyarrays
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance: the larger abs_tol, the more likely that "true" 
%                    will be returned.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement (1 if P==Q, 0 otherwise)
%
% see also NE, LE, GE
%

% Copyright is with the following author(s):
%
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

if ~isa(P, 'polytope') | ~isa(Q, 'polytope')
  error('EQ: Argument MUST be a polytope object');
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3
    Options=[];
end

if ~isfield(Options,'rel_tol')
    Options.rel_tol=mptOptions.rel_tol;    % relative tolerance
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'elementwise')
    Options.elementwise=0;
end

lenP=length(P.Array);
lenQ=length(Q.Array);

isfulldimP = isfulldim(P);
isfulldimQ = isfulldim(Q);
if ~isfulldimP & ~isfulldimQ,
    % both polytopes empty, return TRUE
    status = 1;
    return
end

if ~isfulldimP | ~isfulldimQ,
    % one of the polytopes empty, return FALSE
    status = 0;
    return
end

abs_tol = Options.abs_tol;

if Options.elementwise,
    % we compare P and Q elementwise, i.e. P(1)==Q(1), P(2)==Q(1), ..., P(n)==Q(1);
    % P(1)==Q(2), P(2)==Q(2), ..., P(n)==Q(2); ... P(1)==Q(m), ..., P(n)==Q(m)
    if lenP>0,
        if lenQ==0,
            lenQ=1;
        end
        status=zeros(lenP,lenQ);
        for ii=1:lenP,
            status(ii,:)=(eq(P.Array{ii},Q,Options));    % recursive call
        end
        return
    end
    if lenQ>0,
        status=zeros(1,lenQ);
        for ii=1:lenQ,
            status(ii)=any(eq(P,Q.Array{ii},Options));   % recursive call
        end
        return
    end
else
    % we compare if P==Q, i.e. if P is fully covered with Q and vice versa
    if lenP>0 | lenQ>0,
        % check dimensions
        dimP = dimension(P);
        dimQ = dimension(Q);
        if dimP==0 & dimQ==0,
            % both polyarrays are 0-dimensional, thus equal
            status = 1;
            return
        elseif dimP==0 | dimQ==0,
            % only one polyarray is 0-dimensional, thus not equal
            status = 0;
            return
        elseif dimP ~= dimQ,
            error('EQ: Only polytopes of equal dimensionality can be compared');
        end

        % try to rule out some cases based on bounding boxes
        bboxOpt.noPolyOutput = 1;    % tell bounding_box() not to create a polytope object
        [R, Plow, Pup] = bounding_box(P, bboxOpt);
        [R, Qlow, Qup] = bounding_box(Q, bboxOpt);
        
        if ~isempty(Plow) & ~isempty(Qlow),
            bbox_tol = 1e3*abs_tol;
            if any(abs(Plow - Qlow) > bbox_tol) | any(abs(Pup - Qup) > bbox_tol),
                % bounding boxes differ by more than abs_tol => polytopes cannot be equal
                status = 0;
                return
            end
            % we cannot reach any conclusion based solely on the fact that bounding
            % boxes are identical, therefore we continue...
        end
        
        Options.simplecheck=1;   % to allow premature break of recursion in mldivide
        if lenP>0,
            myOnes=ones(size(P.Array{1}.H,2),1);
            for i=1:lenP
                P.Array{i}.K=P.Array{i}.K+abs(P.Array{i}.H*myOnes*abs_tol);
            end
        end
        status = ~isfulldim(mldivide(Q,P,Options));
        if(~status)
            return
        end
        for i=1:length(P.Array)
            P.Array{i}.K=P.Array{i}.K-abs(P.Array{i}.H*myOnes*abs_tol);
        end
        
        if lenQ>0,
            myOnes=ones(size(Q.Array{1}.H,2),1);
            for i=1:lenQ
                Q.Array{i}.K=Q.Array{i}.K+abs(Q.Array{i}.H*myOnes*abs_tol);
            end
        end
        status = ~isfulldim(mldivide(P,Q,Options));
        return
    end
end

if ~P.minrep
    P=reduce(P,Options);
end
if ~Q.minrep
    Q=reduce(Q,Options);
end
[ncP,nxP]=size(P.H);
[ncQ,nxQ]=size(Q.H);
if nxP~=nxQ
    error('EQ: Only polytopes of equal dimensionality can be compared');
end
if ncP ~= ncQ,
    status = 0;
    return
end

Pbbox = P.bbox;
Qbbox = Q.bbox;
if ~isempty(Pbbox) & ~isempty(Qbbox),
    bbox_tol = 1e4*abs_tol;
    if any(abs(Pbbox(:,1) - Qbbox(:,1)) > bbox_tol) | any(abs(Pbbox(:,2) - Qbbox(:,2)) > bbox_tol),
        % bounding boxes differ by more than abs_tol => polytopes cannot be equal
        status = 0;
        return
    end
    % we cannot reach any conclusion based solely on the fact that bounding
    % boxes are identical, therefore we continue...
end

status=1;
PAB=[P.H P.K];
QAB=[Q.H Q.K];

for ii=1:ncP
    %if all(sum(abs(QAB-repmat(PAB(ii,:),ncQ,1)),2)>abs_tol)
    %if all(sum(abs(QAB-PAB(ones(ncQ,1)*ii,:)),2)>abs_tol)
    Z = sum(abs(QAB-PAB(ones(ncQ,1)*ii,:)),2);
    if all(Z>abs_tol)
        status=0;
        return;
    end
end

for ii=1:ncQ
    %if all(sum(abs(PAB-repmat(QAB(ii,:),ncP,1)),2)>abs_tol)
    %if all(sum(abs(PAB-QAB(ones(ncP,1)*ii,:)),2)>abs_tol)
    Z = sum(abs(PAB-QAB(ones(ncP,1)*ii,:)),2);
    if all(Z>abs_tol)
        status=0;
        return;
    end
end
