function [ctrlStruct,vol]=mpt_simplexContr(sysStruct,probStruct,Options)
%MPT_SIMPLEXCONTR Computes a piecewise affine feedback law defined over simplices
%
% function [ctrlStruct,vol]=mpt_simplexContr(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Partitions the feasible set into simplical cells and subsequently computes 
% the input sequence for each vertex of each simplex. By linear interpolation, 
% we obtain a piecewise affine feedback law defined for each simplex.
%
%               This function is not based on multi-parametric programming !
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% sysStruct       - System structure in the sysStruct format (see Manual)
% probStruct      - Problem structure in the probStruct format (see Manual)
% Options.verbose - Level of verbosity
% Options.maxCtr  - Maximum number of iterations (default is 1000)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%
% ctrlStruct    - controller structure with the following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - open-loop cost associated to each region (x'Aix + x'Bi + Ci)
%   Pfinal      - Defines the feasible state space partition (i.e. union of
%                 all regions) as Phard.H*x<=Phard.K
%   dynamics    - dynamics active in region Pn(i)
%   details     - contains additional information:
%                 +Volume: volume of polySet
%
%  vol          - the volume of the controller partition
%               

% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% This function was used in "Two Level Model Predictive Control for the Maximum
% Control Invariant Set" by P. Grieder, Z. Wan, M. Kothare and M. Morari published 
% at the American Control Conference 2004, Boston, USA
%
% For details on triangulation see K. Fukudas FAQ page:
% http://www.cs.mcgill.ca/~fukuda/soft/polyfaq/polyfaq.html

% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

error(nargchk(2,3,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<3,
    Options=[];
end
if nargin<2,
    error('mpt_simplexContr: Wrong number of input arguments!');
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'maxCtr')
    Options.maxCtr=1000;
end

if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end
if ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end

if isinf(probStruct.N),
    error('mpt_simplexContr: Prediction horizon must be finite! Set probStruct.N and/or probStruct.subopt_lev = 0');
end

tic; %start timer
nx=size(sysStruct.A,2);
nu=size(sysStruct.B,2);


loopCtr=0;
PfinalOld = polytope;
notconverged=1;


disp('Computing feasible Set...')
[G1,W1,E1]=mpt_constructMatrices(sysStruct,probStruct,Options); 

% compute feasible set via projection
P=polytope([-E1 G1],W1);
Pfinal=projection(P,(1:size(E1,2)));
polySet = Pfinal;

%divide polySet into Simplices
disp('Computing simplices...')
[hhPoly,adj,vvPoly,vol,polySet]=triangulate(polySet);

disp('Computing interpolation...')
%see if each simplex is really full dimensional (i.e. islinear interpolation possible)
Fi={}; Gi={}; hPoly=polytope; minIndex={}; vPoly={}; %initialize to empty
Ffi = cell(1,length(vvPoly));
Gfi = cell(1,length(vvPoly));
for i=1:length(vvPoly)
    %see if Simplex is flat; If so, try subpartitioning
    minEig=-Inf;
    for ind=1:nx+1
        tmpX=repmat(vvPoly{i}(ind,:)',1,nx);
        tmpX=vvPoly{i}([1:ind-1 ind+1:end],:)'-tmpX;
        if(min(abs(eig(tmpX)))>minEig)
            minTmpX=tmpX;
            minEig=min(abs(eig(tmpX)));
            minInd=ind;
        end
    end
    tmpX=minTmpX;
    invX=inv(tmpX);
    Ffi{i}=invX;
    Gfi{i}=-invX*vvPoly{i}(minInd,:)';
    %check if everything is correct
    x0=sum(vvPoly{i},1)'/(nx+1);
    alpha=Ffi{i}*x0+Gfi{i};
    if(min(alpha)>=0 & norm(vvPoly{i}'*[alpha(1:minInd-1);1-sum(alpha);alpha(minInd:end)]-x0)>Options.abs_tol)
        %recompute inverse 
        invX=eye(nx)/tmpX;
        Ffi{i}=invX;
        Gfi{i}=-invX*vvPoly{i}(minInd,:)';
        %check if everything is correct
        alpha=Ffi{i}*x0+Gfi{i};
        if(min(alpha)>=0 & norm(vvPoly{i}'*[alpha(1:minInd-1);1-sum(alpha);alpha(minInd:end)]-x0)>Options.abs_tol)
            %subpartiton simplex
            [hNew,vNew,Fnew,Gnew]=sub_partitionSimplex(vvPoly{i},x0,Options);
        end
    else
        hNew=hhPoly(i);
        vNew{1}=vvPoly{i};
        Fnew{1}=Ffi{i};
        Gnew{1}=Gfi{i};
    end
    for j=1:length(vNew)
        vPoly{end+1}=vNew{j};
        hPoly=[hPoly hNew(j)];
        Fi{end+1}=Fnew{j};
        Gi{end+1}=Gnew{j};
        minIndex{end+1}=minInd;
    end    
end
disp('Computing input sequences...')

% this is important! do not include LQR invariant set
Options.includeLQRset = 0;

%go through each simplex and compute the associated input sequence
[G,W,E,H,F,Y,Cf,Cx,Cc,symmetric,bndA,bndb]=mpt_constructMatrices(sysStruct,probStruct,Options);
U = cell(1, length(vPoly));
for i=1:length(vPoly)
    if(mod(i,round(length(vPoly)/10))==0)
         disp(['Computing input sequences:   ' num2str(i) ' / ' num2str(length(vPoly))]);
    end
    U{i}=[];
    for j=1:(nx+1)
        if(probStruct.norm==2)
            [utmp,ll,hh,exitflag]=mpt_solveQP(H,vPoly{i}(j,:)*F,G,W+E*vPoly{i}(j,:)');
        else
            [utmp,ll,hh,exitflag]=mpt_solveLP(H,G,W+E*vPoly{i}(j,:)',[],[],[],mptOptions.lpsolver);
        end
        
        if(exitflag<=0)
            if(max(G*utmp-W-E*vPoly{i}(j,:)')>=Options.abs_tol)
                plot(polySet);
                hold on
                if(nx>=3)
                    plot3(vPoly{i}(j,1),vPoly{i}(j,2),vPoly{i}(j,3),'b*','LineWidth',3);
                else
                    plot(vPoly{i}(j,1),vPoly{i}(j,2),'b*','LineWidth',3);
                end
                error(['No feasible input associated to vertex ' num2str(vPoly{i}(j,:))]);
            end
        end
        U{i}=[U{i} utmp];%store only first input
    end
end

%go through each simplex and compute the associated feedback law
for i=1:length(vPoly)
    tmpU=repmat(U{i}(:,minIndex{i}),1,nx);
    tmpU=U{i}(:,[1:minIndex{i}-1 minIndex{i}+1:end])-tmpU;
    
    Fi{i}=tmpU*Fi{i};
    Gi{i}=tmpU*Gi{i}+U{i}(:,minIndex{i});
    
    %check if everything is correct
    x0=sum(vPoly{i},1)'/(nx+1);
    u0=Fi{i}*x0+Gi{i};
    if(max(G*u0-W-E*x0)>Options.abs_tol)
        save inter_prob
        error('Error in input interpolation')
    end
end

ctrlStruct.sysStruct=sysStruct;
ctrlStruct.probStruct=probStruct;
ctrlStruct.Pn=hPoly;
ctrlStruct.Fi=Fi;
ctrlStruct.Gi=Gi;
ctrlStruct.Ai=cell(length(hPoly),1);
ctrlStruct.Bi=cell(length(hPoly),1);
ctrlStruct.Ci=cell(length(hPoly),1);
for ii=1:length(Fi),
    ctrlStruct.Ai{ii} = zeros(nx);
    ctrlStruct.Bi{ii} = zeros(1,nx);
    ctrlStruct.Ci{ii} = 0;
end
ctrlStruct.overlaps=0;
ctrlStruct.Pfinal=polySet;
ctrlStruct.details.runTime=toc;
ctrlStruct.details.volume=vol;
ctrlStruct.details.origSysStruct = sysStruct;
ctrlStruct.details.origProbStruct = probStruct;
ctrlStruct.dynamics=ones(length(hPoly),1);
ctrlStruct = mptctrl(ctrlStruct);

return




%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function  [hhPoly,vvPoly,Ffi,Gfi]=sub_partitionSimplex(vPoly,xnew,Options)
minRad=1e-6;

nx=size(vPoly,2);
vert=[vPoly;xnew']; %add new center point
adj = delaunayn(vert);  %compute delaunay triangulation using qhull

%compute V representation of simplices
clear vPoly
for i=1:size(adj,1)
    vPoly{i}=vert(adj(i,:)',:);
end 


%compute H representation of simplices
vvPoly={}; Ffi={}; Gfi={};
hhPoly=polytope;
index=[];
for i=1:size(adj,1)
    tmpPoly=hull(vPoly{i},Options);
    if(isfulldim(tmpPoly))
        [H,K]=double(tmpPoly);
        addpoly=1;
        for j=1:size(vPoly{i},1)
            tmp=find(abs(H*vPoly{i}(j,:)'-K)<=Options.abs_tol); %find intersections with hyperplanes
            if length(tmp)~=nx
                addpoly=0;
                %check if each point is really a vertex
            end
        end
        [x0,rad]=chebyball(tmpPoly); %only consider simplices with "reasonable" size
        if(addpoly & rad>=minRad)
            %see if interpolation is now possible
            minEig=-Inf;
            for ind=1:nx+1
                tmpX=repmat(vPoly{i}(ind,:)',1,nx);
                tmpX=vPoly{i}([1:ind-1 ind+1:end],:)'-tmpX;
                if(min(abs(eig(tmpX)))>minEig)
                    minTmpX=tmpX;
                    minEig=min(abs(eig(tmpX)));
                    minInd=ind;
                end
            end
            tmpX=minTmpX;
            
            Ffi{i}=inv(tmpX);
            Gfi{i}=-inv(tmpX)*vPoly{i}(minInd,:)';
            
            %check if everything is correct
            alpha=Ffi{i}*x0+Gfi{i};
            if(min(alpha)>=0 & norm(vPoly{i}'*[alpha(1:minInd-1);1-sum(alpha);alpha(minInd:end)]-x0)>Options.abs_tol)
            else        
                index=[index i];
                hhPoly=[hhPoly tmpPoly]; %build array of simplices
                vvPoly{length(index)}=vPoly{i};
            end
        end
    end
end