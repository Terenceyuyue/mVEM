function [Pn]=mpt_delaunay(P,Options)
%MPT_DELAUNAY Computes the delaunay triangulation of a polytope 
%             (in H or V representation)
%
% [Pn]=mpt_delaunay(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% The delaunay trianglulation partitions a polytopes into simplices.
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P             -   Optional input:
%                   Either a polytope object P or a p times nx matrix which 
%                   contains p points in nx dimensions;
%                   If no argument is passed (P empty), then "mousepoly" is called 
%                   automatically
% Options.plot  -   if set to 1, plots the delanuay triangulation (0 is default)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% Pn            -   Delaunay triangulation

% Copyright is with the following author(s):
%
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%     grieder@control.ee.ethz.ch

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
    mpt_error;
end

if(nargin==0 | isempty(P))
    P=mousepoly;
end


if(isa(P,'polytope'))
    if length(P)>1,
        error('mpt_delaunay: Arrays of polytopes not supported!');
    end
    V=extreme(P);
else
    V = P;
end
if(nargin<2 | ~isfield(Options,'plot'))
    Options.plot=0;
end
if(~isfield(Options,'abs_tol'))
    Options.abs_tol = mptOptions.abs_tol;
end


nx=size(V,2);
nV=size(V,1);

%peturb vertices to obtain general position
V=V+(rand(nV,nx)-0.5)*Options.abs_tol;

G=[];
W=[];
E=[];

%bound lambda_i >0 sum lambda_i =1
G=-eye(nV-1);
G=[G;ones(1,nV-1)];
G=[G zeros(nV,1)]; %add extra variable theta
W=[zeros(nV-1,1);1];
E=zeros(nV,nx);

%Theta >= \sum \lamda_i v_i^2
Gt=[];
for i=1:(nV-1)
    Gt(i)=V(i,:)*V(i,:)'-V(end,:)*V(end,:)';
end
Gt(nV)=-1;
G=[G;Gt];
W=[W;-V(end,:)*V(end,:)'];
E=[E;zeros(1,nx)];

if(0)
    %x=sum_i lambda_i v_i
    slack=1e-5;
    G=[G; V(1:end-1,:)'-repmat(V(end,:)',1,nV-1) zeros(nx,1);-V(1:end-1,:)'+repmat(V(end,:)',1,nV-1) zeros(nx,1)];
    W=[W;-V(end,:)'+slack;V(end,:)'+slack];
    E=[E;eye(nx);-eye(nx)];
else
    %Transform equality and inequality system into system of inequalities
    %in lower dimension via SVD
    Geq = [V(1:end-1,:)'-repmat(V(end,:)',1,nV-1) zeros(nx,1)];
    Weq = -V(end,:)';
    Eeq = eye(nx);
    
    [Uv,Ev,Vv]=svd(Geq);
    nsvd=length(find(abs(diag(Ev))>1e-7));
    
    Ev=Ev(1:nsvd,1:nsvd);
    
    %uhat=inv(Ev)*Uv'*(Weq+Eeq*x0);
    W=W-G*Vv(:,1:nsvd)*inv(Ev)*Uv'*Weq;
    E=E-G*Vv(:,1:nsvd)*inv(Ev)*Uv'*Eeq;
    G=G*Vv(:,nsvd+1:end);
end

Matrices.G=G;
Matrices.E=E;	     
Matrices.W=W;
Matrices.H=zeros(1,size(G,2));
Matrices.H(end)=1;
[Matrices.bndA, Matrices.bndb]=double(unitbox(nx,max(max(abs(V)))*1.5));

[Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options);

if((nx==2 | nx==3) & Options.plot==1)
    plot(P)
    figure
    plot(Pn)
    hold on
    for i=1:size(V,1)
        if(nx==2)
            plot(V(i,1),V(i,2),'k*','LineWidth',2);
        else
            %plot3(V(i,1),V(i,2),V(i,3),'k*','LineWidth',2);
        end
	end
end
